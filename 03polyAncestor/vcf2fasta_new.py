import vcf
import os
import sys
import re
from collections import Counter

# --- Helper Functions & Constants ---
# IUPAC ambiguity codes mapping: tuple of sorted unique bases to IUPAC character
# This map is used by convert_ambigous_base
SIMPLE_IUPAC_MAP = {
    ('A',): 'A', ('C',): 'C', ('G',): 'G', ('T',): 'T',
    ('A', 'G'): 'R', ('C', 'T'): 'Y', ('A', 'C'): 'M',
    ('G', 'T'): 'K', ('C', 'G'): 'S', ('A', 'T'): 'W',
    ('C', 'G', 'T'): 'B', ('A', 'G', 'T'): 'D',
    ('A', 'C', 'T'): 'H', ('A', 'C', 'G'): 'V',
    ('A', 'C', 'G', 'T'): 'N'
}

def convert_ambigous_base(gt_bases_str):
    if gt_bases_str is None or \
       gt_bases_str in {'*/*', '*|*', '*/*/*/*', 'N/N', 'N/N/N/N', './.'}:
        return 'N'
    cleaned_gt_str = gt_bases_str.replace('*', '')
    alleles_raw = re.split(r'/|\|', cleaned_gt_str)
    valid_dna_bases = {'A', 'C', 'G', 'T'}
    alleles_cleaned = [allele for allele in alleles_raw if allele and allele in valid_dna_bases]
    if not alleles_cleaned:
        return 'N'
    unique_alleles = tuple(sorted(list(set(alleles_cleaned))))
    return SIMPLE_IUPAC_MAP.get(unique_alleles, 'N')

def record_region(out_file_handle, chrom, start_pos, end_pos):
    out_file_handle.write(f"{chrom}.{start_pos}.{end_pos}\n")


def write_fasta_segment(output_dir, segment_chrom, segment_start_pos, segment_end_pos, 
                        seq_data, sample_for_len_check, min_len=1000):
    if not seq_data or not sample_for_len_check or sample_for_len_check not in seq_data \
       or not seq_data[sample_for_len_check]:
        return

    if len(seq_data[sample_for_len_check]) >= min_len:
        fasta_filename = os.path.join(output_dir, f"{segment_chrom}.{segment_start_pos}.{segment_end_pos}.fa")
        with open(fasta_filename, "w") as t_ouf:
            for sample_name, sequence in seq_data.items():
                if sequence:
                    t_ouf.write(f">{sample_name}\n{sequence}\n")

# --- Main Processing Function ---

def process_vcf_to_regions_and_fasta(vcf_reader, regions_output_filepath, fasta_output_dir,
                                     gap_threshold=500, fasta_min_len=1000):
    region_current_chrom = None
    region_current_start_pos = 0
    region_prev_pos = 0
    first_locus_for_regions = True

    fasta_segment_current_chrom = None
    fasta_segment_start_genomic_pos = 0
    fasta_prev_pos = 0  # This will serve as the END_POS for the segment being written
    seq_dic = {sample_name: "" for sample_name in vcf_reader.samples}
    first_locus_for_fasta = True

    sample_name_for_length_check = None
    if vcf_reader.samples:
        sample_name_for_length_check = vcf_reader.samples[0]
    else:
        print("Warning: No samples found in VCF. FASTA file generation will be skipped.", file=sys.stderr)

    with open(regions_output_filepath, 'w') as regions_out_file_handle:
        regions_out_file_handle.write("CHROM.START.END\n")

        for loci in vcf_reader:
            raw_chrom_name = loci.CHROM
            if '_' in raw_chrom_name:
                try:
                    current_chrom_name = raw_chrom_name.split('_')[1]
                except IndexError:
                    current_chrom_name = raw_chrom_name.split('_')[-1]
            else:
                current_chrom_name = raw_chrom_name
            
            current_pos = loci.POS

            # --- Region Segmentation Logic ---
            if first_locus_for_regions:
                region_current_chrom = current_chrom_name
                region_current_start_pos = current_pos
                region_prev_pos = current_pos
                first_locus_for_regions = False
            else:
                is_new_region = False
                if current_chrom_name != region_current_chrom:
                    is_new_region = True
                elif (current_pos - region_prev_pos) >= gap_threshold:
                    is_new_region = True
                
                if is_new_region:
                    record_region(regions_out_file_handle, region_current_chrom, region_current_start_pos, region_prev_pos)
                    region_current_chrom = current_chrom_name
                    region_current_start_pos = current_pos
            region_prev_pos = current_pos

            # --- FASTA Segmentation Logic ---
            if not vcf_reader.samples:
                if first_locus_for_fasta: first_locus_for_fasta = False
                # fasta_prev_pos still needs to be updated for region logic if used elsewhere,
                # but here it's mainly for FASTA. If no samples, current_pos won't be used
                # as fasta_prev_pos for writing a fasta segment.
                # However, ensure fasta_prev_pos is updated at the end of the loop iteration
                # if it were to be used for other logic.
                # For now, just skip the FASTA specific part.
                # The fasta_prev_pos = current_pos at the end of loop will handle this for next iteration.
                pass # Skip FASTA specific processing
            elif first_locus_for_fasta: # This 'elif' ensures this block only runs if vcf_reader.samples is true
                fasta_segment_current_chrom = current_chrom_name
                fasta_segment_start_genomic_pos = current_pos
                # fasta_prev_pos will be updated at the end of this iteration
                for sample_obj in loci.samples:
                    seq_dic[sample_obj.sample] = convert_ambigous_base(sample_obj.gt_bases)
                first_locus_for_fasta = False
            else: # Not the first locus, and samples exist
                is_new_fasta_segment = False
                if current_chrom_name != fasta_segment_current_chrom:
                    is_new_fasta_segment = True
                elif (current_pos - fasta_prev_pos) >= gap_threshold: # fasta_prev_pos is from previous locus
                    is_new_fasta_segment = True

                if is_new_fasta_segment:
                    # Write the completed FASTA segment.
                    # The segment ended at fasta_prev_pos (position of the previous locus).
                    # MODIFICATION: Pass fasta_prev_pos as segment_end_pos
                    write_fasta_segment(fasta_output_dir, fasta_segment_current_chrom, 
                                        fasta_segment_start_genomic_pos, fasta_prev_pos, 
                                        seq_dic, 
                                        sample_name_for_length_check, fasta_min_len)
                    
                    fasta_segment_current_chrom = current_chrom_name
                    fasta_segment_start_genomic_pos = current_pos
                    seq_dic = {sample_name: "" for sample_name in vcf_reader.samples} 
                    for sample_obj in loci.samples:
                        seq_dic[sample_obj.sample] = convert_ambigous_base(sample_obj.gt_bases)
                else:
                    for sample_obj in loci.samples:
                        if sample_obj.sample not in seq_dic:
                             seq_dic[sample_obj.sample] = "" 
                        seq_dic[sample_obj.sample] += convert_ambigous_base(sample_obj.gt_bases)
            
            # This assignment is crucial: fasta_prev_pos always stores the position
            # of the locus whose bases have just been processed (or attempted).
            fasta_prev_pos = current_pos
            
        # --- After the loop: Write the last pending region and FASTA segment ---
        if not first_locus_for_regions:
            record_region(regions_out_file_handle, region_current_chrom, region_current_start_pos, region_prev_pos)

        if not first_locus_for_fasta and vcf_reader.samples:
            # The last segment ends at fasta_prev_pos (position of the very last locus processed)
            # MODIFICATION: Pass fasta_prev_pos as segment_end_pos
            write_fasta_segment(fasta_output_dir, fasta_segment_current_chrom, 
                                fasta_segment_start_genomic_pos, fasta_prev_pos,
                                seq_dic, 
                                sample_name_for_length_check, fasta_min_len)

# --- Main Execution Block ---
if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python combined_script.py <vcf_filepath> <regions_output_filepath> <fasta_output_directory>", file=sys.stderr)
        sys.exit(1)

    vcf_filepath = sys.argv[1]
    regions_output_filepath = sys.argv[2]
    fasta_output_dir = sys.argv[3]

    if not os.path.exists(fasta_output_dir):
        try:
            os.makedirs(fasta_output_dir)
        except OSError as e:
            print(f"Error: Could not create FASTA output directory {fasta_output_dir}: {e}", file=sys.stderr)
            sys.exit(1)
    
    vcf_reader_obj = None
    try:
        vcf_reader_obj = vcf.Reader(filename=vcf_filepath)
    except FileNotFoundError:
        print(f"Error: VCF file not found at {vcf_filepath}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error opening or parsing VCF file {vcf_filepath}: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Processing VCF file: {vcf_filepath}")
    process_vcf_to_regions_and_fasta(
        vcf_reader_obj,
        regions_output_filepath,
        fasta_output_dir
    )

    print(f"Regions extracted and saved to: {regions_output_filepath}")
    print(f"FASTA files generated in directory: {fasta_output_dir}")