# The paradoxical prevalence of autopolyploidy among highly hybridizing birches 

A bioinformatics pipeline for studying the evolution and polyploid origins of birch species (*Betula* sect. *Costatae*) in East Asia.

## Overview

This repository contains scripts and analysis pipelines for investigating:
- **Population genetics** of multiple *Betula* species (diploid and tetraploid) using ScanTools
- **Polyploid Polarization and Phylogenetic Analysis** for tetraploid species using PolyAncestor pipeline
- **Demographic history** using fastsimcoal
- **Niche comparison and species distribution** 

## Pipeline Structure

```
├── 00CallingSNP/          # Variant calling pipeline (GATK)
├── 01Diversity/           # Genetic diversity analysis
├── 02Entropy/             # Population structure analysis
├── 03polyAncestor/        # PolyAncestor pipeline and ABC-SimAnnealing
├── 04polyAncestor-Model_evaluation/  # Model validation
├── 05Demographic_modeling_fsc/       # Demographic modeling (fastsimcoal2)
├── 06ENM_Niche_Comparison/           # Ecological niche modeling and niche comparison
└── 07plot/                # Visualization scripts
```

## Module Descriptions

### 1. SNP Calling (`00CallingSNP/`)
Complete variant calling workflow using GATK practices:
- **Quality control**: Trimmomatic for adapter trimming and quality filtering
- **Alignment**: BWA-MEM against *Betula pendula* reference genome
- **Duplicate marking**: GATK MarkDuplicates
- **Variant calling**: GATK HaplotypeCaller (supports both diploid and tetraploid)
- **Joint genotyping**: GenomicsDBImport + GenotypeGVCFs

### 2. Genetic Diversity (`01Diversity/`)
Population-level diversity statistics using ScanTools (https://github.com/pmonnahan/ScanTools):
- Nucleotide diversity (π)
- Fixation index (Fst)
- Absolute divergence (Dxy)
- ρ statistic
- Linkage disequilibrium (LD)

### 3. Population Structure (`02Entropy/`)
Ancestry inference using Entropy software:
- Bayesian clustering for mixed-ploidy populations
- Admixture proportion estimation
- Visualization of population structure
- Randomly subsampled two alleles per each tetraploid locus

### 4. Polyploid Ancestry Analysis (`03polyAncestor/`)
Novel ABC-Simulated Annealing approach to infer polyploid origins:
- Generation of MSAs
- Polyploid polarization using PolyAncestor pipeline (https://github.com/LLN273/PolyAncestor)
- Polyploid model testing following simulated annealing and approximate Bayesian computation (ABC) algorithm following (https://github.com/LLN273/Complex-Polyploids)

### 5. Model Evaluation (`04polyAncestor-Model_evaluation/`)
Validation of polyploid ancestry models:
- L2-norm distance computation
- Generated posterior distributions for each parameter

### 6. Demographic Modeling (`05Demographic_modeling_fsc/`)
Population demographic history using fastsimcoal2:
- Fastsimcoal2 pipeline
- Bootstrap confidence intervals
- Goodness-of-fit testing

### 7. Ecological Niche Modeling (`06ENM_Niche_Comparison/`)
Species distribution modeling and niche comparison:
- GBIF occurrence data retrieval
- Environmental variable selection (World bioclim + soil)
- SDM model
- Niche overlap and equivalency tests

### 8. Visualization (`07plot/`)
- Visualization

## License

This project is released under the MIT License.

