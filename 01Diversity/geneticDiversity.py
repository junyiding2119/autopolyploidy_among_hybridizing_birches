cd ./ScanTools
conda activate gatk
python3
import ScanTools

albo = ScanTools.scantools('/data/home/ljliu/data_analysis/DingJunYi/costataeEvolution/Result/10genetic_diversity/albo')

#check
albo.pops
albo.samps

albo.splitVCFs("/data/home/ljliu/data_analysis/DingJunYi/costataeEvolution/Result/10genetic_diversity/all_species_vcf", ref_path="/data/home/ljliu/data_analysis/reference/Betula_pendula/35080Betula_pendula.fa",repolarization_key="-99",overwrite=True)

albo.calcwpm(albo.split_dirs[0],window_size=10000,min_snps=10,pops = ['BA','BMX','BYS','DTX','HHG','HP','LHS','LS','LYC','MYL','PQG','QLL','RTX','SND','SWP','TBS','TBX','TTH','TZZ','XYB'],sampind=5,use_repol=False)

albo.calcPairwisebpm(albo.split_dirs[0],window_size=10000,min_snps=10,pops=["DTX","HHG","HP","LYC","QLL","SWP","TBX","TZZ","XYB"],use_repol=False)


#ash
ash = ScanTools.scantools('/data/home/ljliu/data_analysis/DingJunYi/costataeEvolution/Result/10genetic_diversity/ash')
#check
ash.pops
ash.samps

ash.splitVCFs("/data/home/ljliu/data_analysis/DingJunYi/costataeEvolution/Result/10genetic_diversity/all_species_vcf", ref_path="/data/home/ljliu/data_analysis/reference/Betula_pendula/35080Betula_pendula.fa",repolarization_key="-99",overwrite=True)

ash.calcwpm(ash.split_dirs[0],window_size=10000,min_snps=10,pops = ["BX","BYS","CKX","DBS","DTX","EDL","FC","HHG","HP","MYL","NSX","SND","SNJ","TBS","TZZ","WLP","XMX"],sampind=5,use_repol=False)

ash.calcPairwisebpm(ash.split_dirs[0],window_size=10000,min_snps=10,pops=["DBS","DTX","EDL","FC","HHG","HP","SND","TZZ","WLP","XMX"],use_repol=False)


#bugg
bug = ScanTools.scantools('/data/home/ljliu/data_analysis/DingJunYi/costataeEvolution/Result/10genetic_diversity/bug')
#check
bug.pops
bug.samps

bug.splitVCFs("/data/home/ljliu/data_analysis/DingJunYi/costataeEvolution/Result/10genetic_diversity/all_species_vcf", ref_path="/data/home/ljliu/data_analysis/reference/Betula_pendula/35080Betula_pendula.fa",repolarization_key="-99",overwrite=True)

bug.calcwpm(bug.split_dirs[0],window_size=10000,min_snps=10,pops = ['CKX', 'EDL', 'NSX', 'SNJ', 'WLP'],sampind=5,use_repol=False)

bug.calcPairwisebpm(bug.split_dirs[0],window_size=10000,min_snps=10,pops=['CKX', 'EDL', 'NSX', 'SNJ'],use_repol=False)


#cos
cos = ScanTools.scantools('/data/home/ljliu/data_analysis/DingJunYi/costataeEvolution/Result/10genetic_diversity/cos')
#check
cos.pops
cos.samps

cos.splitVCFs("/data/home/ljliu/data_analysis/DingJunYi/costataeEvolution/Result/10genetic_diversity/all_species_vcf", ref_path="/data/home/ljliu/data_analysis/reference/Betula_pendula/35080Betula_pendula.fa",repolarization_key="-99",overwrite=True)

cos.calcwpm(cos.split_dirs[0],window_size=10000,min_snps=10,pops = ['BJF', 'EDG', 'FHS', 'HR', 'HTY', 'JYX', 'KCL', 'LLZ', 'LYL', 'QF', 'QSLC', 'WTG', 'XF', 'YCS'],sampind=5,use_repol=False)

cos.calcPairwisebpm(cos.split_dirs[0],window_size=10000,min_snps=10,pops=["BJF","EDG","HR","JYX","KCL","LLZ","LYL","QF","QSLC","WTG","XF","YCS"],use_repol=False)


#erman
erman = ScanTools.scantools('/data/home/ljliu/data_analysis/DingJunYi/costataeEvolution/Result/10genetic_diversity/erman')
#check
erman.pops
erman.samps

erman.splitVCFs("/data/home/ljliu/data_analysis/DingJunYi/costataeEvolution/Result/10genetic_diversity/all_species_vcf", ref_path="/data/home/ljliu/data_analysis/reference/Betula_pendula/35080Betula_pendula.fa",repolarization_key="-99",overwrite=True)

erman.calcwpm(erman.split_dirs[0],window_size=10000,min_snps=10,pops = ['BJF', 'CBS', 'FHSJ', 'LLZ', 'LYL', 'XLA'],sampind=5,use_repol=False)

erman.calcPairwisebpm(erman.split_dirs[0],window_size=10000,min_snps=10,pops=['BJF', 'CBS', 'FHSJ', 'LLZ', 'LYL'],use_repol=False)

#uti
uti = ScanTools.scantools('/data/home/ljliu/data_analysis/DingJunYi/costataeEvolution/Result/10genetic_diversity/uti')
#check
uti.pops
uti.samps

uti.splitVCFs("/data/home/ljliu/data_analysis/DingJunYi/costataeEvolution/Result/10genetic_diversity/all_species_vcf", ref_path="/data/home/ljliu/data_analysis/reference/Betula_pendula/35080Betula_pendula.fa",repolarization_key="-99",overwrite=True)

uti.calcwpm(uti.split_dirs[0],window_size=10000,min_snps=10,pops = ['BX', 'DCX', 'DQ', 'EBY', 'JD', 'JLX', 'JS', 'LFG', 'LJS', 'LP', 'MEK', 'SDX', 'SLS', 'WZ', 'YJA'],sampind=5,use_repol=False)

uti.calcPairwisebpm(uti.split_dirs[0],window_size=10000,min_snps=10,pops=['DCX', 'DQ', 'WZ'],use_repol=False)


all = ScanTools.scantools('/data/home/ljliu/data_analysis/DingJunYi/costataeEvolution/Result/10genetic_diversity/all_species_vcf')
#check
all.pops
all.samps

all.splitVCFs("/data/home/ljliu/data_analysis/DingJunYi/costataeEvolution/Result/10genetic_diversity/all_species_vcf", ref_path="/data/home/ljliu/data_analysis/reference/Betula_pendula/35080Betula_pendula.fa",repolarization_key="-99",overwrite=True)

all.calcPairwisebpm(all.split_dirs[0],window_size=10000,min_snps=10,pops=['albo', 'ash', 'bug', 'cos', 'erman', 'uti'],use_repol=False)