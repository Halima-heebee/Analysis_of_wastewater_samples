#VARSCAN
import os

## Exeter ONT NT001
# directory = "/home/mbxha18/bigzipfile/exeter_results_NT001/results/10945/20230208_1539_1A_PAM81850_9bee2fb9/ill_match/"
# path = [d for d in os.listdir(f'{directory}') if os.path.isdir(d)]
# for e in path:
#     print(e)
#     os.chdir(f"{directory}/{e}")
#     os.system(f"varscan mpileup2snp exeter.{e}.ont.pileup --min-var-freq 0.01  --p-value 1 --min-coverage 1 --min-reads2 1 > {e}.ont.mincov1.mpileup2snp.varscan.tsv")
#     os.system(f"varscan mpileup2snp exeter.{e}.ont.pileup --min-var-freq 0.01  --p-value 1 --min-coverage 1 --min-reads2 1 --strand-filter 0 > {e}.ont.mincov1.mpileup2snp.nostrandbiasfilter.varscan.tsv")
#     # os.system(f"varscan mpileup2snp {e}.ont.pileup --min-var-freq 0.01 --min-coverage 8 >{e}.ont.mincov8.minvar0.01.mpileup2snp.varscan.tsv")
    # os.system(f"varscan mpileup2snp {e}.ont.pileup --min-var-freq 0.02 --min-coverage 8 > {e}.ont.mincov8.minvar0.02.mpileup2snp.varscan.tsv")
    # os.system(f"varscan mpileup2snp {e}.ont.pileup --min-var-freq 0.03 --min-coverage 8 {e}.ont.mincov8.minvar0.03.mpileup2snp.varscan.tsv")
    # os.system(f"varscan pileup2snp {e}.ont.pileup --min-var-freq 0.04 --min-coverage 8 > {e}.ont.mincov8.minvar0.04.mpileup2snp.varscan.tsv")
    # os.system(f"varscan mpileup2snp {e}.ont.pileup --min-var-freq 0.05 --min-coverage 8 > {e}.ont.mincov8.minvar0.05.mpileup2snp.varscan.tsv")
    # os.system(f"varscan mpileup2snp {e}.ont.pileup --min-var-freq 0.03 --min-coverage 8 --min-reads2 6 > {e}.ont.mincov8.minvar0.03.minreads6.mpileup2snp.varscan.tsv")


#Exeter ONT NT002 and 3
# directory = "/home/mbxha18/allcombinations/New_illumina/exeter/exeter_ill_vs_exeter_ont_NT002/"
# path = [d for d in os.listdir(f'{directory}') if os.path.isdir(d)]
# for e in path:
#     print(e)
#     os.chdir(f"{directory}/{e}")
#     file = [d for d in os.listdir('.') if d.endswith(".ont.pileup")]
#     for x in file:
#         print(x)

#         os.system(f"varscan mpileup2snp {x} --min-var-freq 0.01   --min-coverage 16 > {e}.ont.mpileup.mincov16.minvar0.01.varscan.tsv")
#         os.system(f"varscan mpileup2snp {x} --min-var-freq 0.02   --min-coverage 16 > {e}.ont.mpileup.mincov16.minvar0.02.varscan.tsv")
#         os.system(f"varscan mpileup2snp {x} --min-var-freq 0.03   --min-coverage 16 > {e}.ont.mpileup.mincov16.minvar0.03.varscan.tsv")
#         os.system(f"varscan mpileup2snp {x} --min-var-freq 0.04   --min-coverage 16 > {e}.ont.mpileup.mincov16.minvar0.04.varscan.tsv")
#         os.system(f"varscan mpileup2snp {x} --min-var-freq 0.05   --min-coverage 16 > {e}.ont.mpileup.mincov16.minvar0.05.varscan.tsv")
#         os.system(f"varscan mpileup2snp {x} --min-var-freq 0.06   --min-coverage 16 > {e}.ont.mpileup.mincov16.minvar0.06.varscan.tsv")
#         os.system(f"varscan mpileup2snp {x} --min-var-freq 0.07   --min-coverage 16 > {e}.ont.mpileup.mincov16.minvar0.07.varscan.tsv")
        
        
# directory_ill = "/home/mbxha18/allcombinations/New_illumina/exeter/exeter_ill_vs_exeter_ont_NT002/"
# path_ill = [d for d in os.listdir(f'{directory_ill}') if os.path.isdir(d)]
# for e in path_ill:
#     print(e)
#     os.chdir(f"{directory_ill}/{e}")
#     file = [d for d in os.listdir('.') if d.endswith(".ill.pileup")]
#     for x in file:
#         print(x)

#         os.system(f"varscan mpileup2snp {x} --min-var-freq 0.01 --p-value 1 --min-coverage 1 --min-reads2 1 --strand-filter 0 > {e}.ill.mpileup.mincov1.minvar0.01.varscan.tsv")
#         os.system(f"varscan mpileup2snp {x} --min-var-freq 0.02 --p-value 1 --min-coverage 1 --min-reads2 1 --strand-filter 0 > {e}.ill.mpileup.mincov1.minvar0.02.varscan.tsv")
#         os.system(f"varscan mpileup2snp {x} --min-var-freq 0.03 --p-value 1 --min-coverage 1 --min-reads2 1 --strand-filter 0 > {e}.ill.mpileup.mincov1.minvar0.03.varscan.tsv")
#         os.system(f"varscan mpileup2snp {x} --min-var-freq 0.04 --p-value 1 --min-coverage 1 --min-reads2 1 --strand-filter 0 > {e}.ill.mpileup.mincov1.minvar0.04.varscan.tsv")
#         os.system(f"varscan mpileup2snp {x} --min-var-freq 0.05 --p-value 1 --min-coverage 1 --min-reads2 1 --strand-filter 0 > {e}.ill.mpileup.mincov1.minvar0.05.varscan.tsv")
#         os.system(f"varscan mpileup2snp {x} --min-var-freq 0.06 --p-value 1 --min-coverage 1 --min-reads2 1 --strand-filter 0 > {e}.ill.mpileup.mincov1.minvar0.06.varscan.tsv")
#         os.system(f"varscan mpileup2snp {x} --min-var-freq 0.07 --p-value 1 --min-coverage 1 --min-reads2 1 --strand-filter 0 > {e}.ill.mpileup.mincov1.minvar0.07.varscan.tsv")
    
    

    

    
# ##Exeter ILL NT001, 2 and 3
# directory = "/home/mbxha18/allcombinations/New_illumina/exeter/exeter_NT003/"
# path = [d for d in os.listdir(f'{directory}') if os.path.isdir(d)]
# for e in path:
#     print(e)
#     os.chdir(f"{directory}/{e}")
#     os.system(f"varscan mpileup2snp {e}.exeter.ill.pileup --min-var-freq 0.01  --p-value 1 --min-coverage 1 --min-reads2 1 > {e}.ill.mincov1.mpileup2snp.varscan.tsv")
#     os.system(f"varscan mpileup2snp {e}.exeter.ill.pileup --min-var-freq 0.01  --p-value 1 --min-coverage 1 --min-reads2 1 --strand-filter 0 > {e}.ill.mincov1.mpileup2snp.nostrandbiasfilter.varscan.tsv")
    
    
    
# ##Nottingham ILL NT001, 2 and 3
# directory = "/home/mbxha18/allcombinations/New_illumina/nottingham/"
# path = [d for d in os.listdir(f'{directory}') if os.path.isdir(d)]
# for e in path:
#     print(e)
#     os.chdir(f"{directory}/{e}")
#     os.system(f"varscan mpileup2snp {e}.nottingham.ill.pileup --min-var-freq 0.01  --p-value 1 --min-coverage 1 --min-reads2 1 > {e}.ill.mincov1.mpileup2snp.varscan.tsv")
#     os.system(f"varscan mpileup2snp {e}.nottingham.ill.pileup --min-var-freq 0.01  --p-value 1 --min-coverage 1 --min-reads2 1 --strand-filter 0 > {e}.ill.mincov1.mpileup2snp.nostrandbiasfilter.varscan.tsv")



# ##Nottingham ONT HAC_GRID
# directory = "/home/mbxha18/allcombinations/kit12/10.4_HAC_GRID_VARSCAN/ill_match/"
# path = [d for d in os.listdir(f'{directory}') if os.path.isdir(d)]
# for e in path:
#     print(e)
#     os.chdir(f"{directory}/{e}")
#     os.system(f"varscan mpileup2snp {e}.ont.pileup --min-var-freq 0.01  --p-value 1 --min-coverage 1 --min-reads2 1 > {e}.ont.mincov1.mpileup2snp.varscan.tsv")
#     os.system(f"varscan mpileup2snp {e}.ont.pileup --min-var-freq 0.01  --p-value 1 --min-coverage 1 --min-reads2 1 --strand-filter 0 > {e}.ont.mincov1.mpileup2snp.nostrandbiasfilter.varscan.tsv")

 ##Nottingham ONT SUP_GRID
# directory = "/home/mbxha18/allcombinations/kit12/10.4_SUP_GRID/ill_match/"
# path = [d for d in os.listdir(f'{directory}') if os.path.isdir(d)]
# for e in path:
#     print(e)
#     os.chdir(f"{directory}/{e}")
#     os.system(f"mosdepth --no-per-base {e}.ont {e}.primertrimmed.rg.sorted.bam")
#     os.system(f"samtools depth -aa {e}.primertrimmed.rg.sorted.bam > {e}.ont.samtools.cov ")
#     os.system(f"samtools mpileup -f /home/mbxha18/artic/fieldbioinformatics/test-data/primer-schemes/nCoV-2019/V9/nCoV-2019.reference.fasta -q 10 -d 1000000 {e}.primertrimmed.rg.sorted.bam > {e}.ont.pileup")
#     os.system(f"varscan mpileup2snp {e}.ont.pileup --min-var-freq 0.01  --p-value 1 --min-coverage 1 --min-reads2 1 > {e}.ont.mincov1.mpileup2snp.varscan.tsv")
#     os.system(f"varscan mpileup2snp {e}.ont.pileup --min-var-freq 0.01  --p-value 1 --min-coverage 1 --min-reads2 1 --strand-filter 0 > {e}.ont.mincov1.mpileup2snp.nostrandbiasfilter.varscan.tsv")


# ##Liverpool ILL NT001
# directory = "/home/mbxha18/allcombinations/New_illumina/liverpool/"
# path = [d for d in os.listdir(f'{directory}') if os.path.isdir(d)]
# for e in path:
#     print(e)
#     os.chdir(f"{directory}/{e}")
#     os.system(f"varscan mpileup2snp {e}.liverpool.ill.pileup --min-var-freq 0.01  --p-value 1 --min-coverage 1 --min-reads2 1 > {e}.ill.mincov1.mpileup2snp.varscan.tsv")
#     os.system(f"varscan mpileup2snp {e}.liverpool.ill.pileup --min-var-freq 0.01  --p-value 1 --min-coverage 1 --min-reads2 1 --strand-filter 0 > {e}.ill.mincov1.mpileup2snp.nostrandbiasfilter.varscan.tsv")

# ##WW080 real sample ONT 
# directory = "/home/mbxha18/WW080/ill_match/"
# path = [d for d in os.listdir(f'{directory}') if os.path.isdir(d)]
# for e in path:
#     print(e)
#     os.chdir(f"{directory}/{e}")
#     os.system(f"varscan mpileup2snp artic.WW080.{e}.WW080.ont.pileup --min-var-freq 0.01  --p-value 1 --min-coverage 1 --min-reads2 1 > artic.WW080.{e}.ont.mincov1.mpileup2snp.varscan.tsv")
#     os.system(f"varscan mpileup2snp artic.WW080.{e}.WW080.ont.pileup --min-var-freq 0.01  --p-value 1 --min-coverage 1 --min-reads2 1 --strand-filter 0 > artic.WW080.{e}.ont.mincov1.mpileup2snp.nostrandbiasfilter.varscan.tsv")


# ##10.4_SUP_PROM_VARSCAN/KIT14 ONT (old Exeter)
# directory = "/home/mbxha18/allcombinations/kit12/9.4_HAC_GRID_VARSCAN/ill_match/"
# path = [d for d in os.listdir(f'{directory}') if os.path.isdir(d)]
# for e in path:
#     print(e)
#     os.chdir(f"{directory}/{e}")
#     os.system(f"varscan mpileup2snp {e}.ont.pileup --min-var-freq 0.01  --p-value 1 --min-coverage 1 --min-reads2 1 > {e}.ont.mincov1.mpileup2snp.varscan.tsv")
#     os.system(f"varscan mpileup2snp {e}.ont.pileup --min-var-freq 0.01  --p-value 1 --min-coverage 1 --min-reads2 1 --strand-filter 0 > {e}.ont.mincov1.mpileup2snp.nostrandbiasfilter.varscan.tsv")
    # os.system(f"samtools depth -aa artic.9.4_HAC_GRID.{e}.primertrimmed.rg.sorted.bam > {e}.ont.samtools.cov")


#run freyja
import os
directory = "/home/mbxha18/allcombinations/New_illumina/nottingham/nottingham_ill_NT001/"
ill_dir = [d for d in os.listdir(f'{directory}') if os.path.isdir(d)]
for e in ill_dir:
    print(e)
    os.chdir(f"{directory}/{e}")
    # os.system(f"freyja variants {e}.nottingham.ill.manualivar.trimmed.bam  --variants ill.{e}.variants.tsv --depths ill.{e}.depths.tsv --ref /home/mbxha18/waste_water/nextflow2/ww_nf_minimal/static/NC_045512.2.fa")
    # os.system(f"freyja demix ill.{e}.variants.tsv ill.{e}.depths.tsv --output ill{e}.demixedtest.tsv")
    os.system(f"mv {directory}/{e}/ill{e}.demixedtest.tsv {directory}/demixed")
    os.system(f"freyja aggregate {directory}/demixed/ --output ill{e}.demixed.tsv")
    
    
    