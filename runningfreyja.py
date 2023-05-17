import os
ill_dir = [d for d in os.listdir('/home/mbxha18/allcombinations/New_illumina/exeter/exeter_ill_vs_exeter_ont_NT003') if os.path.isdir(d)]
for e in ill_dir:
    print(e)
    os.chdir(f"/home/mbxha18/allcombinations/New_illumina/exeter/exeter_ill_vs_exeter_ont_NT003/{e}")
    os.system(f"freyja variants {e}.exeter.ill*.bam --variants ill.{e}.variants.tsv --depths ill.{e}.depths.tsv --ref /home/mbxha18/waste_water/nextflow2/ww_nf_minimal/static/NC_045512.2.fa")
    os.system(f"freyja demix ill.{e}.variants.tsv ill.{e}.depths.tsv --output ill{e}.demixedtest.tsv")
    #os.mkdir(f'/home/mbxha18/allcombinations/New_illumina/exeter/exeter_ill_vs_exeter_ont_NT003/demixed')
    os.system(f"mv /home/mbxha18/allcombinations/New_illumina/exeter/exeter_ill_vs_exeter_ont_NT003/{e}/ill{e}.demixedtest.tsv /home/mbxha18/allcombinations/New_illumina/exeter/exeter_ill_vs_exeter_ont_NT003/demixed")
    os.system(f"freyja aggregate /home/mbxha18/allcombinations/New_illumina/exeter/exeter_ill_vs_exeter_ont_NT003/demixed/ --output ill-exeterNT003-ww-demixed.tsv")
    #os.system(f"freyja plot ill{e}-ww-demixed.tsv --output ill{e}-ww-demixed.png")
