#Run Alex's nextflow (idk how), and move the sorted bams into the correct directory for varscanning




# # # Run varscan for ill
# ill_dir = [d for d in os.listdir('/home/mbxha18/bigrun/mattrebasecalled+normalise+guppyplex/ill/') if os.path.isdir(d)]

# for e in ill_dir:
#     print(e)
#     os.chdir(f"/home/mbxha18/bigrun/mattrebasecalled+normalise+guppyplex/ill/{e}")
#     os.system(f"samtools mpileup -f /home/mbxha18/waste_water/nextflow2/ww_nf_minimal/static/NC_045512.2.fa -q 10 -d 1000000 {e}.ill.manualivar.trimmed.bam > {e}.ill.pileup")
#     os.system(f"varscan pileup2snp {e}.ill.pileup --min-var-freq 0.01  --p-value 1 --min-coverage 1 --min-reads2 1 > {e}.ill.mincov1.varscan.tsv")
#     os.system(f"varscan pileup2snp {e}.ill.pileup --min-var-freq 0.01 --min-coverage 8 > {e}.ill.mincov8.minvar0.01.varscan.tsv")
#     os.system(f"varscan pileup2snp {e}.ill.pileup --min-var-freq 0.02 --min-coverage 8 > {e}.ill.mincov8.minvar0.02.varscan.tsv")
#     os.system(f"varscan pileup2snp {e}.ill.pileup --min-var-freq 0.03 --min-coverage 8 > {e}.ill.mincov8.minvar0.03.varscan.tsv")
#     os.system(f"mosdepth --no-per-base {e}.ill {e}.ill.manualivar.trimmed.bam")
#     os.system(f"samtools depth -aa {e}.ill.manualivar.trimmed.bam > {e}.ill.samtools.cov ")
