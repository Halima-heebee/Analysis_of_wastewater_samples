#Loose Lab, University of Nottingham. 17.05.2023.
#Takes demultiplexed ONT fastq files barcoded using Nimagen Indexing Kits (completed using X.py), matches them to their corresponding well plate using a sample sheet, and inputs them to the first two stages of the artic pipeline (guppuplex and minion).

import os
import pandas as pd

# RUN guppyplex on each Nimagen folder in ont
ont_dir = [d for d in os.listdir('/home/mbxha18/allcombinations/kit14/10.4_SUP_PROM_VARSCAN/NT001_super_barcoded') if os.path.isdir(d)] #provide path to demultiplexed and barcoded fastq files

for e in ont_dir:
    print(e)
    os.chdir(f"/home/mbxha18/allcombinations/kit14/10.4_SUP_PROM_VARSCAN/NT001_super_barcoded/{e}/") #replace this example path with your ont_dir path
    os.system(f"artic guppyplex --min-length 80 --max-length 480 --directory . --prefix guppyplexed.{e}") #run guppyplex filtering. Edit length of amplicons, for our analysis we selected 80 and 480, as average amplicon size ranged between 230-280bp and artic suggests 200bp intervals either side. 



# #Read in the bigsamplesheet containing Nimagen indexes merged
csv = pd.read_csv("/home/mbxha18/OvsI_prom10.4sup_onebarcode/UDI06samplesheet.txt", sep="\t") #read in relevant Nimagen Index Kit Sample Sheet
csv = csv.rename(columns={"Index combination name": "Index_combination_name"})
csv = csv.rename(columns={"Plate location": "Plate_location"})
csv["index_number"] = csv.Index_combination_name.str[-4:7]
csv

# # match the ont nimagen index folder to the csv file to get the corresponding plate well
plates = []
all_subdirs = [d for d in os.listdir('/home/mbxha18/allcombinations/kit14/10.4_SUP_PROM_VARSCAN/NT001_super_barcoded/') if d.startswith("NIMAGEN")]    #lists all directories beginning with NIMAGEN in our folder. May need changing dependent on data structure.
for name in all_subdirs:
    Nimagen_index = name[7:12]
    for index, row in csv.iterrows():
        if (row["index_number"]  == Nimagen_index):
            plate = row["Plate_location"]
            plates.append(f"{plate}-{Nimagen_index}")
# print (plates)           
            
# # create a list of the NT001 illumina fastq and then a set of their well plate. Change dependent on your data structure.
file_lists= []
illumina_files = [f for f in os.listdir('/home/liam/wastewater/processed/NT001/')]
for x in illumina_files:
    if x[11:42][0:5] == "NT001":
        file_lists.append(x)
    if x[12:42][0:5] == "NT001":
        file_lists.append(x)

directory = []    
file_list_id = ([x.split('_')[0] for x in file_lists])
file_set = set(file_list_id)


# # make a directory that matches both the well plate and the nimagen index
for string in file_set:
    for e in plates:
        if string.split('-')[0] == e.split('-')[0]:
            os.mkdir(f'/home/mbxha18/allcombinations/kit14/10.4_SUP_PROM_VARSCAN/ill_match/{string}-{e[-4:(len(e))]}')


# Move the guppyplex files in each ill_match (ill) folder

ont = [x for x in os.listdir('/home/mbxha18/allcombinations/kit14/10.4_SUP_PROM_VARSCAN/NT001_super_barcoded/')] #path to Nimagen-barcoded ONT data
ill = [x for x in os.listdir('/home/mbxha18/allcombinations/kit14/10.4_SUP_PROM_VARSCAN/ill_match/')] #paht to Illumina data
for b in ont:
    for c in ill:
         try:
             if b[-4:] == c.split('-')[4]:
                 os.system(f"cp /home/mbxha18/allcombinations/kit14/10.4_SUP_PROM_VARSCAN/NT001_super_barcoded/{b}/guppyplexed.{b}_..fastq /home/mbxha18/allcombinations/kit14/10.4_SUP_PROM_VARSCAN/ill_match/{c}")
         except IndexError:
             pass



# # Run artic minion
ill_dir = [d for d in os.listdir('/home/mbxha18/allcombinations/kit14/10.4_SUP_PROM_VARSCAN/ill_match/') if os.path.isdir(d)]
for e in ill_dir:
    print(e)
    os.chdir(f"/home/mbxha18/allcombinations/kit14/10.4_SUP_PROM_VARSCAN/ill_match/{e}")
    os.system(f"artic minion --normalise 1000000 --threads 4 --medaka --medaka-model r1041_e82_400bps_sup_g615 --scheme-directory /home/mbxha18/artic/fieldbioinformatics/test-data/primer-schemes --read-file guppyplexed.NIMAGEN{e[-4:]}_..fastq nCoV-2019/V9 artic.10.4_HAC_PROM.{e}") #Run Artic Minion. Use a very high normalize value to restrict depth normalization so observed variant frequencies are correct. Provide relevant base-calling model and the primer scheme here.




# # Run varscan for ont

import os
ill_dir = [d for d in os.listdir('/home/mbxha18/allcombinations/kit14/10.4_SUP_PROM_VARSCAN/ill_match/') if os.path.isdir(d)]

for e in ill_dir:
    print(e)
    os.chdir(f'/home/mbxha18/allcombinations/kit14/10.4_SUP_PROM_VARSCAN/ill_match/{e}')
    os.system(f"samtools mpileup -f /home/mbxha18/artic/fieldbioinformatics/test-data/primer-schemes/nCoV-2019/V9/nCoV-2019.reference.fasta -q 10 -d 1000000 Kit14.artic.10.4_SUP_PROM.{e}.primertrimmed.rg.sorted.bam > {e}.ont.pileup") #Run samtools mpileup, again with extremely high normalization value.
    os.system(f"varscan mpileup2snp {e}.ont.pileup --min-var-freq 0.01  --p-value 1 --min-coverage 1 --min-reads2 1 --strand-filter 0 > {e}.ont.mincov1.varscan.tsv")
    os.system(f"varscan mpileup2snp {e}.ont.pileup --min-var-freq 0.01 --min-coverage 8 > {e}.ont.mincov8.minvar0.01.varscan.tsv")
    os.system(f"varscan mpileup2snp {e}.ont.pileup --min-var-freq 0.02 --min-coverage 8 > {e}.ont.mincov8.minvar0.02.varscan.tsv")
    os.system(f"varscan mpileup2snp {e}.ont.pileup --min-var-freq 0.03 --min-coverage 8 > {e}.ont.mincov8.minvar0.03.varscan.tsv")
    os.system(f"mosdepth --no-per-base {e}.ont Kit14.artic.10.4_SUP_PROM.{e}.primertrimmed.rg.sorted.bam") #Generate Mosdepth file
    os.system(f"samtools depth -aa Kit14.artic.10.4_SUP_PROM.{e}.primertrimmed.rg.sorted.bam > {e}.ont.samtools.cov ") #Generate Coverage File



