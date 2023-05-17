import os
import subprocess
ill_dir = [d for d in os.listdir('/home/mbxha18/bigzipfile/exeter_results_NT002') if os.path.isdir(d)]
downsample_percentages = [0.01,0.05, 0.1,0.25,0.5,0.75,1]
for dir_path in ill_dir:
	for file_name in os.listdir(dir_path):
		if file_name.startswith("guppyplexed") and not file_name.endswith("downsampled.fastq"):
			print(os.path.join(dir_path, file_name))
			for percentage in downsample_percentages:
	 			os.system(f"seqkit sample -p {percentage} /home/mbxha18/bigzipfile/exeter_results_NT002/{dir_path}/{file_name} > /home/mbxha18/bigzipfile/exeter_results_NT002/{dir_path}/{file_name}.{percentage}.downsampled.fastq")
