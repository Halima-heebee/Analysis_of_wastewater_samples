import os
import subprocess

# Define the directory containing the files
ill_dir = [d for d in os.listdir('/home/mbxha18/bigzipfile/exeter_results_NT002') if os.path.isdir(d)]

# Define the percentages for downsampling
downsample_percentages = [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 1]

# Iterate over each directory in the ill_dir list
for dir_path in ill_dir:
    # Iterate over each file in the directory
    for file_name in os.listdir(dir_path):
        # Check if the file starts with "guppyplexed" and does not end with "downsampled.fastq"
        if file_name.startswith("guppyplexed") and not file_name.endswith("downsampled.fastq"):
            # Print the file path
            print(os.path.join(dir_path, file_name))
            
            # Iterate over each percentage for downsampling
            for percentage in downsample_percentages:
                # Create a downsampling command using seqkit
                # and execute it using os.system()
                os.system(f"seqkit sample -p {percentage} /home/mbxha18/bigzipfile/exeter_results_NT002/{dir_path}/{file_name} > /home/mbxha18/bigzipfile/exeter_results_NT002/{dir_path}/{file_name}.{percentage}.downsampled.fastq")
