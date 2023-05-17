import os
from collections import Counter
# Define the parent directory containing subdirectories with varscan files
parent_dir = '/home/mbxha18/allcombinations/New_illumina/exeter/exeter_ill_vs_exeter_ont_NT002'

all_subdirs = [d for d in os.listdir(parent_dir) if os.path.isdir(os.path.join(parent_dir, d))]

# Create a list to store all unique positions from each subdirectory
all_unique_positions = []

for subdir in all_subdirs:
    subdir_path = os.path.join(parent_dir, subdir)

    # Create a list of varscan files in the current subdirectory
    varscan_files = [file for file in os.listdir(subdir_path) if file.endswith(('.ont.mincov1.mpileup2snp.nostrandbiasfilter.varscan.tsv', '.ont.mincov1.mpileup2snp.varscan.tsv'))]

    # Check that both types of varscan files are present
    if len(varscan_files) == 2:

        # Load the contents of each varscan file into a set
        with open(os.path.join(subdir_path, varscan_files[0]), 'r') as file1, open(os.path.join(subdir_path, varscan_files[1]), 'r') as file2:
            varscan1 = set(line.strip() for line in file1.readlines()[1:])
            varscan2 = set(line.strip() for line in file2.readlines()[1:])

        # Find the positions unique to the first varscan file
        unique_positions = varscan1 - varscan2

        # Save the unique positions to a file in the current subdirectory
        output_filename = os.path.join(subdir_path, 'StrandBiasPositions.tsv')
        with open(output_filename, 'w') as outfile:
            outfile.write('Chrom\tPosition\tRef\tVar\tCons:Cov:Reads1:Reads2:Freq:P-value\tStrandFilter:R1+:R1-:R2+:R2-:pval\tSamplesRef\tSamplesHet\tSamplesHom\tSamplesNC\tCons:Cov:Reads1:Reads2:Freq:P-value\n')
            outfile.writelines(f"{line}\n" for line in sorted(unique_positions))

        print(f"Unique positions in {subdir} have been saved to {output_filename}")
        
        # Append the unique positions to the list of all unique positions
        all_unique_positions.extend(unique_positions)
    else:
        print(f"Both types of varscan files not found in {subdir}")

# Save all unique positions from all subdirectories to a file called "biasedSNPs.tsv"
output_filename = os.path.join(parent_dir, 'biasedSNPs.tsv')
with open(output_filename, 'w') as outfile:
    outfile.write('Position\n')
    # Iterate through all files ending in "StrandBiasPositions.tsv" in each subdirectory
    for subdir in all_subdirs:
        subdir_path = os.path.join(parent_dir, subdir)
        bias_file = os.path.join(subdir_path, 'StrandBiasPositions.tsv')
        # Check that the current subdirectory contains a "StrandBiasPositions.tsv" file
        if os.path.isfile(bias_file):
            # Read all lines except the header and keep only the "Position" column
            with open(bias_file, 'r') as infile:
                positions = [line.strip().split('\t')[1] for line in infile.readlines()[1:]]
            # Write the positions to the output file
            outfile.write('\n'.join(positions) + '\n')

# Remove duplicate positions from the output file

with open('biasedSNPs.tsv', 'r') as infile:
    positions = [line.strip().split('\t')[0] for line in infile.readlines()]

# Count the number of times each value of "Position" is observed
position_counts = Counter(positions)

# Write the counts to a new file called "biasedSNPs_counts.tsv"
with open('biasedSNPs_counts.tsv', 'w') as outfile:
    outfile.write('Position\tCount\n')
    outfile.writelines(f"{position}\t{count}\n" for position, count in position_counts.items())
