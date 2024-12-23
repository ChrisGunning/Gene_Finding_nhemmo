"""
Script for converting .gtf files into the exon file format required by GlimmerHMM.
"""
from dotenv import load_dotenv
import os

# Load environment variables
load_dotenv()
GTF_FILEPATH = os.getenv('GTF_FILEPATH')
OUTPUT_EXON_FILEPATH = os.getenv('OUTPUT_EXON_FILEPATH')
SUMMARY_FILEPATH = os.getenv('SUMMARY_FILEPATH')

# Keys are chromosome sequence numbers
chromosome_sequences = {}

# Open the GTF file
with open(GTF_FILEPATH, "rt") as gtf:
    
    for line in gtf:
        if line.startswith("#") or line.strip() == "":
            continue

        fields = line.strip().split("\t")
        seq_id, feature_type = fields[0], fields[2]
        start, end, attributes = int(fields[3]), int(fields[4]), fields[8]

        # Stop if not longer considering chromosome sequences
        if not seq_id.isdigit():
            break
        
        # Only process exons
        if feature_type != "exon":
            continue
        
        # Extract gene ID from attributes
        gene_id = None
        for attr in attributes.split(";"):
            if "gene_id" in attr:
                gene_id = attr.split("\"")[1]
                break
        
        if seq_id not in chromosome_sequences:
            chromosome_sequences[seq_id] = []
        chromosome_sequences[seq_id].append((gene_id, start, end))


# Summary statistics
count_dict = {}
total_num_sequences = 0
total_num_genes = 0
total_num_exons = 0



# Write data to output file
with open(OUTPUT_EXON_FILEPATH, "w") as out:
    for seq, data in chromosome_sequences.items():
        total_num_sequences += 1
        count_dict[seq] = []

        # Group exons by gene
        genes = {}
        for gene_id, start, end in data:
            if gene_id not in genes:
                genes[gene_id] = []
            genes[gene_id].append((start, end))
        
        # Write exon data for each gene
        for gene_id, exons in genes.items():
            total_num_genes += 1
            gene_exon_count = 0
            for start, end in exons:
                total_num_exons += 1
                gene_exon_count += 1
                out.write(f"seq{seq} {start} {end}\n")
            out.write("\n")  # Blank line between genes
            count_dict[seq].append((gene_id, gene_exon_count))


with open(OUTPUT_EXON_FILEPATH, "r") as infile:
    data = infile.read()

# Split groups by blank lines
groups = data.strip().split("\n\n")

# Sort each group by start coordinate
sorted_groups = []
for group in groups:
    lines = group.strip().split("\n")
    sorted_group = sorted(
        (line.split() for line in lines),  # Split each line into fields
        key=lambda x: int(x[1])  # Sort by the second field (start coordinate)
    )
    sorted_groups.append("\n".join(" ".join(line) for line in sorted_group))

# Write sorted groups back to the output file
with open(OUTPUT_EXON_FILEPATH, "w") as outfile:
    outfile.write("\n\n".join(sorted_groups))


# Write statistics to summary file
with open(SUMMARY_FILEPATH, "w") as out:
    out.write(f"Total number of sequences: {total_num_sequences}\n")
    out.write(f"Total number of genes: {total_num_genes}\n")
    out.write(f"Total number of exons: {total_num_exons}\n\n")

    for seq, data in count_dict.items():
        out.write(f"Sequence {seq}:\n")
        for gene_id, gene_exon_count in data:
            out.write(f"\tGene {gene_id}: {gene_exon_count} exons\n")
        out.write("\n")
        

                
      
