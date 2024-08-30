from Bio import Entrez
import subprocess

# Function to retrieve DNA sequence from NCBI using Entrez
def fetch_sequence(gene_id, output_file):
    Entrez.email = "your.email@example.com"
    handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype="fasta", retmode="text")
    sequence_data = handle.read()
    handle.close()
    
    with open(output_file, "w") as file:
        file.write(sequence_data)
    print(f"Sequence for {gene_id} saved to {output_file}")

# Function to run BLAST using subprocess
def run_blast(input_file, output_file):
    attempts = 0
    while attempts < 5:
        try:
            subprocess.run(["blastn", "-query", input_file, "-db", "nt", "-out", output_file, "-outfmt", "6", "-remote"])
            print(f"BLAST results saved to {output_file}")
            break
        except Exception as e:
            attempts += 1
            if attempts == 5:
                print("All attempts to run BLAST failed.")
            print(f"Error running BLAST: {e}")

def add_headers_with_format(output_file):
    headers = [
        "Query ID", "Subject ID", "% Identity", "Alignment Length", "Mismatch",
        "Gap Opens", "Q. Start", "Q. End", "S. Start", "S. End", "E-value", "Bit Score"
    ]
    # Define a format string with specified column widths
    header_format = "{:<12} {:<12} {:<12} {:<18} {:<10} {:<10} {:<8} {:<8} {:<8} {:<8} {:<10} {:<10}\n"

    with open(output_file, 'r') as original:
        data = original.readlines()
    
    # Write the headers with formatting
    with open(output_file, 'w') as modified:
        modified.write(header_format.format(*headers))
        for line in data:
            formatted_line = header_format.format(*line.strip().split('\t'))
            modified.write(formatted_line)

# Example usage
output_file = "BRCA1_blast_results.txt"
add_headers_with_format(output_file)


# Function to calcualte GC content of Fasta file
def GC_content(output_file):
    gc_total = 0
    length = len(output_file)
    for nuc in output_file:
        if nuc == 'G' or nuc == 'C':
            gc_total += 1
    gc_percentage = (gc_total / length) * 100
    print(f"GC Percentage = {gc_percentage}%")

def AT_content(output_file):
    at_total = 0
    length = len(output_file)
    for nuc in output_file:
        if nuc == 'A' or nuc == 'T':
            at_total += 1
    at_percentage = (at_total / length) * 100
    print(f"AT Percentage = {at_percentage}%")
            
        
# Main script
if __name__ == "__main__":
    gene_id = "NM_007294"  # Example gene ID for BRCA1
    fasta_file = "BRCA1_sequence.fasta"
    blast_output = "BRCA1_blast_results.txt"

    # Fetch the sequence
    fetch_sequence(gene_id, fasta_file)

    # Run BLAST on the fetched sequence
    run_blast(fasta_file, blast_output)

    # Calculate GC and AT %
    GC_content(fasta_file)
    AT_content(fasta_file)

    # Add headers to output file
    add_headers_with_format(output_file)

#To create this code, I first imported Entrez using biopython. This gives me access to the NCBI database.
#Eventually we want to run blast on our dataset so we need to import subprocess for that. In fetching our DNA
#sequence I created a function for that takes 2 parameters: gene_id and output_file. The gene_id is the specific 
#gene you want to obtain and the output file is where your sequence will be stored. Then we specify with queries the 
#information needed. We are searching the nucleotide database from the gene id we provide, then the returning it as a fasta file 
#and as text rather than binary. Once we get all the information we close the NCBI server. Then we just print a message
#saying we the fetched sequence was saved in the output file. 
#    Our Blast function takes 2 parameters: An input file and an output file. Subprocess.run is used to run the blastn command. 
#Just like the previous function all the paramaters are being listed inside subpress, and then the results are saved to an output file.
#If there is an error running blast it will display that there is an error.
#    Our Main script is testing this on the given gene id: NM_007294. The fasta file is the name where the sequence will be saved.
#The the blast output is where the blast results will be saved. We finall call our previous 2 functions and obtain our sequence, then run
#blast on it and we will obtain 2 new files with the outputs.

