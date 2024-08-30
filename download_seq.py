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
    try:
        subprocess.run(["blastn", "-query", input_file, "-out", output_file])
        print(f"BLAST results saved to {output_file}")
    except Exception as e:
        print(f"Error running BLAST: {e}")


# Main script
if __name__ == "__main__":
    gene_id = "NM_007294"  # Example gene ID for BRCA1
    fasta_file = "BRCA1_sequence.fasta"
    blast_output = "BRCA1_blast_results.txt"
    
    # Fetch the sequence
    fetch_sequence(gene_id, fasta_file)
    
    # Run BLAST on the fetched sequence
    run_blast(fasta_file, blast_output)
