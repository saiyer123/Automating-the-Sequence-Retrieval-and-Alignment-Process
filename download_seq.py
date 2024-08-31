from Bio import Entrez
import subprocess
import sys

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


# Function to calcualte GC content of Fasta file and total nucleotide length
def GC_content(file_path):
    gc_total = 0
    atcg_length = 0
    with open(file_path, 'r') as file:
        sequence = file.read().replace('\n', '')  # Read the file and remove any newlines
    for nuc in sequence:
        if nuc in 'ATCG':  # Count all nucleotides
            atcg_length += 1
        if nuc in 'GC':  # Count GC nucleotides
            gc_total += 1
    gc_percentage = (gc_total / atcg_length) * 100
    print(f"Sequence Length: {atcg_length} Nucleotides")
    print(f"GC Percentage = {gc_percentage:.2f}%")
    return gc_percentage

def AT_content(file_path):
    at_total = 0
    atcg_length = 0
    with open(file_path, 'r') as file:
        sequence = file.read().replace('\n', '')  # Read the file and remove any newlines
    for nuc in sequence:
        if nuc in 'ATCG':  # Count all nucleotides
            atcg_length += 1
        if nuc in 'AT':  # Count AT nucleotides
            at_total += 1
    at_percentage = (at_total / atcg_length) * 100
    print(f"AT Percentage = {at_percentage:.2f}%")
    return at_percentage

def Motif_search(motif, output_file):
    with open(output_file, 'r') as file:
        sequence = file.read().replace('\n', '')
    motif_count = 0
    motif_length = len(motif)

    for i in range(len(sequence) - motif_length + 1):
        if sequence[i:i + motif_length] == motif:
            motif_count += 1
    print(f"Motif '{motif}' found {motif_count} times.")               

def translate_seq(output_file):
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }

    amino_acid_dict = {
        'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS',
        'I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN', 
        'G': 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP', 
        '*': 'TER', 'A': 'ALA', 'V': 'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET', 'X': 'XAA'
    }

    try:
        with open(input_file, 'r') as file:
            sequence = file.read().replace('\n', '')
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found.")
        return None
    
    print(f"Read sequence length: {len(sequence)}")
    
    trimmed_seq = sequence[:len(sequence) - (len(sequence) % 3)]
    print(f"Trimmed sequence length: {len(trimmed_seq)}")
    
    protein = ""
    full_protein = []
    
    for i in range(0, len(trimmed_seq), 3):
        codon = trimmed_seq[i:i+3]
        amino_acid = codon_table.get(codon, 'X')
        protein += amino_acid
        full_protein.append(amino_acid_dict.get(amino_acid, 'XAA'))
    
    print(f"Translated protein length: {len(protein)}")
    
    try:
        with open('translated_protein.txt', "w") as output_file:
            output_file.write(protein)
        print("Protein sequence written to translated_protein.txt")
    except Exception as e:
        print(f"Error writing to file: {e}")
    
    return protein


            
            
        
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
    add_headers_with_format(blast_output)

    motif = input("Please enter a motif to search for (or type 'na' to skip): ").strip()

    if motif.lower() != 'na':
        Motif_search(motif, fasta_file)
    else:
        print("Motif search skipped.")

    translate_seq('BRCA1_sequence.fasta')




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

