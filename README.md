# BRCA1 Gene Sequence Analysis

## Project Overview

This project is designed to analyze a gene sequence, specifically the BRCA1 gene. The analysis includes calculating GC content, AT content, motif searching, and translating nucleotide sequences into protein sequences. The project aims to offer an automated pipeline for sequence analysis and provides users with tools to explore various aspects of genetic sequences.

## Features

1. **Fetch Gene Sequence**: 
   - The script fetches the BRCA1 gene sequence from the NCBI database using the gene ID `NM_007294`.
   - The sequence is saved in a FASTA file format (`BRCA1_sequence.fasta`).

2. **Run BLAST**: 
   - The script executes a BLAST search against the NCBI nucleotide database (`nt`) to find regions of similarity between the BRCA1 sequence and other nucleotide sequences.
   - The BLAST results are saved in `BRCA1_blast_results.txt`.

3. **GC and AT Content Calculation**: 
   - The script calculates the GC and AT content of the BRCA1 sequence.
   - The results are displayed as percentages in the terminal.

4. **Motif Search**: 
   - Users can search for specific motifs (short nucleotide sequences) within the BRCA1 sequence.
   - If the user opts out by typing `na`, the motif search is skipped.

5. **Sequence Translation**: 
   - The script translates the BRCA1 nucleotide sequence into a protein sequence.
   - The translated protein sequence is saved in `translated_protein.txt`.

## Requirements

- Python 3.x
- Biopython library
- Internet connection for fetching data from NCBI and running BLAST remotely.

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/your-username/BRCA1-Gene-Sequence-Analysis.git
   cd BRCA1-Gene-Sequence-Analysis
   ```

2. Install the required Python packages:
   ```bash
   pip install biopython
   ```

## Usage

To run the sequence analysis, simply execute the `download_seq.py` script:

```bash
python download_seq.py
```

### Script Workflow:

1. **Fetching the Gene Sequence**:
   - The sequence for the gene ID `NM_007294` is fetched from NCBI and saved as `BRCA1_sequence.fasta`.

2. **Running BLAST**:
   - The fetched sequence is used to perform a BLAST search.
   - The BLAST results are saved in `BRCA1_blast_results.txt`.

3. **Calculating GC and AT Content**:
   - The GC and AT content percentages of the sequence are calculated and displayed.

4. **Motif Search**:
   - The script will prompt you to enter a motif for searching within the BRCA1 sequence.
   - If you do not wish to perform a motif search, type `na`.

5. **Sequence Translation**:
   - The nucleotide sequence is translated into a protein sequence and saved as `translated_protein.txt`.

### Example Output:

```bash
Sequence for NM_007294 saved to BRCA1_sequence.fasta
BLAST results saved to BRCA1_blast_results.txt
GC Percentage = 41.77%
AT Percentage = 58.23%
Please enter a motif to search for (or type 'na' to skip): GCA
Motif 'GCA' found 115 times.
Translated protein saved to translated_protein.txt
```

### Output Files:

- `BRCA1_sequence.fasta`: The FASTA file containing the BRCA1 gene sequence.
- `BRCA1_blast_results.txt`: The file containing the BLAST search results.
- `translated_protein.txt`: The file containing the translated protein sequence.

## Future Work

- **Sequence Classification**: Extend the project by implementing a sequence classification model to predict the function of the translated protein.
- **Protein Structure Prediction**: Incorporate tools like AlphaFold to predict the 3D structure of the translated protein.
- **Interactive User Interface**: Develop a GUI for easier use and visualization of results.

## Contributing

Contributions are welcome! Please submit a pull request or open an issue if you have suggestions for improvements.

## License

This project is licensed under the MIT License.

