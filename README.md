# primer_cheq

Tool for checking primer coverage and drift.

## Description

`primer_cheq.py` is a Python script designed to check the coverage and drift of primer sequences against viral and bacterial genomes. It downloads genome data, processes it, and uses BLAST to align primer sequences to the genomes, providing detailed information on mismatches and potential issues.

## Requirements

- Python 3.x
- BLAST+ (installed and available in the system PATH)
- NCBI Datasets command-line tool (installed and available in the system PATH)

## Installation

1. Clone the repository:
    ```sh
    git clone https://github.com/yourusername/primer_cheq.git
    cd primer_cheq
    ```

2. Install the required Python packages:
    ```sh
    pip install -r requirements.txt
    ```

## Usage

```sh
python primer_cheq.py -p primer.fasta -o sample_name -w working_dir [options]
```

### Arguments
#### Required
- `-p`, `--primers`: FASTA file of primers to check.
- `-o`, `--prefix`: Prefix for output files in the working directory.
- `-w`, `--working_directory`: Path for intermediate and output files.
#### One or more of the following
- `-v`, `--ncbi_virus`: Download a database of viral references from NCBI using a taxid.
- `-b`, `--ncbi_bacteria`: Download a database of microbial references from NCBI using a taxid.
- `-d`, `--directory_db`: A directory of FASTA files, one for each reference.
- `-f`, `--fasta_db`: A single FASTA file, each entry will be treated as its own reference.
- `-g`, `--glob_db`: A glob (i.e. "/path/to/data/*.fasta") used to identify FASTA files each FASTA will be treated as its own reference.

n.b. arguments can be used multiple times i.e. -f genomes.fasta -f another_genome.fasta -v 10244 will create a database using both fastas and download additional FASTAS from NCBI.

### Example

```sh
python primer_cheq.py -p primers.fasta -o my_sample -w ./workdir -v 10244
```

This command will check the primers in `primers.fasta` against viral genomes with taxid `10244` downloaded from NCBI, using `my_sample` as the prefix for output files and `./workdir` as the working directory.

## License

This project is licensed under the MIT License - see the `LICENSE` file for details.
```