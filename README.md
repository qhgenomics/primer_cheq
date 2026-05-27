# primer_cheq

Tool for checking primer coverage and drift.

## Description

`primer_cheq.py` is a Python script designed to check the coverage and drift of primer sequences against viral and bacterial genomes. It downloads genome data, processes it, and uses BLAST to align primer sequences to the genomes, providing detailed information on mismatches and potential issues.

## Requirements

- Python 3.x
- BLAST+ (installed and available in the system PATH)
- NCBI Datasets command-line tool (installed and available in the system PATH)
- unzip bash tool (installed and available in the system PATH)

## Installation

1. Clone the repository:
    ```sh
    git clone https://github.com/yourusername/primer_cheq.git
    cd primer_cheq
    ```

2. Create the conda environment:
    ```sh
    conda env create -f environment.yml
    ```
3. Ensure unzip and datasets are in your path

## Usage

```sh
conda activate primer_cheq
python primer_cheq.py -p primer.fasta -o sample_name -w working_dir [options]
```

### Arguments

#### Required Arguments
- `-p`, `--primers`: FASTA file of primers to check (mutually exclusive with `-T`).
- `-T`, `--primer_table`: Table of primers to check in TSV format (mutually exclusive with `-p`).
- `-s`, `--prefix`: Prefix for output files in the working directory.
- `-w`, `--working_directory`: Path for intermediate and output files.

#### Database Source (one or more required)
- `-V`, `--ncbi_virus`: Download a database of viral references from NCBI using a taxid (e.g., `10244` for influenza).
- `-b`, `--ncbi_bacteria`: Download a database of microbial references from NCBI using a taxid (e.g., `590` for Salmonella).
- `-d`, `--directory_db`: A directory of FASTA files, one for each reference.
- `-f`, `--fasta_db`: A single FASTA file, each entry will be treated as its own reference.
- `-g`, `--glob_db`: A glob pattern (e.g., `/path/to/data/*.fasta`) to identify FASTA files; each FASTA will be treated as its own reference.

**Note:** Arguments can be used multiple times. For example: `-f genomes.fasta -f another_genome.fasta -V 10244` will create a database using both FASTA files and download additional viral sequences from NCBI.

#### Optional Arguments
- `-Y`, `--year`: Year of viral or bacterial collection to download when using `--ncbi_virus` or `--ncbi_bacteria`. Default is `all`.
- `-i`, `--max_indel`: Maximum number of indels allowed in reported hits. Default is `2`.
- `-m`, `--max_mismatch`: Maximum mismatches allowed in reported hits (includes substitutions and indels). Default is `4`.
- `-I`, `--indel_mult`: Multiplier for indel penalties in alignment scoring. Default is `2`. This value is only used to calculate the max_mismatch score.
- `-t`, `--threads`: Number of threads to use for sassy alignment. Default is `1`.
- `-P`, `--max_primer_dist`: Maximum distance between primers or primer and probe to produce a product. Default is `5000` bp.
- `--sassy_loc`: Path to sassy binary if not in PATH. Default is `sassy`.
- `-D`, `--dataset_loc`: Path to NCBI datasets binary if not in PATH. Default is `datasets`.
- `--predict_amb_in_ref`: Threshold for predicting whether primer site is missing in reference due to ambiguous bases. Default is `0.05`.
- `--only_product_forming`: If set, only report alignments that are predicted to form a product based on max_primer_dist.

### Example

```sh
python primer_cheq.py -p primers.fasta -s my_sample -w ./workdir -V 10244
```

This command will check the primers in `primers.fasta` against viral genomes with taxid `10244` (influenza) downloaded from NCBI, using `my_sample` as the prefix for output files and `./workdir` as the working directory.

## Primer Input Formats

### FASTA Format (using `-p` or `--primers`)

A standard FASTA file with primer sequences:

```fasta
>Primer_F1
ACGTACGTACGTACGT
>Primer_R1
TGCATGCATGCATGCA
>Probe_P1
ACGTACGTACGTACGT
```

### Table Format (using `-T` or `--primer_table`)

A TSV file with primer information. Columns should be:
- `Primer`: Name of the primer
- `Primer_Set`: Name of the primer set or assay (used for product size detection)
- `Type`: Type of primer (`F` for forward, `R` for reverse, `P` for probe)
- `Sequence`: DNA sequence of the primer
- `Expected_Product_Size`: Expected amplicon size in bp (integer, used only for F/R pairs)
- `Comments` (optional): Additional notes about the primer

Example primer table:

| Primer | Primer_Set | Type | Sequence | Expected_Product_Size | Comments |
|--------|-----------|------|----------|----------------------|----------|
| InfluA_F1 | InfluA_Set1 | F | ACGTACGTACGTACGT | 150 | Forward primer |
| InfluA_R1 | InfluA_Set1 | R | TGCATGCATGCATGCA | 150 | Reverse primer |
| InfluA_Probe1 | InfluA_Set1 | P | ACGTACGTACGTACGT | 150 | TaqMan probe |
| InfluB_F1 | InfluB_Set1 | F | CGACGACGACGACGAC | 200 | Forward primer |
| InfluB_R1 | InfluB_Set1 | R | GCTGCTGCTGCTGCTG | 200 | Reverse primer |

## Output

The script will generate up to five output files in the working directory, depending on the input:

### \<prefix\>_primer_cheq.tsv

List of all primer hits found in references. This file contains detailed alignment information for every primer match.

**Columns:**

| column | variable       | description                                           |
|--------|----------------|-------------------------------------------------------|
| 1      | `Reference`    | Name of the reference sequence.                       |
| 2      | `Contig`       | Name of the contig where the primer was found.        |
| 3      | `Primer`       | Name of the primer sequence.                          |
| 4      | `Primer_Seq`   | Sequence of the primer.                               |
| 5      | `Alignment`    | Alignment of the primer to the reference contig (`.` = match, letter = mismatch). |
| 6      | `Position`     | Start position of the alignment on the contig.        |
| 7      | `Strand`       | Strand where primer aligns (`+` or `-`).              |
| 8      | `Alert`        | Primer match alert (LOW, MEDIUM, or HIGH).            |
| 9      | `Substitution` | Number of substitution mismatches.                    |
| 10     | `Insertion`    | Number of insertion mismatches.                       |
| 11     | `Deletion`     | Number of deletion mismatches.                        |
| 12     | `Last_3`       | Whether there is a mutation in last 3 bases (0 or 1). |
| 13     | `Last_1`       | Whether there is a mutation in last base (0 or 1).    |

**Alert Levels:**
- **LOW**: < 2 substitutions, no indels, no mutations in last 3 bases
- **MEDIUM**: ≥ 3 substitutions OR substitution in last 3 bases (but no indels, no mutations in last base)
- **HIGH**: Indel present OR mutation in last base

**Example output:**

```
Reference	Contig	Primer	Primer_Seq	Alignment	Position	Strand	Alert	Substitution	Insertion	Deletion	Last_3	Last_1
Influenza_A_H1N1	segment_1	InfluA_F1	ACGTACGTACGTACGT	................ 	1250	+	LOW	0	0	0	0	0
Influenza_A_H3N2	segment_1	InfluA_F1	ACGTACGTACGTACGT	......x......... 	1250	+	LOW	1	0	0	0	0
Influenza_B	segment_1	InfluA_F1	ACGTACGTACGTACGT	....x.......x... 	1250	+	MEDIUM	2	0	0	1	0
```

### \<prefix\>_alignment_summary.tsv

Summary of unique primer alignments, grouping identical alignment patterns across references.

**Columns:**

| column | variable       | description                                           |
|--------|----------------|-------------------------------------------------------|
| 1      | `Primer`       | Name of the primer.                                   |
| 2      | `Matched_Seq`  | The aligned sequence with mismatches shown.           |
| 3      | `Count`        | Number of references where this exact alignment occurred. |
| 4      | `Alert`        | Alert level for this alignment pattern.               |
| 5      | `References`   | Comma-separated list of reference names.              |

**Example output:**

```
Primer	Matched_Seq	Count	Alert	References
InfluA_F1	................	10	LOW	Influenza_A_H1N1,Influenza_A_H1N1_PR8,Influenza_A_H1N1_CA09
InfluA_F1	......x.........	5	LOW	Influenza_A_H3N2,Influenza_A_H3N2_A_Perth
InfluA_F1	....x.......x...	2	MEDIUM	Influenza_B_Lee,Influenza_B_Marin
```

### \<prefix\>_summary.tsv

Summary statistics for each primer across all references.

**Columns:**

| column | variable           | description                                             |
|--------|--------------------|---------------------------------------------------------|
| 1      | `Primer`           | Name of the primer sequence.                            |
| 2      | `Single_target`    | Number of references with exactly one primer hit.       |
| 3      | `Multiple_target`  | Number of references with multiple primer hits.         |
| 4      | `Missing`          | Number of references with no primer hits.               |
| 5      | `amb_in_ref`       | Number of references where primer site may be missing due to ambiguous bases. |
| 6      | `Low`              | Number of references where best hit was a LOW alert.    |
| 7      | `Medium`           | Number of references where best hit was a MEDIUM alert. |
| 8      | `High`             | Number of references where best hit was a HIGH alert.   |

**Example output:**

```
Primer	Single_target	Multiple_target	Missing	amb_in_ref	Low	Medium	High
InfluA_F1	15	2	1	0	17	0	0
InfluA_R1	14	3	1	0	16	1	0
InfluA_Probe1	12	4	2	0	15	1	0
InfluB_F1	10	1	7	0	9	2	0
```

### \<prefix\>_product_size.tsv

Predicted amplicon products based on primer pair positions (only generated when using primer table with F/R pairs).

**Columns:**

| column | variable              | description                                             |
|--------|----------------------|---------------------------------------------------------|
| 1      | `Primer_Set`         | Name of the primer set.                                 |
| 2      | `Max_Alert`          | Maximum alert level among primers in the pair.          |
| 3      | `Product_Size`       | Predicted amplicon size in bp.                          |
| 4      | `Expected_Product_Size` | Expected amplicon size from primer table.             |
| 5      | `product_size_alert` | YES if product size is outside ±20% of expected, NO otherwise. |
| 6      | `Reference`          | Name of the reference sequence.                         |
| 7      | `Contig`             | Name of the contig.                                     |
| 8      | `Start`              | Start position of the amplicon.                         |
| 9      | `Strand_Info`        | Details of primer positions and alerts (colon-separated). |

**Example output:**

```
Primer_Set	Max_Alert	Product_Size	Expected_Product_Size	product_size_alert	Reference	Contig	Start	Strand_Info
InfluA_Set1	LOW	150	150	NO	Influenza_A_H1N1	segment_1	1250	InfluA_F1:F:+:1250:LOW:InfluA_R1:R:-:1400:LOW
InfluA_Set1	LOW	145	150	NO	Influenza_A_H3N2	segment_1	1250	InfluA_F1:F:+:1250:LOW:InfluA_R1:R:-:1395:LOW
InfluA_Set1	MEDIUM	160	150	YES	Influenza_B	segment_1	1250	InfluA_F1:F:+:1250:MEDIUM:InfluA_R1:R:-:1410:LOW
InfluB_Set1	LOW	200	200	NO	Influenza_B_Lee	segment_3	2100	InfluB_F1:F:+:2100:LOW:InfluB_R1:R:-:2300:LOW
```

### \<prefix\>_product_size_summary.tsv

Summary of primer sets and their predicted amplicon products across all references.

**Columns:**

| column | variable                      | description                                             |
|--------|------------------------------|---------------------------------------------------------|
| 1      | `Primer_Set`                 | Name of the primer set.                                 |
| 2      | `single_product`             | Number of references with exactly one product detected. |
| 3      | `two_products`               | Number of references with two products detected.        |
| 4      | `three_or_more_products`     | Number of references with three or more products.       |
| 5      | `missing`                    | Number of references with no products detected.         |
| 6      | `low_alert`                  | Count of products with LOW alert level.                 |
| 7      | `medium_alert`               | Count of products with MEDIUM alert level.              |
| 8      | `high_alert`                 | Count of products with HIGH alert level.                |
| 9      | `average_product_size_all`   | Average size of all detected products.                  |
| 10     | `average_smallest_product_size` | Average size of the smallest product per reference.   |
| 11     | `average_secondary_product_size` | Average size of secondary products (when multiple detected). |

**Example output:**

```
Primer_Set	single_product	two_products	three_or_more_products	missing	low_alert	medium_alert	high_alert	average_product_size_all	average_smallest_product_size	average_secondary_product_size
InfluA_Set1	14	2	0	2	16	0	0	149.8	149.8	155.5
InfluB_Set1	10	1	0	7	10	1	0	200.2	200.2	205.0
```

## License

This project is licensed under the MIT License - see the `LICENSE` file for details.
