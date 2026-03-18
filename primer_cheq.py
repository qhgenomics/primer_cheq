import sys
import subprocess
import glob
import argparse
import os
from collections import defaultdict
import re

# downloads viruses into a single file
def download_virus(taxnum, working_dir, prefix, datasets="datasets"):
    datasets_filename = os.path.join(working_dir, prefix + "_ncbi_dataset.zip")
    subprocess.Popen("{} download virus genome taxon {} --complete-only --filename {}".format(
        datasets, taxnum,datasets_filename
    ), shell=True).wait()
    datasets_unzip_folder = os.path.join(working_dir, prefix + "_downloads")
    subprocess.Popen("unzip -o {f_datasets_filename} -d {f_datasets_unzip_folder} && rm {f_datasets_filename}".format(
        f_datasets_filename=datasets_filename, f_datasets_unzip_folder=datasets_unzip_folder
    ), shell=True).wait()
    fasta_file = os.path.join(datasets_unzip_folder, "ncbi_dataset" , "data", "genomic.fna")
    if not os.path.exists(fasta_file):
        sys.stderr.write("Something went wrong downloading using datasets, please check above for error messages.\n")
        sys.exit(0)
    return(fasta_file)



# downloads bacteria into multiple files
def download_bac(taxnum, working_dir, prefix, datasets="datasets"):
    datasets_filename = os.path.join(working_dir, prefix + "_ncbi_dataset.zip")
    subprocess.Popen("{} download genome taxon {} --assembly-source RefSeq --filename {}".format(
        datasets, taxnum, datasets_filename
    ), shell=True).wait()
    datasets_unzip_folder = os.path.join(working_dir, prefix + "_downloads")
    subprocess.Popen("unzip -o {f_datasets_filename} -d {f_datasets_unzip_folder} && rm {f_datasets_filename}".format(
        f_datasets_filename=datasets_filename, f_datasets_unzip_folder=datasets_unzip_folder
    ), shell=True).wait()
    fasta_files = glob.glob(os.path.join(datasets_unzip_folder, "ncbi_dataset", "*", "*.fna"))
    if len(fasta_files) < 1:
        sys.stderr.write("Something went wrong downloading using datasets, please check above for error messages.\n")
        sys.exit(0)
    return(fasta_files)


def get_db_folder(db_folder):
    fasta_files = []
    for i in os.listdir(db_folder):
        fasta_files.append(os.path.join(db_folder, i))
    if len(fasta_files) < 1:
        sys.stderr.write("Something went wrong finding fasta files, please check your directory.\n")
        sys.exit(0)
    return(fasta_files)


def get_db_glob(theglob):
    fasta_files = glob.glob(theglob)
    if len(fasta_files) < 1:
        sys.stderr.write("Something went wrong finding fasta files, please check your glob.\n")
        sys.exit(0)
    return(fasta_files)


def get_db_fastas(fasta_list, working_dir, prefix):
    reference_count = 0
    with open(os.path.join(working_dir, "{}_db.fasta".format(prefix)), 'a') as o:
        for fasta in fasta_list:
            filename = os.path.basename(fasta)
            filename = os.path.splitext(filename)[0]
            filename.replace("|", "_")
            reference_count += 1
            with open(fasta) as f:
                for line in f:
                    if line.startswith(">"):
                        contig_name = line.split()[0][1:]
                        contig_name.replace("|", "_")
                        o.write(">{}|{}\n".format(filename, contig_name))
                    else:
                        o.write(line.lower())
    return(reference_count)

def get_db_fasta(fasta_file, working_dir, prefix):
    reference_count = 0
    with open(os.path.join(working_dir, "{}_db.fasta".format(prefix)), 'a') as o, open(fasta_file) as f:
        for line in f:
            if line.startswith(">"):
                o.write(line.split()[0] + "|na\n")
                reference_count += 1
            else:
                o.write(line.lower())
    return(reference_count)


def parse_cigar(cigar_str):
    # Find all digits and all characters
    digits = re.findall(r'[0-9]+', cigar_str)
    chars = re.findall(r'[A-Z=]+', cigar_str)
    
    # Pair them together as (length, operation)
    return list(zip(map(int, digits), chars))

def align_primers(primer_dict, database, working_dir, prefix, max_indel=2, max_mismatch=5, threads=1, sassy_loc="sassy"):
    primer_fasta = os.path.join(working_dir, prefix + "_primers.fasta")
    with open(primer_fasta, 'w') as o:
        for name, seq in primer_dict.items():
            o.write(">{}\n{}\n".format(name, seq))
    sassy_output = os.path.join(working_dir, prefix + ".sassy.tsv")
    subprocess.Popen("{} search -f {} --threads {} -k {} {} > {}".format(
        sassy_loc, primer_fasta, threads, max([max_indel, max_mismatch]),  database, sassy_output), shell=True).wait()
    outlist = []
    with open(sassy_output) as f:
        f.readline()
        for line in f:
            primer, ref, cost, strand, start, end, match_region, cigar = line.rstrip().split("\t")
            cigar = parse_cigar(cigar)
            mismatch_count = sum([x[0] for x in cigar if x[1] == "X"])
            insertion_count = sum([x[0] for x in cigar if x[1] == "I"])
            deletion_count = sum([x[0] for x in cigar if x[1] == "D"])
            indel_count = insertion_count + deletion_count
            if mismatch_count <= max_mismatch and indel_count <= max_indel:
                ref, contig = ref.split("|")
                match_region = match_region.lower()
                if cigar[-1][1] == "=":
                    last_1 = False
                else:
                    last_1 = True
                if cigar[-1][0] >= 3 and cigar[-1][1] == "=":
                    last_3 = False
                else:
                    last_3 = True
                if last_1:
                    alert = "HIGH"
                elif indel_count >= 1:
                    alert = "HIGH"
                elif last_3:
                    alert = "MEDIUM"
                elif mismatch_count >= 2:
                    alert = "MEDIUM"
                else:
                    alert = "LOW"
                outlist.append([ref, contig, primer, primer_dict[primer], match_region, alert, 
                str(mismatch_count), str(insertion_count), str(deletion_count), str(last_3), str(last_1)])
    return(outlist)


def create_output(outlist, working_dir, prefix, reference_count):    
    stats_dict = defaultdict(lambda: defaultdict(lambda: "MISSING"))
    single = defaultdict(lambda: set())
    double = defaultdict(lambda: set())
    with open(os.path.join(working_dir, prefix +"_primer_cheq.tsv"), 'w') as report:
        report.write("Reference\tContig\tPrimer\tPrimer_Seq\tMatched_Seq\tAlert\tSubstitution\tInsertion\tDeletion\tLast_3\tLast_1\n")
        for i in outlist:
            ref = i[0]
            primer = i[2]
            alert = i[5]
            if ref in single[primer]:
                double[primer].add(ref)
            else:
                single[primer].add(ref)
            if alert == "LOW":
                stats_dict[primer][ref] = "LOW"
            elif alert == "MEDIUM" and stats_dict[primer][ref] != "LOW":
                stats_dict[primer][ref] = "MEDIUM"
            elif alert == "HIGH" and stats_dict[primer][ref] != "LOW" and stats_dict[primer][ref] != "MEDIUM":
                stats_dict[primer][ref] = "HIGH"
            report.write("\t".join(i) + "\n")
    with open(os.path.join(working_dir, prefix + "_summary.tsv"), 'w') as summary:
        summary.write("Primer\tSingle_target\tMultiple_target\tMissing\tLow\tMedium\tHigh\n")
        for i in stats_dict:
            summary.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(i, len(single[i]) - len(double[i]), len(double[i]),
                                                                reference_count - len(single[i]),
                                                                len([x for x in stats_dict[i].values() if x == "LOW"]),
                                                                len([x for x in stats_dict[i].values() if x == "MEDIUM"]),
                                                                len([x for x in stats_dict[i].values() if x == "HIGH"])))


def get_primer_sequences(primer_fasta):
    primer_dict = {}
    with open(primer_fasta) as f:
        for line in f:
            if line.startswith(">"):
                primer_name = line[1:].rstrip()
                primer_dict[primer_name] = ""
            else:
                primer_seq = line.rstrip()
                primer_dict[primer_name] += primer_seq
    return(primer_dict)

def get_primer_table(primer_table):
    primer_dict = {}
    with open(primer_table) as f:
        for line in f:
            if line.startswith("#") or line.startswith("Primer"):
                continue
            else:
                primer_name = line.split("\t")[0]
                primer_seq = line.split("\t")[1]
                primer_dict[primer_name] = primer_seq
    return(primer_dict)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--primers", help="FASTA file of primers to check.", metavar="primers.fasta")
    parser.add_argument("-t", "--primer_table", help="Table of primers to check.", metavar="primers.tsv")
    parser.add_argument("-s", "--prefix", help="Prefix for output files in workign directory.", required=True, metavar="sample_name")
    parser.add_argument("-w", '--working_directory', help="Path for intermediate and output files.", required=True, metavar="working_dir")
    parser.add_argument("-v", "--ncbi_virus", action="append", help="Download a database of viral references from NCBI using a taxid.", metavar="10244")
    parser.add_argument("-b", "--ncbi_bacteria", action="append", help="Download a database of microbial references from NCBI using a taxid", metavar="590")
    parser.add_argument("-d", "--directory_db", action="append", help="A directory of FASTA files, one for each reference.")
    parser.add_argument("-f", "--fasta_db", action="append", help="A single fasta file, each entry will be treated as its own reference.")
    parser.add_argument("-g", "--glob_db", action="append", help="A glob used to identify FASTA files, each FASTA will be treated as its own reference.")
    parser.add_argument("-i", "--max_indel", type=int, default=2, help="Maximum indels allowed in reported hits.")
    parser.add_argument("-m", "--max_mismatch", type=int, default=5, help="Maximum mismatches allowed in reported hits.")



    args = parser.parse_args()

    if args.primers is None and args.primer_table is None:
        sys.stderr.write("You must provide a FASTA file or a table of primers to check.\n")
        sys.exit(0)
    elif args.primers is not None and args.primer_table is not None:
        sys.stderr.write("You cannot provide both a FASTA file and a table of primers to check.\n")
        sys.exit(0)
    elif args.primers is not None:
        primer_dict = get_primer_sequences(args.primers)
    elif args.primer_table is not None:
        primer_dict = get_primer_table(args.primer_table)

    if not os.path.exists(args.working_directory):
        os.makedirs(args.working_directory)
    elif os.path.exists(args.working_directory) and not os.path.isdir(args.working_directory):
        sys.stderr.write("The working directory you provided is not a directory.\n")
        sys.exit(0)

    if args.ncbi_virus is None and args.ncbi_bacteria is None and args.directory_db is None and args.fasta_db is None and args.glob_db is None:
        sys.stderr.write("You must provide a database of references to check against.\n")
        sys.exit(0)

    fasta_file = os.path.join(args.working_directory, args.prefix + "_db.fasta")
    with open(fasta_file, 'w') as o:
        pass
    reference_count = 0
    if args.ncbi_virus:
        for i in args.ncbi_virus:
            fasta_file = download_virus(i, args.working_directory, args.prefix)
            reference_count += get_db_fasta(fasta_file, args.working_directory, args.prefix)
    if args.ncbi_bacteria:
        for i in args.ncbi_bacteria:
            fasta_files = download_bac(i, args.working_directory, args.prefix)
            reference_count += get_db_fastas(fasta_files, args.working_directory, args.prefix)
    if args.directory_db:
        for i in args.directory_db:
            fasta_files = get_db_folder(i)
            reference_count += get_db_fastas(fasta_files, args.working_directory, args.prefix)
    if args.fasta_db:
        for i in args.fasta_db:
            reference_count += get_db_fasta(i, args.working_directory, args.prefix)
    if args.glob_db:
        for i in args.glob_db:
            fasta_files = get_db_glob(i)
            reference_count += get_db_fastas(fasta_files, args.working_directory, args.prefix)
    database_fasta = os.path.join(args.working_directory, args.prefix + "_db.fasta")
    outlist = align_primers(primer_dict, database_fasta, args.working_directory, args.prefix, args.max_indel, args.max_mismatch)
    create_output(outlist, args.working_directory, args.prefix, reference_count)