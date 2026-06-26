#!/usr/bin/env python3
__version__ = "0.2.1"
__author__ = "Mitchell Sullivan"

import sys
import subprocess
import glob
import argparse
import os
from collections import defaultdict
import re

# downloads viruses into a single file
def download_virus(taxnum, working_dir, prefix, date="all", datasets="datasets", batchnum=5000):
    metadata = subprocess.check_output("{} summary virus genome taxon {} | jq | jq '.reports[] | \"\\(.accession) \\(.isolate.collection_date)\"'".format(datasets, taxnum), shell=True).decode()
    accession_list = []
    fasta_files = []
    for line in metadata.split("\n"):
        if line != "":
            accession, collection_date = line.strip('"').split()
            if date == "all" or collection_date.startswith(date):
                accession_list.append(accession)
    for num in range(0, len(accession_list), batchnum):
        accession_filename = os.path.join(working_dir, "{}_accessions_{}.txt".format(prefix, num//batchnum))
        with open(accession_filename, 'w') as f:
            f.write('\n'.join(accession_list[num:num+batchnum]))
        datasets_filename = os.path.join(working_dir, "{}_ncbi_dataset_{}.zip".format(prefix, num//batchnum))
        subprocess.Popen("{} download virus genome accession --inputfile {} --complete-only --filename {}".format(
            datasets, accession_filename, datasets_filename
        ), shell=True).wait()
        datasets_unzip_folder = os.path.join(working_dir,  "{}_downloads_{}".format(prefix, num//batchnum))
        subprocess.Popen("unzip -o {f_datasets_filename} -d {f_datasets_unzip_folder} && rm {f_datasets_filename}".format(
            f_datasets_filename=datasets_filename, f_datasets_unzip_folder=datasets_unzip_folder), shell=True).wait()
        fasta_file = os.path.join(datasets_unzip_folder, "ncbi_dataset" , "data", "genomic.fna")
        fasta_files.append(fasta_file)
        if not os.path.exists(fasta_file):
            sys.stderr.write("Something went wrong downloading using datasets, please check above for error messages.\n")
            sys.exit(0)
    return(fasta_files)



# downloads bacteria into multiple files
def download_bac(taxnum, working_dir, prefix, date="all", datasets="datasets", batchnum=5000):
    metadata = subprocess.check_output("{} summary genome taxon {} | jq | jq '.reports[] | \"\\(.accession) \\(.assembly_info.biosample.collection_date)\"'".format(datasets, taxnum), shell=True).decode()
    accession_list = []
    fasta_files = []
    for line in metadata.split("\n"):
        if line != "":
            accession, collection_date = line.strip('"').split()
            if date == "all" or collection_date.startswith(date):
                accession_list.append(accession)
    for num in range(0, len(accession_list), batchnum):
        accession_filename = os.path.join(working_dir, "{}_accessions_{}.txt".format(prefix, num//batchnum))
        with open(accession_filename, 'w') as f:
            f.write('\n'.join(accession_list[num:num+batchnum]))
        datasets_filename = os.path.join(working_dir, "{}_ncbi_dataset_{}.zip".format(prefix, num//batchnum))
        subprocess.Popen("{} download genome accession --inputfile {} --filename {}".format(
            datasets, accession_filename, datasets_filename
        ), shell=True).wait()
        datasets_unzip_folder = os.path.join(working_dir,  "{}_downloads_{}".format(prefix, num//batchnum))
        subprocess.Popen("unzip -o {f_datasets_filename} -d {f_datasets_unzip_folder} && rm {f_datasets_filename}".format(
            f_datasets_filename=datasets_filename, f_datasets_unzip_folder=datasets_unzip_folder), shell=True).wait()
        new_fasta_files = glob.glob(os.path.join(datasets_unzip_folder, "ncbi_dataset" , "data", "*", "*.fna"))
        fasta_files += new_fasta_files
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


def get_db_fastas(fasta_list, working_dir, prefix, amb_bases):
    reference_count = 0
    with open(os.path.join(working_dir, "{}_db.fasta".format(prefix)), 'a') as o:
        for fasta in fasta_list:
            filename = os.path.basename(fasta)
            filename = os.path.splitext(filename)[0]
            filename.replace("|", "_")
            reference_count += 1
            with open(fasta) as f:
                seq = None
                for line in f:
                    if line.startswith(">"):
                        contig_name = line.split()[0][1:]
                        contig_name.replace("|", "_")
                        o.write(">{}|{}\n".format(filename, contig_name))
                        if not seq is None:
                            amb_bases[filename] = amb_bases[filename].union(set([i for i, char in enumerate(seq) if char == 'n']))
                        seq = ''
                    else:
                        o.write(line.lower())
                        seq += line.lower().rstrip()
                amb_bases[filename] = amb_bases[filename].union(set([i for i, char in enumerate(seq) if char == 'n']))
    return amb_bases

def get_db_fasta(fasta_file, working_dir, prefix, amb_bases):
    reference_count = 0
    with open(os.path.join(working_dir, "{}_db.fasta".format(prefix)), 'a') as o, open(fasta_file) as f:
        seq = None
        for line in f:
            if line.startswith(">"):
                o.write(line.split()[0] + "|na\n")
                reference_count += 1
                if not seq is None:
                    amb_bases[refname] = set([i for i, char in enumerate(seq) if char == 'n'])
                refname = line.split()[0][1:]
                seq = ''
            else:
                o.write(line.lower())
                seq += line.lower().rstrip()
        amb_bases[refname] = set([i for i, char in enumerate(seq) if char == 'n'])
    return amb_bases


def parse_cigar(cigar_str):
    # Find all digits and all characters
    digits = re.findall(r'[0-9]+', cigar_str)
    chars = re.findall(r'[A-Z=]+', cigar_str)
    
    # Pair them together as (length, operation)
    return list(zip(map(int, digits), chars))

def align_primers(primer_dict, database, working_dir, prefix, max_indel, max_mismatch, indel_mult, threads, sassy_loc="sassy"):
    primer_fasta = os.path.join(working_dir, prefix + "_primers.fasta")
    with open(primer_fasta, 'w') as o:
        for name, seq in primer_dict.items():
            o.write(">{}\n{}\n".format(name, seq))
    sassy_output = os.path.join(working_dir, prefix + ".sassy.tsv")
    subprocess.Popen("{} search -f {} --threads {} -k {} {} > {}".format(
        sassy_loc, primer_fasta, threads, max_mismatch,  database, sassy_output), shell=True).wait()
    outlist = []
    start_locations = defaultdict(lambda: defaultdict(lambda: 0))
    with open(sassy_output) as f:
        f.readline()
        for line in f:
            primer, ref_temp, cost, strand, start, end, match_region, cigar = line.rstrip().split("\t")
            cigar = parse_cigar(cigar)
            alignment = ""
            pos = 0
            for i in cigar:
                if i[1] == "=":
                    alignment += '.' * i[0]
                    pos += i[0]
                elif i[1] == "X":
                    for num in range(i[0]):
                        alignment += match_region[pos].lower()
                        pos += 1
                elif i[1] == "I":
                    alignment += '-' * i[0]
                elif i[1] == "D":
                    for num in range(i[0]):
                        alignment += match_region[pos].upper()
                        pos += 1
            mismatch_count = sum([x[0] for x in cigar if x[1] == "X"])
            insertion_count = sum([x[0] for x in cigar if x[1] == "I"])
            deletion_count = sum([x[0] for x in cigar if x[1] == "D"])
            indel_count = insertion_count + deletion_count
            if mismatch_count + indel_count * indel_mult <= max_mismatch and indel_count <= max_indel:
                ref = "|".join(ref_temp.split("|")[:-1])
                contig = ref_temp.split("|")[-1]
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
                outlist.append([ref, contig, start, strand, primer, primer_dict[primer], alignment, alert, 
                str(mismatch_count), str(insertion_count), str(deletion_count), str(last_3), str(last_1)])
                start_locations[primer][int(start)] += 1
    return outlist, start_locations



def create_output(outlist, working_dir, prefix, reference_count, amb_references):    
    stats_dict = defaultdict(lambda: defaultdict(lambda: "MISSING"))
    single = defaultdict(lambda: set())
    double = defaultdict(lambda: set())
    alignment_dict = defaultdict(lambda: defaultdict(lambda: [0, [], None]))
    with open(os.path.join(working_dir, prefix +"_primer_cheq.tsv"), 'w') as report:
        report.write("Reference\tContig\tPosition\tstrand\tPrimer\tPrimer_Seq\tMatched_Seq\tAlert\tSubstitution\tInsertion\tDeletion\tLast_3\tLast_1\n")
        for i in outlist:
            ref = i[0]
            primer = i[4]
            alignment = i[6]
            alert = i[7]
            alignment_dict[primer][alignment][0] += 1
            alignment_dict[primer][alignment][1].append(ref)
            alignment_dict[primer][alignment][2] = alert
            if ref in amb_references[primer]:
                amb_references[primer].remove(ref)
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
    with open(os.path.join(working_dir, prefix + "_alignment_summary.tsv"), 'w') as alignment_summary:
        alignment_summary.write("Primer\tMatched_Seq\tCount\tAlert\tReferences\n")
        for primer in alignment_dict:
            for alignment in alignment_dict[primer]:
                count = alignment_dict[primer][alignment][0]
                alert = alignment_dict[primer][alignment][2]
                refs = ",".join(alignment_dict[primer][alignment][1])
                alignment_summary.write("{}\t{}\t{}\t{}\t{}\n".format(primer, alignment, count, alert, refs))
    with open(os.path.join(working_dir, prefix + "_summary.tsv"), 'w') as summary:
        summary.write("Primer\tSingle_target\tMultiple_target\tMissing\tamb_in_ref\tLow\tMedium\tHigh\n")
        for i in stats_dict:
            summary.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(i, len(single[i]) - len(double[i]), len(double[i]),
                                                                reference_count - len(single[i]) - len(amb_references[i]),
                                                                len(amb_references[i]),
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
    return primer_dict

def get_primer_table(primer_table):
    primer_dict = {}
    primer_info = {}
    with open(primer_table) as f:
        for line in f:
            if line.startswith("#") or line.startswith("Primer"):
                continue
            else:
                if line.count("\t") == 4:
                    primer_name, primer_set, primer_type, primer_seq, expected_product_size = line.rstrip("\n").split("\t")
                    comments = ""
                elif line.count("\t") == 5:
                    primer_name, primer_set, primer_type, primer_seq, expected_product_size, comments = line.rstrip("\n").split("\t")
                else:
                    sys.stderr.write("Your primer table is not formatted correctly, please check the documentation for formatting instructions.\n")
                    sys.exit(0)
                primer_dict[primer_name] = primer_seq
                primer_info[primer_name] = (primer_set, primer_type, expected_product_size, comments)
    return primer_dict, primer_info


def create_product_size_output(align_list, primer_info, working_dir, prefix, primer_dict, max_primer_dist, reference_count):
    binding_sites = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: [])))
    for i in align_list:
        ref = i[0]
        contig = i[1]
        pos = int(i[2])
        strand = i[3]
        primer = i[4]
        alert = i[7]
        primerset = primer_info[primer][0]
        binding_sites[primerset][ref][contig].append((pos, strand, primer, alert))
    outlist = []
    filtered_set = set()
    for primerset in binding_sites:
        for ref in binding_sites[primerset]:
            for contig in binding_sites[primerset][ref]:
                bs = binding_sites[primerset][ref][contig]
                bs.sort()
                primer_strand = []
                templist = []
                for i in bs:
                    pos, strand, primer, alert = i
                    primer_type = primer_info[primer][1]
                    expected_product_size = int(primer_info[primer][2])
                    if strand == "+" and primer_type == "F":
                        primer_strand = [(primer, primer_type, strand, pos, alert)]
                    elif strand == '+' and primer_type == "R":
                         primer_strand = [(primer, primer_type, strand, pos, alert)]
                    elif strand == "-" and primer_type == "R" and primer_strand != [] and primer_strand[0][1] == "F" \
                      and primer_strand[0][2] == "+" and pos - primer_strand[-1][3] <= max_primer_dist:
                        primer_strand.append((primer, primer_type, strand, pos, alert))
                        templist.append(primer_strand)
                        primer_strand = []
                    elif strand == "-" and primer_type == "F" and primer_strand != [] and primer_strand[0][1] == "R" \
                      and primer_strand[0][2] == "-" and pos - primer_strand[-1][3] <= max_primer_dist:
                        primer_strand.append((primer, primer_type, strand, pos, alert))
                        templist.append(primer_strand)
                        primer_strand = []
                    elif primer_type == "P" and primer_strand != [] and pos - primer_strand[-1][3] <= max_primer_dist:
                        primer_strand.append((primer, primer_type, strand, pos, alert))
                    else:
                        primer_strand = []
                for i in templist:
                    temp = []
                    maxalert = "LOW"
                    for j in i:
                        filtered_set.add((j[0], ref, contig, j[2], j[3]))
                        temp.append(":".join(map(str, j)))
                        if j[4] == "HIGH":
                            maxalert = "HIGH"
                        elif j[4] == "MEDIUM" and maxalert != "HIGH":
                            maxalert = "MEDIUM"                    
                    start = i[0][3]
                    size = abs(i[0][3] - i[-1][3]) + len(primer_dict[i[0][0]])
                    if size < expected_product_size * 0.8 or size > expected_product_size * 1.2:
                        product_size_alert = "YES"
                    else:
                        product_size_alert = "NO"
                    strandstr = ";".join(temp)
                    outlist.append((primerset, maxalert, size, expected_product_size, product_size_alert, ref, contig, start, strandstr))
    filtered_list = []
    for i in align_list:
        ref = i[0]
        contig = i[1]
        pos = int(i[2])
        strand = i[3]
        primer = i[4]
        if (primer, ref, contig, strand, pos) in filtered_set:
            filtered_list.append(i)  
    summary_dict = defaultdict(lambda: defaultdict(lambda: [[], 0, 0, 0]))
    with open(os.path.join(working_dir, prefix + "_product_size.tsv"), 'w') as product_size_report:
        product_size_report.write("Primer_Set\tMax_Alert\tProduct_Size\tExpected_Product_Size\tproudct_size_alert\tReference\tContig\tStart\tStrand_Info\n")
        for i in outlist:
            primer_set, maxalert, size, expected_product_size, product_size_alert, ref, contig, start, strandstr = i
            product_size_report.write("\t".join(map(str, i)) + "\n")
            summary_dict[primer_set][ref][0].append(size)
            if maxalert == "HIGH":
                summary_dict[primer_set][ref][3] += 1
            elif maxalert == "MEDIUM":
                summary_dict[primer_set][ref][2] += 1
            elif maxalert == "LOW":
                summary_dict[primer_set][ref][1] += 1
    with open(os.path.join(working_dir, prefix + "_product_size_summary.tsv"), 'w') as product_size_summary:
        product_size_summary.write("Primer_Set\tsingle_product\ttwo_products\tthree_or_more_products\tmissing\tlow_alert\tmedium_alert\thigh_alert\taverage_product_size_all\taverage_smallest_product_size\taverage_secondary_product_size\n")
        for ps in summary_dict:
            sizes_all = []
            sizes_smallest = []
            sizes_large = []
            alert_counts = [0, 0, 0]
            count_counts = [0, 0, 0]
            for ref in summary_dict[ps]:
                sizes = summary_dict[ps][ref][0]
                sizes.sort()
                sizes_all += sizes
                sizes_smallest.append(sizes[0])
                sizes_large += sizes[1:]
                alert_counts[0] += summary_dict[ps][ref][1]
                alert_counts[1] += summary_dict[ps][ref][2]
                alert_counts[2] += summary_dict[ps][ref][3]
                if len(sizes) == 1:
                    count_counts[0] += 1
                elif len(sizes) == 2:
                    count_counts[1] += 1
                elif len(sizes) >= 3:
                    count_counts[2] += 1
            product_size_summary.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                ps, count_counts[0], count_counts[1], count_counts[2], reference_count - sum(count_counts), 
                alert_counts[0], alert_counts[1], alert_counts[2], sum(sizes_all)/len(sizes_all), 
                sum(sizes_smallest)/len(sizes_smallest), 
                sum(sizes_large)/len(sizes_large) if len(sizes_large) > 0 else "NA"))
    return filtered_list

        

def get_amb_in_ref(primer_dict, amb_bases, start_locations, fractions_threshold=0.05):
    amb_refs = {}
    for primer in primer_dict:
        amb_refs[primer] = set()
        total_starts = sum(start_locations.get(primer, {}).values())
        for start in start_locations[primer]:
            if start_locations[primer][start] / total_starts >= fractions_threshold:
                start_bases = set(range(start, start + len(primer_dict[primer])))
                for i in amb_bases:
                    if len(start_bases.intersection(amb_bases[i])) > 3:
                        amb_refs[primer].add(i)
    return amb_refs

def create_xls(working_dir, prefix, xls_file):
    try:
        import pandas as pd
    except ImportError:
        sys.stderr.write("pandas is not installed, please install it to use the xls output option.\n")
        sys.exit(0)
    the_files = []
    for i in ["_primer_cheq.tsv", "_alignment_summary.tsv", "_summary.tsv", "_product_size.tsv", "_product_size_summary.tsv"]:
        path = os.path.join(working_dir, prefix + i)
        if os.path.exists(path):
            the_files.append(path)
    writer = pd.ExcelWriter(xls_file, engine='xlsxwriter')
    for f in the_files:
        df = pd.read_csv(f, delimiter="\t", index_col=False)
        df.to_excel(writer,sheet_name=f[len(os.path.join(working_dir, prefix))+1:-4],index=False)
    writer.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='primer_cheq', description='Check primer sequences against a database of reference genomes to predict potential mismatches and their impacts on primer binding.',
                                     epilog='Thankyou for using primer_cheq! If you have any issues or feature requests, please submit them at https://github.com/qhgenomics/primer_cheq/issues')
    parser.add_argument("-p", "--primers", help="FASTA file of primers to check.", metavar="primers.fasta")
    parser.add_argument("-T", "--primer_table", help="Table of primers to check.", metavar="primers.tsv")
    parser.add_argument("-s", "--prefix", help="Prefix for output files in workign directory.", required=True, metavar="sample_name")
    parser.add_argument("-w", '--working_directory', help="Path for intermediate and output files.", required=True, metavar="working_dir")
    parser.add_argument("-V", "--ncbi_virus", action="append", help="Download a database of viral references from NCBI using a taxid.", metavar="10244")
    parser.add_argument("-b", "--ncbi_bacteria", action="append", help="Download a database of microbial references from NCBI using a taxid", metavar="590")
    parser.add_argument("-Y", "--year", default="all", help="Year of viral or bacterial collection to download when using --ncbi_virus or --ncbi_bacteria, default is all.")
    parser.add_argument("-d", "--directory_db", action="append", help="A directory of FASTA files, one for each reference.")
    parser.add_argument("-f", "--fasta_db", action="append", help="A single fasta file, each entry will be treated as its own reference.")
    parser.add_argument("-g", "--glob_db", action="append", help="A glob used to identify FASTA files, each FASTA will be treated as its own reference.")
    parser.add_argument("-i", "--max_indel", type=int, default=2, help="Maximum indels allowed in reported hits. Default value is %(default)s.")
    parser.add_argument("-m", "--max_mismatch", type=int, default=4, help="Maximum mismatches allowed in reported hits, this includes substitutions and indels. Default value is %(default)s.")
    parser.add_argument("-I", "--indel_mult", type=int, default=2, help="Multiplier for indel penalties in alignment scoring. Default value is %(default)s. This value is only used to calculate max_mismatch score.")
    parser.add_argument("--sassy_loc", default="sassy", help="Path to sassy binary if not in PATH.")
    parser.add_argument("-t", "--threads", default=1, type=int, help="Threads to use for sassy alignment, default is %(default)s.")
    parser.add_argument("-D", "--dataset_loc", default="dataset", help="Path to NCBI datasets binary if not in PATH.")
    parser.add_argument("-P", "--max_primer_dist", type=int, default=5000, help="Maximum distance between primers or primer and probe to produce a product.")
    parser.add_argument("--predict_amb_in_ref", default=0.05, type=float, help="Predict whether primer site is missing in reference due to ambiguous bases.")
    parser.add_argument("--only_product_forming", action="store_true", help="Only report alignments that are predicted to form a product based on max_primer_dist.")
    parser.add_argument("-x", "--xlsx", help="Create an excel table with all output in tabs. Requires openpyxl to be installed.")
    parser.add_argument("-v", "--version", action="version", version='%(prog)s ' + __version__)

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
        primer_dict, primer_info = get_primer_table(args.primer_table)

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
    amb_bases = defaultdict(lambda:set())
    reference_count = 0
    if args.ncbi_virus:
        for i in args.ncbi_virus:
            fasta_files = download_virus(i, args.working_directory, args.prefix)
            for j in fasta_files:
                amb_bases = get_db_fasta(j, args.working_directory, args.prefix, amb_bases)
    if args.ncbi_bacteria:
        for i in args.ncbi_bacteria:
            fasta_files = download_bac(i, args.working_directory, args.prefix)
            amb_bases = get_db_fastas(fasta_files, args.working_directory, args.prefix, amb_bases)
    if args.directory_db:
        for i in args.directory_db:
            fasta_files = get_db_folder(i)
            amb_bases = get_db_fastas(fasta_files, args.working_directory, args.prefix, amb_bases)
    if args.fasta_db:
        for i in args.fasta_db:
            amb_bases = get_db_fasta(i, args.working_directory, args.prefix, amb_bases)
    if args.glob_db:
        for i in args.glob_db:
            fasta_files = get_db_glob(i)
            amb_bases = get_db_fastas(fasta_files, args.working_directory, args.prefix, amb_bases)
    reference_count = len(amb_bases)
    database_fasta = os.path.join(args.working_directory, args.prefix + "_db.fasta")
    outlist, start_locations = align_primers(primer_dict, database_fasta, args.working_directory, args.prefix, args.max_indel, args.max_mismatch, args.indel_mult, args.threads, args.sassy_loc)
    if not args.primer_table is None:
        filtered_list = create_product_size_output(outlist, primer_info, args.working_directory, args.prefix, primer_dict, args.max_primer_dist, reference_count)
        if args.only_product_forming:
            outlist = filtered_list
    amb_references = get_amb_in_ref(primer_dict, amb_bases, start_locations, args.predict_amb_in_ref)
    create_output(outlist, args.working_directory, args.prefix, reference_count, amb_references)
    if not args.xlsx is None:
        create_xls(args.working_directory, args.prefix, args.xls)