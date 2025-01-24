import sys
import subprocess
import glob
import argparse
import os
from collections import defaultdict
# read in primer sequences
def get_primer_sequences(primers):
    with open(primers) as f:
        primer_dict = {}
        for line in f:
            if line.startswith(">"):
                name = line.rstrip()[1:]
                primer_dict[name] = ""
            else:
                primer_dict[name] += line.rstrip().lower()
    return primer_dict


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
    fasta_file = os.path.join(datasets_unzip_folder, "ncbi_dataset" , "data" "genomic.fna")
    if not os.path.exists(fasta_file):
        sys.stderr.write("Something went wrong downloading using datasets, please check above for error messages.\n")
        sys.exit(0)
    return(fasta_file)



# downloads bacteria into multiple files and then puts them into a single file with correct names
# returns a dictionary of names that will
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
    with open(os.path.join(working_dir, "{}_db.fasta".format(prefix)), 'a') as o:
        for fasta in fasta_list:
            filename = os.path.basename(fasta)
            filename = os.path.splitext(filename)[0]
            with open(fasta) as f:
                for line in f:
                    if line.startswith(">"):
                        contig_name = line.split()[0][1:]
                        o.write(">{}|{}\n".format(filename, contig_name))
                    else:
                        o.write(line)

def get_db_fasta(fasta_file, working_dir, prefix):
    with open(os.path.join(working_dir, "{}_db.fasta".format(prefix)), 'a') as o, open(fasta_file) as f:
        for line in f:
            if line.startswith(">"):
                o.write(line.split()[0] + "|na\n")
            else:
                o.write(line)



def blast_primers(database, primer_dict, working_dir, prefix, blastn_loc="blastn"):
    num = 0
    subject_names = {}
    outlist = []
    reference_count = set()
    with open(os.path.join(working_dir, "{}_db.fasta".format(prefix))) as f:
        for line in f:
            if line.startswith(">"):
                num += 1
                subject_names["Subject_{}".format(num)] = (line.split("|")[0][1:], "|".join(line.rstrip().split("|")[1:]))
                reference_count.add(line.split("|")[0])
    reference_count = len(reference_count)
    for primer, primer_seq in primer_dict.items():
        sys.stderr.write(primer + "\n")
        single_primer_file = os.path.join(working_dir, prefix + ".tmp.fa")
        tmp_alignment = os.path.join(working_dir, prefix + '.tmp.aln')
        with open(single_primer_file, "w") as o:
            o.write(">{}\n{}".format(primer, primer_seq))
        subprocess.Popen(
            "{} -max_target_seqs 1000000 -query {} -subject {} -outfmt 3 -task blastn-short > {}".format(
                blastn_loc, single_primer_file, database, tmp_alignment), shell=True).wait()
        qdict, mutdict = {}, {}
        # parse the alignment
        with open(tmp_alignment) as f:
            while not line.startswith("Query_"):
                line = f.readline()
            # read the Query line
            refseq = line.split()[2]
            refstart = line.find(refseq)
            refend = refstart + len(refseq)
            lastseq = None
            actual_bases = 0
            last3pos = set()
            sampledict = {}
            qdict = {}
            for num, i in enumerate(refseq[::-1]):
                if i != '-':
                    actual_bases += 1
                if actual_bases == 1:
                    lastbase = len(refseq) - 1 - num
                elif actual_bases <= 3:
                    last3pos.add(len(refseq) - 1 - num)
            # read the alignment lines
            for line in f:
                if line.rstrip() == "":
                    break
                # only get the top hit for each contig
                if line.split()[0] == lastseq:
                    continue
                lastseq = line.split()[0]
                seq = line[refstart:refend]
                thename = subject_names[line.split()[0]]
                sampledict[thename] = seq
                if seq in qdict:
                    qdict[seq] += 1
                    continue
                insert, deletion, substitution = 0, 0, 0
                last3, last1 = False, False
                for num, (i, j) in enumerate(zip(refseq, seq)):
                    if i == '-' and j == '-':
                        pass
                    elif i == "-" or i == " ":
                        insert += 1
                    elif j == '-' or j == " ":
                        deletion += 1
                    elif j != '.':
                        substitution += 1
                    if j != '.' and num == lastbase:
                        last1 = True
                    elif j != '.' and num in last3pos:
                        last3 = True
                if deletion >= 3:
                    continue
                qdict[seq] = 1
                mutdict[seq] = (substitution, insert, deletion, last3, last1)

        seqlist = list(qdict)
        seqlist.sort(key=lambda x: qdict[x], reverse=True)
        for i in seqlist:
            j = qdict[i]
            substitution, insert, deletion, last3, last1 = mutdict[i]
            if last1:
                alert = "HIGH"
            elif insert + deletion > 0:
                alert = "HIGH"
            elif last3:
                alert = "MEDIUM"
            elif substitution > 2:
                alert = "MEDIUM"
            else:
                alert = "LOW"
            for qq in sampledict:
                if sampledict[qq] == i:
                    outlist.append(list(map(str,
                                        [qq[0], qq[1], primer, primer_seq, i, j, alert, substitution, insert, deletion,
                                         last3, last1])))

    stats_dict = defaultdict(lambda: defaultdict(lambda: "MISSING"))
    single = defaultdict(lambda: set())
    double = defaultdict(lambda: set())
    with open(os.path.join(working_dir, prefix +"_primer_cheq.tsv"), 'w') as report:
        report.write("Reference\tContig\tPrimer\tPrimer_Seq\tAlignment\tCount\tAlert\tSubstitution\tInsertion\tDeletion\tLast_3\tLast_1\n")
        for i in outlist:
            ref = i[0]
            primer = i[2]
            alert = i[6]
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




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--primers", help="FASTA file of primers to check.", required=True, metavar="primer.fasta")
    parser.add_argument("-s", "--prefix", help="Prefix for output files in workign directory.", required=True, metavar="sample_name")
    parser.add_argument("-w", '--working_directory', help="Path for intermediate and output files.", required=True, metavar="working_dir")
    parser.add_argument("-v", "--ncbi_virus", action="append", help="Download a database of viral references from NCBI using a taxid.", metavar="10244")
    parser.add_argument("-b", "--ncbi_bacteria", action="append", help="Download a database of microbial references from NCBI using a taxid", metavar="590")
    parser.add_argument("-d", "--directory_db", action="append", help="A directory of FASTA files, one for each reference.")
    parser.add_argument("-f", "--fasta_db", action="append", help="A single fasta file, each entry will be treated as its own reference.")
    parser.add_argument("-g", "--glob_db", action="append", help="A glob used to identify FASTA files, each FASTA will be treated as its own reference.")




    args = parser.parse_args()


    primer_dict = get_primer_sequences(args.primers)

    if not os.path.exists(args.working_directory):
        os.makedirs(args.working_directory)
    elif os.path.exists(args.working_directory) and not os.path.isdir(args.working_directory):
        sys.stderr.write("The working directory you provided is not a directory.\n")
        sys.exit(0)

    if args.ncbi_virus is None and args.ncbi_bacteria is None and args.directory_db is None and args.fasta_db is None and args.glob_db is None:
        sys.stderr.write("You must provide a database of references to check against.\n")
        sys.exit(0)

    primer_file = os.path.join(args.working_directory, args.prefix + "_db.fasta")
    open(primer_file, 'w').close()

    if args.ncbi_virus:
        for i in args.ncbi_virus:
            fasta_file = download_virus(i, args.working_directory, args.prefix)
            get_db_fasta(fasta_file, args.working_directory, args.prefix)
    if args.ncbi_bacteria:
        for i in args.ncbi_bacteria:
            fasta_files = download_bac(i, args.working_directory, args.prefix)
            get_db_fastas(fasta_files, args.working_directory, args.prefix)
    if args.directory_db:
        for i in args.directory_db:
            fasta_files = get_db_folder(i)
            get_db_fastas(fasta_files, args.working_directory, args.prefix)
    if args.fasta_db:
        for i in args.fasta_db:
            get_db_fasta(i, args.working_directory, args.prefix)
    if args.glob_db:
        for i in args.glob_db:
            fasta_files = get_db_glob(i)
            get_db_fastas(fasta_files, args.working_directory, args.prefix)

    blast_primers(primer_file, primer_dict, args.working_directory, args.prefix)