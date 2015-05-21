
import random
import itertools
import os
from cStringIO import StringIO
import argparse
import sys
import shutil

path_human = "/home/ashwini/ash/testing_tools/genomes/hg18.fa"
path_phix = "/home/ashwini/ash/testing_tools/genomes/phix174.fasta"
path_bacteria = "/home/ashwini/ash/testing_tools/genomes/bacteria/all.fna"
path_virus = "/home/ashwini/ash/testing_tools/genomes/phage/all.fna"


def parse_stdin_args():
    """This function reads the input
    arguments passed through standard in"""
    
    #Initialize parser
    parser = argparse.ArgumentParser(description = "Creates synthetic paired end DNA from given fraction of Human, Bacteria, Phix174 and Virus/Phage. \
                                     Note: List of Human, Bacteria, Phix174 and Virus/Phage fractions should be space separated. Example: -hu 0.5 0.2 \
                                     -b 0.3 0.01 -x 0.1 0.01. The Virus/Phage fraction is taken 1 - Human - Bacteria - Phix174.")
    
    #Required command line arguments
    parser.add_argument("-s", help="Number of total reads", required = True, type = int)
    parser.add_argument("-p", help="Path to directory for output files", required = True)
    parser.add_argument("-hu", help="Human DNA percentage. Default: [0.5, 0.1, 0.01, 0.001]", required = False, type= float, nargs = "+", default = [0.5, 0.1, 0.01, 0.001])
    parser.add_argument("-b", help="Bacterial DNA percentage. Default: [0.4, 0.25, 0.1, 0.5]", required = False, type= float, nargs = "+", default = [0.4, 0.25, 0.1, 0.5])
    parser.add_argument("-x", help="Phix174 DNA percentage. Default: [0.01, 0.001]", required = False, type= float, nargs = "+", default = [0.01, 0.001])
    parser.add_argument("-isfastq", help="Specify output format 1=FASTQ 0=FASTA (Default: 1)", required = False, type= int, default = 1)
    parser.add_argument("-n", help="Read length. Default: 250", required = False, type= int, default = 250)
    parser.add_argument("-err", help="Error rate. Default: 0", required = False, type= int, default = 0)

    
    
    #Read the command line arguments
    args = parser.parse_args()
    
       
    #check parameters
    if args.s < 100:
        print "Number of total reads must be greater than 100."
        sys.exit(0)
    
    if args.n < 35:
        print "Read length should be greater than 35"
        sys.exit(0)

    return args

def load_fasta(filepath, trim_desc=True):
    """Load all sequences from a FASTA file
        
        Parameters
        ----------
        fasta_fp : Input filepath, FASTA format.
        
        Returns
        -------
        A dictionary mapping sequence identifiers to sequences.
        """
    
    with open(filepath) as f:
        seqs = parse_fasta(f, trim_desc=trim_desc)
        return dict(seqs)


def parse_fasta(f, trim_desc=False):
    """Parse a FASTA format file.

    Parameters
    ----------
    f : File object or iterator returning lines in FASTA format.

    Returns
    -------
    An iterator of tuples containing two strings
        First string is the sequence description, second is the
        sequence.

    Notes
    -----
    This function removes whitespace in the sequence and translates
    "U" to "T", in order to accommodate FASTA files downloaded from
    SILVA and the Living Tree Project.
    """
    
    f = iter(f)
    desc = next(f).strip()[1:]
    if trim_desc:
        desc = desc.split()[0]
    seq = StringIO()
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            yield desc, seq.getvalue()
            desc = line[1:]
            if trim_desc:
                desc = desc.split()[0]
            seq = StringIO()
        else:
            seq.write(line.replace(" ", "").replace("U", "T"))
    yield desc, seq.getvalue()


def get_chromosome_length(genome):
    """
    Function to get chromosome list and the size and returns a dictionary
    of chromosome as keys and size as values
    """
    
    chr_list = {}
    
    for key in genome:
        chr_list[key] = len(genome[key])

    return chr_list


def get_random_sequence(genome):
    """
    Function to get random sequences from the genomes by selecting chromosome,
    sequence length and start position at random and returning the extracted 
    random sequence from get_fragment function
    """
    
    chr_list = get_chromosome_length(genome)
    
    random_seq = {}
    chr = random.sample(chr_list.keys(),1)  #select chromosome
    slen = random.randint(300,1000) #select sequence length
    if chr_list[chr[0]] - slen > 0:
        spos = random.randint(1,chr_list[chr[0]] - slen)    #select start position
   
        seq = get_fragment(genome, chr[0], slen, spos)
        if seq.count("N") > 0.1 * slen:
            seq = get_random_sequence(genome)
    else:
        seq = get_random_sequence(genome)
    
    return seq
   
def get_fragment(genome, chr, slen, spos):
    """
    Function to extract sequence from genome using chromosome, length of sequence
    and start position
    """
    
    return genome[chr][spos:spos+slen]


def make_paired_end_reads(sequence):
    """
    Function to make paired end reads by taking 250 bases from start and nd of the
    selected sequence and taking a reverse complement of one of the sequence at random
    """
    
    R1 = sequence[0:n]
    R2 = sequence[len(sequence) - n:len(sequence)]


    #one reads are reverse complement, so make reverse complement of R2
    R2 = make_reverse_complement(R2)

    return [R1, R2]

'''
def introduce_errors(per_error, pair):

    """
    This function introduces insertions in the sequence
    """
    en = int(round((per_error * n)/100, 2))
    error_bp = ''.join(random.choice(["a","t","g","c"]) for i in range(0,en))
    
    s = random.randint(0,n - en)

    pair[0] = ''.join(pair[0][0:s] + error_bp + pair[0][s + en:n])
    pair[1] = ''.join(pair[1][0:s] + error_bp + pair[1][s + en:n])
    
    
    return pair

'''
def introduce_errors(per_error, pair):
    
    """
    This function introduces errors in the sequence
    """
    en = int(round((per_error * n)/100, 2))
    
    
    for i in range(0,en):
        s = random.randint(0,n)
    
        pair[0] = ''.join(pair[0][0:s] + random.choice(["a","t","g","c"])  + pair[0][s+1:n])
        pair[1] = ''.join(pair[1][0:s] + random.choice(["a","t","g","c"])  + pair[1][s+1:n])
    
    
    return pair

def make_reverse_complement(seq):
    """
    Function to make reverse complement of a sequence
    """
    
    comp ={"A":"T", "T":"A", "G":"C", "C":"G", "a":"t", "t":"a", "g":"c", "c":"g", "N":"N", \
            "R":"R", "K":"K", "M":"M", "B":"B", "V":"V", "S":"S", "W":"W", "D":"D", "-":"-",\
            "U":"U", "Y":"Y", "H":"H", " ":""}

    rev = ""
    for i in range(0,len(seq)):
        rev += comp[seq[i]]

    return rev[::-1]


def make_fastq(pair, filename, id):
    """
    Function to write sequence in FASTQ files for both reads
    """
    
    fname = filename + "-R1.fastq"
    with open(fname, "w") as r1:
        r1.write("@" + id + "\n")
        r1.write(pair[0])
        r1.write("\n+\n")
        r1.write("E" * len(pair[0]))

    fname = filename + "-R2.fastq"
    with open(fname, "w") as r2:
        r2.write("@" + id + "\n")
        r2.write(pair[1])
        r2.write("\n+\n")
        r2.write("E" * len(pair[1]))
    

def make_fasta(pair, filename, id):
    """
    Function to write sequence in FASTQ files for both reads
    """
    
    fname = filename + "-R1.fasta"
    with open(fname,"w") as r1:
        r1.write(">" + id + "\n")
        r1.write(pair[0])
        r1.write("\n")
    
    fname = filename + "-R2.fasta"
    with open(fname,"w") as r2:
        r2.write(">" + id + "\n")
        r2.write(pair[1])
        r2.write("\n")


def concatenate_fastq(path, isfastq, sample_name):
    """
    Function to concatenate all the R1 and R2 FASTQ/FASTA files.
    """
    
    r1 = []
    r2 = []
    filenames = get_filesnames_in_dir(path)
    
    for i in filenames:
        if "fake_genome" in i:
            continue
        elif "R1" in i:
            r1.append(i)
        elif "R2" in i:
            r2.append(i)
    if isfastq:
        nameR1 = sample_name + "-R1.fastq"
        nameR2 = sample_name + "-R2.fastq"
    else:
        nameR1 = sample_name + "-R1.fasta"
        nameR2 = sample_name + "-R2.fasta"

    #concatinate R1
    with open(path +  nameR1, 'w') as outfile:
        for fname in sorted(r1):
            with open(path + fname) as infile:
                outfile.write(infile.read())
                outfile.write("\n")

    #concatinate R2
    with open(path +  nameR2, 'w') as outfile:
        for fname in sorted(r2):
            with open(path + fname) as infile:
                outfile.write(infile.read())
                outfile.write("\n")

    
    for i in r1 + r2:
        os.remove(path + i)

def write_samples(path, sample_number, isfastq):

    with open(path + "samples.txt", "w") as outfile:
        for i in range(0, sample_number):
            if isfastq:
                outfile.write("fake_genome" + str(i+1) + "\t" + path + "fake_genome" + str(i+1) + "-R1.fastq" +  "\t" + path + "fake_genome" + str(i+1) + "-R2.fastq"  + "\n")
            else:
                outfile.write("fake_genome" + str(i+1) + "\t" + path + "fake_genome" + str(i+1) + "-R1.fasta" +  "\t" + path + "fake_genome" + str(i+1) + "-R2.fasta"  + "\n")

def permutate_genome_percent(human, phix, bacteria):
    """
    Function to permutate to get all combinations of the compositions of fake genomes
    """
    
    per = list(itertools.product(human, phix, bacteria))
    sum_per = [sum(i) for i in zip(*per)]
    
    #check percentage sum < 1
    if all(i > 1 for i in sum_per):
        print "Some combinations of human, phix and bacteria greater than 1"
        sys.exit(0)
    
    return per


def make_synthetic_genome(human, phix, bacteria, size, dir, isfastq):
    """
    Function to make synthetic genomes from human, phix, bacteria and viruses.
    """
    
    # generate human reads
    get_human_reads(human, size, dir, isfastq)
    
    # generate phix reads
    get_phix_reads(phix, size, dir, isfastq)
    
    # generate bacteria reads
    get_bacteria_reads(bacteria, size, dir, isfastq)
    
    # generate virus reads
    get_virus_reads(1 - human - phix - bacteria, size, dir, isfastq)


def get_human_reads(percent, size, dir, isfastq):
    """
    Function to get random human reads
    """
    
    for i in range(0,int(size * percent)):
        seq = get_random_sequence(human_genome)
        pair = make_paired_end_reads(seq)
        
        global errr
  
        if errr:
            pair = introduce_errors(errr, pair)
            #errr = 0
        
        if isfastq:
             make_fastq(pair, dir + "human" + str(i+1), "human" + str(i+1))
        else:
             make_fasta(pair, dir + "human" + str(i+1), "human" + str(i+1))
    
def get_phix_reads(percent, size, dir, isfastq):
    """
    Function to get random phix174 reads
    """
    
    genome = load_fasta(path_phix)
   
    for i in range(0,int(size * percent)):
        seq = get_random_sequence(genome)
        pair = make_paired_end_reads(seq)
        
        if errr:
            pair = introduce_errors(errr, pair)
            #errr = 0
        
        if isfastq:
             make_fastq(pair, dir + "phix" + str(i+1), "phix" + str(i+1))
        else:
             make_fasta(pair, dir + "phix" + str(i+1), "phix" + str(i+1))
    
    
def get_bacteria_reads(percent, size, dir, isfastq):
    """
    Function to get random bacteria reads
    """
    bac_select = random.sample(get_filesnames_in_dir(path_bacteria), 1)
    gen = random.sample(get_filesnames_in_dir(path_bacteria + "/" + bac_select[0]), 1)
    path = path_bacteria + "/" + bac_select[0] + "/" + gen[0]

    genome = load_fasta(path)
   
    for i in range(0,int(size * percent)):
        seq = get_random_sequence(genome)
        pair = make_paired_end_reads(seq)
        
        if errr:
            pair = introduce_errors(errr, pair)
            #errr = 0
            
        if isfastq:
             make_fastq(pair, dir +  "bacteria" + str(i+1), "bacteria" + str(i+1))
        else:
             make_fasta(pair, dir +  "bacteria" + str(i+1), "bacteria" + str(i+1))
       
def get_virus_reads(percent, size, dir, isfastq):
    """
    Function to get random virus/phage reads
    """
    bac_select = random.sample(get_filesnames_in_dir(path_virus), 1)
    gen = random.sample(get_filesnames_in_dir(path_virus + "/" + bac_select[0]), 1)
    path = path_virus + "/" + bac_select[0] + "/" + gen[0]

    genome = load_fasta(path)
  
    for i in range(0,int(size * round(percent,2))):
        seq = get_random_sequence(genome)
        pair = make_paired_end_reads(seq)
        if errr:
            pair = introduce_errors(errr, pair)
            #errr = 0
        
        if isfastq:
            make_fastq(pair, dir + "phage" + str(i+1), "phage/virus" + str(i+1))
        else:
            make_fasta(pair, dir + "phage" + str(i+1), "phage/virus" + str(i+1))
        

def get_filesnames_in_dir(path):
    """
    This function gets file names in directory and returns a list of file names after removing .xxx system files
    """
    
    file_list = os.listdir(path)
    if ".DS_Store" in file_list:
        del file_list[file_list.index(".DS_Store")]
    return file_list

'''
def make_readme(human, phix, bacteria, dir, sample_name):
    """
    Function creates readme file listing the composition of fake genome
    """
    
    with open(dir + "Readme.txt", "w") as r:
        r.write("sample_name\thuman\tphix\tbacteria\tvirus")
        r.write("\n" + str(human) +"\t" + str(phix) +"\t" + str(bacteria) +"\t" + str(round(1 - human - phix - bacteria, 2)))

'''
def main():
    args = parse_stdin_args()
    global n
    n = args.n
    global errr
    errr = args.err
    
    list = permutate_genome_percent(args.hu, args.x, args.b)
    
    global human_genome
    human_genome = load_fasta(path_human)

    for i in range(0,len(list)):
        make_synthetic_genome(list[i][0], list[i][1], list[i][2], args.s, args.p, args.isfastq)  #creates FASTQ files from randomly selected sequences from different organisms
        concatenate_fastq(args.p, args.isfastq, "fake_genome" + str(i+1))  #concatenates all FASTQ files to make R1 and R2 FASTQ files and removes other files
        #make_readme(list[i][0], list[i][1], list[i][2], args.p, "fake_genome" + str(i+1))    #creates readme file denoting percentage of DNA for each organism in Synthetic genome
        print "Finished making synthetic genome " + str(i+1)

    write_samples(args.p, len(list), args.isfastq)
    return

if __name__ == '__main__': main()


 
