
import random
import itertools
import os
from cStringIO import StringIO
import argparse
import sys

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
    parser.add_argument("-hu", help="Human DNA percentage. Default: [0.5, 0.1, 0.01, 0.001]", 
        required = False, type= int, nargs = "+", default = [0.5, 0.1, 0.01, 0.001])
    parser.add_argument("-b", help="Bacterial DNA percentage. Default: [0.4, 0.25, 0.1, 0.5]", 
        required = False, type= int, nargs = "+", default = [0.4, 0.25, 0.1, 0.5])
    parser.add_argument("-x", help="Phix174 DNA percentage. Default: [0.01, 0.001]", 
        required = False, type= int, nargs = "+", default = [0.01, 0.001])
    
    #Read the command line arguments
    args = parser.parse_args()
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
    '''
    Function to get chromosome list and the size and returns a dictionary
    of chromosome as keys and size as values
    '''
    
    chr_list = {}
    
    for key in genome:
        chr_list[key] = len(genome[key])

    return chr_list


def get_random_sequence(genome):
    '''
    Function to get random sequences from the genomes by selecting chromosome,
    sequence length and start position at random and returning the extracted 
    random sequence from get_fragment function
    '''
    
    chr_list = get_chromosome_length(genome)
    
    random_seq = {}
    # please dont overwrite reserved python words
    chr = random.sample(chr_list.keys(),1)  #select chromosome
    slen = random.randint(300,1000) #select sequence length
    spos = random.randint(1,chr_list[chr[0]] - slen)    #select start position
   
    seq = get_fragment(genome, chr[0], slen, spos)
    if seq.count("N") > 0.1 * slen:
        seq = get_random_sequence(genome)

    return seq
   
def get_fragment(genome, chr, slen, spos):
    '''
    Function to extract sequence from genome using chromosome, length of sequence
    and start position
    '''
    
    return genome[chr][spos:spos+slen]


def make_paired_end_reads(sequence):
    '''
    Function to make paired end reads by taking 250 bases from start and nd of the
    selected sequence and taking a reverse complement of one of the sequence at random
    '''
    read_length = 250 
    R1 = sequence[0:read_length]
    R2 = sequence[len(sequence) - read_length:len(sequence)]

    if random.randint(0,1):
        R1 = make_reverse_complement(R1)
    else:
        R2 = make_reverse_complement(R2)

    return [R1, R2]



def make_reverse_complement(seq):
    '''
    Function to make reverse complement of a sequence
    '''
    
    comp ={"A":"T","T":"A","G":"C","C":"G","a":"t","t":"a","g":"c","c":"g", "N":"N"}

    rev = ""
    for i in range(0,len(seq)):
        rev += comp[seq[i]]

    return rev[::-1]


def make_fastq(pair, filename, id):
    '''
    Function to write sequence in FASTQ files for both reads
    '''
    
    fname = filename + "-R1.fastq"
    r1 = open(fname,"w")
    r1.write("@" + id + "\n")
    r1.write(pair[0])
    r1.write("\n+\n")
    r1.write("E" * len(pair[0]))
    r1.close()

    fname = filename + "-R2.fastq"
    r2 = open(fname,"w")
    r2.write("@" + id + "\n")
    r2.write(pair[1])
    r2.write("\n+\n")
    r2.write("E" * len(pair[1]))
    r2.close()


def concatenate_fastq(path):
    '''
    Function to concatenate all the R1 and R2 FASTQ files.
    '''
    
    r1 = []
    r2 = []
    filenames = os.listdir(path)
    
    for i in filenames:
        if "R1" in i:
            r1.append(i)
        elif "R2" in i:
            r2.append(i)
    
    #concatinate R1
    with open(path +  "fake_genome-R1.fastq", 'w') as outfile:
        for fname in r1:
            with open(path + fname) as infile:
                outfile.write(infile.read())
                outfile.write("\n")

    infile.close()
    outfile.close()

    #concatinate R2
    with open(path +  "fake_genome-R2.fastq", 'w') as outfile:
        for fname in r2:
            with open(path + fname) as infile:
                outfile.write(infile.read())
                outfile.write("\n")


    infile.close()
    outfile.close()

    for i in filenames:
        if i == "Readme.txt":
            continue
        os.remove(path + i)

def permutate_genome_percent(human, phix, bacteria):
    '''
    Function to permutate to get all combinations of the compositions of fake genomes
    '''
    
    return list(itertools.product(human, phix, bacteria))


def make_synthetic_genome(human, phix, bacteria, size, dir):
    '''
    Function to make synthetic genomes from human, phix, bacteria and viruses.
    '''
    
    # generate human reads
    get_human_reads(human, size, dir)
    
    # generate phix reads
    get_phix_reads(phix, size, dir)
    
    # generate bacteria reads
    get_bacteria_reads(bacteria, size, dir)
    
    # generate virus reads
    get_virus_reads(1 - human - phix - bacteria, size, dir)


def get_human_reads(percent, size, dir):
    '''
    Function to get random human reads
    '''
    
    genome = load_fasta("/home/ashwini/ash/testing_tools/genomes/hg18.fa")
       
    for i in range(0,int(size * percent)):
        seq = get_random_sequence(genome)
       
        pair = make_paired_end_reads(seq)
        make_fastq(pair, dir + "human" + str(i+1), "human" + str(i+1))

    
def get_phix_reads(percent, size, dir):
    '''
    Function to get random phix174 reads
    '''
    
    # should not be hard-coded
    genome = load_fasta("/home/ashwini/ash/testing_tools/genomes/phix174.fasta")
   
    for i in range(0,int(size * percent)):
        seq = get_random_sequence(genome)
        pair = make_paired_end_reads(seq)
        make_fastq(pair, dir + "phix" + str(i+1), "phix" + str(i+1))


def get_bacteria_reads(percent, size, dir):
    '''
    Function to get random bacteria reads
    '''
    
    bac_select = random.sample(os.listdir("/home/ashwini/ash/testing_tools/genomes/bacteria/all.fna"), 1)
    gen = random.sample(os.listdir("/home/ashwini/ash/testing_tools/genomes/bacteria/all.fna/" + bac_select[0]), 1)
    path = "/home/ashwini/ash/testing_tools/genomes/bacteria/all.fna/" + bac_select[0] + "/" + gen[0]

    genome = load_fasta(path)
   
    for i in range(0,int(size * percent)):
        seq = get_random_sequence(genome)
        pair = make_paired_end_reads(seq)
        make_fastq(pair, dir +  "bacteria" + str(i+1), "bacteria" + str(i+1))
    

def get_virus_reads(percent, size, dir):
    '''
    Function to get random virus/phage reads
    '''
    
    bac_select = random.sample(os.listdir("/home/ashwini/ash/testing_tools/genomes/phage/all.fna"), 1)
    gen = random.sample(os.listdir("/home/ashwini/ash/testing_tools/genomes/phage/all.fna/" + bac_select[0]), 1)
    path = "/home/ashwini/ash/testing_tools/genomes/phage/all.fna/" + bac_select[0] + "/" + gen[0]

    genome = load_fasta(path)
  
    for i in range(0,int(size * round(percent,2))):
        seq = get_random_sequence(genome)
        pair = make_paired_end_reads(seq)
        make_fastq(pair, dir + "phage" + str(i+1), "phage/virus" + str(i+1))
    
def make_readme(human, phix, bacteria, dir):
    '''
    Function creates readme file listing the composition of fake genome
    '''
    
    r = open(dir + "Readme.txt","w")
    r.write("Human :" + str(human))
    r.write("\nPhix174 :" + str(phix))
    r.write("\nBacteria :" + str(bacteria))
    r.write("\nVirus/Phage :" + str(round(1 - human - phix - bacteria, 2)))
    r.close()



def main(argv):
    args = parse_stdin_args()
    
    list = permutate_genome_percent(args.hu, args.x, args.b)
    
    for i in range(0,len(list)):
        genome_dir = args.p +  "/"  + "fake_genome" + str(i+1) + "/"
        if not os.path.exists(genome_dir):
            os.makedirs(genome_dir)
        
        (human, phiX, bacteria) = list[i]
        make_synthetic_genome(human, phiX, bacteria, args.s, genome_dir)  #creates FASTQ files from randomly selected sequences from different organisms
        make_readme(human, phiX, bacteria, genome_dir)    #creates readme file denoting percentage of DNA for each organism in Synthetic genome
        concatenate_fastq(genome_dir)  #concatenates all FASTQ files to make R1 and R2 FASTQ files and removes other files
        print "Finished making synthetic genome " + str(i+1)

if __name__ == '__main__': 
    main(sys.argv)


 
