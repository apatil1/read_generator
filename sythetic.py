
import random
import itertools
import os
from cStringIO import StringIO #? where do you use it?
import argparse

def parse_stdin_args(): # this function is not used
    """This function reads the input
    arguments passed through standard in"""
    #Initialize parser
    parser = argparse.ArgumentParser()
    #Required command line arguments
    parser.add_argument( "-s", help="Size of synthetic genome", required=True)
    parser.add_argument("-p", help="Path to directory for output files", required=True)
    
    #Read the command line arguments
    args = parser.parse_args()
    
                        
    
    return args



def parse_fasta(f, trim_desc=False): # BioPython can do it for free
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


def write_fasta(f, seqs):
    for desc, seq in seqs:
        f.write(">%s\n%s\n" % (desc, seq))


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
        # you used generators but in the next line convert it to dictionary  - that get rid of generators
        return dict(seqs) # use BioPython SeqIO.index


# please remove unused function
def parse_greengenes_accessions(f):
    for line in f:
        if line.startswith("#"):
            continue
        line = line.strip()
        yield line.split("\t")




def get_chromosome_length(genome):
    # doc string is missing here
    
    #Function to get chromosome list and the size and returnds a dictionary of chromosome as keys and size as values
    chr_list = {}
    
    for key in genome:
        chr_list[key] = len(genome[key])

    return chr_list


def get_random_sequence(genome):
    
    #Function to get random sequences from the genomes by selecting chromosome, sequence length and start position at random and returning the extracted random sequence from get_fragment function
    chr_list = get_chromosome_length(genome)
    
    random_seq = {}
    # soorry at being pedantic but there is a space after comma
    chr = random.sample(chr_list.keys(),1)  #select chromosome
    #dont hard-code the values here
    slen = random.randint(300,1000) #select sequence length
    # the distribution here is uniform 
    spos = random.randint(1,chr_list[chr[0]] - slen)    #select start position
   
    seq = get_fragment(genome, chr[0], slen, spos)
    if seq.count("N") > 0.1 * slen:
        # you can potentially end up in an infinity loop here please bail after several attmpts
        seq = get_random_sequence(genome)

    return seq
   
def get_fragment(genome, chr, slen, spos): # chr ? 
    # please protect from overflow
    # what is chromosome when you are working with bacteria or virus?
    
    #Function to extract sequence from genome using chromosome, length of sequence and start position
    return genome[chr][spos:spos+slen]


def make_paired_end_reads(sequence):
    
    # here we have an assumtion that all reads are 250...
    
    #Function to make paired end reads by taking 250 bases from start and nd of the selected sequence and taking a reverse complement of one of the sequence at random
    R1 = sequence[0:250]  # 250 should not be hard-coded
    R2 = sequence[len(sequence) - 250:len(sequence)]

    if random.randint(0,1): # you can be more explicit here: random.choice([True, False])
        R1 = make_reverse_complement(R1)
    else:
        R2 = make_reverse_complement(R2)

    return [R1, R2]



def make_reverse_complement(seq):  # BioPyhton has API to write FASTQ
    
    #Function to make reverse complement of a sequence
    comp ={"A":"T","T":"A","G":"C","C":"G","a":"t","t":"a","g":"c","c":"g", "N":"N"}

    rev = ""
    for i in range(0,len(seq)):
        rev += comp[seq[i]]

    return rev[::-1]


def make_fastq(pair, filename):  # BioPyhton has API to write FASTQ
    #Function to write sequence in FASTQ files for both reads

    # note here code duplication for R1 and R2 please create function
    fname = filename + "-R1.fastq"
    r1 = open(fname,"w")
    r1.write("@" + filename + "\n")
    r1.write(pair[0])
    r1.write("\n+\n")
    r1.write("E" * len(pair[0]))
    r1.close()

    fname = filename + "-R2.fastq"
    r2 = open(fname,"w")
    r2.write("@" + filename + "\n")
    r2.write(pair[1])
    r2.write("\n+\n")
    r2.write("E" * len(pair[1]))
    r2.close()


def permutate_genome_percent(human, phix, bacteria):
    
    #Function to permutate to get all combinations of the compositions of fake genomes
    return list(itertools.product(human, phix, bacteria))


def make_synthetic_genome(human, phix, bacteria, size, dir):
    
    #Function to make synthetic genomes from human, phix, bacteria and viruses.
    
    # generate human reads    # this comment is redundant
    get_human_reads(human, size, dir)  # please rename dir to something else
    
    # generate phix reads
    get_phix_reads(phix, size, dir)
    
    # generate bacteria reads
    get_bacteria_reads(bacteria, size, dir)
    
    # generate virus reads
    get_virus_reads(1 - human - phix - bacteria, size, dir)
    # note that this is actually can be different from what you are writing to the make_readme due to round off errors


def get_human_reads(percent, size, dir):
    # please rename size it to read_count , read_number of something similar
    
    #Function to get random human reads
    genome = load_fasta("genomes/hg18.fa") # we don't need to hard-code data
       
    for i in range(0,int(size * percent)):
        seq = get_random_sequence(genome)
       
        pair = make_paired_end_reads(seq)
        make_fastq(pair, dir + "human" + str(i+1))

    
def get_phix_reads(percent, size, dir):

    #Function to get random phix174 reads
    genome = load_fasta("genomes/phix174.fasta")
   
    for i in range(0,int(size * percent)):
        seq = get_random_sequence(genome)
        pair = make_paired_end_reads(seq)
        make_fastq(pair, dir + "phix" + str(i+1))


def get_bacteria_reads(percent, size, dir):
 
    #Function to get random bacteria reads
    bac_select = random.sample(os.listdir("genomes/bacteria/all.fna"), 1) # dont hard-code the path have it as an argument
    gen = random.sample(os.listdir("genomes/bacteria/all.fna/" + bac_select[0]), 1)
    path = "genomes/bacteria/all.fna/" + bac_select[0] + "/" + gen[0]

    genome = load_fasta(path)
   
    for i in range(0,int(size * percent)):
        seq = get_random_sequence(genome)
        pair = make_paired_end_reads(seq)
        make_fastq(pair, dir +  "bacteria" + str(i+1), "bacteria")
    

def get_virus_reads(percent, size, dir):
    
    #Function to get random virus/phage reads
    bac_select = random.sample(os.listdir("genomes/phage/all.fna"), 1)
    gen = random.sample(os.listdir("genomes/phage/all.fna/" + bac_select[0]), 1)
    path = "genomes/phage/all.fna/" + bac_select[0] + "/" + gen[0]

    genome = load_fasta(path)
  
    for i in range(0,int(size * round(percent,2))):
        seq = get_random_sequence(genome)
        pair = make_paired_end_reads(seq)
        make_fastq(pair, dir + "phage" + str(i+1))
    
def make_readme(human, phix, bacteria, size, dir):
    # size variable is not used in the function
   
    #Function creates readme file listing the composition of fake genome
    r = open(dir + "Readme.txt","w") # please use functions from os.path for working with paths
    r.write("Human :" + str(human))
    r.write("\nPhix174 :" + str(phix))
    r.write("\nBacteria :" + str(bacteria))
    r.write("\nVirus/Phage :" + str(round(1 - human - phix - bacteria, 2)))
    r.close()



def main(argv):
    args = parse_stdin_args() # this function is not used
    
    human = [0.5, 0.1, 0.01, 0.001]
    phix = [0.01, 0.001]
    bacteria = [0.4, 0.25, 0.1, 0.5]
    # you need to have a check that total is always 1 or less
    list = permutate_genome_percent(human, phix, bacteria) # please dont use list as aname of variable

    main_dir = "//"
    for i in range(0,len(list)):
    
        # please call it simulated not fake genome :)
        dir = main_dir + "fake_genome" + str(i+1) + "/"
        if not os.path.exists(dir):
            os.makedirs(dir)

        make_synthetic_genome(list[i][0], list[i][1], list[i][2], 100, dir) # what is 100
        make_readme(list[i][0], list[i][1], list[i][2], 100, dir) # what is 100
    

    return

if __name__ == '__main__': main(sys.argv)


 
