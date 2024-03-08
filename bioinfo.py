# Author: Emily Bartlett ebart@uoregon.edu

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "0.6"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning
#as of 20240304
DNA_bases = "ATCG"
RNA_bases = "AUCG"

def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    return(ord(letter)) - 33

"""Check that convert_phred returns the correct value for several different inputs"""
assert convert_phred("I") == 40, "wrong phred score for 'I'"
assert convert_phred("C") == 34, "wrong phred score for 'C'"
assert convert_phred("2") == 17, "wrong phred score for '2'"
assert convert_phred("@") == 31, "wrong phred score for '@'"
assert convert_phred("$") == 3, "wrong phred score for '$'"


def qual_score(phred_score: str) -> float:
    '''Calculates the average quality score of the whole phred string'''
    total_score = 0
    count = 0
    for letter in phred_score:
        pscores = convert_phred(letter)
        total_score += pscores
        count += 1
    averageQ = total_score/count
    return averageQ

"""Check that qual_score returns the correct value for several different inputs"""
assert qual_score("EEE") == 36
assert qual_score("#I") == 21
assert qual_score("EJ") == 38.5

def calc_median(lst):
    """Calculates and returns the median of a sorted one-dimensional list."""
    length = len(lst)
    if length % 2 == 0:
        # If the list has an even length, calculate the average of the middle two elements
        mid1 = length // 2 - 1
        mid2 = length // 2
        #print(lst, mid1, mid2)
        median = (lst[mid1] + lst[mid2]) / 2
    else:
        # If the list has an odd length, the median is the middle element
        mid = length // 2
        median = lst[mid]
    return median

    # tests for calc_median
    assert calc_median([10, 20, 30]) == 20, "function not working :C"
    assert calc_median([2, 4, 6, 8]) == 5, "function not working :C"
    assert calc_median([798, 58, 12589]) == 58, "function not working :C"

def init_list(lst: list, value: float=0.0) -> list:
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0.'''
    i = 0
    for x in range(101):
        lst.append(value)
    return lst

def populate_list(file: str) -> tuple[list, int]:
    '''Creates and empty list using init_list(), opens Fastq file, loops through every individual sequence record and converts the phred score to a number, then adds the qual scores to an cumulative sum in the list, keeps track of total number of lines in the file, and returns the list and number of lines as a tuple.'''
    num_lines = 0 #Initialize line counter
    my_list: list = []
    my_list = init_list(my_list)
    with open(file, "r") as fastq: #open fastq file
        for line in fastq: #iterate over each line in the file
            num_lines += 1 #increment the line counter
            if num_lines%4==0: #use only the 4th line of the record that contains quality scores
                phred_scores = line.strip() #strip whitespaces,new lines, etc.
                for i, phred in enumerate(phred_scores): #iterate over each quality score character
                    my_list[i] += bioinfo.convert_phred(phred) #convert quality scores to phred scores and add the quality scores to my_list

    return my_list, num_lines # return list of quality scores (my_list) and line count (num_lines)



def validate_base_seq(seq, RNAflag=False):
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return len(seq) == seq.count("A") + seq.count("U" if RNAflag else "T") + seq.count("G") + seq.count("C")

def gc_content(str):
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    seq = seq.upper()
    return (seq.count('G') + seq.count('C'))/ len(seq)

def kmer_count(filename):
    '''Returns dictionary with kmers as keys and numer of times the corresponding k-mer occurs as value'''
    kmer_counts = {}
    k = 10 #SET
    read_length = 10 #SET
    kmers_per_record = read_length - k + 1
    with open(filename,"r") as fh:
        for line_num, line in enumerate(fh):
            line = line.strip()
            if line_num % 4 == 1:  # Process sequence lines only
                for i in range(kmers_per_record):
                    kmer = line[i:i+k] #kut out kmers
                    if kmer not in kmer_counts:#increment kmer count by one if kmer already in dict, if not start the counter at 1
                        kmer_counts[kmer] = 1
                    else:
                        kmer_counts[kmer] += 1



def oneline_fasta(filename):
    '''This function takes a fasta file, goes line by line, find the sequence line and stips them of their new line character to buid a fasta that only has header and one line of sequence '''
    with open(filename, "r") as fh1, open(f'tempfile.fa', 'w') as fh2:
        fh2.write(fh1.readline()) #creates and opens a tempoary file, open and read through input file line by line
        for line in fh1:
            if line.startswith('>'): #find lines that start with > (header)
                fh2.write(f'\n{line}') #write these lines to fh2 beginning with a new line character to separate reads
            else:
                line = line.strip() #if the line doesnt begin with > strip the newcharacter off the end and write to fh2
                fh2.write(line)
    return 'tempfile.fa'


#find the longest gene in your one line fasta
def find_longest_gene_length(filename):
    longest_gene_length = 0
    with open(filename, 'r') as fh:
        current_gene_length = 0
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                # If this is a header line, update the longest_gene_length if needed
                longest_gene_length = max(longest_gene_length, current_gene_length)
                current_gene_length = 0  # Reset gene length for the next sequence
            else:
                # If this is a sequence line, add its length to current_gene_length
                current_gene_length += len(line)
        # Check the length of the last gene
        longest_gene_length = max(longest_gene_length, current_gene_length)
    return longest_gene_length


#define genetic bases
dnabases = set('ATGCNatcgn')
rnabases = set('AUGCNaucgn')
iupacbases = {
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': '[TU]',
    'U': '[TU]',
    'R': '[AG]',
    'Y': '[CTU]',
    'S': '[GC]',
    'W': '[ATU]',
    'K': '[GTU]',
    'M': '[AC]',
    'B': '[CGTU]',
    'D': '[AGTU]',
    'H': '[ACTU]',
    'V': '[ACG]',
    'N': '[ACGTU]'
}


if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")