#!/usr/bin/env python

import cairo
import argparse
import re
import os
import itertools
from bioinfo import oneline_fasta
from bioinfo import find_longest_gene_length

#assumptions
#-exons in sequence are uppercase and introns are lowercase
#-motfs are set up as one motif per line

#   ./motifmark_oop.py -f /Users/emilybartlett/pycairo/motif-mark/Figure_1.fasta -m /Users/emilybartlett/pycairo/motif-mark/Fig_1_motifs.txt

#-----------------------------------------------------------------------------------------------------------------------------------
#                             Implement Argparse
#-----------------------------------------------------------------------------------------------------------------------------------

def get_args():
    parser = argparse.ArgumentParser(description='Generate an image using pycairo showing a genetic sequence the location of introns exons and motifs.')
    parser.add_argument("-f", "--fasta", help='filename of input fasta file', required=True, default='test.fa')
    parser.add_argument("-m", "--motifs", help='filename of input file containing one motif per line', required=True, default='test_motif.txt')
    parser.add_argument("-o", "--output", help="name of output file", required=False)
    return parser.parse_args()

args = get_args()

motifs_file = args.motifs
output_file = args.output
fasta_file = args.fasta
name = fasta_file.split(".")[0]
motif_f = motifs_file.split(".")[0]
svg_filename = name + '.svg'
png_filename = name + '.png'
image_name_fasta = name.split("/")[5]
image_name_motif = motif_f.split("/")[5]

#-----------------------------------------------------------------------------------------------------------------------------------
#                             Set Colors and iupac bases
#-----------------------------------------------------------------------------------------------------------------------------------

#Set up a list of colors
colors = [
    (1.0, 0.5, 0.5),#light red
    (0.5, 1.0, 0.5),#ight green
    (0.0, 1.0, 1.0),#cyan
    (0.7, 0.5, 0.9),#purple
]

#set up a dictionary with iupac bases
iupacbases = {
    'A': ['A'],
    'C': ['C'],
    'G': ['G'],
    'T': ['T', 'U'],
    'U': ['T'],
    'R': ['A', 'G'],
    'Y': ['C', 'T', 'U'],
    'S': ['G', 'C'],
    'W': ['A', 'T', 'U'],
    'K': ['G', 'T', 'U'],
    'M': ['A', 'C'],
    'B': ['C', 'G', 'T', 'U'],
    'D': ['A', 'G', 'T', 'U'],
    'H': ['A', 'C', 'T', 'U'],
    'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T', 'U']
}
#-----------------------------------------------------------------------------------------------------------------------------------
#                             Set Up Functions
#-----------------------------------------------------------------------------------------------------------------------------------

def create_motif_list(motifs_file):
    """Creates list from motif text file."""
    motif_list = []#initiate list to hold motifs
    with open(motifs_file, "r") as fh:
        for line in fh:
            motif = line.strip()
            if motif:# Skip empty lines
                motif_list.append(motif)
    return motif_list
motif_list = create_motif_list(motifs_file)

def expand_motif_list(motif, iupacbases):
    """Creates a list containing the original motif and all its iupac options."""
    expanded_motif_list = [motif]#initialize a list with the original motif
    #create a list: of possible bases 
    possible_bases = [iupacbases[base] if base in iupacbases else [base] for base in motif]
    for combination in itertools.product(*possible_bases):
        new_motif = ''.join(combination)
        # Replace 'U' with 'T'
        new_motif = new_motif.replace('U', 'T')
        expanded_motif_list.append(new_motif)
    return expanded_motif_list

def draw_legend(MOTIFCONTEXT, motif_list, longest_gene_length):
    keylength = 45#set up a variable to position each motif color box, adjusted as motifs added to legend
    for index, motif in enumerate(motif_list):
        MOTIFCONTEXT.set_source_rgba(0, 0, 0, 1)#black black black
        MOTIFCONTEXT.select_font_face("Purisa", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
        MOTIFCONTEXT.set_font_size(10)
        MOTIFCONTEXT.move_to(longest_gene_length + 70, keylength + 7)
        MOTIFCONTEXT.show_text(motif)
        MOTIFCONTEXT.set_source_rgba(*colors[index % len(colors)])#assign color based on index
        #draw a square, starting on the y axis based on the longest gene lenngth, and on the x axis based on how many items in the legend
        MOTIFCONTEXT.rectangle(longest_gene_length + 35, keylength, 10, 10)
        MOTIFCONTEXT.fill()

        keylength += 15

    return keylength

# Using one-line fasta function to create "tempfile.fa", function from bioinfo.py
oneline_fasta(args.fasta)
#find the longest gene length and store it as a variable, this will be used for y axis orientation of drawings
#function from bioinfo.py
longest_gene_length = find_longest_gene_length('tempfile.fa')

#-----------------------------------------------------------------------------------------------------------------------------------
#                             Set Up pycairo
#-----------------------------------------------------------------------------------------------------------------------------------

#Open SVG surface
#Set the width to accomodate the longest gene so drawings dont look weird
surface = cairo.SVGSurface(svg_filename, longest_gene_length + 200, 500)
MOTIFCONTEXT = cairo.Context(surface)
MOTIFCONTEXT.set_source_rgb(1, 1, 1)
MOTIFCONTEXT.paint()

#-----------------------------------------------------------------------------------------------------------------------------------
#                             Set Up Classes
#-----------------------------------------------------------------------------------------------------------------------------------

class Gene:
    def __init__(self, length, xaxis_loc):
        self.length = length
        self.xaxis_loc = xaxis_loc

    def draw(self, header, gene_stop, xaxis_loc):
        '''Draw a gene that is represented by a black line.'''
        MOTIFCONTEXT.set_source_rgba(0, 0, 0, 1)# set color to black
        MOTIFCONTEXT.move_to(30, xaxis_loc - 15)# moves cursor to x=25, leaving a margin of 15 from the top of the gene
        MOTIFCONTEXT.show_text(header)
        MOTIFCONTEXT.set_line_width(2)
        MOTIFCONTEXT.move_to(30, xaxis_loc)
        MOTIFCONTEXT.line_to(30 + gene_stop, xaxis_loc)#draws the line, corresponding to the length of that gene
        MOTIFCONTEXT.stroke()

class Exon:
    def __init__(self, start, stop, xaxis_loc):
        self.start = start
        self.stop = stop
        self.xaxis_loc = xaxis_loc

    def draw(self, start, stop, xaxis_loc):
        '''Draw an exon that is represented by a black box.'''
        MOTIFCONTEXT.set_source_rgba(0, 0, 0, 1)#black color
        MOTIFCONTEXT.set_line_width(9)
        MOTIFCONTEXT.move_to(30 + start, xaxis_loc)#place cursor at 25 pixels plus beginning location of exon, and the xaxis_loc
        MOTIFCONTEXT.line_to(30 + start + stop, xaxis_loc)#draw the box to 25 pixels plus end location of exon, and the xaxis_loc
        MOTIFCONTEXT.stroke()

class Motif:
    def __init__(self, start, stop, xaxis_loc, color):
        self.start = start
        self.stop = stop
        self.xaxis_loc = xaxis_loc
        self.color = color  # Assign color attribute during initialization

    def draw(self, start, stop, xaxis_loc):
        '''Draw a motif that is represented by a colored box.'''
        MOTIFCONTEXT.set_source_rgba(*self.color)
        MOTIFCONTEXT.set_line_width(10)
        MOTIFCONTEXT.move_to(30 + start, xaxis_loc)#place cursor at 25 pixels plus beginning location of motif, and the xaxis_loc
        MOTIFCONTEXT.line_to(30 + start + stop, xaxis_loc)#draw the box to 25 pixels plus end location of motif, and the xaxis_loc
        MOTIFCONTEXT.stroke()

#-----------------------------------------------------------------------------------------------------------------------------------
#                             Going through the fasta file, creating class instances
#-----------------------------------------------------------------------------------------------------------------------------------

with open('tempfile.fa', 'r') as fh1:
    xaxis_loc = 55  # set variable to set where on the y axis we are going to draw our objects
    for line in fh1:
        if line.startswith('>'):
            header = line.strip()[1:]  # If the line is a header line, strip /n and assign to variable 'header'
        else:
            seq = line.strip()  # assign sequence variable
            gene_stop = len(seq)  # assign the stopping point of the gene by calculating its total length from the isolated sequence
            xaxis_loc += 50  # add a value of 50 for each gene sequence to space drawings
            # Find exons (any amount of capital letters) and return a list of all matches, selecting the first item in the list
            exon_seq = re.findall("[A-Z]+", seq)[0]
            exon_stop = len(exon_seq)
            # Find all occurrences of a lowercase sequence followed by an uppercase letter, return a list and select the first item in the list
            exon_start = len(re.findall("([a-z]+)[A-Z]", seq)[0])

            # Create an instance of the Gene class
            gene_instance = Gene(gene_stop, xaxis_loc)
            # Call the draw method of the Gene instance
            gene_instance.draw(header, gene_stop, xaxis_loc)

            # Create an instance of the Exon class
            exon_instance = Exon(exon_start, exon_stop, xaxis_loc)
            # Call the draw method of the Exon instance
            exon_instance.draw(exon_start, exon_stop, xaxis_loc)

            seq_u = seq.upper()
            #print("Sequence:", seq_u)
            print(f'debug 747 -------------')
            for index, motif in enumerate(motif_list):
                motif = motif.upper()# Ensure motif is in uppercase
                motif_stop = len(motif)
                extended_motif_list = expand_motif_list(motif, iupacbases)
                print("Motif:", motif)

                for extended_motif in extended_motif_list:

                    motif_positions = [m.start() for m in re.finditer(extended_motif, seq_u)]  # Find motif positions for the current motif
                    #print(motif_positions)
                    for motif_start in motif_positions:
                        motif_instance = Motif(motif_start, motif_stop, xaxis_loc, colors[index % len(colors)])
                        motif_instance.draw(motif_start, motif_stop, xaxis_loc)

#-----------------------------------------------------------------------------------------------------------------------------------
#                             Set Up Legend
#-----------------------------------------------------------------------------------------------------------------------------------
legend_xaxis_loc = draw_legend(MOTIFCONTEXT, motif_list, longest_gene_length)


# Set title based on the input file name
MOTIFCONTEXT.set_source_rgba(0, 0, 0, 1)
MOTIFCONTEXT.select_font_face("Purisa", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_BOLD)
MOTIFCONTEXT.set_font_size(20)
MOTIFCONTEXT.move_to(25, 25)
MOTIFCONTEXT.show_text("Visualization of Motifs")
MOTIFCONTEXT.set_font_size(12)
MOTIFCONTEXT.move_to(25, 42)
MOTIFCONTEXT.show_text("Sequences from:  " + image_name_fasta)
MOTIFCONTEXT.set_font_size(12)
MOTIFCONTEXT.move_to(25, 57)
MOTIFCONTEXT.show_text("Motifs from:  " + image_name_motif)

# Write out to a PNG for preview
surface.write_to_png(png_filename)

surface.finish()

print("Image Created Successfully and Saved")
