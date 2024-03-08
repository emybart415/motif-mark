
# motif-mark
----------------------------------------------------------------------------------------------------------------
A python script that uses object-orientated code to visualize the locations of motifs and exons on a gene.
Bellow is an example: 
----------------------------------------------------------------------------------------------------------------

![Figure 1](https://github.com/emybart415/motif-mark/blob/main/Figure_1.png)


It requires the input of a fasta file and a text tile that contains the motifs of interest. The text file mus be set up as a mone motif per line file.

### Input Files

• fasta.fasta file

• motif.txt file

## Files required to run the script

• 'motif-mark-oop.py' - main script

• 'bioinfo.py' - script that inclueds some functions that are imported and used by the main script

## Required Packages (scripts themselves will handle importation)

motif-mark-oop.py

• re

• cairo

• math

• os

• itertools

• argparse

• bioinfo

bioinfo.py

• re

• matplotlib.pylab as plt

• numpy as np

• argparse

• math


## Usage and Output
Argparse is being used in the main script, to run you enter the following in your terminal once you are in your working directory:

$ ./motif-mark-oop.py -f <fata_file_name_and_location>.fasta -m <motif_file_name_and_location>.txt

The output file will be of both .png and .svg file format. It will be named based on the name of the fasta input file.
The script will output the phrase "Image Created Successfully and Saved" when the script is run successfully. 
