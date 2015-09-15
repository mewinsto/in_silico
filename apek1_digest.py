##IN SILICO DIGEST OF GENOME WITH ApeK1
#EXAMPLE IS SOLENOPSIS GENOME "Sinv_1.0_scaffolds.fa"
#WRITES TO FILE 'frag_lengths_Sinv'

from Bio import SeqIO
import re
import numpy as np

reader = SeqIO.parse("/home/mwinston/ant_genomes/Sinv_1.0_scaffolds.fa", format = 'fasta')

frag_lengths_out = open('frag_lengths_Sinv','w')
#frag_lengths_scaf_out = open('FLS','w')

seqs = []
fragment_lengths = []
frag_lengths_scafs = []

for rec in reader:
    seqs.append(str(rec.seq))

for i in range(len(seqs)):
    frags = re.split('CGACG|CGTCG',seqs[i])
    frag_lengths_scafs.append(frags)
    for x in frags:
        z = len(x)
        fragment_lengths.append(z)

for n in range(len(fragment_lengths)):
    print >>frag_lengths_out, fragment_lengths[n]
