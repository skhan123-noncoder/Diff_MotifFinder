#/home/usr/bin/python

from Bio import Align
from Bio import SeqIO
import os
import sys
import argparse

def main():
     
     func1=parse_arguments()
     func2=quality_check()
     func3=run_alignment()
     func4=run_mast()
     #func4=analyze_mast_output()

def parse_arguments():

    parser=argparse.ArgumentParser(prog='Differential Motif finder', description="Provide fasta file of containing a pair of simialr DNA sequences of equal lengths", epilog='none')
    parser.add_argument('-i', '--input', action='store', dest='infile', type=str, required=True, help='The input fasta file to the script.')
    parser.add_argument('-d', '--database', action='store', dest='TFdatabase', type=str, required=True, help="path to PWM file for TFs in MEME format")
    parser.add_argument('-o', '--output', action='store', dest='outfile', type=str, help="name of the output file", default="out")
    args=parser.parse_args()

    infile=str(os.path.basename(args.infile))
    TFdatabase=str(os.path.basename(args.TFdatabase))
    outfile=str(os.path.basename(args.outfile))
    data=[infile, TFdatabase, outfile]
    return(data)

    #os.system('mast -o '+ outfile+ ' -nostatus -minseqs 1 -remcorr -ev 10.0 ' +TFdatabase+ ' '+infile)

def quality_check():
    data=parse_arguments()
    seqID=[]
    fasta=[]

    for record in SeqIO.parse(data[0], "fasta"):
        seqID.append(record.id)
        fasta.append(record.seq)
    #print(seqID)
    #print(fasta)
    if len(fasta[0])!=len(fasta[1]):
        raise Exception("Sequences need to be of the same length")
    else:
        print("Sequences are of same length. Performing pairwise alignment")

    return(fasta)

def run_alignment():

    func=parse_arguments()

    sys.stdout=open(func[2], 'w')

    fasta_list=quality_check()

    aligner=Align.PairwiseAligner()
    alignment= aligner.align(fasta_list[0], fasta_list[1])
    for alignments in alignment:
        print("Score = %.1f" %alignments.score)
        print(alignments)

def run_mast():

    argument_list=parse_arguments()
    os.system('module load meme')
    os.system('mast -o '+ argument_list[2]+ ' -nostatus -minseqs 1 -remcorr -ev 10.0 ' +argument_list[1]+ ' '+argument_list[0])
    #output needs to be stored here and then checked 


#def analyze_mast_output():
 #   pass

if __name__=='__main__':
    main()

