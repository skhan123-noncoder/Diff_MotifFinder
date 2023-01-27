#/home/usr/bin/python

from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import SeqIO
import os
import sys
import argparse

def main():
     
    user_provided_inputs=parse_arguments()
    fasta_length_check=quality_check(user_provided_inputs.infile1, user_provided_inputs.infile2)
    alignment_check=run_alignment(fasta_length_check[0][0],fasta_length_check[0][1])
    primary_outfile=run_mast(user_provided_inputs.infile1, user_provided_inputs.infile2, user_provided_inputs.TFdatabase, user_provided_inputs.outfile)
     #func4=analyze_mast_output()

def parse_arguments():
    parser=argparse.ArgumentParser(prog='Differential Motif finder', description="Provide fasta file of containing a pair of simialr DNA sequences of equal lengths", epilog='none')
    parser.add_argument('-i', '--fasta_file1', action='store', dest='infile1', type=str, required=True, help='The first fasta file for comparison')
    parser.add_argument('-j', '--fasta_file2', action='store', dest='infile2', type=str, required=True, help='The second fasta file for comparison')
    parser.add_argument('-d', '--database', action='store', dest='TFdatabase', type=str, required=True, help="path to the PWM file for TFs in MEME format")
    parser.add_argument('-o', '--output', action='store', dest='outfile', type=str, help="name of the output file", default="out_file")
    args=parser.parse_args()
    path=str(args.outfile)
    return(args)

def quality_check(fasta1, fasta2):
    seqID=[]
    fasta=[]
    for record1 in SeqIO.parse(fasta1, "fasta"):
        seqID.append(record1.id)
        fasta.append(record1.seq)
    for record2 in SeqIO.parse(fasta2, "fasta"):
        seqID.append(record2.id)
        fasta.append(record2.seq)
    if len(fasta[0])!=len(fasta[1]):
        raise Exception("Sequences need to be of the same length")
    else:
        print("Sequences are of same length. Performing pairwise alignment")

    return(fasta,seqID)

def run_alignment(fasta_seq1,fasta_seq2):
    alignment= pairwise2.align.globalxx(fasta_seq1, fasta_seq2)
    print(format_alignment(*alignment[0]))

#This function calls the meme module and runs the mast function
##Unable to load the meme module while running this. Need to have a closer look at this
def run_mast(fasta1,fasta2,TF,outfile):

    subprocess.run('mast -o '+outfile+ ' -nostatus -minseqs 1 -remcorr -ev 10.0 ' +TF+ ' '+fasta1)
    subprocess.run('mast -o '+outfile+ ' -nostatus -minseqs 1 -remcorr -ev 10.0 ' +TF+ ' '+fasta2)
    subprocess.run('mast -hit_list -nostatus -minseqs 1 -remcorr -ev 10.0 HOCOMOCOv11_core_HUMAN_mono_meme_format.meme fasta_test_2seq.txt >'+ argument_list[2]+'/hit_list_temp.txt') 

def analyze_mast_output():
    seqID=[]
    outfile_location=parse_arguments()
    sys.stdout=open("Differential_analysis_file.txt", 'w')

    for record in SeqIO.parse(outfile_location[0], "fasta"):
        seqID.append(record.id)

    with open(outfile_location[2]+'/hit_list_temp.txt', 'r') as file:
        list1=[]
        list2=[]
        for each in file:
            if re.search('# sequence_name', each):
                header=each
                header=re.sub(r'# sequence_name',"sequence_name", each)
            
            if not re.search('#', each):
                if str(each.split()[0])==str(seqID[0]):
                    list1.append(re.sub(str(seqID[0]), "", each))
                elif str(each.split()[0])==str(seqID[1]):
                    list2.append(re.sub(str(seqID[1]), "", each))
        i=len(list1)
        j=len(list2)
        if i == j:
            print("The two sequences have equal number of TF binding motifs, check further for detailed analysis")
            diff=0
            while i>=0:
                var1=str(list1[i-1])
                var2=str(list2[i-1])
                if var1 != var2:
                    diff=diff+1
                print('\n')
                #print(f"################  Difference number {diff} #######################")
                sys.stdout.write(header)
                sys.stdout.write(str(seqID[0])+var1)
                print(str(seqID[1])+var2)
                if var1.split()[1] != var2.split()[1]:
                    print("Analysis: There is a difference in binding of the TF {var1.split()[1]} vs {var2.split()[1]} between the two sequences")
                    print("Analysis: The position where this difference occurs is at the {var1.split()[0]}")
                    
                elif var1.split()[1] == var2.split()[1]:
                    if var1.split()[0] != var2.split()[0]:
                        print("Analysis: There is no difference in binding of the TF, but the position of binding has shifted from {var1.split()[1]} to {var2.split()[1]} between the two sequences")
                    if var1.split()[-1] != var2.split()[-1]:
                        print("Analysis: There is no difference in binding of the TF, but the p-value of binding has changed from {var1.split()[-1]} to {var2.split()[-2]} between the two sequences")
                    if var1.split()[3] != var2.split()[3] or var1.split()[4] != var2.split()[4]:
                        print("Analysis: There is no difference in binding of the TF, but the position of binding has chnaged from {var1.split()[3]}&{var1.split()[4]} to {var2.split()[3]}&{var2.split()[4]}  between the two sequences")
                    
            i=i-1

        if diff == 0:
            sys.stdout.write("Sorry, not much difference in your DNA sequences in terms of Motif binding")
   
                          

if __name__=='__main__':
    main()

