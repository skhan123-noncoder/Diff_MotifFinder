#/home/usr/bin/python
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import SeqIO
import os
import subprocess
import argparse
import tempfile
import re

def main():
     
    user_provided_inputs=parse_arguments()
    fasta_length_check=quality_check(user_provided_inputs.infile1, user_provided_inputs.infile2)
    alignment_check=run_alignment(fasta_length_check[0][0],fasta_length_check[0][1],fasta_length_check[3])
    primary_outfile=run_mast(fasta_length_check[2], user_provided_inputs.TFdatabase, user_provided_inputs.outfile,fasta_length_check[3])
    func4=analyze_mast_output(primary_outfile,fasta_length_check[1][0],fasta_length_check[1][1],user_provided_inputs.outfile)

def parse_arguments():
    parser=argparse.ArgumentParser(prog='Differential Motif finder', description="Provide fasta file of containing a pair of simialr DNA sequences of equal lengths", epilog='none')
    parser.add_argument('-i', '--fasta_file1', action='store', dest='infile1', type=str, required=True, help='The first fasta file for comparison')
    parser.add_argument('-j', '--fasta_file2', action='store', dest='infile2', type=str, required=True, help='The second fasta file for comparison')
    parser.add_argument('-d', '--database', action='store', dest='TFdatabase', type=str, required=True, help="path to the PWM file for TFs in MEME format")
    parser.add_argument('-o', '--output', action='store', dest='outfile', type=str, help="name of the output file", default="out_file")
    args=parser.parse_args()
    path=str(args.outfile)

    return(args)

#Checking for the length of the fasta files. If unequal, the analysis might be different. Needs to have more robustness here
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
        raise Exception("Sequences need to be of the same length. Unequal sequences will have inconsistent results!")
    else:
        print("Sequences are of same length. Performing pairwise alignment")
    
    #Creating a temp directory where you store your combined fasta file and store the 
    temp_dir=tempfile.mkdtemp(dir=os.path.dirname(fasta1))
    mast_fa_input=open(os.path.abspath(temp_dir)+'/mast_fasta_file','w')
    with open(fasta1,'r') as f1:
        for each in f1.readlines():
            mast_fa_input.write(each)
        mast_fa_input.write('\n')
    f1.close()
    with open(fasta2,'r') as f2:
        for each in f2.readlines():
            mast_fa_input.write(each)
        mast_fa_input.write('\n')
    f2.close()
    mast_fa_input.close()
    
    return(fasta,seqID,mast_fa_input.name,temp_dir)

def run_alignment(fasta_seq1,fasta_seq2,temp_dir):

    aligned_file=open(os.path.abspath(temp_dir)+'/alignment.txt', 'w')
    alignment= pairwise2.align.globalxx(fasta_seq1, fasta_seq2)
    aligned_file.write(format_alignment(*alignment[0]))


def run_mast(mast_fa_input,TF,outfile,temp_dir):
    subprocess.run(['mast','-nostatus','-minseqs','1','-remcorr','-ev','10.0','-o',f"{outfile}",TF,mast_fa_input])
    hit=open(temp_dir+'/hit_list.txt', 'w')
    subprocess.run(['mast','-hit_list','-nostatus','-minseqs','1','-remcorr','-ev','10.0',TF,mast_fa_input], stdout=hit)
    hit.close()
    return(hit.name)

def analyze_mast_output(hit_list,seq1,seq2,outfile):
    with open(hit_list,'r') as hit_file:
        comp_list1=[]
        comp_list2=[]
        for each in hit_file:
            if re.search('# sequence_name', each):
                header=re.sub(r'# sequence_name',"sequence_name", each)
            if each.split()[0]==seq1:
                comp_list1.append(each)
            elif each.split()[0]==seq2:
                comp_list2.append(each)
    hit_file.close()

    analysis=open(os.path.abspath(outfile)+'/analysis_file.txt', 'w')

    if len(comp_list1)==len(comp_list2):
        analysis.write("The two sequences have equal number of TF binding motifs, check further for detailed analysis")
    analysis.write('\n')
            
    TF_entries=len(comp_list1)
    while TF_entries>=0:
        if comp_list1[TF_entries-1].split()[1:]!=comp_list2[TF_entries-1].split()[1:]:
            if comp_list1[TF_entries-1].split()[2]!=comp_list2[TF_entries-1].split()[2]:
                analysis.write(header)
                #analysis.write('\n')
                analysis.write('%s'%(comp_list1[TF_entries-1]))
                analysis.write('%s'%(comp_list2[TF_entries-1]))
                analysis.write('\n')
                analysis.write('Analysis: There is a difference in binding of the TF %s vs %s between the sequences' %({comp_list1[TF_entries-1].split()[2]},{comp_list2[TF_entries-1].split()[2]}))
                #analysis.write('\n')
                analysis.write(f'Analysis: The position where this difference occurs is at the position: %s' %({comp_list1[TF_entries-1].split()[1]}))
                analysis.write('\n')
            analysis.write('\n')        
            if comp_list1[TF_entries-1].split()[2] == comp_list2[TF_entries-1].split()[2]:
                analysis.write(header)
                analysis.write('%s'%(comp_list1[TF_entries-1]))
                analysis.write('%s'%(comp_list2[TF_entries-1]))
                analysis.write('\n')
                if comp_list1[TF_entries-1].split()[1] != comp_list2[TF_entries-1].split()[1]:
                    analysis.write("Analysis: There is no difference in binding of the TF, but the position of binding has shifted from %s to %s between the two sequences" %({comp_list1[TF_entries-1].split()[1]},{comp_list2[TF_entries-1].split()[1]}))
                    analysis.write('\n')
                if comp_list1[TF_entries-1].split()[-1] != comp_list2[TF_entries-1].split()[-1]:
                    analysis.write("Analysis: There is no difference in binding of the TF, but the p-value of binding has changed from %s to %s between the two sequences" %({comp_list1[TF_entries-1].split()[-1]},{comp_list2[TF_entries-1].split()[-2]}))
                    analysis.write('\n')
                if comp_list1[TF_entries-1].split()[4] != comp_list2[TF_entries-1].split()[4] or comp_list1[TF_entries-1].split()[5] != comp_list2[TF_entries-1].split()[5]:
                    analysis.write("Analysis: There is no difference in binding of the TF, but the position of binding has chnaged from %s & %s to %s and %s between the two sequences" %({comp_list1[TF_entries-1].split()[4]},{comp_list1[TF_entries-1].split()[5]},{comp_list2[TF_entries-1].split()[4]},{comp_list2[TF_entries-1].split()[5]}))
                    analysis.write('\n')
            analysis.write('\n')
        TF_entries=TF_entries-1

    analysis.close()                

if __name__=='__main__':
    main()

