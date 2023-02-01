#/home/usr/bin/python
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import SeqIO
import os
import subprocess
import argparse
import tempfile
import re
import shutil

def main():
     
    usr_input=parse_args()
    fa_len_check=qual_check(usr_input.infile1, usr_input.infile2)
    align=run_align(fa_len_check[0][0],fa_len_check[0][1],fa_len_check[3])
    mast_out=run_mast(fa_len_check[2], usr_input.TFdatabase, usr_input.outfile,fa_len_check[3])
    analysis=compare_mast(fa_len_check[3],mast_out,fa_len_check[1][0],fa_len_check[1][1],usr_input.outfile,align)

#Ask fasta files, TF database and name of the output file from the user
def parse_args():
    parser=argparse.ArgumentParser(prog='Differential Motif finder', description="Provide fasta file of containing a pair of simialr DNA sequences of equal lengths", epilog='none')
    parser.add_argument('-i', '--fasta_file1', action='store', dest='infile1', type=str, required=True, help='The first fasta file for comparison')
    parser.add_argument('-j', '--fasta_file2', action='store', dest='infile2', type=str, required=True, help='The second fasta file for comparison')
    parser.add_argument('-d', '--database', action='store', dest='TFdatabase', type=str, required=True, help="path to the PWM file for tf in MEME format")
    parser.add_argument('-o', '--output', action='store', dest='outfile', type=str, help="name of the output file", default="out_file")
    args=parser.parse_args()

    return(args)

#Checking for the length of the fasta files. If unequal, the analysis might be different. Needs to have more robustness here
def qual_check(fasta1, fasta2):
    seqID=[]
    seq=[]
    for record1 in SeqIO.parse(fasta1, "fasta"):
        seqID.append(record1.id)
        seq.append(record1.seq)        
    for record2 in SeqIO.parse(fasta2, "fasta"):
        seqID.append(record2.id)
        seq.append(record2.seq)
    if len(seq[0])!=len(seq[1]):
        raise Exception("Sequences need to be of the same length. Unequal sequences will have inconsistent results!")
    else:
        print("Sequences are of same length. Performing pairwise alignment")
    
    #Creating a temp directory where you store your combined fasta file and store the alignemnt file as well. This directory will be deleted at the end
    temp_dir=tempfile.mkdtemp(dir=os.path.dirname(fasta1))
    multi_fa=open(os.path.abspath(temp_dir)+'/mast_fasta_file','w')
    with open(fasta1,'r') as f1:
        for each in f1.readlines():
            multi_fa.write(each)
        multi_fa.write('\n')
    f1.close()
    with open(fasta2,'r') as f2:
        for each in f2.readlines():
            multi_fa.write(each)
        multi_fa.write('\n')
    f2.close()
    multi_fa.close()
    
    return(seq,seqID,multi_fa.name,temp_dir)

#Program to align the fasta files provided by the user. This can provide some hint into the differences in the two sequences. This will be present in the out folder
def run_align(seq1,seq2,temp_dir):

    aligned=open(os.path.abspath(temp_dir)+'/alignment.txt', 'w')
    alignment= pairwise2.align.globalxx(seq1, seq2)
    aligned.write(format_alignment(*alignment[0]))
    aligned.close()

    return(aligned)

#The meme module needs to be loaded here. The mast output file is sent to the outfile described by the user. The hitlist is sent to the temp dir
def run_mast(multi_fa,TF,outfile,temp_dir):
    subprocess.run(['mast','-nostatus','-minseqs','1','-remcorr','-ev','10.0','-o',f"{outfile}",TF,multi_fa])
    hits=open(temp_dir+'/TF_list.txt', 'w')
    subprocess.run(['mast','-hit_list','-nostatus','-minseqs','1','-remcorr','-ev','10.0',TF,multi_fa], stdout=hits)
    hits.close()
    return(hits.name)


#Analyze the hit_list file and point out the differences. The alignment file is moved from the temporary folder to the output folder and temp dir is deleted
def compare_mast(temp_dir,hits,seqID_1,seqID_2,outfile,aligned):
    with open(hits,'r') as hit_file:
        hits_1=[]
        hits_2=[]
        for each in hit_file:
            if re.search('# sequence_name', each):
                header=re.sub(r'# sequence_name',"sequence_name", each)
            if each.split()[0]==seqID_1:
                hits_1.append(each)
            elif each.split()[0]==seqID_2:
                hits_2.append(each)
    hit_file.close()

    analysis=open(os.path.abspath(outfile)+'/analysis_file.txt', 'w')

    if len(hits_1)==len(hits_2):
        analysis.write("The two sequences have equal number of TF binding motifs, check further for detailed analysis\n")
    analysis.write("###################################################################################################\n")
            
    tf=len(hits_1)
    while tf>=0:
        if hits_1[tf-1].split()[1:]!=hits_2[tf-1].split()[1:]:
            if hits_1[tf-1].split()[2]!=hits_2[tf-1].split()[2]:
                analysis.write(f'{header}\n')
                analysis.write(f'{hits_1[tf-1]}')
                analysis.write(f'{hits_2[tf-1]}\n')
                analysis.write('Analysis: There is a difference in binding of the TF %s vs %s between the sequences' %({hits_1[tf-1].split()[2]},{hits_2[tf-1].split()[2]})+'\n')
                analysis.write(f'Analysis: The position where this difference occurs is at the position: %s' %({hits_1[tf-1].split()[1]})+'\n')
                analysis.write("#######################################################################################################################################\n")        
            if hits_1[tf-1].split()[2] == hits_2[tf-1].split()[2]:
                analysis.write(f'{header}\n')
                analysis.write(f'{hits_1[tf-1]}')
                analysis.write(f'{hits_2[tf-1]}\n')
                if hits_1[tf-1].split()[1] != hits_2[tf-1].split()[1]:
                    analysis.write("Analysis: There is no difference in binding of the TF, but the position of binding has shifted from %s to %s between the two sequences" %({hits_1[tf-1].split()[1]},{hits_2[tf-1].split()[1]})+'\n')
                    analysis.write("#######################################################################################################################################\n")
                if hits_1[tf-1].split()[-1] != hits_2[tf-1].split()[-1]:
                    analysis.write("Analysis: There is no difference in binding of the TF, but the p-value of binding has changed from %s to %s between the two sequences" %({hits_1[tf-1].split()[-1]},{hits_2[tf-1].split()[-2]})+'\n')
                    analysis.write("#######################################################################################################################################\n")
                if hits_1[tf-1].split()[4] != hits_2[tf-1].split()[4] or hits_1[tf-1].split()[5] != hits_2[tf-1].split()[5]:
                    analysis.write("Analysis: There is no difference in binding of the TF, but the position of binding has chnaged from %s & %s to %s and %s between the two sequences" %({hits_1[tf-1].split()[4]},{hits_1[tf-1].split()[5]},{hits_2[tf-1].split()[4]},{hits_2[tf-1].split()[5]})+'\n')
                    analysis.write("#######################################################################################################################################\n")
        tf=tf-1

    analysis.close()
    subprocess.run(['mv',str(aligned.name),str(os.path.abspath(outfile))])
    shutil.rmtree(temp_dir)

if __name__=='__main__':
    main()

