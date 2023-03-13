#/home/usr/bin/python
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import SeqIO
import os
import subprocess
import argparse
import tempfile
import shutil
import logging
import pandas as pd


logger=logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter= logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
file_handler= logging.FileHandler('DiffMotifFinder.log')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)


def main():
    usr_input=parse_args()
    temp_dir= make_temp_dir(usr_input.outdir)

    try:
        fa_len_check=qual_check(usr_input.infile1, usr_input.infile2, usr_input.TFdatabase,usr_input.outdir,temp_dir)
        mast_out=run_mast(usr_input.infile1, usr_input.infile2, usr_input.TFdatabase,usr_input.outdir,temp_dir,fa_len_check)
        analysis=compare_mast(usr_input.infile1, usr_input.infile2, usr_input.TFdatabase,usr_input.outdir,temp_dir,fa_len_check,mast_out)

    except:
       print("Something wrong. Cleaning up temp files")

    finally:
      clean_temp=clean_temp_dir(temp_dir)


def parse_args():
    
    parser=argparse.ArgumentParser(prog='Differential Motif finder', description="Provide fasta file of containing a pair of simialr DNA sequences of equal lengths", epilog='none')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('-i', '--fasta_file1', action='store', dest='infile1', type=str, required=True, help='The first fasta file for comparison')
    parser.add_argument('-j', '--fasta_file2', action='store', dest='infile2', type=str, required=True, help='The second fasta file for comparison')
    parser.add_argument('-d', '--database', action='store', dest='TFdatabase', type=str, required=True, help="path to the PWM file for tf in MEME format")
    parser.add_argument('-o', '--output', action='store', dest='outdir', type=str, help="name of the output directory", default="out_dir")
    args=parser.parse_args()

    logger.debug('Imported arguments from user')

    return(args)

def make_temp_dir(outdir):
    os.mkdir(os.path.join(outdir)+os.path.dirname(outdir))
    temp_dir=tempfile.mkdtemp(dir=os.path.abspath(outdir))
    return(temp_dir)


def clean_temp_dir(temp_dir):
    shutil.rmtree(temp_dir)


def qual_check(fasta1, fasta2, TFdatabase, outdir, temp_dir):

    with open (os.path.abspath(temp_dir)+'/mast_fasta_file','w') as multi_fa:
        for record1 in SeqIO.parse(fasta1, "fasta"):
            multi_fa.writelines('>'+str(record1.id)+'\n')
            multi_fa.writelines(str(record1.seq)+'\n')
            
        for record2 in SeqIO.parse(fasta2, "fasta"):
            multi_fa.writelines('>'+str(record2.id)+'\n')
            multi_fa.writelines(str(record2.seq)+'\n')
    
    if len(str(record1.seq))!=len(str(record2.seq)):
        raise Exception("Sequences need to be of the same length. Unequal sequences will have inconsistent results! Nonetheless, performing Alignemnt")
        logger.debug("Sequences are not of same length")
    else:
        logger.debug("Sequences are of same length. Performing pairwise alignment")
        
    with open (os.path.abspath(outdir)+'/alignment.txt', 'w') as aligned:
        alignment= pairwise2.align.globalxx(str(record1.seq),str(record2.seq))
        aligned.write(format_alignment(*alignment[0]))
    logger.debug("Alignment completed and written to file")

    return(multi_fa.name)


def run_mast(fasta1, fasta2, TFdatabase, outdir, temp_dir, multi_fa):
    subprocess.run(['mast','-nostatus','-minseqs','1','-remcorr','-ev','10.0','-o',f"{temp_dir+'/mastout'}",TFdatabase,multi_fa])
    logger.debug("MAST run complete")
    with open(temp_dir+'/TF_list.txt', 'w') as hits:
        subprocess.run(['mast','-hit_list','-nostatus','-minseqs','1','-remcorr','-ev','10.0',TFdatabase,multi_fa], stdout=hits)
    return(hits.name)


def compare_mast(fasta1, fasta2, TFdatabase, outdir, temp_dir, multi_fa, hits):
    hits_df= pd.read_table (hits, skiprows=2, skipfooter=1, sep=r'\s{1,}', engine='python', header=None,
                   names=["sequence_name","strand_sign_motif", "TF_id", "alt_id", "hit_start", "hit_end", "score", "pvalue"],
                   index_col=False)

    hits_differences= hits_df.drop_duplicates(subset=["strand_sign_motif", "TF_id", "alt_id", "hit_start", "hit_end", "score", "pvalue"], keep=False, ignore_index=True)
    seq_IDs=hits_differences.sequence_name.unique()

    hits_diff_set1 = hits_differences[hits_differences.sequence_name==seq_IDs[0]]
    hits_diff_set1 = hits_diff_set1.reset_index()
    hits_diff_set1= hits_diff_set1.drop(columns="index")
    hits_diff_set2 = hits_differences[hits_differences.sequence_name==seq_IDs[1]]
    hits_diff_set2 = hits_diff_set2.reset_index()
    hits_diff_set2= hits_diff_set2.drop(columns="index")


    if len(hits_diff_set1) >= len(hits_diff_set2):
        hits_bigger =hits_diff_set1
        hits_smaller=hits_diff_set2
    else:
        hits_bigger =hits_diff_set2
        hits_smaller=hits_diff_set1

    i=len(hits_bigger)
    j=len(hits_smaller)

    diff_category=[]
    analysis=[]
    sequence_id=[]
    result_dict={'sequence_name2':sequence_id, 'differential_category':diff_category,'analysis':analysis}

    while i>0:
        while j>0:
            if (hits_bigger.TF_id[i-1]==hits_smaller.TF_id[j-1]) & (hits_bigger.hit_start[i-1]==hits_smaller.hit_start[j-1]) & (hits_bigger.hit_end[i-1]==hits_smaller.hit_end[j-1]):
                sequence_id.append(hits_bigger.sequence_name[i-1])
                diff_category.extend(["hit_p-value"])
                analysis.extend([f'p-value of binding changed from {hits_bigger.pvalue[i-1]} to {hits_smaller.pvalue[i-1]}'])
                j=j-1
        
            elif (hits_bigger.TF_id[i-1]==hits_smaller.TF_id[j-1]) & (hits_bigger.hit_start[i-1]!=hits_smaller.hit_start[j-1]) or (hits_bigger.hit_end[i-1]==hits_smaller.hit_end[j-1]):
                sequence_id.append(hits_bigger.sequence_name[i-1])
                diff_category.extend(["hit_start hit_end"])
                analysis.extend([f'TF binding position changed from {hits_bigger.hit_start[i-1]} and {hits_bigger.hit_end[i-1]} to {hits_smaller.hit_start[j-1]} and {hits_smaller.hit_end[j-1]}'])
                j=j-1

            elif (hits_bigger.TF_id[i-1]!=hits_smaller.TF_id[j-1]):
                sequence_id.append(hits_bigger.sequence_name[i-1])
                diff_category.append("TF_id")
                analysis.append(f'The Transcription factor {hits_bigger.TF_id[i-1]} is present only in {hits_bigger.sequence_name[i-1]}. Check the adjoining position for more info')
                j=j-1
            
            i=i-1

    result_dict_df=pd.DataFrame(result_dict)
    final_out=pd.concat([hits_bigger, result_dict_df], axis=1)
    final_out=final_out.drop(columns="sequence_name2")
    final_out.to_csv(os.path.join(outdir, 'analysis.txt'), sep='\t')


if __name__=='__main__':
    main()




