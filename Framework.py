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

def analyze_mast_output():
    seqID=[]
    outfile_location=parse_arguments()
    mast_hit_file= outfile_location[2]
    path_out=outfile_location[3]
    sys.stdout=open("Differential_analysis_file.txt", 'w')

    for record in SeqIO.parse(path_out+mast_hit_file, "fasta"):
        seqID.append(record.id)

    with open(path_out+mast_hit_file, 'r') as file:
        list1=[]
        list2=[]
        for each in file:
            if re.search('# sequence_name', each):
                  header=each
                  header=re.sub(r'# sequence_name',"sequence_name", each)
            
            if not re.search('#', each):
                  re.sub(r'\n',"", each)
                  if str(each.split()[0])==str(seqID[0]):
                        list1.append(each)
                  elif str(each.split()[0])==str(seqID[1]):
                        list2.append(each)

        i=len(list1)
        j=len(list2)
        if i == j:
              print("The two sequences have equal number of TF binding motifs, check further for details regard")
              while i>=0:
                    var1=list1[i-1]
                    var2=list2[i-1]
                    i=i-1
                    if str(var1) != str(var2):
                          #this is for difference in TF binding motif itself
                          if var1.split()[2] != var2.split()[2]:
                                print(var1)
                                print (var2, '\n')

                                print(f"There is a difference in binding of the TF {var1.split()[2]} vs {var2.split()[2]} between the two sequences")
                                print(f"The position where this difference occurs is at the {var1.split()[1]}")
                                #print(f"There is a difference in binding of the TF {var1.split()[2]} vs {var2.split()[2]} between the two sequences")
                          elif var1.split()[2] == var2.split()[2]:
                                if var1.split()[1] != var2.split()[1]:
                                      
                                      print(f"There is no difference in binding of the TF, but the position of binding has shifted from {var1.split()[1]} to {var2.split()[1]} between the two sequences")
                                if var1.split()[-1] != var2.split()[-1]:
                                      print(f"There is no difference in binding of the TF, but the p-value of binding (which measures how strongly a TF can bind to a region) has changed from {var1.split()[-1]} to {var2.split()[-2]} between the two sequences")

                    else:
                          print("Sorry, not much difference in your DNA sequences in terms of Motif binding")


if __name__=='__main__':
    main()

