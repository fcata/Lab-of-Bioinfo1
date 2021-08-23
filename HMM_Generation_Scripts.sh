less ~/Downloads/rcsb_pdb_custom_report_20210811004133.csv  #to view the csv file report

awk -F "," '{print $(NF-1), $5, $3, $4}' ~/Downloads/rcsb_pdb_custom_report_20210811004133.csv | tr -d '"' |tail -n +2 >clean_pdb.txt #to remove the quotes and specify columns, then place into a new file named clean_pdb.txt

awk '{if ($4>=40 && $4<=80) print $0}' clean_pdb.txt >clean_pdb_len40_80.txt #to print only those that are within the range we specified and redirect to a new filed named as such

#Go to pdbefold or rupee and create a csv of the target/reference protein 3TGI, export to csv

#Open it with
less ~/Downloads/data.csv  #output:
    # Chain Id,RMSD,TM-Score,TM-Align
    # 1,3tgiI,0.00,1.00,text 3d pdb
    # 2,1f7zI,0.11,1.00,text 3d pdb

cut -d "," -f 2-4 ~/Downloads/data.csv # to specify columns and separate
# output:
    # Chain Id,RMSD,TM-Score
    # 3tgiI,0.00,1.00
    # 1f7zI,0.11,1.00

awk -F ',' '{print $2,$3,$4}' ~/Downloads/data.csv >3tgi_aln.list # to specify those columns and direct into a new file named as such

awk -F ',' '{print $2,$3,$4}' ~/Downloads/data.csv |tail -n +2 >3tgi_aln.list # remove the first row with the column name
# output:
    # 3TGII 0.00 1.00
    # 1F7ZI 0.11 1.00
    # 1FY8I 0.12 1.00

awk -F ',' '{print toupper($2),$3,$4}' ~/Downloads/data.csv |tail -n +2 >3tgi_aln.list # to make upper case like above!

#Now we have the list of 3TGI IDs and their scores….next we will get the list of those from the PDB of len40-80 and their IDs sorted
awk '{print $1$2}' clean_pdb_len40_80.txt >list1 #puts only the IDs into a new file named list1
awk '{print $1$2}' clean_pdb_len40_80.txt |sort -u >list1 #to sort the list
awk '{print $1}' 3tgi_aln.list  |sort -u >list2 # do the same for the 3TGI IDs and put into a new list called list2
comm list1 list2 # Now use comm to compare the two sorted files:

#The above output is three columns where the first column is separated by zero tab and contains names only present in file1.txt
# the second column contains names only present in file2.txt and separated by one tab
# the third column contains names common to both the files and is separated by two tabs from the beginning of the line

comm <(sort list1) <(sort list2 ) |less  # to sort the two lists and comm
comm -12 <(sort list1) <(sort list2) |less # suppress printing of column 1 and 2
comm -12 <(sort list1) <(sort list2) >comm_list.txt # put into a new file only the ones that overlap in column 3, named as such comm_list

#Make a python script that will allow you to merge the two lists together by selecting just the sequence from the 3TGI and the PDB IDs from the list

#!/usr/bin/env python3
import sys


def get_dic(filename):
    d = {}
    f = open(filename)
    for line in f:
        v = line.rstrip().split()
        d[v[0]] = d.get(v[0],[])
        d[v[0]].append(line)
    return d

def get_common(d1,d2):
    s1 = set(list(d1.keys()))
    s2 = set(list(d2.keys()))
    c = list(s1.intersection(s2))
    # print(len(list(c)))
    for i in c:
        print(d1[i][0].rstrip())


if __name__ == "__main__":
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    d1 = get_dic(file1)
    d2 = get_dic(file2)
    get_common(d1,d2)

python compare.py list1 list2 # Call it in shell to get the list of common sequences between the two
python compare.py 3tgi_aln.list list1 |less # you can do the same with the sequences  and their IDS from the reference protein and the list from list1 which is the cleaned exhaustive PDB files
python3 compare.py clean_pdb.seq list2 |awk '{print ","$1;print $2}' |less # taking the cleaned PDB file  and comparing them to the seq of the PDB reference template protein IDs which gives you a fasta file of all the sequences that are in common from these two lists!
python3 compare.py clean_pdb.seq list2 |awk '{print ","$1;print $2}' >comm_seq.fasta # Put that into a .fasta file!
python3 compare.py 3tgi_aln.list list1 |sort -nk 2 |less # Double check those with a RMSD score that is over 1.54 and see if you want to remove them or not as the lower the RMSD the better as it indicates a better fit and can predict the model data more accurately
#Now with the fasta file we want to perform a cluster on the, using blastclust !OR! CDHIT:

blastclust -i comm_seq.fasta -o comm_seq.cluster -S 95 -L 0.8
blastclust -i comm_seq.fasta -o comm_seq.cluster -S 95 -L 0.8 # Start clustering on the queries - We can do the same on blastclust using this command
head -n 1 comm_seq.cluster |awk '{split($0,a," "); for (i=1; i<=length(a);i++) {print substr(a[i],1,4)}}' |sort -u >cluster_1
python compare.py resolu.idx cluster_1 #to find the common IDs between cluster 1 and the resolu.idx

wget "http://weizhong-lab.ucsd.edu/cdhit-web-server/cgi-bin/result.cgi?JOBID=1628857817" -O cd-hit.cluster #to pull the output from CDHIT anf place into a file cd-hit.cluster

awk '{if (substr($0,1,1)==">") {print ""} else {printf "%s ", substr($3,2,5)}}' cd-hit.cluster  |tail -n +2 >cdhit-seq.list #leaves you only with the PDB IDs without anything else, I.e a list of the CD-HIT

awk '{print substr($1,1,4)}' ../cdhit-seq.list |sort -u >list_pdb.txt #to select only the FIRST IDs of each cluster to be titled as the REPRESENTATIVE protein for that cluster

# Now to start the multiple sequence alignment of the structures. First by downloading the structures from PDB and extracting the chains of interest
for i in 'cat list_pdb.txt'; do wget -q https://files.rcsb.org/view/$i.pdb; done #this will extract the chains from the PDB based off the PDB IDS in the text
awk '{print substr($1,1,4),substr($1,5,1)}' ../cdhit-seq.list |sort -u >list_chain.txt #extract that ifo and place into a new text file for the chains

#create a file to run from the shell
#!/bin/bash
pdbfile=$1
chain=$2
awk -v c=$chain '{if ((substr($0,1,4)=="ATOM" || substr($0,1,3)=="TER") && substr($0,22,1)==c) print $0}' $pdbfile
pdbfile=$1
chain=$2
---------------
chmod a+x selch.sh # if permission is needed
awk '{print "./selch.sh",$1.pdb,$2" >chains"$1$2".pdb"}' list_chain.txt >run.sh #run the file from shell getting the chains of interest

head -n 179 positives.rids >set_pos1.ids
tail -n 180 positives.rids >set_pos2.ids
head -n 277821 negatives.rids >set_neg1.ids
tail -n 277822 negatives.rids >set_neg2.ids
#THIS IS THE SPLITTING^

#next is mTmAlign, #submit a zipped file only of the PDB IDs and their CHAINS
wget https://yanglab.nankai.edu.cn/mTM-align/output/mTM015676/seq.fasta -O tm-ali.fasta # downloading the output of the multiple sequence alignment from mTM-align (file tm-ali.fasta) and cleaning to remove misalignments
awk '{if (substr($0,1,1)==">") {printf "\n%s ",$0} else {printf "%s",$0}}' tm-ali.fasta |awk '{print substr($1,1,6);print substr($2,28,62)}' |tail -n +3 >bpti-kunitz.ali #make bpti-kuniiz alignment file that will be submitted to HMMER


#hmmbuild run to generate the model, call it bpti-kunitz.hmm, perfomance testing
hmmbuild bpti-kunitz.hmm bpti-kunitz.ali

#download from uniprot reviewed and unreviewwed to create the two lsits of positives and negatives
zcat ~/Downloads/uniprot-reviewed_yes+database+\(type_pfam+pf00014\)+length_\{40+TO+_%--.list.gz |less
zcat ~/Downloads/uniprot-database%3A%28type%3Apfam+pf00014%29+length%3A%5B40+TO+\*%5D+reviewed%3Ay--.list |less
gunzip -c ~/Downloads/uniprot-database_\(type_pfam+pf00014\)+length_\[40+TO+_\]+reviewed_y--.list.gz >positives.txt
mv ~/Downloads/uniprot-database_\(type_pfam+pf00014\)+length_\[40+TO+_\]+reviewed_y--.fasta.gz positives.fasta.gz
gunzip positives.fasta.gz
cat positives.fasta |awk '{if (substr($0,1,1)==">") {split($0,a,"|");print a[2]} }' |less
less positives.fasta >positives.ids

wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
du -h ~/Downloads/uniprot-NOT+database_\(type_pfam+pf00014\)+length_\[40+TO+_\]+reviewed--.fasta.gz
mv ~/Downloads/uniprot-NOT+database_\(type_pfam+pf00014\)+length_\[40+TO+_\]+reviewed--.fasta.gz negatives.fasta.gz
gunzip -c negatives.fasta.gz |grep "^>"
gunzip -c negatives.fasta.gz |awk '{if (substr($0,1,1)==">") {split($0,a,"|");print a[2]} }' >negatives.ids


# 2-fold cross-validation: splitting the list of identifiers in 2 equal parts both for the positive and the negative set to perform the optimization on one set, identifying the optimal threshold and testing on the second set
cat positives.fasta |awk '{if (sub-str($0,1,1)==">") {split($0,a,"|");print a [2]} }' >positives.ids
cat negatives.fasta |awk '{if (sub-str($0,1,1)==">") {split($0,a,"|");print a [2]} }' >negatives.ids
sort -R positives.ids >positives.rids
sort -R negatives.ids >negatives.rids
head -n 180 positives.rids >set_pos1.ids
tail -n +181 positives.rids >set_pos2.ids
head -n 277821 negatives.rids >set_neg1.ids
tail -n +277822 negatives.rids >set_neg2.ids
#THIS IS THE SPLITTING^

# RUN THE PYTHON SCRIPT  ON THE NEW SPLITTED IDS
#!/usr/bin/env python3
'''
import sys

def get_ids(fileids):
    d = {}
    listid = open(fileids).read().rstrip().split('\n')
    d = dict([(i, True) for i in listid])
    return d

def get_sequences(fileseq, listid):
    c=0
    f = open(fileseq)
    for line in f:
        if line.find('>')==0:
            pid=line.split('|')[1]
            if listid.get(pid):
                c=1
                print(">"+pid)
                continue
            else:
                c=0

        if c==1: print(line.rstrip())
    return

''''
if __name__ == "__main__":
    fileids= sys.argv[1]
    fileseq = sys.argv[2]
    listid = get_ids(fileids)
    get_sequences(fileseq, listid)'

    ''''

#this will select the sequences from each set alog with the IDs, placing into a fasta file
python3 select-seqs.py pdb/set_neg1.ids pdb/negatives.fasta >pdb/set_neg1.fasta
python3 select-seqs.py pdb/set_neg2.ids pdb/negatives.fasta >pdb/set_neg2.fasta
python3 select-seqs.py pdb/set_pos1.ids pdb/positives.fasta >pdb/set_pos1.fasta
python3 select-seqs.py pdb/set_pos2.ids pdb/positives.fasta >pdb/set_pos2.fasta


#Perform a HMMSEARCH on the FASTA splits to get some statistical values like E and P..etc
hmmsearch -Z 1 --max --noali -o set_pos1.hmmout bpti-kunitz.hmm set_pos1.fasta
hmmsearch -Z 1 --max --noali -o set_pos2.hmmout bpti-kunitz.hmm set_pos2.fasta
hmmsearch -Z 1 --max --noali -o set_neg2.hmmout bpti-kunitz.hmm set_neg2.fasta
hmmsearch -Z 1 --max --noali -o set_neg1.hmmout bpti-kunitz.hmm set_neg1.fasta

hmmsearch -Z 1 --max  -o set_pos1.hmmout bpti-kunitz.hmm set_pos1.fasta
hmmsearch -Z 1 --max  -o set_pos2.hmmout bpti-kunitz.hmm set_pos2.fasta
hmmsearch -Z 1 --max  -o set_neg2.hmmout bpti-kunitz.hmm set_neg2.fasta
hmmsearch -Z 1 --max  -o set_neg1.hmmout bpti-kunitz.hmm set_neg1.fasta


#To make the output of the sets a little clearer as a tabular output and removing the alignments
hmmsearch -Z 1 --noali --max --tblout set_pos1.out bpti-kunitz.hmm set_pos1.fasta
hmmsearch -Z 1 --noali --max --tblout set_pos2.out bpti-kunitz.hmm set_pos2.fasta
hmmsearch -Z 1 --noali --max --tblout set_pos1.out bpti-kunitz.hmm set_pos1.fasta
hmmsearch -Z 1 --noali --max --tblout set_neg2.out bpti-kunitz.hmm set_neg2.fasta


#taking the first hlaf and removing the tail in order to shape the dataset, and get only the E values from the best domains, 1 is for the positives and 0 is for the negatives
head -n 126988 set_neg1.hmmout |tail -n +19 |awk '{print $NF, $4, 0}' > set_neg1.res
head -n 126997 set_neg2.hmmout |tail -n +19 |awk '{print $NF, $4, 0}' >set_neg2.res
tail -n +19 set_pos1.hmmout |head -n 179 |awk '{print $NF, $4, 1}' > set_pos1.res
tail -n +19 set_pos2.hmmout |head -n 180 |awk '{print $NF, $4, 1}' > set_pos2.res

#to remove extra symbols:
grep -v "^#" set_pos1.out |awk '{print $1, $8, 1}' >set_pos1.res
grep -v "^#" set_pos2.out |awk '{print $1, $8, 1}' >set_pos2.res
grep -v "^#" set_neg2.out |awk '{print $1, $8, 0}' >set_neg2.res
grep -v "^#" set_neg1.out |awk '{print $1, $8, 0}' >set_neg1.res

#perform the reinclusion commands on the hmmsearch set that removed thousands of models
#Give a default of E value higher than 10, since we know it must be greater than 10. Using the comm command to match the initial list with the final list and eject only the elements in the final list.
comm -23 <(sort set_neg1.ids) <(awk '{print $1}'  set_neg1.res |sort) > set_neg1.add
comm -23 <(sort set_neg2.ids) <(awk '{print $1}'  set_neg2.res |sort) > set_neg2.add
#Then to add info about the e value we can do this instead:
comm -23 <(sort set_neg1.ids) <(awk '{print $1}'  set_neg1.res |sort) |awk '{print $1,10,0}' > set_neg1.add
comm -23 <(sort set_neg2.ids) <(awk '{print $1}'  set_neg2.res |sort) |awk '{print $1,10,0}' > set_neg2.add

#concatenate all the sets from 1 and 2 into two different files
cat set_neg1.res set_neg1.add set_pos1.res >set_all1.res
cat set_neg2.res set_neg2.add set_pos2.res >set_all2.res

#create a python script to calcualte the performance: Accuracy and MCC

'''
#!/usr/bin/env python3

import sys
import numpy as np



def get_preds(filename, sp=-2, rc=-1):
    pred_list = []
    f = open(filename)
    for line in f:
        v = line.rstrip().split()
        pred_list.append([v[0], float(v[sp]), int(v[rc])])
    return pred_list


def get_cm(data, th=0.1):
    cm = [[0, 0], [0, 0]]
    for i in data:
        # print(i)
        rc = i[-1]
        # if (i[-2] <=th):
        if (i[-2] >= th): #when evaluating all the sets
            pc = 1
        else:
            pc = 0
        cm[pc][rc] = cm[pc][rc] + 1
        # print(cm)
    return cm


def calculate_performance(cm):
    n = float(cm[0][0] + cm[1][1] + cm[0][1] + cm[1][0])
    d = np.sqrt((cm[0][0] + cm[0][1]) * (cm[0][0] + cm[1][0]) * (cm[1][1] + cm[1][0]) * (cm[1][1] + cm[0][1]))
    mcc = (cm[0][0] * cm[1][1] - cm[0][1] * cm[1][0]) / d
    acc = (cm[0][0] + cm[1][1]) / n
    return acc, mcc


def opt_th(pred_list):
    for i in range(20):
        cm = get_cm(pred_list, 10 ** -i)
        acc, mcc = calculate_performance(cm)
        print('TH:', 10 ** -i, 'ACC:', acc, 'MCC:', mcc)


if __name__ == "__main__":
    predfile = sys.argv[1]
    pred_list = get_preds(predfile)
    # # th = float(sys.argv[2])
    # opt_th(pred_list)
    if len(sys.argv) > 2:
        th = float(sys.argv[2])
        cm = get_cm(pred_list, th)
        acc, mcc = calculate_performance(cm)
        print('TH:', th, 'ACC:', acc, 'MCC:', mcc)
    else:
        opt_th(pred_list)
    # cm=get_cm(pred_list,th)
    # calculate_performance(cm)
    '''

#use 6 to just select the column for MCC value, the output will give the best E value, which can be the new threshold to comapre against the opposing set
python3 performance.py set_all1.res |sort -nrk 6 >set_all1.out
python3 performance.py set_all2.res |sort -nrk 6 >set_all2.out



#check for wrong predictions:
awk '{if ($2!=$3) print $0}' set_all.res > wrong_pred.res

# OUTPUT:
# P40500 1 0 -false +
# P78746 1 0 -false +
# C5H8E7 1 0 -false +
# O62247 0 1 -false -
# Q11101 0 1 -false -
# D3GGZ8 0 1 -false -



#Calc the alignment with these sequences specifically or look at the alignment of them
awk '{p=0;if ($2<1e-5) {p=1};print $1,p, $3}' set_all1.res > set_all.res
awk '{p=0;if ($2<1e-7) {p=1};print $1,p, $3}' set_all2.res >> set_all.res

#check for the best performance
python3 performance_all.py pdb/set_all.res 0.5 >best_set_all.out

#Look for wrong cases in set 1
awk '{if ($2<=1e-5 && $3==0 || ($2>1e-5 && $3==1))print $0}' set_all1.res
P40500 5.5e-06 0
P78746 9.2e-06 0
C5H8E7 6.2e-06 0
O62247 0.0011 1



#Wrong cases in SET 2
awk '{if ($2<=1e-9 && $3==0 || ($2>1e-9 && $3==1))print $0}' set_all2.res
P86963 8.2e-09 1
Q11101 1.9e-07 1
D3GGZ8 0.0028 1