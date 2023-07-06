######################################################################################
# CLBtope is developed for predicting, desigining and scanning the B-cell epitopes.  #
# It is developed by Prof G. P. S. Raghava's group.       #
# Please cite: https://webs.iiitd.edu.in/raghava/clbtope/                            #
######################################################################################
import argparse
import warnings
import subprocess
import pkg_resources
import os
import sys
import numpy as np
import pandas as pd
import math
import itertools
from collections import Counter
import pickle
import re
import glob
import time
import uuid
from time import sleep
from tqdm import tqdm
from sklearn.ensemble import RandomForestClassifier
import urllib.request
import zipfile
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser(description='Please provide following arguments') 

## Read Arguments from command
parser.add_argument("-i", "--input", type=str, required=True, help="Input: protein or peptide sequence(s) in FASTA format or single sequence per line in single letter code")
parser.add_argument("-o", "--output",type=str, help="Output: File for saving results by default outfile.csv")
parser.add_argument("-j", "--job",type=int, choices = [1,2,3], help="Job Type: 1:Predict, 2: Design, 3:Scan, by default 1")
parser.add_argument("-t","--threshold", type=float, help="Threshold: Value between 0 to 1 by default 0.53")
parser.add_argument("-w","--winleng", type=int, choices =range(8,31), help="Window Length: 8 to 30 (scan mode only), by default 9")
parser.add_argument("-d","--display", type=int, choices = [1,2], help="Display: 1:B-cell epitope, 2: All peptides, by default 1")
args = parser.parse_args()

# Function for generating all possible mutants
def mutants(file1,file2):
    std = list("ACDEFGHIKLMNPQRSTVWY")
    cc = []
    dd = []
    ee = []
    df2 = file2
    df2.columns = ['Name']
    df1 = file1
    df1.columns = ['Seq']
    for k in range(len(df1)):
        cc.append(df1['Seq'][k])
        dd.append('Original_'+'Seq'+str(k+1))
        ee.append(df2['Name'][k])
        for i in range(0,len(df1['Seq'][k])):
            for j in std:
                if df1['Seq'][k][i]!=j:
                    dd.append('Mutant_'+df1['Seq'][k][i]+str(i+1)+j+'_Seq'+str(k+1))
                    cc.append(df1['Seq'][k][:i] + j + df1['Seq'][k][i + 1:])
                    ee.append(df2['Name'][k])
    xx = pd.concat([pd.DataFrame(ee),pd.DataFrame(dd),pd.DataFrame(cc)],axis=1)
    xx.columns = ['Seq_ID','Mutant_ID','Seq']
    return xx
# Function for generating pattern of a given length
def seq_pattern(file1,file2,num):
    df1 = file1
    df1.columns = ['Seq']
    df2 = file2
    df2.columns = ['Name']
    cc = []
    dd = []
    ee = []
    for i in range(len(df1)):
        for j in range(len(df1['Seq'][i])):
            xx = df1['Seq'][i][j:j+num]
            if len(xx) == num:
                cc.append(df2['Name'][i])
                dd.append('Pattern_'+str(j+1)+'_Seq'+str(i+1))
                ee.append(xx)
    df3 = pd.concat([pd.DataFrame(cc),pd.DataFrame(dd),pd.DataFrame(ee)],axis=1)
    df3.columns= ['Seq_ID','Pattern_ID','Seq']
    return df3
# Function to check the seqeunce
def readseq(file):
    with open(file) as f:
        records = f.read()
    records = records.split('>')[1:]
    seqid = []
    seq = []
    for fasta in records:
        array = fasta.split('\n')
        name, sequence = array[0].split()[0], re.sub('[^ACDEFGHIKLMNPQRSTVWY-]', '', ''.join(array[1:]).upper())
        seqid.append('>'+name)
        seq.append(sequence)
    if len(seqid) == 0:
        f=open(file,"r")
        data1 = f.readlines()
        for each in data1:
            seq.append(each.replace('\n',''))
        for i in range (1,len(seq)+1):
            seqid.append(">Seq_"+str(i))
    df1 = pd.DataFrame(seqid)
    df2 = pd.DataFrame(seq)
    return df1,df2
# Function to check the length of seqeunces
def lenchk(file1):
    cc = []
    df1 = file1
    df1.columns = ['seq']
    for i in range(len(df1)):
        if len(df1['seq'][i])>30:
            cc.append(df1['seq'][i][0:30])
        else:
            cc.append(df1['seq'][i])
    df2 = pd.DataFrame(cc)
    df2.columns = ['Seq']
    return df2
# Function to generate the features out of seqeunces
def feature_gen(file,q=1):
    std = list("ACDEFGHIKLMNPQRSTVWY")
    df1 = file
    df1.columns = ['Seq']
    zz = df1.Seq
    dd = []
    for i in range(0,len(zz)):
        cc = []
        for j in std:
            for k in std:
                count = 0
                temp = j+k
                for m3 in range(0,len(zz[i])-q):
                    b = zz[i][m3:m3+q+1:q]
                    b.upper()
                    if b == temp:
                        count += 1
                    composition = (count/(len(zz[i])-(q)))*100
                cc.append(composition)
        dd.append(cc)
    df2 = pd.DataFrame(dd)
    head = []
    for s in std:
        for u in std:
            head.append("DPC"+str(q)+"_"+s+u)
    df2.columns = head
    return df2
# Function to process the blast output
def BLAST_processor(blast_result,name1,ml_results,thresh):
    if os.stat(blast_result).st_size != 0:
        df1 = pd.read_csv(blast_result, sep="\t", names=['name','hit','identity','r1','r2','r3','r4','r5','r6','r7','r8','r9'])
        df__2 = name1
        df2 = pd.DataFrame()
        df2 = df2.append(df__2.values.tolist())
        df3 = ml_results
        cc = []
        for i in df2[0]:
            kk = i.replace('>','')
            if len(df1.loc[df1.name==kk])>0:
                df4 = df1[['name','hit']].loc[df1['name']==kk].reset_index(drop=True)
                if df4['hit'][0].split('_')[0]=='P':
                    cc.append(0.5)
                if df4['hit'][0].split('_')[0]=='N':
                    cc.append(-0.5)
            else:
                cc.append(0)
        df6 = pd.DataFrame()
        df6['Seq_ID'] = [i.replace('>','') for i in df2[0]]
        df6['ML_Score'] = df3['ML_score']
        df6['BLAST_Score'] = cc
        df6['Total_Score'] = df6['ML_Score']+df6['BLAST_Score']
        df6['Prediction'] = ['B-cell epitope' if df6['Total_Score'][i]>thresh else 'Non B-cell epitope ' for i in range(0,len(df6))]
    else:
        df__2 = name1
        df3 = ml_results
        df2 = pd.DataFrame()
        df2 = df2.append(df__2.values.tolist())
        ss = []
        vv = []
        for j in df2[0]:
            ss.append(j.replace('>',''))
            vv.append(0)
        df6 = pd.DataFrame()
        df6['Seq_ID'] = ss
        df6['ML_Score'] = df3['ML_score']
        df6['BLAST_Score'] = vv
        df6['Total_Score'] = df6['ML_Score']+df6['BLAST_Score']
        df6['Prediction'] = ['B-cell epitope' if df6['Total_Score'][i]>thresh else 'Non B-cell epitope' for i in range(0,len(df6))]
    return df6
# Function to read and implement the model
def model_run(file1,file2):
    a = []
    data_test = file1
    clf = pickle.load(open(file2,'rb'))
    y_p_score1=clf.predict_proba(data_test)
    y_p_s1=y_p_score1.tolist()
    a.extend(y_p_s1)
    df = pd.DataFrame(a)
    df1 = df.iloc[:,-1].round(2)
    df2 = pd.DataFrame(df1)
    df2.columns = ['ML_score']
    return df2
    
('############################################################################################')
print('# This program CLBTope is developed for predicting, desigining and scanning B-cell epitopes #')
print('# developed by Prof G. P. S. Raghava group.               #')
print('# Please cite: CLBTope; available at https://webs.iiitd.edu.in/raghava/clbtope/  #')
print('############################################################################################')

# Parameter initialization or assigning variable for command level arguments

Sequence= args.input        # Input variable 
 
# Output file 
 
if args.output == None:
    result_filename= "outfile.csv" 
else:
    result_filename = args.output
         
# Threshold 
if args.threshold == None:
        Threshold = 0.53
else:
        Threshold= float(args.threshold)
# Job Type 
if args.job == None:
        Job = int(1)
else:
        Job = int(args.job)
# Window Length 
if args.winleng == None:
        Win_len = int(9)
else:
        Win_len = int(args.winleng)

# Display
if args.display == None:
        dplay = int(1)
else:
        dplay = int(args.display)

###################################Model Calling##########################################
filepath = os.path.dirname(os.path.abspath(__file__))
if os.path.isdir('model') == False:
    url = "https://webs.iiitd.edu.in/raghava/clbtope/model.zip"
    file_name_dn = filepath+"/model.zip"
    urllib.request.urlretrieve(url, file_name_dn)
    with zipfile.ZipFile(filepath+'/model.zip', 'r') as zip_ref:
        zip_ref.extractall('.')
else:
    pass
#####################################BLAST Path############################################
if os.path.exists('envfile'):
    with open('envfile', 'r') as file:
        data = file.readlines()
    output = []
    for line in data:
        if not "#" in line:
            output.append(line)
    if len(output)==2:
        paths = []
        for i in range (0,len(output)):
            paths.append(output[i].split(':')[1].replace('\n',''))
        blastp = paths[0]
        blastdb = paths[1]
    else:
        print("####################################################################################")
        print("Error: Please provide paths for BLAST, and required files", file=sys.stderr)
        print("####################################################################################")
        sys.exit()

else:
    print("####################################################################################")
    print("Error: Please provide the '{}', which comprises paths for BLAST".format('envfile'), file=sys.stderr)
    print("####################################################################################")
    sys.exit()
###########################################################################################

if Job==2:
    print("\n");
    print('##############################################################################')
    print('Summary of Parameters:')
    print('Input File: ',Sequence,'; Threshold: ', Threshold,'; Job Type: ',Job)
    print('Output File: ',result_filename,'; Window Length: ',Win_len,'; Display: ',dplay)
    print('##############################################################################')
else:
    print("\n");
    print('##############################################################################')
    print('Summary of Parameters:')
    print('Input File: ',Sequence,'; Threshold: ', Threshold,'; Job Type: ',Job)
    print('Output File: ',result_filename,'; Display: ',dplay)
    print('# ############################################################################')
#======================= Prediction Module start from here =====================
if Job == 1:
    print('\n======= Thanks for using Predict module of CLBTope. Your results will be stored in file :',result_filename,' =====\n')
    df_2,dfseq = readseq(Sequence)
    df1 = lenchk(dfseq)
    X = feature_gen(df1)
    mlres = model_run(X,'model/clbtope_dpc400model.sav')
    filename = str(uuid.uuid4())
    df11 = pd.concat([df_2,df1],axis=1)
    df11.to_csv(filename,index=None,header=False,sep="\n")
    mlres = mlres.round(3)
    os.system(blastp + " -task blastp-short -db " + blastdb + " -query " + filename  + " -out RES_1_6_6.out -outfmt 6 -evalue 1")
    df44 = BLAST_processor('RES_1_6_6.out',df_2,mlres,Threshold)
    df44['Sequence'] = df1.Seq
    df44 = df44[['Seq_ID','Sequence','ML_Score','BLAST_Score','Total_Score','Prediction']]
    if dplay == 1:
        df44 = df44.loc[df44.Prediction=="B-cell epitope"]
    else:
        df44 = df44
    df44 = round(df44,3)
    df44.to_csv(result_filename, index=None)
    os.remove('RES_1_6_6.out')
    os.remove(filename)
    print("\n=========Process Completed. Have an awesome day ahead.=============\n")    
#===================== Design Model Start from Here ======================
elif Job == 2:
    print('\n======= Thanks for using Design module of CLBTope. Your results will be stored in file :',result_filename,' =====\n')
    print('==== Designing Peptides: Processing sequences please wait ...')
    df_2,dfseq = readseq(Sequence)
    df1 = lenchk(dfseq)
    df_1 = mutants(df1,df_2)
    dfseq = df_1[['Seq']]
    X = feature_gen(dfseq)
    mlres = model_run(X,'model/clbtope_dpc400model.sav')
    filename = str(uuid.uuid4())
    df_1['Mutant'] = ['>'+df_1['Mutant_ID'][i] for i in range(len(df_1))]
    df11 = df_1[['Mutant','Seq']] 
    df11.to_csv(filename,index=None,header=False,sep="\n")
    mlres = mlres.round(3)
    os.system(blastp + " -task blastp-short -db " + blastdb + " -query " + filename  + " -out RES_1_6_6.out -outfmt 6 -evalue 1")
    df44 = BLAST_processor('RES_1_6_6.out',df11[['Mutant']],mlres,Threshold)
    df44['Mutant_ID'] = ['_'.join(df44['Seq_ID'][i].split('_')[:-1]) for i in range(len(df44))]
    df44.drop(columns=['Seq_ID'],inplace=True)
    df44['Seq_ID'] = [i.replace('>','') for i in df_1['Seq_ID']]
    df44['Sequence'] = df_1.Seq
    df44 = df44[['Seq_ID','Mutant_ID','Sequence','ML_Score','BLAST_Score','Total_Score','Prediction']]
    if dplay == 1:
        df44 = df44.loc[df44.Prediction=="B-cell epitope"]
    else:
        df44 = df44
    df44 = round(df44,3)
    df44.to_csv(result_filename, index=None)
    os.remove('RES_1_6_6.out')
    os.remove(filename)
    print("\n=========Process Completed. Have an awesome day ahead.=============\n")
#=============== Scan Model start from here ==================
elif Job==3:
    print('\n======= Thanks for using Scan module of CLBTope. Your results will be stored in file :',result_filename,' =====\n')
    print('==== Scanning Peptides: Processing sequences please wait ...')
    df_2,dfseq = readseq(Sequence)
    df_1 = seq_pattern(dfseq,df_2,Win_len)
    dfseq = df_1[['Seq']]
    X = feature_gen(dfseq)
    mlres = model_run(X,'model/clbtope_dpc400model.sav')
    filename = str(uuid.uuid4())
    df_1['Pattern'] = ['>'+df_1['Pattern_ID'][i] for i in range(len(df_1))]
    df11 = df_1[['Pattern','Seq']]
    df11.to_csv(filename,index=None,header=False,sep="\n")
    mlres = mlres.round(3)
    os.system(blastp + " -task blastp-short -db " + blastdb + " -query " + filename  + " -out RES_1_6_6.out -outfmt 6 -evalue 1")
    df44 = BLAST_processor('RES_1_6_6.out',df11[['Pattern']],mlres,Threshold)
    df44['Pattern_ID'] = ['_'.join(df44['Seq_ID'][i].split('_')[:-1]) for i in range(len(df44))]
    df44.drop(columns=['Seq_ID'],inplace=True)
    df44['Seq_ID'] = [i.replace('>','') for i in df_1['Seq_ID']]
    df44['Sequence'] = df_1.Seq
    df44 = df44[['Seq_ID','Pattern_ID','Sequence','ML_Score','BLAST_Score','Total_Score','Prediction']]
    if dplay == 1:
        df44 = df44.loc[df44.Prediction=="B-cell epitope"]
    else:
        df44 = df44
    df44 = round(df44,3)
    df44.to_csv(result_filename, index=None)
    os.remove('RES_1_6_6.out')
    os.remove(filename)
    print("\n=========Process Completed. Have an awesome day ahead.=============\n")
print('\n======= Thanks for using CLBTope. Your results are stored in file :',result_filename,' =====\n\n')
