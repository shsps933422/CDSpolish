from Bio import SeqIO
import numpy as np
import pandas as pd
import sys
from collections import Counter

np.seterr(divide='ignore',invalid='ignore')

def tocsv(read, db_np, truth_np):
    record = SeqIO.read(read,"fasta")
    genome_size = len(record)
    db_arr = np.load(db_np)  
    truth_arr = np.load(truth_np)

    db_arr, db_base, db_Ins = db_arr['arr_0'], db_arr['arr_1'], db_arr['arr_2']
    truth_arr, truth_base, truth_Ins = truth_arr['arr_0'], truth_arr['arr_1'], truth_arr['arr_2']

    homo_arr = np.zeros(genome_size, dtype=np.int)
    
    position = []
    draft = []
    A, T, C , G = [], [], [], []
    gap = []
    ins_A, ins_T, ins_C, ins_G = [], [], [], []
    homo = []
    coverage = []
    label = []
    #deletion_length=[]
    #insertion_length=[]
    indel_length=[]
    
    
    for i in range(genome_size):
        count = 0
        count2=1
        #count3=0
        if truth_base[i] != 0 and db_base[i] != 0:
            position.append(i)
            draft.append(record[i])
            A.append(db_arr[i][0])
            T.append(db_arr[i][1])
            C.append(db_arr[i][2])
            G.append(db_arr[i][3])
            gap.append(db_arr[i][4])
            ins_A.append(0)
            ins_T.append(0)
            ins_C.append(0)
            ins_G.append(0)
            coverage.append(db_base[i])
            

            index = i
            nuc = record[i]
            while i + 1 < genome_size and record[i+1] == nuc:                               
                count += 1
                i += 1                
            i = index
            
            while record[i] and record[i] == nuc:
                count += 1
                i -= 1
            i = index
            
            homo_arr[i] = count
            
            if db_arr[i][4]==0:
                indel_length.append(0)
                
            else:
                index2=i
                while db_arr[i][4]/db_base[i]>0.5 and db_arr[i+1][4]/db_base[i+1]>0.5:
                    count2 += 1
                    i += 1
                i=index2
            
                while db_arr[i][4]/db_base[i]>0.5 and db_arr[i-1][4]/db_base[i-1]>0.5:
                    count2 += 1
                    i -= 1
                i = index2
                indel_length.append(count2)
                
            
            
                
            if truth_arr[i][4] > 0:
                label.append(4)
            else:
                label.append(5)
 
            for j in range(4):
                if db_Ins[i][j][0] != 0 or db_Ins[i][j][1] != 0 or db_Ins[i][j][2] != 0 or db_Ins[i][j][3] != 0:
                    position.append(i)
                    draft.append('-')
                    A.append(0)
                    T.append(0)
                    C.append(0)
                    G.append(0)
                    gap.append(0)
       
                    
                    if j==0:
                        count3=1
                        if (db_Ins[i][j+1][0]/db_base[i]>0.5 or db_Ins[i][j+1][1]/db_base[i]>0.5 or db_Ins[i][j+1][2]/db_base[i]>0.5 or db_Ins[i][j+1][3]/db_base[i]>0.5) :
                            count3+=1
                            if (db_Ins[i][j+2][0]/db_base[i]>0.5 or db_Ins[i][j+2][1]/db_base[i]>0.50 or db_Ins[i][j+2][2]/db_base[i]>0.5 or db_Ins[i][j+2][3]/db_base[i]>0.5):
                                count3+=1
                                if (db_Ins[i][j+3][0]/db_base[i]>0.5 or db_Ins[i][j+3][1]/db_base[i]>0.5 or db_Ins[i][j+3][2]/db_base[i]>0.5 or db_Ins[i][j+3][3]/db_base[i]>0.5):
                                    count3+=1
                                    indel_length.append(count3)
                                    
                                else:
                                    indel_length.append(count3)
                                    
                            else:
                                indel_length.append(count3)
                                
                                
                        else:
                            indel_length.append(count3)
                            
                    elif j==1:
                        count3=2
                        if (db_Ins[i][j+1][0]/db_base[i]>0.5 or db_Ins[i][j+1][1]/db_base[i]>0.5 or db_Ins[i][j+1][2]/db_base[i]>0.5 or db_Ins[i][j+1][3]/db_base[i]>0.5) :
                            count3+=1
                            if (db_Ins[i][j+2][0]/db_base[i]>0.5 or db_Ins[i][j+2][1]/db_base[i]>0.5 or db_Ins[i][j+2][2]/db_base[i]>0.5 or db_Ins[i][j+2][3]/db_base[i]>0.5):
                                count3+=1
                                indel_length.append(count3)
                                
                            else:
                                indel_length.append(count3)
                                
                                
                        else:
                            indel_length.append(count3)
                            
                    elif j==2:
                        count3=3
                        if (db_Ins[i][j+1][0]/db_base[i]>0.5 or db_Ins[i][j+1][1]/db_base[i]>0.5 or db_Ins[i][j+1][2]/db_base[i]>0.5 or db_Ins[i][j+1][3]/db_base[i]>0.5) :
                            count3+=1
                            indel_length.append(count3)
                            
                        else:
                            indel_length.append(count3)
                            
                    elif j==3:
                        count3=4
                        indel_length.append(count3)
                        
                        
                        
                    
                    
                    ins_A.append(db_Ins[i][j][0])
                    ins_T.append(db_Ins[i][j][1])
                    ins_C.append(db_Ins[i][j][2])
                    ins_G.append(db_Ins[i][j][3])
                    coverage.append(db_base[i])
                    
                   
                    if truth_Ins[i][j][0] == 0 and truth_Ins[i][j][1] == 0 and truth_Ins[i][j][2] == 0 and truth_Ins[i][j][3] == 0 :
                        label.append(5)              
                    else:
                        for k in range(4):
                            if truth_Ins[i][j][k] != 0:
                                label.append(k)
                                break
    for i in range(genome_size):
        if truth_base[i] != 0 and db_base[i] != 0:
            homo.append(homo_arr[i])
            for j in range(4):
                if db_Ins[i][j][0] != 0 or db_Ins[i][j][1] != 0 or db_Ins[i][j][2] != 0 or db_Ins[i][j][3] != 0:                
                    if max(db_Ins[i][j][0],db_Ins[i][j][1],db_Ins[i][j][2],db_Ins[i][j][3]) == db_Ins[i][j][0] and record[i+1]=='A' :
                        homo.append(homo_arr[i+1])
                    elif max(db_Ins[i][j][0],db_Ins[i][j][1],db_Ins[i][j][2],db_Ins[i][j][3]) == db_Ins[i][j][1] and record[i+1]=='T' :
                        homo.append(homo_arr[i+1])
                    elif max(db_Ins[i][j][0],db_Ins[i][j][1],db_Ins[i][j][2],db_Ins[i][j][3]) == db_Ins[i][j][2] and record[i+1]=='C' :
                        homo.append(homo_arr[i+1])
                    elif max(db_Ins[i][j][0],db_Ins[i][j][1],db_Ins[i][j][2],db_Ins[i][j][3]) == db_Ins[i][j][3] and record[i+1]=='G' :
                        homo.append(homo_arr[i+1])
                    else:
                        homo.append(0)     
                        
    dict = {"position": position,
            "draft": draft,
            "A": A,
            "T": T,
            "C": C,
            "G": G,
            "gap": gap,
            "Ins_A": ins_A,
            "Ins_T": ins_T,
            "Ins_C": ins_C,
            "Ins_G": ins_G,
            "coverage": coverage,
            "homopolymer": homo,
            #"deletion_length":deletion_length,
            #"insertion_length":insertion_length,
            "indel_length":indel_length,
            "label": label
        }
    # print(len(position))
    # print(len(A))
    # print(len(gap))
    # print(len(ins_C))
    #print(len(homo))
    #print(len(deletion_length))
    #print(len(insertion_length))
    #print(len(label))
    df = pd.DataFrame(dict)
    # df.to_csv("alignment.csv",index=False)
    #df.to_hdf('alignment.h5', key='df', mode='w')
    df.to_feather('alignment.feather')
#tocsv(sys.argv[1], sys.argv[2], sys.argv[3])