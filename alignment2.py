from Bio import SeqIO
import numpy as np
import sys

#np.set_printoptions(threshold=np.inf)

LONG_DELETION_LENGTH = 25

def aligment(filename, Genome_size):
    """[summary]
    
    Arguments:
        filename {paf file name} -- db alignment result
        Genome_size {integer} -- TGS assembly only genome
    
    Returns:
        arr1: genome each base [0-3]: ATCG, [4]: deletion
        total: genome each base total alignment result
        Ins: genome each base [0-3]: 1st Ins, 2nd Ins, 3rd Ins; [0-3]: ATCG
    """
    arr1 = np.zeros((Genome_size,5), dtype=np.int)
    #arr1_i = np.zeros((Genome_size,5), dtype=np.int)
    total = np.zeros(Genome_size, dtype=np.int)
    Ins = np.zeros((Genome_size, 4, 4), dtype=np.int)
    #Ins_i = np.zeros((Genome_size, 4, 4), dtype=np.int)
    
    #total_gen = np.zeros((Genome_size,5), dtype=np.int)
    line_count = 0  
    over_ins = []  
    with open(filename, 'r') as f:      
        for line in f:
            line_count += 1 #paf: the number of lines
            line = line.split()      
            if line[11] != '0':#and length_percent>=85: #mapping quality != 0
                cigar = line[len(line)-1]
                cigar = cigar[5:]
                start_pos = int(line[7])                
                flag = 0
                longdel_count = 0
                longdel_status = 0
                for i in cigar: # ex. =ATC*c
                    
                    if i == 'A' or i == 'a':  # A:0, T:1, C:2, G:3
                        base = 0
                    elif i == 'T' or i == 't':
                        base = 1
                    elif i == 'C' or i == 'c':
                        base = 2
                    elif i == 'G' or i == 'g':
                        base = 3
                    
                    if i == '=': #match:1
                        flag = 1
                    elif i == '*': #mismatch:2
                        flag = 2
                        mismatch = 0
                    elif i == '+': #insertion
                        flag = 3
                        Ins_pos = 0
                        over_ins = []
                    elif i == '-': #deletion
                        flag = 4
                        longdel_status = 0
                    
                    elif flag == 1:
                        longdel_count = 0
                        arr1[start_pos][base] += 1
                        #arr1_i[start_pos][base] += weight 
                        total[start_pos] += 1
                        start_pos += 1 
                    elif flag == 2: 
                        # *gc
                        # -01
                        # 01
                        longdel_count = 0
                        if mismatch != 1:
                            mismatch += 1
                        else:
                            arr1[start_pos][base] += 1                            # 
                            total[start_pos] += 1 
                            start_pos += 1
                        
                    elif flag == 3:  #insertion
                        #+AAAAA
                        #-0123
                        #01234
                        longdel_count = 0
                        if Ins_pos < 4:
                            Ins[start_pos-1][Ins_pos][base] += 1
                            #Ins_i[start_pos-1][Ins_pos][base] += weight
                            Ins_pos += 1
                            over_ins.append(base)
                        elif Ins_pos == 4:
                            # for i in range(4):
                            #     for j in over:
                            #         Ins[start_pos-1][i][j] = 0
                        
                            # Ins_pos += 1
                            for x,y in zip(range(4), over_ins):
                                Ins[start_pos-1][x][y] -= 1
                                #Ins_i[start_pos-1][x][y] -= weight
                            over_ins = []
                            
                            
                    elif flag == 4:    #deletion  
                       
                        if longdel_status == 0:
                            longdel_count += 1
                        if longdel_count > LONG_DELETION_LENGTH and longdel_status == 0:
                            for i in range(1,LONG_DELETION_LENGTH + 1):
                                arr1[start_pos-i][4] -= 1
                                #arr1_i[start_pos-i][4] -= weight
                                total[start_pos-i] -= 1                                                  
                            longdel_status = 1
                            longdel_count = 0
                        elif longdel_status != 1:
                            arr1[start_pos][4] += 1
                            #arr1_i[start_pos][4] += weight
                            total[start_pos] += 1
                        start_pos+=1
        return arr1, total, Ins
def align(read, db_paf, truth_paf):
    record = SeqIO.read(read,"fasta")
    Genome_size = len(record)
    #db_arr, db_base, db_Ins,db_arri,db_Insi = aligment(db_paf, Genome_size)
    truth_arr, truth_base, truth_Ins = aligment(truth_paf,Genome_size)
    #print("align ok")
    #np.savez('db.npz', db_arr, db_base, db_Ins,db_arri,db_Insi)
    np.savez('truth.npz', truth_arr, truth_base, truth_Ins)
    #print("npz saved")
    return 'truth.npz'
    
align(sys.argv[1], sys.argv[2], sys.argv[3])

