from Bio import SeqIO
import numpy as np
import sys
import math

LONG_DELETION_LENGTH = 25
DB_COVERAGE=10

def aligment(filename, Genome_size,genome_id):
    """[summary]
    
    Arguments:
        filename {paf file name} -- db alignment result
        Genome_size {integer} -- TGS assembly only genome
    
    Returns:
        arr1: genome each base [0-3]: ATCG, [4]: deletion
        total: genome each base total alignment result
        Ins: genome each base [0-3]: 1st Ins, 2nd Ins, 3rd Ins; [0-3]: ATCG
    """
    arr1 = np.zeros((Genome_size,5), dtype=np.float)
    total = np.zeros(Genome_size, dtype=np.int)
    #total2 = np.zeros(Genome_size, dtype=np.int)
    Ins = np.zeros((Genome_size, 4, 4), dtype=np.float)
    
    total_weight = np.zeros(Genome_size, dtype=np.float)
    
    #total_gen = np.zeros((Genome_size,5), dtype=np.int)
    line_count = 0  
    over_ins = []  
    db_dict={}
    
    with open(filename, 'r') as f:      
        for line in f:
            line_count += 1 #paf: the number of lines
            line = line.split()      
            #mapping quality == 60 , tp:A:P, 
            db_id=line[0]
            identity=float(line[13])
            query_cover=float(line[6])
            query_length=int(line[3])
            cigar = line[14]
            cigar = cigar[5:]
            start_pos = int(line[9])
            end_pos = int(line[10])
            middle_pos=int((end_pos+start_pos)/2)
            #if query_length%3 != 0 or query_cover!=100:
            #    continue  
            #if identity<95:
             #   continue  
            if total[middle_pos]>=DB_COVERAGE:
                continue 
            flag = 0
            longdel_count = 0
            longdel_status = 0
            rank=0
            for i in cigar: # ex. =ATC*c
                db_dict.setdefault(start_pos,{})[db_id]=identity
                
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
                    rank = list(db_dict[start_pos].keys()).index(db_id)
                    if rank <=2 :                       
                        arr1[start_pos][base] += DB_COVERAGE/math.pow(2,rank) 
                        total_weight[start_pos] += DB_COVERAGE/math.pow(2,rank) 
                    else:
                        arr1[start_pos][base] += DB_COVERAGE/math.pow(2,3) 
                        total_weight[start_pos] += DB_COVERAGE/math.pow(2,3)              
                    total[start_pos] += 1
                    start_pos += 1 
                elif flag == 2: 
                    longdel_count = 0
                    if mismatch != 1:
                        mismatch += 1
                    else:
                        rank = list(db_dict[start_pos].keys()).index(db_id)
                        if rank <=2 :  
                            arr1[start_pos][base] +=  DB_COVERAGE/math.pow(2,rank) 
                            total_weight[start_pos] += DB_COVERAGE/math.pow(2,rank) 
                        else:
                            arr1[start_pos][base] +=  DB_COVERAGE/math.pow(2,3)
                            total_weight[start_pos] += DB_COVERAGE/math.pow(2,3) 
                        total[start_pos] += 1 
                        start_pos += 1
                        
                elif flag == 3:  #insertion
                    #+AAAAA
                    #-0123
                    #01234
                    longdel_count = 0
                    if Ins_pos < 4:
                        rank = list(db_dict[start_pos-1].keys()).index(db_id)
                        if rank <=2 : 
                            Ins[start_pos-1][Ins_pos][base] += DB_COVERAGE/math.pow(2,rank)
                            #total_weight[start_pos-1] += DB_COVERAGE/math.pow(2,rank) 
                        else:
                            Ins[start_pos-1][Ins_pos][base] += DB_COVERAGE/math.pow(2,3)
                            #total_weight[start_pos-1] += DB_COVERAGE/math.pow(2,3) 
                        Ins_pos += 1
                        over_ins.append(base)
                    elif Ins_pos == 4:
                        # for i in range(4):
                        #     for j in over:
                        #         Ins[start_pos-1][i][j] = 0
                        
                        # Ins_pos += 1
                        for x,y in zip(range(4), over_ins):
                            rank = list(db_dict[start_pos-1].keys()).index(db_id)
                            if rank <=2 : 
                                Ins[start_pos-1][x][y] -= DB_COVERAGE/math.pow(2,rank) 
                                #total_weight[start_pos-1] -= DB_COVERAGE/math.pow(2,rank) 
                            else:
                                Ins[start_pos-1][x][y] -= DB_COVERAGE/math.pow(2,3)
                                #total_weight[start_pos-1] -= DB_COVERAGE/math.pow(2,3) 
                        over_ins = []
                            
                            
                elif flag == 4:    #deletion  
                       
                    if longdel_status == 0:
                        longdel_count += 1
                    if longdel_count > LONG_DELETION_LENGTH and longdel_status == 0:
                        for i in range(1,LONG_DELETION_LENGTH + 1):
                            rank = list(db_dict[start_pos-i].keys()).index(db_id)
                            if rank <=2 : 
                                arr1[start_pos-i][4] -= DB_COVERAGE/math.pow(2,rank)
                                total_weight[start_pos-i] -= DB_COVERAGE/math.pow(2,rank) 
                            else:
                                arr1[start_pos-i][4] -= DB_COVERAGE/math.pow(2,3) 
                                total_weight[start_pos-i] -= DB_COVERAGE/math.pow(2,3)                             
                            total[start_pos-i] -= 1                                                  
                        longdel_status = 1
                        longdel_count = 0
                    elif longdel_status != 1:
                        rank = list(db_dict[start_pos].keys()).index(db_id)
                        if rank <=2 : 
                            arr1[start_pos][4] += DB_COVERAGE/math.pow(2,rank) 
                            total_weight[start_pos] += DB_COVERAGE/math.pow(2,rank)   
                        else:
                            arr1[start_pos][4] += DB_COVERAGE/math.pow(2,3) 
                            total_weight[start_pos] += DB_COVERAGE/math.pow(2,3)   
                                           
                        total[start_pos] += 1
                    start_pos+=1
        fp=open("identity.txt",'w')   
        for k,v in sorted(db_dict.items()):
	          fp.write(str(k)+' '+str(v)+'\n')
        #fp.write(db_dict)
        fp.close()
        return arr1, total_weight, Ins#, cds
    
def align(read, db_paf, truth_paf):
    record = SeqIO.read(read,"fasta")
    Genome_size = len(record)
    genome_id=record.id
    #db_arr, db_base, db_Ins,db_cds = aligment(db_paf, Genome_size,gff,genome_id)
    db_arr, db_base, db_Ins = aligment(db_paf, Genome_size,genome_id)
    
    #np.savez('db.npz', db_arr, db_base, db_Ins, db_cds)
    np.savez('db.npz', db_arr, db_base, db_Ins)
    
    return 'db.npz'
    
align(sys.argv[1], sys.argv[2], sys.argv[3])

