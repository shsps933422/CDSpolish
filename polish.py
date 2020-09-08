import sys
import time
import numpy as np
import pandas as pd
import more_itertools as mit

from Bio import SeqIO

def get_nuc(label):
    if label == 0:
        return 'A'
    elif label == 1:
        return 'T'
    elif label == 2:
        return 'C'
    elif label == 3:
        return 'G'
    
start = time.time()    
result = sys.argv[1]
draft = sys.argv[2]

d_position = []
df = pd.read_feather(result)
record = SeqIO.read(draft,"fasta")
deletion = df[df['predict'] == 4].position.values #deletion
insertion = df[(df['predict'] != 5) & (df['predict'] != 4)] #insertion

d_temp = [list(group) for group in mit.consecutive_groups(deletion)]  
for key in d_temp:
    if len(key) <= 6:
        d_position = np.concatenate((d_position, key), axis=0)

seq = []
seq.append('>test\n')

for i in range(len(record)): #insertion/deletion
    if i not in d_position: #deletion -> label = 4 
        if i in insertion['position'].values: #insertion
            match = insertion[insertion['position'] == i]     
            seq.append(record[i])  
            for index, row in match.iterrows():                       
                seq.append(get_nuc(row.predict))
        else: #
            seq.append(record[i])
polished = open('polished.fasta','w')
polished.write(''.join(seq)) 

end = time.time()   
print(end-start) 
