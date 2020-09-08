import pandas as pd
import numpy as np

def preprocessing(data):
    new_draft = []
    draft = data['draft']
    for i in draft:
        if i == 'A':
            new_draft.append(0)
        elif i == 'T':
            new_draft.append(1)
        elif i == 'C':
            new_draft.append(2)
        elif i == 'G':
            new_draft.append(3)
        elif i == '-':
            new_draft.append(4)

    data['draft'] = new_draft
    '''
    data['A'] = data['A'].astype(int) / data['coverage'].astype(int)
    data['T'] = data['T'].astype(int) / data['coverage'].astype(int)
    data['C'] = data['C'].astype(int) / data['coverage'].astype(int)
    data['G'] = data['G'].astype(int) / data['coverage'].astype(int)
    data['gap'] = data['gap'].astype(int) / data['coverage'].astype(int)
    data['Ins_A'] = data['Ins_A'].astype(int) / data['coverage'].astype(int)
    data['Ins_T'] = data['Ins_T'].astype(int) / data['coverage'].astype(int)
    data['Ins_C'] = data['Ins_C'].astype(int) / data['coverage'].astype(int)
    data['Ins_G'] = data['Ins_G'].astype(int) / data['coverage'].astype(int)
    '''
    #data.loc[data['homopolymer'] > 10, 'homopolymer'] = 10   
    return data
    
def double_ins_chaos(data):
    chaos = np.zeros((data.shape[0],), dtype=int)
    ins_index = data[data.draft == '-'].index.values
    side = 10
    dis = np.diff(ins_index)
    idxs = np.where((dis < 10) & (dis > 1))[0]
    

    for i in idxs:        
        head = data[data.index == ins_index[i]]
        head_total = head[['Ins_A','Ins_T','Ins_C','Ins_G']].sum(axis=1).values[0]
        tail = data[data.index == ins_index[i+1]]
        tail_total = tail[['Ins_A','Ins_T','Ins_C','Ins_G']].sum(axis=1).values[0]
        coverage = head.coverage.values[0]
        if (head_total + tail_total) == coverage:   
            if head_total > tail_total:
                chaos[ins_index[i]] = 1
                chaos[ins_index[i+1]] = 2
            elif head_total < tail_total:
                chaos[ins_index[i]] = 2
                chaos[ins_index[i+1]] = 1
            else:
                chaos[ins_index[i]] = 1
                chaos[ins_index[i+1]] = 1
    data['chaos'] = chaos
    return data