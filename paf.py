import os
import time

#DRAFTS=["Bacillus","ccu063","ecoli","Ecoli","Enterococcus","goku","hide","Listeria","Pseudomonas","Salmonella","sawa","Staphylococcus","sue","vgh117"]
#GENUS=["Bacillus","Proteus","Escherichia","Escherichia","Enterococcus","Proteus","Shewanella","Listeria","Pseudomonas","Salmonella","Klebsiella","Staphylococcus","Elizabethkingia","Shewanella"]
DRAFTS=["ecoli"]
GENUS=["Escherichia"]


db_path='/big6_disk/shiuanrung107/cdspolish/seqkit_cdsdb/'
truth_path='/bip7_disk/shiuanrung107/20200817/ecoli/truth/'


draft_path='/bip7_disk/shiuanrung107/20200817/ecoli/draft/'    # draft

if not os.path.exists('paf'):
    os.mkdir('paf')
if not os.path.exists('final_paf'):
    os.mkdir('final_paf')

os.chdir('paf')

for i,j in zip(DRAFTS,GENUS):
    print("======")
    start_time = time.time()
    os.system('minimap2 -cx asm20 -A2 --cs=long -t 32  {path}{draft}.fasta {path2}{genus}_CDS.fasta > {draft}_db.paf'.format(draft=i,genus=j,path=draft_path,path2=db_path))
    end_time = time.time()
    print('{draft}_db.paf time:'.format(draft=i))
    print(end_time-start_time)
    print("")
    start_time2 = time.time()
    os.system('minimap2 -cx asm20 -A2 --cs=long -t 32  {path}{draft}.fasta {path2}{draft}.fasta > ../final_paf/{draft}_truth.paf'.format(draft=i,path=draft_path,path2=truth_path))
    end_time2 = time.time()
    print('{draft}_truth.paf time:'.format(draft=i))
    print(end_time2-start_time2)
    print("======")



for i in DRAFTS:
    print("======")
    start_time = time.time()
    os.system('sh ../test_length.sh {draft}_db.paf'.format(draft=i))
    end_time = time.time()
    print('final_{draft}_db.paf time:'.format(draft=i))
    print(end_time-start_time)
    print("======")
       


    