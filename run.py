import os
import time
import multiprocessing

start_time = time.time()

DRAFTS=["Bacillus","ccu063","ecoli","Ecoli","Enterococcus","goku","hide","Listeria","Pseudomonas","Salmonella","sawa","Staphylococcus","sue","vgh117"]

draft_path='/big6_disk/shiuanrung107/cdspolish/draft_medaka/'
paf_path='/bip7_disk/shiuanrung107/20200817/final_paf/'

#paf
os.system('python paf.py > paf.out')

#alignment
def alignment(DRAFTS):
    if not os.path.exists('{directory}'.format(directory=DRAFTS)):
        os.mkdir('{directory}'.format(directory=DRAFTS))
    os.chdir('{directory}'.format(directory=DRAFTS))
    start_time2 = time.time()
    os.system('python ../main.py {p1}{draft}.fasta {p2}{draft}_db.paf {p2}{draft}_truth.paf'.format(draft=DRAFTS,p1=draft_path,p2=paf_path))
    end_time2 = time.time()
    print("------")
    print('{draft} align ok'.format(draft=DRAFTS))
    print('time:')
    print(end_time2-start_time2)
    print("------")
    os.chdir('../')
pool = multiprocessing.Pool(14)
pool.map(alignment, DRAFTS)
pool.close()
pool.join()

#predict
for i in DRAFTS:
    os.chdir('{directory}'.format(directory=i))
    os.system('python ../predict.py alignment.feather ../RF.model > confusion.txt')
    os.chdir('../')

#polish
def polish(DRAFTS):
    os.chdir('{directory}'.format(directory=DRAFTS))
    os.system('python ../polish.py result.feather {p1}{draft}.fasta'.format(draft=DRAFTS,p1=draft_path))
    os.chdir('../')
pool = multiprocessing.Pool(14)
pool.map(polish, DRAFTS)
pool.close()
pool.join()

os.system('python accuracy.py')

end_time = time.time()
print('time:')
print(end_time-start_time)
