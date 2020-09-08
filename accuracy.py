import os
import multiprocessing

TEST=["Bacillus","ccu063","ecoli","hide","Listeria","Salmonella","sue"]
TRAIN=["Ecoli","Enterococcus","goku","Pseudomonas","sawa","Staphylococcus","vgh117"]

fastmer = '/big6_disk/shiuanrung107/assembly_accuracy/fastmer.py'
truth='/big6_disk/shiuanrung107/cdspolish/truth/'

os.system('touch accuracy_test.txt')
os.system('touch accuracy_train.txt')
    
def qscore(TEST):
    os.chdir('{directory}'.format(directory=TEST))
    os.system('mv polished.fasta {TEST}_polished_RF.fasta'.format(TEST=TEST))
    os.system('python {py} --reference {truth}{TEST}.fasta --assembly {TEST}_polished_RF.fasta >> ../accuracy_test.txt'.format(py=fastmer,truth=truth,TEST=TEST))
    os.chdir('../')
pool = multiprocessing.Pool(7)
pool.map(qscore, TEST)
pool.close()
pool.join()

def qscore(TRAIN):
    os.chdir('{directory}'.format(directory=TRAIN))
    os.system('mv polished.fasta {TRAIN}_polished_RF.fasta'.format(TRAIN=TRAIN))
    os.system('python {py} --reference {truth}{TRAIN}.fasta --assembly {TRAIN}_polished_RF.fasta >> ../accuracy_train.txt'.format(py=fastmer,truth=truth,TRAIN=TRAIN))
    os.chdir('../')
pool = multiprocessing.Pool(7)
pool.map(qscore, TRAIN)
pool.close()
pool.join()

