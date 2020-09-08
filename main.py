import sys
#import mash
#import download as d
import os
import alignment2csv
import alignment
import alignment2

def main():    
    read = sys.argv[1]
    db_paf = sys.argv[2]
    truth_paf = sys.argv[3]
    #dup_paf = sys.argv[4]
    #gff_paf = sys.argv[5]
    db_np = alignment.align(read, db_paf, truth_paf)#,gff_paf)
    truth_np = alignment2.align(read, db_paf, truth_paf)
    alignment2csv.tocsv(read, db_np, truth_np)
    
if __name__ == "__main__":
    main()