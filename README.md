# M.S. Thesis
# Assembly polishing for Nanopore sequencing using random forest with genus-specific coding sequences 

# Abstract
Nanopore sequencing has the advantages of fast sequencing, real-time determination of sequences and ultra-long reads, etc. Ultra-long reads is essential for reconstructing a complete
genome. However, the assembled genome requires further quality improvement as a result of the high error rate of Nanopore sequencing. Existing methods correct errors by training models
from reads or signals; nevertheless, they are unable to correct Nanopore systemic errors dominated by insertions and deletions. We observed coding sequences are highly conservation, where
insertions or deletions are rarely occurred. This thesis aims to train a random forest model using genus-specific coding sequences to correct Nanopore systemic errors. The experimental results
indicated that polishing by coding sequences can achieve higher accuracy in comparison with Racon and Medaka. We showed that this method is particularly suitable for removing errors 
within mobile genetic elements (e.g., insertion sequences).