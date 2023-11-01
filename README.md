# bulkRNAseq
Script repo for bulk RNA seq analyses

Pipeline: 
STAR > featureCounts > filtercounts > DESeq2 OR edgeR
  -DESeq2 & edgeR do the same thing. Try DESeq2 first, then edge R 
  -EdgeR is less stringent
