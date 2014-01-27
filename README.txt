##
## Pipeline for analysing doing viral metagenomics ##
##

Pipeline consists of 4 major steps: 
1. pre-assembly screening against Human Genomic and Human Genomic Transcripts
2. assembly with MIRA
3. search pipline against NT/NR
4. post search analysis


Requirements:
1. sffinfo is expected to be installed and to be in path
2. blastall/megablast is expected to be installed and in path
3. RepeatMasker is needed and should be placed in bin/repeatmasker
4. 454 files should be placed in data/454
5. sanger files should be placed in data/sanger and be NCBI tracedb formatted (see MIRA manual for more information on how to format Sanger input)


Example run (expects you to be in the pipeline folder):

# download all needed databases
./db/blastdb/download-all
./db/taxdb/download

# run pipeline on data placed in data folder
./pipeline test data/454/test.sff

# results can be found in result/test
