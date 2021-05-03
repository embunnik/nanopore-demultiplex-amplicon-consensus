#usage: sh nanopore-demultiplex-amplicon-consensus.sh <data_directory> <minimum sequence length in nt> <maximum sequence length in nt>

subdir=$(echo $1 | sed 's/\/$//')

rm -rf ./TEMP; mkdir ./TEMP
rm -rf ./$subdir/merged_fastq; mkdir ./$subdir/merged_fastq
rm -rf ./$subdir/filtered_fastq; mkdir ./$subdir/filtered_fastq
rm -rf ./$subdir/demultiplexed_fastq; mkdir ./$subdir/demultiplexed_fastq
rm -rf ./$subdir/consensus_fastq; mkdir ./$subdir/consensus_fastq
rm -rf ./$subdir/consensus_seqs; mkdir ./$subdir/consensus_seqs
ls -R ./$subdir > ./TEMP/content.txt

perl ./DO_NOT_REMOVE/parse_dir.pl ./TEMP/content.txt $1 $2 $3

mv parse.sh ./TEMP/parse.sh

source ./TEMP/parse.sh

ls -R ./$subdir/demultiplexed_fastq > ./TEMP/content2.txt

perl ./DO_NOT_REMOVE/parse_dir2.pl ./TEMP/content2.txt $1 $3

mv parse2.sh ./TEMP/parse2.sh

source ./TEMP/parse2.sh

rm -rf ./TEMP
