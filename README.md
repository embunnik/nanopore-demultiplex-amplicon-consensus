# nanopore-demultiplex-amplicon-consensus

This script is a canu wrapper for generating amplicon consensus sequences using the multiplexed fastq files from Oxford Nanopore Technologies DNA sequencing (ONT-seq) as input. 

## Prerequisits
Canu and all of its requirements need to be installed. The directory containing the canu executable needs to be added to PATH. See Canu documentation for installation instructions (https://canu.readthedocs.io/en/latest/index.html).

## Setting up the correct directory structure
Multiplexed ONT-seq fastq files should be copied into a separate data directory within the directory that contains this README file. There is no need to rename any of the fastq files generated by ONT sequencing. In summary: nanopore-demultiplex-amplicon-consensus/\<data>/<files.fastq>. This repository contains a list of barcodes, supplied as the file barcodes.csv in the directory 'DO_NOT_REMOVE'. If using different barcodes than listed in this file, this file should be replaced. It is important that the location and the name of the barcodes.csv file remain the same.

## Running the nanopore-amplicon-consensus wrapper
Use the following command:
```
sh nanopore-demultiplex-amplicon-consensus.sh <name_of_run_directory> <minimum_amplicon_length> <maximum_amplicon_length>
```
This will first remove any existing TEMP or result folders that may have been generated during a previous execution of the script on the same data directory. Next, multiple fastq files will be merged and filtered for minimum and maximum amplicon length. The script will then prompt for input on how many barcodes were used during the preparation of this data set. The supplied list contains 100 barcodes. For example, if only the first 25 barcodes were used, provide '25' as input. In that case, barcodes 26 - 100 will be ignored. This will significantly speed up the execution of the demultiplexing step. The merged and filtered files will be searched for the presence of barcodes and split into separate barcode files upon finding a match. Finally, these demultiplexed fastq files will be used by canu to generate consensus sequences for each barcode. The fasta files of the resulting consensus sequences will be copied to the output folder "consensus_seqs" for convenience. Note that this output directory and directories with merged, filtered, and consensus fastq files created in the process will be removed when the script is executed again on the same data directory. 
