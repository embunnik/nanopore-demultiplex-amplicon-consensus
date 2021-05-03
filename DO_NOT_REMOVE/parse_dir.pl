#!/usr/bin/perl -w
use strict;

# parse_dir.pl 07NOV2019 by EMB to write shell script that will extract consensus sequence from nanopore data
# usage perl parse_dir.pl <content.txt> where content.txt contains the content of subdirectories in the folder of interest

open(IN, "<", $ARGV[0]) || die "cannot open input file: $!";
open(OUT, ">", "parse.sh") || die "cannot open output file: $!";

my $subdir = $ARGV[1];
$subdir =~ s/\/$//;
my $min = $ARGV[2];
my $max = $ARGV[3];
#my %dirs;
#my $flag;
#my $rank = 1;
#my @merged_files; my @filtered_files;
#my $merged_file; my $filtered_file;

# parse content list of data directory to extract fastq files
print OUT "cat ";

while(my $line = <IN>){
  chomp $line;
  if($line =~ m/\.fastq$/){
    print OUT "./$subdir/$line", " ";
  }
}

print OUT "> ./$subdir/merged_fastq/merged.fastq", "\n";
print OUT "printf \"\n\tFastq files merged\n\n\"", "\n";


# print command to filter by minimum and maximum amplicon size
print OUT "perl ./DO_NOT_REMOVE/filter.pl ", "./$subdir/merged_fastq/merged.fastq", " ", $min, " ", $max, "\n";
print OUT "printf \tFastq files filtered by length\n\n";


# print command to move fastq files filtered by size to output folder
print OUT "mv ./*filtered.fastq ./$subdir/filtered_fastq/", "\n";

# print command to demultiplex filtered fastq file
print OUT "perl ./DO_NOT_REMOVE/demultiplex.pl ", "./$subdir/filtered_fastq/filtered.fastq", "\n";

# print command to remove empty demultiplex files
print OUT "find . -empty -type f -delete", "\n";

# print command to move demultiplex files to demultiplex folder
print OUT "mv ./barcode* ./$subdir/demultiplexed_fastq/", "\n";
print OUT "printf \"\tFastq files demultiplexed\n\n\"", "\n";

close IN;
close OUT;
