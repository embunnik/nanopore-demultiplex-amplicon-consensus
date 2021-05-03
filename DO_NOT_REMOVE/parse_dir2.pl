#!/usr/bin/perl -w
use strict;

# parse_dir.pl 07NOV2019 by EMB to write shell script that will extract consensus sequence from nanopore data
# usage perl parse_dir.pl <content.txt> where content.txt contains the content of subdirectories in the folder of interest

open(IN, "<", $ARGV[0]) || die "cannot open input file: $!";
open(OUT, ">", "parse2.sh") || die "cannot open output file: $!";

my $subdir = $ARGV[1];
my $max = $ARGV[2];
$subdir =~ s/\/$//;

my @fastq;

# parse content list of data directory to extract fastq files
while(my $line = <IN>){
  chomp $line;
  if($line =~ m/\.fastq$/){
    push(@fastq, $line);
  }
}

# print canu command and copy result consensus files to separate folder
foreach my $file (@fastq){
  my @data2 = split('\.', $file);
  print OUT "canu -p ", $data2[0], " -d ./$subdir/consensus_fastq/", $data2[0], " genomeSize=", $max/1000, "k corMhapSensitivity=high corMinCoverage=0 corOutCoverage=20000 -nanopore ", "./$subdir/demultiplexed_fastq/", $file, "\n";
  print OUT "cp ./$subdir/consensus_fastq/", $data2[0], "/", $data2[0], ".contigs.fasta", " ./$subdir/consensus_seqs/", $data2[0], ".contigs.fasta", "\n";
}

close IN;
close OUT;
