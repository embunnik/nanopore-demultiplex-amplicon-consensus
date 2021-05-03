#!/usr/bin/perl -w
use strict;

# filter.pl by EMB 03NOV2019 to remove sequences that are too long or too short from fastq sequence file
# usage: perl filter.pl <file.fastq> <minimal length in nt> <maximum length in nt>

open(IN, "<", "$ARGV[0]") || die;
open(OUT, ">", "filtered.fastq") || die;

my $min = $ARGV[1];
my $max = $ARGV[2];

my $line1; my $line2; my $line3; my $line4;
my $counter = 1;

while(my $line = <IN>){
  chomp $line;
  if($counter == 1){
    $line1 = $line;
    $counter++;
  }elsif($counter == 2){
    $line2 = $line;
    $counter++;
  }elsif($counter == 3){
    $line3 = $line;
    $counter++;
  }elsif($counter == 4){
    $line4 = $line;
    $counter = 1;
    if(length($line2) >= $min && length($line2) <= $max){
      print OUT $line1, "\n", $line2, "\n", $line3, "\n", $line4, "\n";
    }else{
      next;
    }
  }
}

close IN;
close OUT;
