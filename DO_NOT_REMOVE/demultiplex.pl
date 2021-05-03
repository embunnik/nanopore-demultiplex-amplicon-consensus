use warnings;
use strict;

# demultiplex.pl 26APR2021 by EMB to demultiplex nanopore sequence data containing forward and reverse 12-nt barcodes
# script searches for perfect matches for forward and reverse barcodes first
# if no perfect match, then searches for up to two mismatches, or one indel

print "\n", "\tWhat is the highest barcode number in the sample?", "\n";
print "\tFor example: 100", "\n";
print "\t>";
my $barcodes = <STDIN>;
my @barcodes;
for(my $i = 1; $i < $barcodes + 1; $i++){
  push(@barcodes, $i);
}

my $start_time = time();

my %fh;
for my $barcode(@barcodes){
  my $filename = "barcode_${barcode}.fastq";
  open($fh{$barcode}, ">", $filename) or die "$filename: $!";
}

open(IN, "<", "./DO_NOT_REMOVE/barcodes.csv") or die $!;
open(SEQ, "<", $ARGV[0]) or die $!;

my %barcodes;

while(my $line = <IN>){
  chomp $line;
  my @data1 = split(',', $line);
  if($data1[0] > $barcodes){last}
  $barcodes{$data1[0]}{forward}{perfect} = $data1[1];
  $barcodes{$data1[0]}{reverse}{perfect} = $data1[2];
  $barcodes{$data1[0]}{forward_RC}{perfect} = revcom($data1[1]);
  $barcodes{$data1[0]}{reverse_RC}{perfect} = revcom($data1[2]);

  my @imperfect_fw;
  mismatch($data1[1], \@imperfect_fw);
  $barcodes{$data1[0]}{forward}{imperfect} = [@imperfect_fw];
  my @imperfect_rv;
  mismatch($data1[2], \@imperfect_rv);
  $barcodes{$data1[0]}{reverse}{imperfect} = [@imperfect_rv];
  my @imperfect_fw_RC;
  mismatch(revcom($data1[1]), \@imperfect_fw_RC);
  $barcodes{$data1[0]}{forward_RC}{imperfect} = [@imperfect_fw_RC];
  my @imperfect_rv_RC;
  mismatch(revcom($data1[2]), \@imperfect_rv_RC);
  $barcodes{$data1[0]}{reverse_RC}{imperfect} = [@imperfect_rv_RC];
}

my $line1; my $line2; my $line3; my $line4;
my $counter = 1;
my $flag = 0;
my $match;
my $total;
my $n_matches;
my $below27 = 0;
my $above27 = 0;

while(my $line = <SEQ>){
  chomp $line;
  if($counter == 1){
    $line1 = $line;
    $total++;
    my $count = substr($total, -3);
    if($count eq "000"){
      print "\tProcessed $total reads\n";
    }
    $counter++;
  }elsif($counter == 2){
    $counter++;
    $line2 = $line;
    my $start = substr($line2, 0, 100);
    my $end = substr($line2, -100);

    foreach my $seq (sort {$a <=> $b} keys %barcodes){
      if($start =~ m/$barcodes{$seq}{forward}{perfect}/ || $start =~ m/$barcodes{$seq}{reverse}{perfect}/ || $end =~ m/$barcodes{$seq}{forward_RC}{perfect}/ || $end =~ m/$barcodes{$seq}{reverse_RC}{perfect}/){
        if($start =~ m/$barcodes{$seq}{forward}{perfect}/ && $end =~ m/$barcodes{$seq}{reverse_RC}{perfect}/){
          $flag = 1;
          $match = $seq;
        }elsif($start =~ m/$barcodes{$seq}{reverse}{perfect}/ && $end =~ m/$barcodes{$seq}{forward_RC}{perfect}/){
          $flag = 1;
          $match = $seq;
        }elsif($start =~ m/$barcodes{$seq}{forward}{perfect}/){
          ($flag, $match) = match(\@{$barcodes{$seq}{reverse_RC}{imperfect}}, $end, $seq);
          if(!$flag){
            ($flag, $match) = match2($barcodes{$seq}{reverse_RC}{perfect}, $end, $seq);
          }
        }elsif($start =~ m/$barcodes{$seq}{reverse}{perfect}/){
          ($flag, $match) = match(\@{$barcodes{$seq}{forward_RC}{imperfect}}, $end, $seq);
          if(!$flag){
            ($flag, $match) = match2($barcodes{$seq}{forward_RC}{perfect}, $end, $seq);
          }
        }elsif($end =~ m/$barcodes{$seq}{forward_RC}{perfect}/){
          ($flag, $match) = match(\@{$barcodes{$seq}{reverse}{imperfect}}, $start, $seq);
          if(!$flag){
            ($flag, $match) = match2($barcodes{$seq}{reverse}{perfect}, $start, $seq);
          }
        }elsif($end =~ m/$barcodes{$seq}{reverse_RC}{perfect}/){
          ($flag, $match) = match(\@{$barcodes{$seq}{forward}{imperfect}}, $start, $seq);
          if(!$flag){
            ($flag, $match) = match2($barcodes{$seq}{forward}{perfect}, $start, $seq);
          }
        }
      }

      if($flag){
        last;
      }else{
        my($partial_flag1, $partial_flag2, $partial_flag3, $partial_flag4);
        my($partial_match1, $partial_match2, $partial_match3, $partial_match4);
        ($partial_flag1, $partial_match1) = match(\@{$barcodes{$seq}{forward}{imperfect}}, $start, $seq);
        if(!$partial_flag1){
          ($partial_flag2, $partial_match2) = match2($barcodes{$seq}{forward}{perfect}, $start, $seq)
        }
        if($partial_flag1 || $partial_flag2){
          ($flag, $match) = match(\@{$barcodes{$seq}{reverse_RC}{imperfect}}, $end, $seq);
          if(!$flag){
            ($flag, $match) = match2($barcodes{$seq}{reverse_RC}{perfect}, $end, $seq);
          }
        }else{
          ($partial_flag3, $partial_match3) = match(\@{$barcodes{$seq}{reverse}{imperfect}}, $start, $seq);
          if(!$partial_flag3){
            ($partial_flag4, $partial_match4) = match2($barcodes{$seq}{reverse}{perfect}, $start, $seq)
          }
          if($partial_flag3 || $partial_flag4){
            ($flag, $match) = match(\@{$barcodes{$seq}{forward_RC}{imperfect}}, $end, $seq);
            if(!$flag){
              ($flag, $match) = match2($barcodes{$seq}{forward_RC}{perfect}, $end, $seq);
            }
          }
        }
      }
      if($flag){
        last;
      }
    }
  }elsif($counter == 3){
    $line3 = $line;
    $counter++;
  }elsif($counter == 4){
    $line4 = $line;
    $counter = 1;
    if($flag){
      print { $fh{$match} } $line1, "\n", $line2, "\n", $line3, "\n", $line4, "\n";
      $flag = 0;
      $n_matches++;
      if($match < 28){
        $below27++;
      }else{
        $above27++;
      }
      undef $match;
    }
  }
}

my $end_time = time();

print "\n\tTotal number of reads: ", $total, "\n\n";
print "\tPercentage of reads matched to barcode: ", $n_matches / ($total) * 100, "%", "\n\n";
print "\tNumber of reads matched to barcodes 1 - 27: ", $below27, "\n\n";
print "\tNumber of reads matched to barcodes 28 - 100: ", $above27, "\n\n";
printf("\tExecution Time: %0.02f s\n", $end_time - $start_time);

sub revcom{
  my ($seq) = @_;
  my $seq2 = reverse $seq;
  $seq2 =~ tr/ATGC/TACG/;
  return $seq2;
}

sub mismatch{
  my ($barcode, $array) = @_;
  for(my $i = 0; $i < 12; $i++){
    @$array[$i] = $barcode;
    substr(@$array[$i], $i, 1, ".{0,2}");
  }
}

sub match{
  my ($bc, $string, $seq) = @_;
  my $flag;
  for(my $i = 0; $i < 12; $i++){
    if($string =~ m/${$bc}[$i]/){
      $flag = 1;
      last;
    }
  }
  return ($flag, $seq);
}

sub match2{
  my ($bc, $string, $seq) = @_;
  my @str1 = split(//, $bc);
  my $flag;
  for(my $i = 0; $i < (length($string) - 12); $i++){
    my $nts2 = substr($string, $i, 12);
    my @str2 = split(//, $nts2);
    my $common;

    for(my $j = 0; $j < 12; $j++){
      if($str1[$j] eq $str2[$j]){
        $common++;
      }
    }
    if($common && $common > 9){
      $flag = 1;
      last;
    }
  }
  return ($flag, $seq);
}

exit;
