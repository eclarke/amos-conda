#!/usr/bin/perl
use lib "/bluearc/data/schatz/mschatz/devel/amos/lib";

use strict;
use Getopt::Long;

my $USAGE = "ctg_motif <options> ctg.fa\n"; 

## Options
my $help  = 0;

my $WINDOW = 50;
my $MIN_CTG_LEN = 500;
my $INTERIOR="I";
my $END="E";
my $MOTIF = undef;

my $res = GetOptions("help"       => \$help,
                     "window=s"   => \$WINDOW,
                     "minctg=s"   => \$MIN_CTG_LEN,
                     "motif=s"    => \$MOTIF,
                     "interior=s" => \$INTERIOR,
                     "end=s"      => \$END,
                     );
 
if ($help || !$res)
{
  print $USAGE;
  print "\n";
  print "  Compute the number of occurrences of a specified motif in the interior versus contig ends\n";
  print "  Output format: <code> \\t <window motif cnt>\n";
  print "\n";
  print "Required\n";
  print "  ctg.fa          : path to draft assembly\n";
  print "  -motif <s>      : motif to count\n";
  print "\n";
  print "Options\n";
  print "  -window <len>  : window length (default: $WINDOW)\n";
  print "  -minctg <len>  : minimum contig length to examine (default: $MIN_CTG_LEN)\n";
  print "  -interior <c>  : Code for ctg interior (default: $INTERIOR)\n";
  print "  -end <c>       : code for ctg end (default: $END)\n";
  exit 0;
}

die "Must specify -motif\n" if !defined $MOTIF;


## Count Occurences of a motif
###############################################################################
sub countmotif
{
  my $str = $_[0];
  my @all = ($str =~ /TCTA/g); 
  my $g = scalar @all;
  return $g+0;
}



## Process a single contig
###############################################################################
sub processContig
{
  my $seq = $_[0];

  return if !defined $seq;

  if (length($seq) >= $MIN_CTG_LEN)
  {
    my $gc = countmotif(substr($seq, 0, $WINDOW));
    print "$END\t$gc\n";

    for (my $pos = $WINDOW; $pos+1.5*$WINDOW < length($seq); $pos+=$WINDOW)
    {
      my $wstr = substr($seq, $pos, $WINDOW);
      $gc = countmotif($wstr);
      print "$INTERIOR\t$gc\n";
    }

    $gc = countmotif(substr($seq, -$WINDOW));
    print "$END\t$gc\n";
  }
}



## Scan the multifasta file
###############################################################################
my $seq = undef;

while (<>)
{
  if (/^>/)
  {
    processContig($seq);
    $seq = "";
  }
  else
  {
    chomp;
    $seq .= $_;
  }
}

processContig($seq);
