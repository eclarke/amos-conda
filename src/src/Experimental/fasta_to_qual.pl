#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $help = 0;
my $QUAL = 60;

my $USAGE = "fasta_to_qual <options> ctg.fa > ctg.qual\n";

my $res = GetOptions("help"  => \$help,
                     "qv=s" => \$QUAL);

if ($help || !$res)
{
  print $USAGE;
  print "\n";
  print "  Create a TIGR/CA qual file from a fasta file\n";
  print "\n";
  print "Options\n";
  print "  -qv <qual> : quality value to use (default: $QUAL)\n";
  exit 0;
}

die $USAGE if (scalar @ARGV == 0);


sub outputQual
{
  my $seqid = shift;
  my $seq = shift;

  return if !defined $seq;

  my $seqlen = length($seq);

  print "$seqid\n";

  for (my $i = 0; $i < $seqlen; $i++)
  {
    print " " if ($i > 0);
    print "$QUAL";
  }
  print "\n";
}


my $seqid = undef;
my $seq = undef;

while (<>)
{
  chomp;

  if (/^>/)
  {
    outputQual($seqid, $seq);

    $seqid = $_;
    $seq = "";
  }
  else
  {
    $seq .= $_;
  }
}

outputQual($seqid, $seq);
