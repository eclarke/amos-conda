#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $help = 0;
my $GENOMELEN = 1000000;
my $READLEN = 100;
my $COV = 10;



my $USAGE = "covanalysis.pl [-genome len] [-read len] [-cov cov]\n";

my $res = GetOptions("help"     => \$help,
                     "genome=s" => \$GENOMELEN,
                     "read=s"   => \$READLEN,
                     "cov=s"    => \$COV);

if ($help || !$res)
{
  print $USAGE;
  print "\n";
  print "  Simulates sequencing of a genome to report the coverage distribution\n";
  print "\n";
  print "Options\n";
  print " -genome <len> : Specify the genome length (Default: $GENOMELEN)\n";
  print " -read <len>   : Specify the read length (Default: $READLEN)\n";
  print " -cov <cov>    : Specify the desired coverage (Default: $COV)\n";

  exit 0;
}


my $desiredbases = $GENOMELEN * $COV;

print STDERR "Simulating $COV X cov for the $GENOMELEN bp genome with $READLEN bp reads $desiredbases total bases\n"; 


## getReadLen
###########################################################

sub getReadLen
{
  my $rl = $READLEN;
  return $rl;
}


## initialize coverage
###########################################################

my @cov;
for (my $i = 0; $i < $GENOMELEN; $i++)
{
  $cov[$i] = 0;
}


## simulate reads
###########################################################

my $sumbases = 0; 

while ($sumbases < $desiredbases)
{
  my $start = rand($GENOMELEN);  ## random start position
  my $rl = getReadLen();

  $sumbases += $rl;

  for (my $i = 0; $i < $rl; $i++)
  {
    $cov[($start + $i) % $GENOMELEN]++;
  }
}


## determine max cov
###########################################################

my $maxcov = 0;
for (my $i = 0; $i < $GENOMELEN; $i++)
{
  if ($cov[$i] > $maxcov)
  {
    $maxcov = $cov[$i];
  }
}

print STDERR "Max coverage: $maxcov\n";


## compute histogram
###########################################################

my @hist;
for (my $c = 0; $c <= $maxcov; $c++)
{
  $hist[$c] = 0;
}

for (my $i = 0; $i < $GENOMELEN; $i++)
{
  $hist[$cov[$i]]++;
}


## print histogram
###########################################################

print "depth\tcnt\n";
for (my $c = 0; $c <= $maxcov; $c++)
{
  my $cnt = $hist[$c];
  print "$c\t$cnt\n";
}


