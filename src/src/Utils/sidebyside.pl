#!/usr/bin/perl -w

use strict;
use Getopt::Long;

my $USAGE = "sidebyside.pl f1 f2\n";
my $help  = 0;
my $width = 50;

my $marker_sam = "==";
my $marker_dif = "!!";

my $res = GetOptions("help" => \$help,
                     "width=s" => \$width,
                     "same=s"  => \$marker_sam,
                     "dif=s"   => \$marker_dif);

if ($help || !$res)
{
  print $USAGE;
  print "\n";
  print "  Print lines from two files side by side, and mark them as the same (==) or different (!!)\n";
  print "\n";
  print "Required\n";
  print " f1 : first file\n";
  print " f2 : second file\n";
  print "\n";
  print "Options\n";
  print "  -width <n>  : line width for each file (default: $width)\n";
  print "  -same <s>   : marker for lines that are the same (default: $marker_sam)\n";
  print "  -dif <s>    : marker for lines that are different (default: $marker_dif)\n";

  exit 0;
}

die $USAGE if (scalar @ARGV != 2);

my $f1 = shift @ARGV;
my $f2 = shift @ARGV;

open F1, "<$f1" or die "Can't open $f1 ($!)\n";
open F2, "<$f2" or die "Can't open $f2 ($!)\n";

while (<F1>)
{
  chomp;

 # if (/^\s*$/){next;}

  my $l2 = <F2>;
  chomp $l2;

  my $marker = ($_ eq $l2) ? $marker_sam : $marker_dif;

  printf("%-${width}s    $marker    %-${width}s\n",
         substr($_,0,$width), 
         substr($l2,0,$width));
}
