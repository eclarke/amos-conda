#!/usr/bin/perl -w
use strict;

my $USAGE = "find-redundant.pl <delta> > redundant.list\n";

my $delta = shift @ARGV or die $USAGE;

open COORDS, "show-coords -rclH $delta |" 
  or die "Can't run show-coords -rclH $delta\n";

my %contigs;

my $MIN_PERC_ID = 95;
my $MIN_PERC_COV = 95;

while (<COORDS>)
{
  chomp;
  my @vals = split /\s+/, $_;

  my $pid  = $vals[10];
  my $lenr = $vals[12];
  my $lenq = $vals[13];

  my $covr = $vals[15];
  my $covq = $vals[16];

  my $rid  = $vals[18];
  my $qid  = $vals[19];

  next if ($rid ge $qid);
  next if ($pid < $MIN_PERC_ID);

  if (($covr >= $MIN_PERC_COV) && ($covq >= $MIN_PERC_COV))
  {
    print "$_\n";
    if ($lenr >= $lenq)
    {
      push @{$contigs{$rid}->{contains}}, $qid;
      push @{$contigs{$qid}->{parent}}, $rid;
    }
    else
    {
      push @{$contigs{$qid}->{contains}}, $rid;
      push @{$contigs{$rid}->{parent}}, $qid;
    }
  }
  elsif ($covq >= $MIN_PERC_COV)
  {
    print "$_\n";
    push @{$contigs{$rid}->{contains}}, $qid;
    push @{$contigs{$qid}->{parent}}, $rid;
  }
  elsif ($covr >= $MIN_PERC_COV)
  {
    print "$_\n";
    push @{$contigs{$qid}->{contains}}, $rid;
    push @{$contigs{$rid}->{parent}}, $qid;
  }
}

my $cnt = 0;
foreach my $ctg (sort keys %contigs)
{
  if (exists $contigs{$ctg}->{parent})
  {
    $cnt++;
    print "$cnt $ctg :";
    foreach my $parent (@{$contigs{$ctg}->{parent}})
    {
      print " $parent";
    }

    print "\n";
  }
}


