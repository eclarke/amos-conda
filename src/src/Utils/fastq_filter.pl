#!/usr/bin/perl 

my $USAGE = "fastq_filter skiplist file.fq > filtered.fq\n";

my $skiplist = shift @ARGV or die $USAGE;
my $fastq    = shift @ARGV or die $USAGE;

open SKIP, "$skiplist" or die "Can't open $skiplist ($!)\n";
open FASTQ, "$fastq"   or die "Can't open $fastq ($!)\n";

my %filter;


## Load Skip list
###############################################################################

print STDERR "Loading reads to ignore from $skiplist...";

my $cnt = 0;

while (<SKIP>)
{
  $cnt++;
  if (($cnt % 100000) == 0) { print STDERR " $cnt"; }

  chomp;

  my ($name) = split /\s/;
  $filter{$name} = 1;
  #print STDERR "saw $name\n";
}

my $cnt = scalar keys %filter;
print STDERR "\n$cnt reads loaded\n";



## Scan fastq file
###############################################################################

my $printed = 0;
my $all = 0;

print STDERR "Processing $fastq ... ";

while (<FASTQ>)
{
  if ($_ =~ /^@(\S+)\/\d/)
  {
    $all++;

    if (($all % 1000000) == 0) { print STDERR " $all"; }

    #print STDERR "Checking $1 of $_";
    my $h1 = $1;
    my $s  = <FASTQ>;
    my $h2 = <FASTQ>;
    my $q  = <FASTQ>;

    if (!exists $filter{$h1})
    {
      $printed++;

      print $_;
      print $s;
      print "+\n";
      print $q;
    }
  }
  else
  {
    die "bad fastq format\n";
  }
}

my $skip = $all - $printed;
print STDERR "\nPrinted $printed of $all (skipped $skip)\n";
