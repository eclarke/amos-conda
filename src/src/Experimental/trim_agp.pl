#!/usr/bin/perl -w 
use strict;

## AGP Codes
my $GAPCODE = "N";
my $CTGCODE = "W";

## Debugging options
#####################################################################

## Minimum contig length to report
my $MINCTGLEN = 200;

## Set to 0 to use old ids, 1 to reindex
my $UPDATEIDX = 1;

## Set to 0 to not slide scaffold because of a leading gap
my $SLIDESCAFFOLD = 1;


my $USAGE = "trim_agp.pl trimlist agp fa\n";

my $trimlist = shift @ARGV or die $USAGE;
my $agpfile  = shift @ARGV or die $USAGE;
my $fafile   = shift @ARGV or die $USAGE;


## Load the list of trims to apply
#####################################################################

open TRIM, "$trimlist" or die "Cant open $trimlist ($!)\n";

my @alltrims;
my $trimcnt = 0;

while (<TRIM>)
{
  next if /^#/;
  chomp;

  my ($contig, $len, $poslist, $source) = split /\s/, $_;

  my $item;
  $item->{ctg} = $contig;
  $item->{len} = $len;
  $item->{poslist} = $poslist;
  $item->{source} = $source;

 # print STDERR "$contig\t[$len]\t$source\t:\t$poslist\n";

  push @alltrims, $item;

  $trimcnt++;
}

print STDERR "Loaded $trimcnt trims\n";


## Load the original AGP file
#####################################################################

open AGP, "$agpfile" or die "Can't open $agpfile ($!)\n";

my %ctg2scf;
my %scaffolds;

my $agplines = 0;
my $ctgcnt = 0;

while (<AGP>)
{
  chomp;
  $agplines++;

  next if /^#/; 

  my ($scf, $sstart, $send, $idx, $code, $ctgid, $cstart, $cend, $oo) = split /\t/, $_;

  ## is a ctg
  if ($code eq "W") 
  { 
    $ctg2scf{$ctgid} = $scf; 
    $ctgcnt++;
  }

  my $item;
  $item->{sstart} = $sstart;
  $item->{send}   = $send;
  $item->{idx}    = $idx;
  $item->{code}   = $code;
  $item->{ctgid}  = $ctgid;
  $item->{cstart} = $cstart;
  $item->{cend}   = $cend; 
  $item->{oo}     = $oo;
  
  push @{$scaffolds{$scf}}, $item;
}

my $scfcnt = scalar keys %scaffolds;

print STDERR "Loaded $ctgcnt contigs; $scfcnt scaffolds\n";


## Create a log of the ctg regions
open CTGEDIT, "> $agpfile.ctgedit" or die "Cant open $agpfile.ctgedit\n";
my %ctgregions;



## Map the trims to the scaffolds
#####################################################################

my %scaffchanges;

for(my $trimidx = 0; $trimidx < scalar @alltrims; $trimidx++)
{
  my $ctg = $alltrims[$trimidx]->{ctg};
  my $scf = $ctg2scf{$ctg};

  push @{$scaffchanges{$scf}}, $trimidx;
}



## Edit the scaffolds
#####################################################################

print STDERR "Editing the scaffolds\n";

foreach my $scfid (sort keys %scaffchanges)
{
  my @trims = @{$scaffchanges{$scfid}};

  my $trimcnt = scalar @trims;

  print STDERR "$scfid: $trimcnt\n";

  my $scf = $scaffolds{$scfid};

  foreach my $t (@trims)
  {
    my $trim    = $alltrims[$t];
    my $ctg     = $trim->{ctg};
    my $clen    = $trim->{len};
    my $poslist = $trim->{poslist};
    my $source  = $trim->{source};

    print STDERR "  $t: $ctg\t[$clen]\t$source\t:\t$poslist\n";

    for(my $i = 0; $i < scalar @{$scf}; $i++)
    {
      if ($scf->[$i]->{ctgid} eq $ctg)
      {
         my $item = $scf->[$i];

         ## remove the old contig
         splice(@{$scf}, $i, 1);

         my $sstart = $item->{sstart};
         my $send   = $item->{send};
         my $code   = $item->{code};
         my $idx    = $item->{idx};
         my $ctgid  = $item->{ctgid};
         my $cstart = $item->{cstart};
         my $cend   = $item->{cend};
         my $oo     = $item->{oo};

         print STDERR "   $scfid\t$sstart\t$send\t$idx\t$code\t$ctgid\t$cstart\t$cend\t$oo\n";

         ## Parse the regions to exclude
         ############################################################

         my @posarr;
         my @raw = split /,/, $poslist;

         foreach my $pair (@raw)
         {
           my ($junkstart, $junkend) = split /\.\./, $pair;

           my $junk;
           $junk->{start} = $junkstart;
           $junk->{end}   = $junkend;

           push @posarr, $junk;
         }

         @posarr = sort {$a->{start} <=> $b->{start}} @posarr;


         ## Now clean up th regions to exclude
         ############################################################

         my $nummerges = 0;

         my @mergedsplits;
         push @mergedsplits, $posarr[0];
          
         ## Check if there is valid contig before the first split
         if ($mergedsplits[0]->{start} < $MINCTGLEN)
         {
           print STDERR "   ***";
           $mergedsplits[0]->{start} = 1;
           $nummerges++;
         }

         ## merge consecutive intervals that are too close together
         for (my $j = 1; $j < scalar @posarr; $j++)
         {
           my $nummerged = scalar @mergedsplits;
           if ($mergedsplits[$nummerged-1] + $MINCTGLEN >= $posarr[$j]->{start})
           {
             ## Merge posarr[j] with the last mergedsplits to join them together
             $mergedsplits[$nummerged-1]->{end} = $posarr[$j]->{end};
             print STDERR " +++";
             $nummerges++;
           }
           else
           {
             ## Start a new interval to exclude
             push @mergedsplits, $posarr[$j];
           }
         }

         ## check if there is a valid contig after the last split
         if ($mergedsplits[scalar @mergedsplits - 1]->{end} + $MINCTGLEN >= $clen)
         {
           $mergedsplits[scalar @mergedsplits - 1]->{end} = $clen;
           print STDERR " \$\$\$";
           $nummerges++;
         }

         if ($nummerges > 0)
         {
           print STDERR "\n";
         }

         ## print the bases to remove after merging
         for (my $j = 0; $j < scalar @mergedsplits; $j++)
         {
           my $js = $mergedsplits[$j]->{start};
           my $je = $mergedsplits[$j]->{end};

           print STDERR "    $js - $je\n";
         }


         ## now apply the edits
         ############################################################

         if ((scalar @mergedsplits == 1) &&
             ($mergedsplits[0]->{start} == 1) && 
             ($mergedsplits[0]->{end} == $clen))
         {
             ## replace entire contig with a gap
             print STDERR "   => gap [1 $clen] | [$sstart $send]\n";

             ## Mark that the entire contig should be excluded
             print CTGEDIT "$ctgid\t1\t-1\t-1\n";
             $ctgregions{$ctgid}->[1]->{start} = -1;
             $ctgregions{$ctgid}->[1]->{end}   = -1;

             my $ngap;
             $ngap->{sstart} = $sstart;
             $ngap->{send}   = $send;
             $ngap->{code}   = $GAPCODE;
             $ngap->{idx}    = -1;
             $ngap->{ctgid}  = $clen;
             $ngap->{cstart} = "scaffold";
             $ngap->{cend}   = "yes";
             $ngap->{oo}     = "paired-ends";

             push @{$scf}, $ngap;
         }
         else
         {
           my $segment = 0;

           for(my $j = 0; $j < scalar @mergedsplits; $j++)
           {
             my $js = $mergedsplits[$j]->{start};
             my $je = $mergedsplits[$j]->{end};

             if ($j == 0)
             {
               if ($js > 1)
               {
                 ## start with a contig prefix
                 $segment++;

                 my $ncend   = $js - 1;
                 my $nclen   = $ncend;
                 my $nsstart = $sstart;
                 my $nsend   = $sstart + $ncend-1;
                 my $nslen   = $nsend - $nsstart + 1;

                 print STDERR "   => $ctgid.$segment [1 $ncend] $nclen | [$nsstart $nsend] $nslen\n";

                 print CTGEDIT "$ctgid\t$segment\t1\t$ncend\n";
                 $ctgregions{$ctgid}->[$segment]->{start} = 1;
                 $ctgregions{$ctgid}->[$segment]->{end}   = $ncend;

                 $ctg2scf{$ctgid} = $scfid;

                 my $nctg;
                 $nctg->{sstart} = $nsstart;
                 $nctg->{send}   = $nsend;
                 $nctg->{code}   = $CTGCODE;
                 $nctg->{idx}    = -1;
                 $nctg->{ctgid}  = "$ctgid.$segment";
                 $nctg->{cstart} = 1;
                 $nctg->{cend}   = $nclen;
                 $nctg->{oo}     = "+";

                 push @{$scf}, $nctg;
               }
             }
             else
             {
               ## interior substring
               $segment++;
               
               my $ncstart = $mergedsplits[$j-1]->{end}+1;
               my $ncend   = $js - 1;
               my $nclen   = $ncend - $ncstart;

               my $nsstart = $sstart + $ncstart - 1;
               my $nsend   = $sstart + $ncend   - 1;

               print STDERR "   => $ctgid.$segment $ncstart $ncend | $nsstart $nsend\n"; 

               print CTGEDIT "$ctgid\t$segment\t$ncstart\t$ncend\n";
               $ctgregions{$ctgid}->[$segment]->{start} = $ncstart;
               $ctgregions{$ctgid}->[$segment]->{end}   = $ncend;

               $ctg2scf{$ctgid} = $scfid;

               die "TODO";

               my $nctg;
               $nctg->{sstart} = $nsstart;
               $nctg->{send}   = $nsend;
               $nctg->{code}   = $CTGCODE;
               $nctg->{idx}    = -1;
               $nctg->{ctgid}  = "$ctgid.$segment";
               $nctg->{cstart} = 1;
               $nctg->{cend}   = $nclen;
               $nctg->{oo}     = "+";

               push @{$scf}, $nctg;
             }

             my $nglen = $je - $js+1;
             my $nsstart = $sstart + $js-1;
             my $nsend   = $sstart + $je-1;
             my $nslen   = $nsend - $nsstart + 1;

             print STDERR "   => gap [$js $je] $nglen | [$nsstart $nsend] $nslen\n";

             my $ngap;
             $ngap->{sstart} = $nsstart;
             $ngap->{send}   = $nsend;
             $ngap->{code}   = $GAPCODE;
             $ngap->{idx}    = -1;
             $ngap->{ctgid}  = $nglen;
             $ngap->{cstart} = "scaffold";
             $ngap->{cend}   = "yes";
             $ngap->{oo}     = "paired-ends";

             push @{$scf}, $ngap;
           }

           my $lje = $mergedsplits[scalar @mergedsplits - 1]->{end};

           if ($lje != $clen)
           {
             ## update final suffix of original contig
             $segment++;

             my $ncstart = $lje+1;
             my $ncend   = $clen;
             my $nclen   = $ncend - $ncstart + 1;

             my $nsstart = $sstart + $ncstart - 1;
             my $nsend   = $sstart + $ncend - 1;
             my $nslen   = $nsend - $nsstart + 1;

             print STDERR "   => $ctgid.$segment [$ncstart $ncend] $nclen | [$nsstart $nsend] $nslen\n"; 

             print CTGEDIT "$ctgid\t$segment\t$ncstart\t$ncend\n";
             $ctgregions{$ctgid}->[$segment]->{start} = $ncstart;
             $ctgregions{$ctgid}->[$segment]->{end}   = $ncend;

             $ctg2scf{$ctgid} = $scfid;

             my $nctg;
             $nctg->{sstart} = $nsstart;
             $nctg->{send}   = $nsend;
             $nctg->{code}   = $CTGCODE;
             $nctg->{idx}    = -1;
             $nctg->{ctgid}  = "$ctgid.$segment";
             $nctg->{cstart} = 1;
             $nctg->{cend}   = $nclen;
             $nctg->{oo}     = "+";

             push @{$scf}, $nctg;
           }
         }
      }
    }
  }

  ## Resort to put the items into order along the scaffold
  @{$scaffolds{$scfid}} = sort {$a->{sstart} <=> $b->{sstart}} @{$scaffolds{$scfid}};

  ## now check for consecutive gaps

  my $i = 0;

  while ($i < (scalar @{$scaffolds{$scfid}})-1)
  {
    if (($scaffolds{$scfid}->[$i]->{code} eq $GAPCODE) &&
        ($scaffolds{$scfid}->[$i+1]->{code} eq $GAPCODE))
    {
      print STDERR "Found consecutive gaps in $scfid\n";
      my $gap1 = $scaffolds{$scfid}->[$i];
      my $gap2 = $scaffolds{$scfid}->[$i+1];

      my $g1sstart = $gap1->{sstart};
      my $g1send   = $gap1->{send};
      my $g1len    = $gap1->{ctgid};

      my $g2sstart = $gap2->{sstart};
      my $g2send   = $gap2->{send};
      my $g2len    = $gap2->{ctgid};

      my $j = $i+1;

      print STDERR " $i: [$g1sstart $g1send] $g1len\n";
      print STDERR " $j: [$g2sstart $g2send] $g2len\n";

      $scaffolds{$scfid}->[$i]->{send}   = $scaffolds{$scfid}->[$i+1]->{send};
      $scaffolds{$scfid}->[$i]->{ctgid} += $scaffolds{$scfid}->[$i+1]->{ctgid};

      $g1sstart = $scaffolds{$scfid}->[$i]->{sstart};
      $g1send   = $scaffolds{$scfid}->[$i]->{send};
      $g1len    = $scaffolds{$scfid}->[$i]->{ctgid};

      print STDERR " => $i: [$g1sstart $g1send] $g1len\n";

      ## delete the second gap
      splice(@{$scaffolds{$scfid}}, $i+1, 1);
    }

    $i++;
  }

  ## finally trim any gaps at the very beginning or end

  if ($scaffolds{$scfid}->[0]->{code} eq $GAPCODE)
  {
    print STDERR "scf $scfid begins with a gap\n";

    if (scalar @{$scaffolds{$scfid}} == 1)
    {
      print STDERR "  => entire scaffold is a gap!\n";
      delete $scaffolds{$scfid};
    }
    else
    {
      ## remove the leading gap
      splice (@{$scaffolds{$scfid}}, 0, 1);

      my $shift = $scaffolds{$scfid}->[0]->{sstart} - 1;

      print STDERR "  => shifting over by $shift\n";

      if ($SLIDESCAFFOLD)
      {
        ## now slide everything over
        for (my $i = 0; $i < scalar @{$scaffolds{$scfid}}; $i++)
        {
          $scaffolds{$scfid}->[$i]->{sstart} -= $shift;
          $scaffolds{$scfid}->[$i]->{send}  -= $shift;
        }
      }
    }
  }

  if (defined $scaffolds{$scfid})
  {
    if ($scaffolds{$scfid}->[scalar @{$scaffolds{$scfid}} - 1]->{code} eq $GAPCODE)
    {
      print STDERR "scf $scfid ends with a gap, deleting\n";
      splice(@{$scaffolds{$scfid}}, scalar @{$scaffolds{$scfid}}-1, 1);
    }
  }

  print STDERR "\n\n";
}

close CTGEDIT;



## Print the edited AGP File
#####################################################################

open AGPEDIT, "> $agpfile.edited" or die "Cant open $agpfile.edited\n";

print AGPEDIT "## AGP-version 2.0\n";
print AGPEDIT "## AGP trimmed to remove contamination\n";

foreach my $scf (sort {substr($a, 9) <=> substr($b, 9)} keys %scaffolds)
{
  my $idx = 0;

  foreach my $item (@{$scaffolds{$scf}})
  {
    $idx++;

    my $sstart = $item->{sstart};
    my $send   = $item->{send};
    my $oidx   = $item->{idx};
    my $code   = $item->{code};
    my $ctgid  = $item->{ctgid};
    my $cstart = $item->{cstart};
    my $cend   = $item->{cend};
    my $oo     = $item->{oo};

    if ($UPDATEIDX)
    {
      $oidx = $idx;
    }

    print AGPEDIT "$scf\t$sstart\t$send\t$oidx\t$code\t$ctgid\t$cstart\t$cend\t$oo\n";
  }
}





## Now output the updated contig file
#####################################################################

print STDERR "\n\n";
print STDERR "Processing $fafile to split contigs\n";

open FASTA, "$fafile" or die "Cant open $fafile ($!)\n";
open OUT, ">$fafile.edited" or die "Cant open $fafile.edited ($!)\n";

sub processSeq
{
  my $seqid = shift;
  my $seq   = shift;

  return if !defined $seqid;

  if (!exists $ctg2scf{$seqid})
  {
    print STDERR "Cant find scaffold for $seqid!\n";
  }

  if (exists $ctgregions{$seqid})
  {
    for(my $i = 0; $i < scalar @{$ctgregions{$seqid}}; $i++)
    {
      next if !defined $ctgregions{$seqid}->[$i];

      my $s = $ctgregions{$seqid}->[$i]->{start};
      my $e = $ctgregions{$seqid}->[$i]->{end};
      my $l = $e - $s + 1;

      next if $s == -1;  ## skip contigs that were entirely removed

      my $substring = substr($seq, $s-1, $l);

      print OUT ">$seqid.$i\n";
      print OUT "$substring\n";
    }
  }
  else
  {
    print OUT ">$seqid\n";
    print OUT "$seq\n";
  }
}


my $seqid = undef;
my $seq;

my $factgcnt = 0;

while (<FASTA>)
{
  if (/^>(\S+)/)
  {
    processSeq($seqid, $seq);
    $factgcnt++;

    $seqid = $1;
    $seq = "";
  }
  else
  {
    chomp;
    $seq .= $_;
  }
}

processSeq($seqid, $seq);
$factgcnt++;

print STDERR "Processed $factgcnt contigs\n";


