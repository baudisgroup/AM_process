package PGX::GenomeIntervals;

use Data::Dumper;

require Exporter;
@ISA        =   qw(Exporter);
@EXPORT     =   qw(
  make_genome_intervals
  get_reference_base_limits
  get_genome_basecount
  segments_add_statusmaps
  interval_cnv_frequencies
);

################################################################################

sub make_genome_intervals {

=pod

Expects:
  - a list reference of genome interval objects, usually representing cytobands:
    [
      {
        no    =>  __integer__,            # not used
        reference_name  =>  __string__,
        start =>  __integer__,
        end   =>  __integer__,
        stain =>  __string__,            # not used
        label =>  __string__,            # not used
      },
      {
      ...
      },
    ]
  - the binning size in bases (optional; defaults to 1000000)

Returns:
  - a list of genome intervals in a similar structure, based on the binning size
  - the intervals are per reference, i.e. each reference starts with a new interval,
    leading to the last interval < binning size
    [
      {
        no    =>  __integer__,          # 1 -> n
        reference_name  =>  __string__,
        start =>  __integer__,
        end   =>  __integer__,
        label =>  __string__,
      },
      {
      ...
      },
    ]

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my (
    $cytobands,
    $intSize
  )             =   @_;

  if ($intSize !~ /^\d{3,9}$/) { $intSize = 1000000 }

  my $refLims   =   get_reference_base_limits($cytobands);
  my $gi        =   [];

  # references are sorted with numerical ones first, then others (e.g. 1 -> 22, X, Y)
  my @refNames  =   ((sort {$a <=> $b } grep{ /^\d\d?$/ } keys %$refLims), (sort grep{ ! /\d/ } keys %$refLims));

  my $intI      =   1;
  for my $i (0..$#refNames) {

    my $refName =   $refNames[$i];
    my $start   =   $refLims->{ $refName }->[0];
    my $end     =   $intSize - 1;

    while ($start < $refLims->{ $refName }->[1]) {

      # adjusting the end of the last interval
      if ($end > $refLims->{ $refName }->[1]) { $end = $refLims->{ $refName }->[1] };
      my $thisSize  =   $end - $start +1;
      push(
        @$gi,
        {
          no    =>  $intI,
          reference_name  =>  $refName,
          start =>  $start,
          end   =>  $end,
          length  =>  $thisSize,
          label =>  $refName.':'.$start.'-'.$end,
        }
      );

      $start    +=  $thisSize;
      $end      +=  $thisSize;
      $intI++;

  }}

  return $gi;

}

################################################################################

sub get_reference_base_limits {

=pod

Expects:
  - a list reference of genome interval objects:
    [
      {
        reference_name  =>  __string__,
        start           =>  __integer__,
        end             =>  __integer__,
        (others)...
      },
      {
        ...
      },
    ]

Returns:
  - a hash reference with each value consisting of a list of 2 integers, representing
    the reference's min and max bases:
    [
      { __string__  =>  [  __integer__,  __integer__ ] },
      { ... },
    ]

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $refLims   =   {};
  my $allRefs   =   $_[0];

  foreach my $ref (map{ $_->{reference_name} } @$allRefs) {
    my @refRefs =   grep{ $_->{reference_name} =~ /^$ref$/i } @$allRefs;
    my @bases   =   sort { $a <=> $b } ((map{ $_->{start} } @refRefs), (map{ $_->{end} } @refRefs));
    $refLims->{$ref}  =   [ $bases[0], $bases[-1] ];
  }

  return $refLims;

}

################################################################################

sub get_genome_basecount {

  my $genBases;
  my ($allRefs, $chr2plot)  =   @_;

  my %refNames  =   map { $_ => 1 } @$chr2plot;

  foreach my $ref (keys %refNames) {
    my @refRefs =   grep{ $_->{reference_name} =~ /^$ref$/i } @$allRefs;
    my @bases   =   sort { $a <=> $b } ((map{ $_->{start} } @refRefs), (map{ $_->{end} } @refRefs));
    $genBases   +=  ($bases[-1] - $bases[0]);
  }

  return $genBases;

}

################################################################################

sub segments_add_statusmaps {

  no warnings 'uninitialized';
  
  # 
  my $pgx       =   shift;
  my $segments  =   shift;
  
  if (ref $segments ne 'ARRAY') {
    $segments   =   $pgx->{segmentdata} }
  
=podmd

### Sub "segments_add_statusmaps"

The subroutine returns an object containing statusvalues (DUP, DEL) and the (min, max)
values of the overlapping variants, foreach of the provided $pgx->{genomeintervals} genome intervals.

Structure:

maps:
  dupmap:
    - null or DUP
    ...
  delmap:
    - null or DEL
    ...
  dupmax:
    - 0 or pos. value
    ...
  delmin:
    - 0 or neg. value
    ...

=cut

  my $maps      =   {
    intervals   =>  scalar(@{ $pgx->{genomeintervals} }),
    binning     =>  $pgx->{genomeintervals}->[0]->{end} - $pgx->{genomeintervals}->[0]->{start} + 1,
  };

  my %intStatLabs   =   (
    DUP         =>  'dupmap',
    DEL         =>  'delmap',
  );
  my %intCoverageLabs   =   (
    DUP         =>  'dupcoverage',
    DEL         =>  'delcoverage',
  );
  my %intValLabs   =   (
    DUP         =>  'dupmax',
    DEL         =>  'delmin',
  );

  my $cnvcoverage   =   {};

  foreach (values %intStatLabs) {
    $maps->{$_} =   [ map{''} 0..$#{ $pgx->{genomeintervals} } ] }
  foreach (values %intCoverageLabs) {
    $maps->{$_} =   [ map{ 0 } 0..$#{ $pgx->{genomeintervals} } ] }
  foreach (values %intValLabs) {
    $maps->{$_} =   [ map{ 0 } 0..$#{ $pgx->{genomeintervals} } ] }

  my $valueMap  =   [ map{[0]} 0..$#{ $pgx->{genomeintervals} } ];

  foreach my $csVar (@{ $segments }) {
  
    if (! grep{ $csVar->{variant_type} eq $_ } keys %intStatLabs) { next }

    # the index of  intervals with a match to the current variant is created and used
    # to assign the status value and collect the segment value (since several variants
    # may overlap the same interval)
    foreach my $ind (grep{
      $csVar->{reference_name} eq $pgx->{genomeintervals}->[ $_ ]->{reference_name}
      &&
      $csVar->{start}->[0] <= $pgx->{genomeintervals}->[ $_ ]->{end}
      &&
      $csVar->{end}->[-1] >=  $pgx->{genomeintervals}->[ $_ ]->{start}
    } 0..$#{ $pgx->{genomeintervals} }) {

      my $ovEnd     =   (sort { $a <=> $b } ($pgx->{genomeintervals}->[ $ind ]->{end},  $csVar->{end}->[-1]) )[0];
      my $ovStart   =   (sort { $b <=> $a } ($pgx->{genomeintervals}->[ $ind ]->{start},  $csVar->{start}->[0]) )[0];
      my $overlap   =   $ovEnd - $ovStart + 1;
    
      $maps->{ $intStatLabs{ $csVar->{variant_type } } }->[$ind] = $csVar->{variant_type};
      $maps->{ $intCoverageLabs{ $csVar->{variant_type } } }->[$ind]  +=  $overlap;
      push(
        @{ $valueMap->[$ind] },
        $csVar->{info}->{value},
      );

  }}
  
  foreach my $cLab (values %intCoverageLabs) {
    foreach my $ind (grep{ $maps->{$cLab}->[$_] > 0 } 0..$#{ $pgx->{genomeintervals} }) {
      $cnvcoverage->{$cLab}   +=  $maps->{$cLab}->[$ind];
      $cnvcoverage->{cnvcoverage} +=  $maps->{$cLab}->[$ind];
      $maps->{$cLab}->[$ind]  =   1 * ( sprintf "%.3f", $maps->{$cLab}->[$ind] / ($pgx->{genomeintervals}->[$ind]->{end} - $pgx->{genomeintervals}->[$ind]->{start} + 1) );
  }}

  # the values for each interval are sorted, to allow extracting the min/max 
  # values by position
  $valueMap     =   [ map{ [ sort { $a <=> $b } @{ $valueMap->[$_] } ] } 0..$#{ $pgx->{genomeintervals} } ];

  # the last of the sorted values is assigned iF > 0
  foreach my $ind (grep{ $valueMap->[$_]->[-1] > 0 } 0..$#{ $pgx->{genomeintervals} }) {
    $maps->{dupmax}->[$ind] =   1 * (sprintf "%.4f", $valueMap->[$ind]->[-1]) }

  # the first of the sorted values is assigned iF < 0
  foreach my $ind (grep{ $valueMap->[$_]->[0] < 0 } 0..$#{ $pgx->{genomeintervals} }) {
    $maps->{delmin}->[$ind] =   1 * (sprintf "%.4f", $valueMap->[$ind]->[0]) }

  $pgx->{statusmaps}    =   $maps;
  $pgx->{cnvstatistics} =   $cnvcoverage;
  
  if ($pgx->{genomesize} > 1) {
    foreach my $covK (keys %$cnvcoverage) {
      my $fracK =   $covK;
      $fracK    =~  s/coverage/fraction/;
      $pgx->{cnvstatistics}->{$fracK} =   1 * (sprintf "%.3f", $pgx->{cnvstatistics}->{$covK} / $pgx->{genomesize});
    } 
  }

  return $pgx;

}


################################################################################

sub interval_cnv_frequencies {

  no warnings 'uninitialized';
  
=pod

Expects:

Returns:

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my $pgx       =   shift;
  my $cnvmaps   =   shift;
  my $name      =   shift;
  my $labels    =   shift;
  my $maps      =   {
    intervals   =>  scalar(@{ $pgx->{genomeintervals} }),
    binning     =>  $pgx->{genomeintervals}->[0]->{end} - $pgx->{genomeintervals}->[0]->{start} + 1,
    name        =>  ($name =~ /\w/ ? $name : q{}),
    labels      =>  (@$labels > 0 ? $labels : []),
  };
  my %intLabs   =   (
    DUP         =>  'dupmap',
    DEL         =>  'delmap',
  );
  my %freqLabs      =   (
    DUP         =>  'dupfrequencies',
    DEL         =>  'delfrequencies',
  );

  # avoiding division by 0 errors if improperly called
  my $fFactor   =   100;
  if (@{ $cnvmaps } > 1) { $fFactor = 100 / @{ $cnvmaps } }

  foreach my $type (keys %intLabs) {
    for my $i (0..$#{ $pgx->{genomeintervals} }) {
      $maps->{ $freqLabs{ $type } }->[$i]   =   sprintf "%.3f", ($fFactor * ( grep{ $_->{ $intLabs{ $type } }->[$i] eq $type } @{ $cnvmaps } ));
    }
  }

  push(@{ $pgx->{frequencymaps} }, $maps);
    
  return $pgx;

}

################################################################################
########    utility subs    ####    ####    ####    ####    ####    ####    ####
################################################################################

1;
