package PGX::ArrayPlotter;

use Data::Dumper;
use GD::Simple;
use MIME::Base64 qw(encode_base64);
use YAML::XS qw(LoadFile DumpFile);
use PGX::CytobandsPlotter;
# use PGX::GenomePlots::PlotParameters;

require Exporter;
@ISA    = qw(Exporter);
@EXPORT = qw(
  return_arrayplot_svg
  get_arrayplot_area
  args_modify_plot_parameters
  hex2rgb);

################################################################################

sub return_arrayplot_svg {

    our $pgx = shift;

    $pgx->{Y} = $pgx->{parameters}->{size_plotmargin_top_px};
    my $plotW = $pgx->{parameters}->{size_plotimage_w_px};
    $pgx->{areastartx} = $pgx->{parameters}->{size_plotmargin_px} +
      $pgx->{parameters}->{size_title_left_px};
    $pgx->{areawidth} = $plotW -
      ( $pgx->{areastartx} + $pgx->{parameters}->{size_plotmargin_px} );
    if ( $pgx->{parameters}->{do_chromosomes_proportional} =~ /y/i
        && @{ $pgx->{parameters}->{chr2plot} } == 1 )
    {
        $pgx->{areawidth} *=
          ( $pgx->{referencebounds}->{ $pgx->{parameters}->{chr2plot}->[0] }
              ->[1] / $pgx->{referencebounds}->{'1'}->[1] );
        $plotW =
          $pgx->{areawidth} + 2 * $pgx->{parameters}->{size_plotmargin_px};
    }
    $pgx->{basepixfrac} = (
        $pgx->{areawidth} - (
            $#{ $pgx->{parameters}->{chr2plot} } *
              $pgx->{parameters}->{size_chromosome_padding_px}
        )
      ) /
      $pgx->{genomesize};
    $pgx               = svg_add_title($pgx);
    $pgx               = svg_add_cytobands($pgx);
    $pgx->{areastarty} = $pgx->{Y};
    $pgx               = get_arrayplot_area($pgx);
    $pgx               = get_fracb_area($pgx);
    $pgx               = svg_add_markers($pgx);
    $pgx               = svg_add_bottom_text($pgx);
    $pgx->{Y} += $pgx->{parameters}->{size_plotmargin_bottom_px};
    my $plotH = sprintf "%.0f", $pgx->{Y};
    $plotW = sprintf "%.0f", $plotW;
    $pgx->{svg} = '<svg
xmlns="http://www.w3.org/2000/svg"
xmlns:xlink="http://www.w3.org/1999/xlink"
version="1.1"
id="' . $pgx->{plotid} . '"
width="' . $plotW . 'px"
height="' . $plotH . 'px"
style="margin: auto; font-family: Helvetica, sans-serif;">

<rect x="0" y="0" width="'
      . $plotW
      . '" height="'
      . $plotH
      . '" style="fill: '
      . $pgx->{parameters}->{color_plotbackground_hex} . '; " />

' . $pgx->{svg} . '
</svg>';

    return $pgx;

}

################################################################################

sub get_arrayplot_area {

=podmd

* Expects
    - the plot object with current Y parameter for placing the plot elements on the SVG

* Returns
    - the  extended SVG and the increased end Y value as start for the next elements

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

    if ( $pgx->{parameters}->{size_plotarea_h_px} < 1 ) { return $pgx }

    $pgx->{svg} .= '
<style type="text/css"><![CDATA[
  .DUP {stroke-width: '
      . $pgx->{parameters}->{size_segments_stroke_px}
      . 'px; stroke: '
      . $pgx->{parameters}->{color_var_dup_hex}
      . '; opacity: 0.8  }
  .DEL {stroke-width: '
      . $pgx->{parameters}->{size_segments_stroke_px}
      . 'px; stroke: '
      . $pgx->{parameters}->{color_var_del_hex}
      . '; opacity: 0.8  }
  .cen {stroke-width: '
      . $pgx->{parameters}->{size_centerline_stroke_px}
      . 'px; stroke: '
      . $pgx->{parameters}->{color_plotgrid_hex}
      . '; opacity: 0.8 ; }
]]></style>';

    $pgx->{Y} += $pgx->{parameters}->{size_plotarea_padding};

    svg_add_labels_y($pgx);

    my $area_x0   = $pgx->{areastartx};
    my $area_y0   = $pgx->{Y};
    my $area_ycen = $pgx->{Y} + $pgx->{parameters}->{size_plotarea_h_px} / 2;
    my $area_yn   = $pgx->{Y} + $pgx->{parameters}->{size_plotarea_h_px};
    my $lowSegY   = $area_yn + $pgx->{parameters}->{size_segments_stroke_px};

    $pgx->{areaendy} = $area_yn;

    my $probeSize = $pgx->{parameters}->{factor_probedots};
    if ( scalar @{ $pgx->{probedata} } < 200000 ) { $probeSize *= 2 }
    if ( scalar @{ $pgx->{probedata} } < 8000 )   { $probeSize *= 2 }

    my $segBLcorr = $pgx->{parameters}->{plot_adjust_baseline};
    if ( $pgx->{parameters}->{segbaseline} ne 0 ) {
        $segBLcorr = 0;
    }

    my $probeBLcorr = $pgx->{parameters}->{plot_adjust_baseline};
    if ( $pgx->{parameters}->{probebaseline} ne 0 ) {
        $probeBLcorr = 0;
    }

  # probe area, probes & segments   ####    ####    ####    ####    ####    ####

    foreach my $refName ( @{ $pgx->{parameters}->{chr2plot} } ) {

        my $areaW = sprintf "%.1f",
          ( $pgx->{referencebounds}->{$refName}->[1] -
              $pgx->{referencebounds}->{$refName}->[0] ) *
          $pgx->{basepixfrac};

        $pgx->{svg} .= '
<rect x="'
          . $area_x0 . '" y="'
          . $pgx->{Y}
          . '" width="'
          . $areaW
          . '" height="'
          . $pgx->{parameters}->{size_plotarea_h_px}
          . '" style="fill: '
          . $pgx->{parameters}->{color_plotarea_hex}
          . '; fill-opacity: 0.8; " />';

    # probes ###    ####    ####    ####    ####    ####    ####    ####    ####
        my $areaProbes =
          [ grep { $_->{reference_name} eq $refName } @{ $pgx->{probedata} } ];
        $areaProbes =
          [ grep { $_->{position} <= $pgx->{referencebounds}->{$refName}->[1] }
              @$areaProbes ];
        $areaProbes =
          [ grep { $_->{position} >= $pgx->{referencebounds}->{$refName}->[0] }
              @$areaProbes ];

        # probes are plotted using GD
        my $probeArea =
          GD::Image->new( $areaW, $pgx->{parameters}->{size_plotarea_h_px}, 1 );
        my $gdDotS    = 1 * $probeSize;
        my $gdAreaCol = $probeArea->colorAllocate(
            @{ hex2rgb( $pgx->{parameters}->{color_plotarea_hex} ) } );
        $probeArea->transparent($gdAreaCol);
        my $gdDotcol = $probeArea->colorAllocateAlpha( 32, 32, 32, 63 );
        $probeArea->filledRectangle( 0, 0, $areaW,
            $pgx->{parameters}->{size_plotarea_h_px}, $gdAreaCol );

        foreach (@$areaProbes) {
            my $dotX = sprintf "%.2f", $pgx->{basepixfrac} *
              ( $_->{position} - $pgx->{referencebounds}->{$refName}->[0] );
            my $dotY = sprintf "%.2f",
              ( $pgx->{parameters}->{size_plotarea_h_px} / 2 -
                  ( $_->{value} + $probeBLcorr ) *
                  $pgx->{parameters}->{pixyfactor} );
            $probeArea->filledEllipse( $dotX, $dotY, $gdDotS, $gdDotS,
                $gdDotcol );
        }

        # the GD object is encoded as base64 and embedded into the svg
        $pgx->{svg} .= '
<image
  x="' . $area_x0 . '"
  y="' . $pgx->{Y} . '"
  width="' . $areaW . '"
  height="' . $pgx->{parameters}->{size_plotarea_h_px} . '"
  xlink:href="data:image/png;base64,' . encode_base64( $probeArea->png ) . '"
/>';

    # / probes #    ####    ####    ####    ####    ####    ####    ####    ####

    # segments #    ####    ####    ####    ####    ####    ####    ####    ####
#         my $areaSegments = [ grep { $_->{reference_name} eq $refName }
#               @{ $pgx->{samples}->[0]->{variants} } ];
        my $areaSegments = [ grep { $_->{reference_name} eq $refName }
              @{ $pgx->{segmentdata} } ];
        $areaSegments = [ grep { $_->{variant_type} =~ /\w/ } @$areaSegments ];
        $areaSegments = [
            grep {
                $_->{start}->[0] <= $pgx->{referencebounds}->{$refName}->[1]
            } @$areaSegments
        ];
        $areaSegments = [
            grep { $_->{end}->[-1] >= $pgx->{referencebounds}->{$refName}->[0] }
              @$areaSegments
        ];
        foreach my $seg (@$areaSegments) {
            my $start = $seg->{start}->[0];
            my $end   = $seg->{end}->[-1];
            if ( $start < $pgx->{referencebounds}->{$refName}->[0] ) {
                $start = $pgx->{referencebounds}->{$refName}->[0];
            }
            if ( $end > $pgx->{referencebounds}->{$refName}->[1] ) {
                $end = $pgx->{referencebounds}->{$refName}->[1];
            }

            # providing a minimum sub-pixel segment plot length
            my $segPixLen = sprintf "%.1f",
              ( $pgx->{basepixfrac} * ( $end - $start ) );
            if ( $segPixLen < 0.2 ) { $segPixLen = 0.2 }

            #print Dumper($seg);
            my $seg_x0 = sprintf "%.1f",
              $area_x0 +
              $pgx->{basepixfrac} *
              ( $start - $pgx->{referencebounds}->{$refName}->[0] ) ;
            my $seg_xn = $seg_x0 + $segPixLen;
            my $seg_y  = sprintf "%.1f",
              $area_ycen -
              ( $seg->{info}->{value} + $segBLcorr ) *
              $pgx->{parameters}->{pixyfactor};
            $pgx->{svg} .= '
<line x1="'
              . $seg_x0
              . '"  y1="'
              . $seg_y
              . '"  x2="'
              . $seg_xn
              . '"  y2="'
              . $seg_y
              . '"  class="'
              . $seg->{variant_type} . '"  />';
            if ( scalar @{ $pgx->{parameters}->{chr2plot} } == 1 ) {
                $pgx->{svg} .= '
<a
xlink:href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='
                  . $pgx->{parameters}->{genome}
                  . '&amp;position=chr'
                  . $refName . '%3A'
                  . $start . '-'
                  . $end . '"
xlink:show="new"
xlink:title="'
                  . ( $seg->{info}->{value} + $segBLcorr ) . ' at '
                  . $refName . ':'
                  . $start . '-'
                  . $end . '">';
            }
            $pgx->{svg} .= '
<line x1="'
              . $seg_x0
              . '"  y1="'
              . $lowSegY
              . '"  x2="'
              . $seg_xn
              . '"  y2="'
              . $lowSegY
              . '"  class="'
              . $seg->{variant_type} . '"  />';
            if ( scalar @{ $pgx->{parameters}->{chr2plot} } == 1 ) {
                $pgx->{svg} .= '</a>';
            }

        }

    # / segments    ####    ####    ####    ####    ####    ####    ####    ####

        # moving x to the next chromosome area
        $area_x0 += $areaW + $pgx->{parameters}->{size_chromosome_padding_px};

    }

    # adding a baseline at 0
    $pgx->{svg} .= '
<line x1="'
      . $pgx->{parameters}->{size_plotmargin_px}
      . '"  y1="'
      . $area_ycen
      . '"  x2="'
      . ( $pgx->{parameters}->{size_plotmargin_px} + $pgx->{areawidth} )
      . '"  y2="'
      . $area_ycen
      . '"  class="cen"  />';

    $pgx->{Y} = $lowSegY + $pgx->{parameters}->{size_segments_stroke_px} / 2;

    return $pgx;

}

################################################################################

sub get_fracb_area {
=pod

Expects:
  - the plot object with current Y parameter for placing the plot elements on
    the SVG

Returns:
  - the  extended SVG and the increased end Y value as start for
    the next elements

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

    if (   @{ $pgx->{probedata_fracb} } < 1
        && @{ $pgx->{segmentdata_fracb} } < 1 )
    {
        return $pgx;
    }

    if ( $pgx->{parameters}->{size_fracbarea_h_px} < 1 ) { return $pgx }

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

    $pgx->{Y} += $pgx->{parameters}->{size_chromosome_padding_px};

    $pgx->{svg} .= '
<style type="text/css"><![CDATA[
  .fb {stroke-width: '
      . $pgx->{parameters}->{size_segments_stroke_px}
      . 'px; stroke: #FF3333; opacity: 0.8  }
]]></style>';

    # area init; here, $area_y0 is the lower margin, corresponding then to
    # a 0-value
    my $area_x0   = $pgx->{areastartx};
    my $area_yn   = $pgx->{Y} + $pgx->{parameters}->{size_fracbarea_h_px};
    my $area_ycen = $pgx->{Y} + $pgx->{parameters}->{size_fracbarea_h_px} / 2;
    my $fbPixYfac = 1 / $pgx->{parameters}->{size_fracbarea_h_px};

    $pgx->{areaendy} = $area_yn;

    my $probeSize = $pgx->{parameters}->{factor_probedots};
    if ( scalar @{ $pgx->{probedata_fracb} } < 200000 ) { $probeSize *= 2 }
    if ( scalar @{ $pgx->{probedata_fracb} } < 8000 )   { $probeSize *= 2 }

    foreach my $refName ( @{ $pgx->{parameters}->{chr2plot} } ) {

        my $areaW = sprintf "%.1f",
          ( $pgx->{referencebounds}->{$refName}->[1] -
              $pgx->{referencebounds}->{$refName}->[0] ) *
          $pgx->{basepixfrac};

        $pgx->{svg} .= '
<rect x="'
          . $area_x0 . '" y="'
          . $pgx->{Y}
          . '" width="'
          . $areaW
          . '" height="'
          . $pgx->{parameters}->{size_fracbarea_h_px}
          . '" style="fill: '
          . $pgx->{parameters}->{color_plotarea_hex}
          . '; fill-opacity: 0.8; " />';

    # probes ###    ####    ####    ####    ####    ####    ####    ####    ####

        my $areaProbes = [ grep { $_->{reference_name} eq $refName }
              @{ $pgx->{probedata_fracb} } ];
        $areaProbes =
          [ grep { $_->{position} <= $pgx->{referencebounds}->{$refName}->[1] }
              @$areaProbes ];
        $areaProbes =
          [ grep { $_->{position} >= $pgx->{referencebounds}->{$refName}->[0] }
              @$areaProbes ];

        # probes are plotted using GD

        my $probeArea =
          GD::Image->new( $areaW, $pgx->{parameters}->{size_fracbarea_h_px},
            1 );
        my $gdDotS    = 1 * $probeSize;
        my $gdAreaCol = $probeArea->colorAllocate(
            @{ hex2rgb( $pgx->{parameters}->{color_plotarea_hex} ) } );
        $probeArea->transparent($gdAreaCol);
        my $gdDotcol = $probeArea->colorAllocateAlpha( 0, 0, 0, 111 );
        $probeArea->filledRectangle( 0, 0, $areaW,
            $pgx->{parameters}->{size_fracbarea_h_px}, $gdAreaCol );

        foreach (@$areaProbes) {
            my $dotX = sprintf "%.2f", $pgx->{basepixfrac} *
              ( $_->{position} - $pgx->{referencebounds}->{$refName}->[0] );
            my $dotY = sprintf "%.2f",
              ( $pgx->{parameters}->{size_fracbarea_h_px} -
                  $_->{value} / $fbPixYfac );
            $probeArea->filledEllipse( $dotX, $dotY, $gdDotS, $gdDotS,
                $gdDotcol );
        }

        # the GD object is encoded as base64 and embedded into the svg
        $pgx->{svg} .= '
<image
  x="' . $area_x0 . '"
  y="' . $pgx->{Y} . '"
  width="' . $areaW . '"
  height="' . $pgx->{parameters}->{size_fracbarea_h_px} . '"
  xlink:href="data:image/png;base64,' . encode_base64( $probeArea->png ) . '"
/>';

    # / probes #    ####    ####    ####    ####    ####    ####    ####    ####

   # fracb segments s###    ####    ####    ####    ####    ####    ####    ####
        my $areaSegments = [ grep { $_->{reference_name} eq $refName }
              @{ $pgx->{segmentdata_fracb} } ];
        $areaSegments =
          [ grep { $_->{start}->[0] <= $pgx->{referencebounds}->{$refName}->[1] }
              @$areaSegments ];
        $areaSegments =
          [ grep { $_->{end}->[-1] >= $pgx->{referencebounds}->{$refName}->[0] }
              @$areaSegments ];
        
        foreach my $seg (@$areaSegments) {
            my $start = $seg->{start}->[0];
            my $end   = $seg->{end}->[-1];
            if ( $start < $pgx->{referencebounds}->{$refName}->[0] ) {
                $start = $pgx->{referencebounds}->{$refName}->[0];
            }
            if ( $end > $pgx->{referencebounds}->{$refName}->[1] ) {
                $end = $pgx->{referencebounds}->{$refName}->[1];
            }

            # providing a minimum sub-pixel segment plot length
            my $segPixLen = sprintf "%.1f",
              ( $pgx->{basepixfrac} * ( $end - $start ) );
            if ( $segPixLen < 0.2 ) { $segPixLen = 0.2 }

            my $seg_x0 = sprintf "%.1f",
              $area_x0 +
              $pgx->{basepixfrac} *
              ( $start - $pgx->{referencebounds}->{$refName}->[0] );
            my $seg_xn = $seg_x0 + $segPixLen;

            my $val   = $seg->{info}->{value} / $fbPixYfac;

            my $seg_y = sprintf "%.2f", ( $area_yn - $val );

            $pgx->{svg} .= '
<line x1="'
              . $seg_x0
              . '" y1="'
              . $seg_y
              . '" x2="'
              . $seg_xn
              . '" y2="'
              . $seg_y
              . '" class="fb" />';

        }

    # / segments    ####    ####    ####    ####    ####    ####    ####    ####

        $area_x0 += $areaW + $pgx->{parameters}->{size_chromosome_padding_px};

    }

    # adding a baseline at 0.5
    $pgx->{svg} .= '
<line x1="'
      . $pgx->{parameters}->{size_plotmargin_px}
      . '"  y1="'
      . $area_ycen
      . '"  x2="'
      . ( $pgx->{parameters}->{size_plotmargin_px} + $pgx->{areawidth} )
      . '"  y2="'
      . $area_ycen
      . '"  class="cen"  />';

    $pgx->{Y} += $pgx->{parameters}->{size_fracbarea_h_px};

    return $pgx;

}

################################################################################

sub args_modify_plot_parameters {

  no warnings 'uninitialized';

=pod

Expects:
  - the dash prefixed input args
  - the plot parameter object

Returns:
  - modified plot parameter object

=cut

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

  my ($plotPars, $args) =   @_;

  # local defaults overwrite general values; but command line parameters
  # take precedence later on
  my $locDefaults       =   {};
  my $defaultsDir       =   '';

  if (defined($args->{'-defaultsfile'})) {
    $defaultsDir        = $args->{'-defaultsfile'};
    $defaultsDir        =~  s/\/[\w\.\,]+?$//;
    if (-f $args->{'-defaultsfile'}) {
      $locDefaults        =   LoadFile($args->{'-defaultsfile'});
      foreach my $par (keys %$plotPars) {
        if ($locDefaults->{$par}) {
          $plotPars->{$par}   =   $locDefaults->{$par};
        }
      }
    }
  }

  # the -plottype | -colorschema specific mappings are processed first & removed
  # thereafter (there are still fallbacks if no parameters given)
  if (defined($args->{'-colorschema'}) && grep{ $args->{'-colorschema'} eq $_ } keys %{ $plotPars->{colorschemas} }) {
    my $colorschema     =   $args->{'-colorschema'};
    foreach (keys %{ $plotPars->{colorschemas}->{ $colorschema } }) {
      if ($plotPars->{colorschemas}->{ $colorschema }->{$_} =~ /^(\#\w{6})$/) {
        $plotPars->{$_} =   $1 }
    }
    delete $plotPars->{colorschemas};
    delete $args->{'-colorschema'};
  }

  # adjusting arguments for the selected plot type
  if (grep{ $args->{'-plottype'} eq $_ } keys %{ $plotPars->{plottype_values} }) {
    foreach (keys %{ $plotPars->{plottype_values}->{ $args->{'-plottype'} } }) {
      $plotPars->{$_} =   $plotPars->{plottype_values}->{ $args->{'-plottype'} }->{$_} }
    delete $plotPars->{plottype_values};
  }

  # arguments to parameters
  foreach my $par (keys %$plotPars) {

    if (! defined($args->{'-'.$par}) || $args->{'-'.$par} !~ /\w/) { next }
    # special evaluation: regions
    if ($par eq 'plotregions') {

      if (ref $args->{'-'.$par} eq 'ARRAY') {      
        $args->{'-'.$par} =   join(',', @{ $args->{'-'.$par} }) }

      foreach my $plotregion (split(',', $args->{'-plotregions'})) {

        if ($plotregion =~ /^(?:chro?)?(\w\d?)\:(\d+?)\-(\d+?)$/) {
          my $plotR =   {
            reference_name  =>  $1,
            start           =>  $2,
            end             =>  $3,
          };
          push(@{ $plotPars->{'plotregions'} }, $plotR);
    }}}

    # special evaluation: markers
    elsif ($par eq 'markers') {

      if (ref $args->{'-markers'} eq 'ARRAY') {      
        $args->{'-markers'} =   join(',', @{ $args->{'-markers'} }) }
      foreach (split(',', $args->{'-markers'})) {

        my @markervals  =   split(':', $_);

        if (
          $markervals[0]  =~ /^(chro?)?([\dxy]\d?)$/i
          &&
          $markervals[1]  =~ /^\d+?\-\d+?$/
        ) {
          my $mark     =   { reference_name  =>  $markervals[0] };
          $mark->{reference_name} =~  s/[^xy\d]//gi;
          ($mark->{start}, $mark->{end})  =   split('-', $markervals[1]);
          if ($markervals[2] =~ /^\w[\w \-\(\)\[\]]+?$/) {
            $mark->{label}  =  $markervals[2] }
          if ($markervals[3] =~ /^\#\w\w\w(\w\w\w)?$/) {
            $mark->{color}  =  $markervals[3] }
          if ($mark->{color} !~ /^\#\w\w\w(?:\w\w\w)?$/) {
            $mark->{color}  =   random_hexcolor() }
          push(@{ $plotPars->{'markers'} }, $mark);
    }}}
    # list style parameters are provided comma concatenated => deparsed
    elsif (grep{ $par eq $_ } qw(chr2plot label_y_m)) {
      if (ref $args->{'-'.$par} eq 'ARRAY') {      
        $args->{'-'.$par} =   join(',', @{ $args->{'-'.$par} }) }

      $plotPars->{$par}   =   [ split(',', $args->{'-'.$par}) ] }
    elsif (
      ($par =~/^color/ || $par =~/color$/)
      &&
      $args->{'-'.$par} =~  /^\w{6}$/
    ) {
      $plotPars->{$par}   =   '#'.$args->{'-'.$par} }
    else {
      $plotPars->{$par}   =   $args->{'-'.$par} }
      
    if ($par eq 'title') {
      $plotPars->{$par} =~  s/arraymap/arrayMap/i;
      $plotPars->{$par} =~  s/progenetix/Progenetix/i;
    }

  }

  # derived
  $plotPars->{pixyfactor}   =   1 * $plotPars->{size_plotarea_h_px} / (2 * $plotPars->{value_plot_y_max});
  foreach my $override (keys %$locDefaults) {
    if (! grep{ $_ eq $override } @{ $plotPars->{local_overrides} }) {
      delete $locDefaults->{$override};
    }
  }

  if (-d $defaultsDir) {
    DumpFile($args->{'-defaultsfile'}, $locDefaults) }

  return $plotPars;

}

################################################################################

sub hex2rgb {

    my ($r, $g, $b)     =   $_[0] =~  m/^\#?(\w{2})(\w{2})(\w{2})$/;

    return [ CORE::hex($r), CORE::hex($g), CORE::hex($b) ];

}
1;
