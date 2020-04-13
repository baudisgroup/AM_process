#!/usr/bin/perl -w

$| = 1;

# CPAN packages
use Data::Dumper;
use File::Basename;
use strict;
use JSON::XS;

# local packages
use lib './PGX';
use PGX;

my $plotargs    =   {};

# command line input

my %args        =   @ARGV;

$args{-f}      									||= q{};
$args{-callset_id}							||=	q{};
$args{-biosample_id}						||=	q{};
$args{-outdir}  								||= q{};
$args{-variants_file_name}			||=	'variants.json';
$args{-callset_file_name}				||=	'callset.json';
$args{-cna_gain_threshold}			||=	q{};
$args{-cna_loss_threshold}			||=	q{};
$args{-plot_adjust_baseline}		||=	q{};
$args{-json_pretty}							||=	0;		# 1 for nice multi-line JSON

foreach (grep{ /^\-\w/ } keys %args) {
	if ($args{$_} =~ /\w/) { $plotargs->{$_} = $args{$_} }
}

check_command(\%args);

if (! -d $args{-outdir}) {
	$args{-outdir}	=	 dirname($args{-f}) }


#mkdir $args{'-outdir'};
my $varFile			=		$args{-outdir}.'/'.$args{-variants_file_name};
my $csFile			=		$args{-outdir}.'/'.$args{-callset_file_name};
my @pathParts		=		split('/', $args{-f});

my $pgx         =   new PGX($plotargs);
$pgx->pgx_add_segments_from_file($args{'-f'});

if ($args{-callset_id} !~ /^\w[\w\-\:]+?\w$/) {
	$args{-callset_id}	=		'pgxcs::'.$pathParts[-3].'::'.$pathParts[-2];
	print "\nusing inferred callset_id $args{-callset_id}\n";
}
if ($args{-biosample_id} !~ /^\w[\w\-\:]+?\w$/) {
	$args{-biosample_id}	=		'PGX_AM_BS_'.$pathParts[-2];
	# $args{-biosample_id}	=		'pgxcs::'.$pathParts[-3].'::'.$pathParts[-2];
	print "\nusing inferred biosample_id $args{-biosample_id}\n";
}

# renaming the ids ...
foreach my $id_type (qw(callset_id biosample_id)) {
	if ($args{'-'.$id_type} =~ /^\w[\w\-\:]+?\w$/) {
		for my $i (0..$#{ $pgx->{segmentdata} }) {
			$pgx->{segmentdata}->[$i]->{$id_type}	=		$args{'-'.$id_type};
		}
	}
}


my $varNo				=		scalar @{$pgx->{segmentdata}};
open (FILE, ">:utf8", $varFile) || warn 'output file '.$varFile.' could not be opened';
print  FILE  JSON::XS->new->pretty( $args{-json_pretty} )->canonical()->allow_blessed->convert_blessed->encode( $pgx->{segmentdata} )."\n";
close FILE;

print "\nWrote $varNo variants to $varFile\n";

$pgx->segments_add_statusmaps();

my $callset			=		{
	callset_id		=>	$args{-callset_id},
	biosample_id	=>	$args{-biosample_id},
	info					=>	{
		cna_gain_threshold	=>	1 * $pgx->{parameters}->{cna_gain_threshold},
		cna_loss_threshold	=>	1 * $pgx->{parameters}->{cna_loss_threshold},
		plot_adjust_baseline		=>	1 * $pgx->{parameters}->{plot_adjust_baseline},
		statusmaps	=>	$pgx->{statusmaps},
		cnvstatistics	=>	$pgx->{cnvstatistics},
	}
};

open (FILE, ">:utf8", $csFile) || warn 'output file '.$csFile.' could not be opened';
print  FILE  JSON::XS->new->pretty( $args{-json_pretty} )->canonical()->allow_blessed->convert_blessed->encode( $callset )."\n";
close FILE;

print "Wrote callset to $csFile\n\n";

################################################################################
################################################################################
################################################################################

sub check_command {

	my $args			=		shift;

	if (! -f $args->{'-f'}) {
 	 print "
No input file was specified:

  -f _path_to_my_file_/segments,cn.tsv

Optional parameters:\n\t".join("\n\t", sort keys %$args)."

More can be found in ./PGX/rsrc/plotdefaults.yaml (prefix with '-')\n\n";

	exit;
	
}
	
	
}


