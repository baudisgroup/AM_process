package PGX;

use File::Basename;
use YAML::XS qw(LoadFile DumpFile);

use lib '.';
use CytobandReader;
use GenomeIntervals;
use IntervalStatistics;
use Genomeplot;
use PlotParameters;
use HistoPlotter;
use ArrayPlotter;
use StripPlotter;
use CytobandsPlotter;
use UtilityLibs;
use PGXfileReader;
use PGXfileWriter;
use PGXdataAggregation;

require Exporter;
@ISA = qw(Exporter);
@EXPORT = qw(new read_plot_defaults);

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

sub new {

	my $class = shift;
	my $args = shift;
	$args = args_modify_plot_parameters(read_plot_defaults(), $args);
	my $self = {
		parameters => $args,
		config => read_config(),
		cytobands => read_cytobands($args->{genome}),
		datasetid => q{},
		plotid => $args->{plotid},
		svg => q{},
		Y => 0,
		errors => [],
	};

	bless $self, $class;

	$self->{genomeintervals} = make_genome_intervals(
		$self->{cytobands},
		$self->{parameters}->{binning},
	);
	$self->{referencebounds} = get_reference_base_limits($self->{cytobands});
	$self->{genomesize} =  get_genome_basecount(
		$self->{cytobands},
		$self->{parameters}->{chr2plot},
	);
	$self->{matrixindex}  = [ 0..$#{ $self->{genomeintervals} } ];
	$self = pgx_get_genome_regions($self);

	return $self;

}

################################################################################

sub read_config {

  my $path_of_this_module = File::Basename::dirname( eval { ( caller() )[1] } );
  my $config = LoadFile($path_of_this_module.'/rsrc/config.yaml');
  return  $config;

}

################################################################################

sub read_plot_defaults {

  my $path_of_this_module = File::Basename::dirname( eval { ( caller() )[1] } );
  my $plotPars = LoadFile($path_of_this_module.'/rsrc/plotdefaults.yaml');
  return  $plotPars;

}

########    ####    ####    ####    ####    ####    ####    ####    ####    ####

1;
