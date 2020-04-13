#!/usr/bin/perl

$| = 1;

# CPAN packages
use Data::Dumper;
use File::Basename;
use MongoDB;
use MongoDB::Database;
use MongoDB::Collection;
use MongoDB::MongoClient;
use POSIX 'strftime';
use strict;
use Term::ProgressBar;
use YAML::XS qw(LoadFile);

=podmd

This script queries the content of the different databases and collections and
stores the values (counts or counts of distinct values, according to the
configuration file `./PGXDB/rsrc/config/config.yaml`) in the
"progenetix.dbstats" collection.

For each of the database entries the current date (YYY-MM-DD) is used as the 
`_id`. Therefore, not more than one entry exists per day; each new run on the 
same day overwrites the previous values.

#### Parameters:

* `-filter`
    - if a string is provided here, it will be matched against the database 
    names and only (partially) matching databases will be processed
    - in that case, only the processed databases are updated
    - if no database entry for the day exists, the "filter" value will be 
    ignored and a complete set will be created

=cut

my $here_path   =   File::Basename::dirname( eval { ( caller() )[1] } );
my $config      =   LoadFile($here_path.'/PGXDB/rsrc/config/config.yaml') or die print 'Content-type: text'."\n\n¡No config.yaml file in this path!";
bless $config;

my %args        =   @ARGV;

my $today				=		strftime("%F", gmtime());
my $client		  =		MongoDB::MongoClient->new();

my @datasets    =   keys %{ $config->{databases} };
my @test_dbs    =   $client->database_names();
@test_dbs       =   grep{ /test/i } @test_dbs;
push(@datasets, @test_dbs);

my $dbstats			=		{ "date" => $today };
my $dbexists		=		$client->get_database('progenetix')->get_collection('dbstats')->find_one( { "date" => $today } );
if ($dbexists->{"date"} eq $today) { 
  if ($args{'-filter'} =~ /^\w\+?/) {
    @datasets   =   grep{ /$args{'-filter'}/ } @datasets }
  $dbstats      =   $dbexists;
}
elsif ($args{'-filter'} =~ /^\w\+?/) {
	print "\nThe filter will be ignored since no statistics exist yet for today.\n" }

################################################################################

foreach my $db (@datasets) {

	my $dbconn		=		$client->get_database($db);

# Distinct variants cannot be counted using MongoDB's `distinct` command, since
# the result would be too big. This has to be done with a dedicated loop.

	my %digests		=		();
  my $varcoll	  =		$dbconn->get_collection( "variants" );
	my $cursor	  =		$varcoll->find()->fields({ digest => 1, _id => 0 });
  my $varNo			=		$varcoll->estimated_document_count();     #estimated_document_count()

  my $progBar   =   Term::ProgressBar->new({name => $varNo.' '.$db.' variants', count => $varNo });
  $progBar->minor(0);
	my $i					=		0;
	
	while (my $var = $cursor->next) {
		$digests{ $var->{digest} }	+=	1;
    $progBar->update($i++);
	}
	
	$varNo        =   $i;
	
	$progBar->update($varNo);
	my $distNo		=		scalar keys %digests;
	print $varNo."\n";
	print $distNo."\n";
	
	if ($distNo < 1) { next }
	
	$dbstats->{$db.'__variants'}->{'count'}	                  =		$varNo;	
	$dbstats->{$db.'__variants'}->{'distincts_count_digest'}	=		$distNo;	
	
	print "\nIdentified $distNo distinct variants out of $varNo total from $db\n";

=podmd

=cut

	foreach my $coll (grep{! /variants/ } keys %{ $config->{databases}->{$db}->{collections} }) {

		$dbstats->{$db.'__'.$coll}->{count}	=		$dbconn->get_collection($coll)->count();
		
		foreach (keys %{ $config->{databases}->{$db}->{collections}->{$coll}->{distincts} }) {
			my $dist 	=		$config->{databases}->{$db}->{collections}->{$coll}->{distincts}->{$_};
			my $distinct  =   $dbconn->run_command([
													"distinct"  =>  $coll,
													"key"       =>  $dist->{attr},
													"query"     =>  {},
												]);
			$dbstats->{$db.'__'.$coll}->{'distincts_count_'.$_} =		grep{ /$dist->{match}/ } @{ $distinct->{values} };
														
}}}

# publication counts fields
my $cursor			=		MongoDB::MongoClient->new()->get_database('progenetix')->get_collection('publications')->find( $config->{databases}->{'progenetix'}->{collections}->{'publications'}->{query} );

while (my $pub = $cursor->next) {
	foreach my $count (keys %{ $pub->{counts} }) {
		$dbstats->{'progenetix__publications'}->{'sum_counts_'.$count}	+=	$pub->{counts}->{$count};
}}

print Dumper($dbstats);

my $statscoll		=		$client->get_database('progenetix')->get_collection('dbstats')->replace_one(
	{ "date" 			=>	$today },
	$dbstats,
	{ upsert			=>	'true' }
);

