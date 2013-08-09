#!/usr/bin/env perl

package Modules::Alignment::BlastRun;

use strict;
use warnings;
use Carp;
use FindBin;
use lib "$FindBin::Bin/../../";
use Parallel::ForkManager;
use Modules::Fasta::FastaFileSplitter;
use File::Copy;
use Log::Log4perl;

sub new {
	my ($class) = shift;
	my $self = {};
	bless( $self, $class );
	$self->_initialize(@_);
	return $self;
} 

sub blastDirectory{
	my $self=shift;
	$self->{'_programDirectory'}=shift // return $self->{'_programDirectory'};
}

sub task{
	my $self=shift;
	$self->{'_task'}=shift // return $self->{'_task'};
}

sub dust{
	my $self=shift;
	$self->{'_dust'}=shift // return $self->{'_dust'};
}

sub query{
	my $self=shift;
	$self->{'_query'}=shift // return $self->{'_query'};
}

sub db{
	my $self=shift;
	$self->{'_db'}=shift // return $self->{'_db'};
}

sub out{
	my $self=shift;
	$self->{'_out'}=shift // return $self->{'_out'};
}

sub evalue{
	my $self=shift;
	$self->{'_evalue'}=shift // return $self->{'_evalue'};
}

sub word_size{
	my $self=shift;
	$self->{'_word_size'}=shift // return $self->{'_word_size'};
}

sub outfmt{
	my $self=shift;
	$self->{'_outfmt'}=shift // return $self->{'_outfmt'};
}

sub num_descriptions{
	my $self=shift;
	$self->{'_num_descriptions'}=shift // return $self->{'_num_descriptions'};
}

sub num_alignments{
	my $self=shift;
	$self->{'_num_alignments'}=shift // return $self->{'_num_alignments'};
}

sub max_target_seqs{
	my $self=shift;
	$self->{'_max_target_seqs'}=shift // return $self->{'_max_target_seqs'};
}

sub num_threads{
	my $self=shift;
	$self->{'_num_threads'}=shift // return $self->{'_num_threads'};
}

sub best_hit_score_edge{
	my $self=shift;
	$self->{'_best_hit_score_edge'}=shift // return $self->{'_best_hit_score_edge'};
}

sub best_hit_overhang{
	my $self=shift;
	$self->{'_best_hit_overhang'}=shift // return $self->{'_best_hit_overhang'};
}

sub numberOfSplits{
	my $self=shift;
	$self->{'_numberOfSplits'}=shift // return $self->{'_numberOfSplits'};
}

sub logger{
	my $self=shift;
	$self->{'_logger'}=shift // return $self->{'_logger'};
}

sub outputXMLfiles{
	my $self=shift;
	$self->{'_outputXMLfiles'}=shift // return $self->{'_outputXMLfiles'};
}

sub splitFileDatabase{
	my $self=shift;
	$self->{'_splitFileDatabase'}=shift // return $self->{'_splitFileDatabase'};
}

sub _filesToRemove{
	my $self=shift;
	$self->{'__filesToRemove'}=shift // return $self->{'__filesToRemove'};
}


#methods
sub _initialize{
	my($self)=shift;
	
	#logging
	$self->logger(Log::Log4perl->get_logger());
	$self->logger->debug("Logger initialized in Modules::Alignment::BlastRun\n");

	#init values
	my %params = @_;

	foreach my $key(keys %params){
		if($self->can($key)){
			$self->$key($params{$key});
		}
		else{
			$self->logger->logconfess("$key is not a valid parameter in Modules::Alignment::BlastRun");
		}
	}

	unless(
		defined $self->db 
		&& defined $self->query
		&& defined $self->blastDirectory
		&& defined $self->task
	){
		$self->logger->logconfess("BlastRun requires:\n
			\t'query'\n
			\t'db'\n
			\t'blastDirectory'\n
			\t'task'\n
			to be defined. Missing one or more required parameters.
		");
	}	

	#init
	$self->outputXMLfiles([]);	
	$self->_filesToRemove([]);

	#defaults
	unless(defined $self->numberOfSplits){
		$self->numberOfSplits(1);
	}
}

=head3 run

Optionally splits a single query file into multiple files of approximately equal size.
Parallel::ForkManager is used, and each fork defines its own $self->query, $self->db, $self->out
parameter. Because these are children, it does nothing to the parent parameters.
$self->out can be a file or directory, and is used as the basis for creating each output file.

=cut

sub run{
	my($self)=shift;		

	my $splitter;
	$splitter = Modules::Fasta::FastaFileSplitter->new(
		'inputFile'=>$self->query,
		'databaseFile'=>$self->splitFileDatabase // undef,
		'numberOfSplits'=>$self->numberOfSplits // undef
	);


	if($self->numberOfSplits > 1){
		$splitter->splitFastaFile();
	}
	
	my $forker = Parallel::ForkManager->new($self->numberOfSplits);

	my $counter=0;
	foreach my $splitFile(@{$splitter->arrayOfSplitFiles}){
		$counter++;
		my $outputXMLname = $self->out . 'blastout' . $counter . '.out';
		push @{$self->outputXMLfiles}, $outputXMLname;
		$forker->start and next;
			$self->query($splitFile);
			$self->out($outputXMLname);
			my $systemLine = $self->_createBlastLine();
			$self->logger->info("Launching " . $self->task . " with: $systemLine");
			system($systemLine);
			exit(1);
			unlink $splitFile;
			$self->_removeTempFiles();
		$forker->finish;
	}
	$forker->wait_all_children();
}


sub _createBlastLine{
	my $self=shift;

	my $systemLine = $self->blastDirectory . $self->task;
	
	$systemLine .= ' -dust ' . $self->dust if (defined $self->dust);
	$systemLine .= ' -query ' . $self->query if (defined $self->query);
	$systemLine .= ' -db ' . $self->db if (defined $self->db);
	$systemLine .= ' -out ' . $self->out if (defined $self->out);
	$systemLine .= ' -evalue ' . $self->evalue if (defined $self->evalue);
	$systemLine .= ' -word_size ' . $self->word_size if (defined $self->word_size);
	$systemLine .= ' -outfmt ' . $self->outfmt if (defined $self->outfmt);
	$systemLine .= ' -num_descriptions ' . $self->num_descriptions if (defined $self->num_descriptions);
	$systemLine .= ' -num_alignments ' . $self->num_alignments if (defined $self->num_alignments);
	$systemLine .= ' -max_target_seqs ' . $self->max_target_seqs if (defined $self->max_target_seqs);
	$systemLine .= ' -num_threads ' . $self->num_threads if (defined $self->num_threads);
	$systemLine .= ' -best_hit_score_edge ' . $self->best_hit_score_edge if (defined $self->best_hit_score_edge);
	$systemLine .= ' -best_hit_overhang ' . $self->best_hit_overhang if (defined $self->best_hit_overhang);

	return $systemLine;
}

=head3 _removeTempFiles

Removes each of the files in $self->_filesToRemove
from the filesystem.

=cut

sub _removeTempFiles{
	my $self=shift;

	foreach my $file(@{$self->_filesToRemove}){
		$self->logger->debug("Removing $file");
		#unlink $file;
	}
}

1;