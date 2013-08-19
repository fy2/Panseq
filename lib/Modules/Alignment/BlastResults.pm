#!/usr/bin/env perl
package Modules::Alignment::BlastResults;

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../../";
use IO::File;
use Log::Log4perl;
use Carp;
use Modules::Fasta::SequenceName;

sub new {
	my ($class) = shift;
	my $self = {};
	bless( $self, $class );
	$self->_initialize(@_);
	return $self;
}

=head3 logger

Stores a logger object for the module.

=cut

sub logger{
	my $self=shift;
	$self->{'_logger'} = shift // return $self->{'_logger'};
}


sub _outFH{
	my $self=shift;
	$self->{'__outFH'} = shift // return $self->{'__outFH'};
}

sub _storedQid{
	my $self=shift;
	$self->{'__storedQid'} = shift // return $self->{'__storedQid'};
}

sub _storedResults{
	my $self=shift;
	$self->{'__storedResults'} = shift // return $self->{'__storedResults'};
}



#methods
sub _initialize {
	my ($self) = shift;

	#logging
	$self->logger(Log::Log4perl->get_logger());
	$self->logger->debug("Logger initialized in Modules::Alignment::BlastResults");
	my $outFile = shift // $self->logger->logdie("No file sent to BlastResults");
	my $cutoff = shift // $self->logger->logdie("percentIdentityCutoff required");

	$self->percentIdentityCutoff($cutoff);
	$self->_outFH(IO::File->new('<' . $outFile)) // $self->logger->logdie("$!");
	$self->_storedResults({});
	$self->_storedQid(undef);
}

sub percentIdentityCutoff{
	my $self=shift;
	$self->{'_percentIdentityCutoff'} = shift // return $self->{'_percentIdentityCutoff'};
}

sub getNextResult{
	my $self=shift;

	my $results=$self->_storedResults();
	my $qid=$self->_storedQid();
	$self->logger->info("getNextResult");

	my $counter=0;
	while(my $line = $self->_outFH->getline()){
		$counter++;
		$line =~ s/\R//g;
		my @la = split("\t",$line);

		if(defined $qid && ($qid ne $la[1])){
			#move back a line in the file
			$self->logger->info("qid defined as $qid la[1] as $la[1]");
			#seek($self->_outFH,$self->_fhPosition,0);
			#last;
			$self->_storedQid($la[1]);
			return $results;
		}

		#'outfmt'=>'"6 
		# [0]sseqid 
		# [1]qseqid 
		# [2]sstart 
		# [3]send 
		# [4]qstart 
		# [5]qend 
		# [6]slen 
		# [7]qlen 
		# [8]pident 
		# [9]length"',
		# [10]sseq,
		# [11]qseq
		my $sName = Modules::Fasta::SequenceName->new($la[0]);
		my $sNameName= $sName->name;

		unless(defined $results->{$sNameName}){
			
			unless($self->_isOverCutoff($la[7],$la[8],$la[4],$la[5])){
				$self->logger->info($sNameName . " is not over cutoff");
				next;
			}
			

			my %lineResults=(
				qseqid=>$la[1],
				sseqid=>$la[0],
				sstart=>$la[2],
				send=>$la[3],
				qstart=>$la[4],
				qend=>$la[5],
				pident=>$la[8],
				qlen=>$la[7],
				slen=>$la[6],
				sseq=>$la[10] // undef,
				qseq=>$la[11] // undef
			);
			$results->{$sNameName}=\%lineResults;
			$self->_storedResults({});
			$self->_storedResults->{$sNameName}=\%lineResults;
		}	
		#$self->logger->info($sName->name . " is defined");	
		$qid = $la[1];
	}
	continue{
		#$self->_fhPosition(tell $self->_outFH);
		 if(eof){
		# 	$self->logger->info("EOF");
			return $results;
		# 	#return \%results;
		}
	}
	$self->logger->info("Counter before return: $counter qid: $qid");
	return 0;
}

sub _isOverCutoff{
	my $self=shift;
	my $qlen =shift // $self->logger->logdie("Missing qlen");
	my $pident =shift // $self->logger->logdie("Missing pident");
	my $qstart = shift // $self->logger->logdie("Missing qstart");
	my $qend =shift // $self->logger->logdie("Missing qend");

	my $queryHitLength = abs($qend - $qstart) +1;
	my $percentId = $queryHitLength / $qlen * $pident;
	 # $self->logger->info("queryHitLength: $queryHitLength");
	 # $self->logger->info("qlen: $qlen");
	 # $self->logger->info("pident: $pident");
	 # $self->logger->info("percentId: $percentId");
	if($percentId >= $self->percentIdentityCutoff){
		return 1;
	}
	else{
		return 0;
	}
}

1;
