#!/usr/bin/env perl

#written by Chad Laing; last updated April 05, 2013

#Pan-genome sequence analysis using Panseq: an online tool for the rapid analysis of core and accessory genomic regions
#Chad Laing, Cody Buchanan, Eduardo Taboada, Yongxiang Zhang, Andrew Kropinski, Andre Villegas, James E Thomas and Victor PJ Gannon
#BMC Bioinformatics 2010, 11:461
#http://www.biomedcentral.com/1471-2105/11/461

package Modules::LociSelector::LociSelector;

#usage
# my $obj = $LociSelector->new();
# $obj->getBestLoci(<file name (mandatory)>,<loci number or 'best' (mandatory)>,<output filehandle (optional, STDOUT defualt)>);

use FindBin;
use lib "$FindBin::Bin/../../";
use IO::File;
use Log::Log4perl;

#object creation
sub new {
    my ($class) = shift;
    my $self = {};
    bless( $self, $class );
    $self->_initialize(@_);
    return $self;
}

sub _initialize{
    my($self)=shift;
    my %params=@_;

    #logging
    $self->logger(Log::Log4perl->get_logger());    
    $self->logger->info("Logger initialized in Modules::LociSelector::LociSelector");  

   #on object construction set all parameters
    foreach my $key(keys %params){
        if($self->can($key)){
            $self->$key($params{$key});
        }
        else{
            #logconfess calls the confess of Carp package, as well as logging to Log4perl
            $self->logger->fatal("$key is not a valid parameter in Modules::LociSelector::LociSelector");
            exit(1);
        }
    }   

    #init data structurs
    $self->_missingChars({});
    $self->_matchedPairs({});
    $self->_data({});
    $self->_selectedLoci({});
    $self->_currentFingerprints([]);
    $self->_dataHeader([]);

    #defaults
    $self->_setMissingChars(['-','?','.']);
    $self->_exhaustedFingerprints(0);
    $self->_numberOfFingerprints(0);
}

sub logger{
    my $self=shift;
    $self->{'_logger'}=shift // return $self->{'_logger'};
}

sub _missingChars{
    my $self=shift;
    $self->{'__missingChars'}=shift // return $self->{'__missingChars'};
}

sub _matchedPairs{
    my $self=shift;
    $self->{'__matchedPairs'}=shift // return $self->{'__matchedPairs'};
}


sub inputFile{
    my $self=shift;
    $self->{'_inputFile'}=shift // return $self->{'_inputFile'};
}

sub lociNumber{
    my $self=shift;
    $self->{'_lociNumber'}=shift // return $self->{'_lociNumber'};
}

sub _fileHandle{
    my $self=shift;
    $self->{'__fileHandle'}=shift // return $self->{'__fileHandle'};
}

sub _data{
    my $self=shift;
    $self->{'__data'}=shift // return $self->{'__data'};
}

sub _currentFingerprints{
    my $self=shift;
    $self->{'__currentFingerprints'}=shift // return $self->{'__currentFingerprints'};
}

sub _selectedLoci{
    my $self=shift;
    $self->{'__selectedLoci'}=shift // return $self->{'__selectedLoci'};
}

sub _exhaustedFingerprints{
    my $self=shift;
    $self->{'__exhaustedFingerprints'}=shift // return $self->{'__exhaustedFingerprints'};
}

sub _dataHeader{
    my $self=shift;
    $self->{'__dataHeader'}=shift // return $self->{'__dataHeader'};
}

sub _numberOfFingerprints{
    my $self=shift;
    $self->{'__numberOfFingerprints'}=shift // return $self->{'__numberOfFingerprints'};
}

sub run{
    my $self =shift;

    $self->_storeDataInMemory($self->inputFile);

    while(my $locus = $self->_getNextLocus()){
        $self->_selectedLoci->{$locus}=1;
        $self->_addLocusToCurrentFingerprint($self->_data->{$locus});
    }
    continue{
        if(($self->lociNumber eq 'best') && ($self->_exhaustedFingerprints == 1) ){
            last;
        }
        elsif(scalar(keys %{$self->_selectedLoci}) == $self->lociNumber){
            last;
        }
    }

    #do the output
    $self->_printResults([sort keys %{$self->_selectedLoci}]);
}


=head2 _addLocusToCurrentFingerprint

Need to add the values from the current locus to the growing fingerprints

=cut

sub _addLocusToCurrentFingerprint{
    my $self=shift;
    my $locus=shift;

    $locus = $self->_substituteMissingChars($locus);
    for my $i(0..(scalar(@{$locus})-1)){
        $self->_currentFingerprints->[$i] = $self->_currentFingerprints->[$i] . $locus->[$i];
    }
}


sub _printResults{
    my $self = shift;
    my $loci = shift;

    $self->logger->info("Printing results " . ref($loci));

    foreach my $header(@{$self->_dataHeader}){
        print "\t" . $header;
    }
    print "\n";

    foreach my $locus(@{$loci}){
        print $locus . "\t" . "@{$self->_data->{$locus}}" . "\n";
    }
}

=head2 _setMissingChars

Sets the characters to be ignored when calculating POD and fingerprints.

=cut

sub _setMissingChars{
    my $self=shift;
    my $chars=shift;

    foreach my $char(@{$chars}){
        $self->logger->debug("Setting missing char " . $char);
         $self->_missingChars->{$char}=1;
    }
}


=head2 _getNextLocus

Iterates through the remaining loci until the lociNumber has been hit.

=cut


sub _getNextLocus{
    my $self=shift;

    my $loci;
    unless($self->_exhaustedFingerprints == 1){
        $loci = $self->_getAllBestFingerprintLoci();
    }  

    my $podLoci;
    if(!defined $loci->[0]){
        $self->_exhaustedFingerprints(1);
        $self->logger->info("Fingerprints exhausted");
        $loci = sort keys %{$self->_data};
    }
    else{
         $self->logger->debug("Sending " . scalar(@{$loci} . " to _getAllBestPodLoci"));
         $podLoci = $self->_getAllBestPodLoci($loci);
    }   

    unless(defined $podLoci->[0]){
        $self->logger->info("Resetting matched pairs");
        $self->_matchedPairs({});
        $podLoci = $self->_getAllBestPodLoci($loci);
    }
    return $podLoci->[0] // undef;
}


=head2 _getAllBestPodLoci

Given a list of input loci, return all that have the maximum POD,
taking into account those that have already been matched in $self->_matchedPairs

=cut

sub _getAllBestPodLoci{
    my $self=shift;
    my $loci = shift;

    my @bestLoci;
    my $topPod=0;
    foreach my $locus(@{$loci}){

        if(defined $self->_selectedLoci->{$locus}){
            next;
        }    
            
        my $pod = $self->_calculatePod($self->_data->{$locus});
       # $self->logger->debug("Pod locus: $locus pod: $pod");
        if($pod > $topPod){
            @bestLoci=($locus);
            $topPod=$pod;
        }
        elsif($pod == $topPod && $pod !=0){
            push @bestLoci, $locus;
        }
    }
    $self->logger->debug("Top pod: $topPod: loci: @bestLoci");
    return \@bestLoci;
}

=head2 _calculatePod

We only need to compare each set once, not bi-directionally.
Thus {i}->{j} is sufficient, we don't need {j}->{i}.
eg. [0][1][2]
The outer loop starts at [0] and goes to [1].
The inner loop starts at [1] and goes to [2].
This gives us:
[0]->[1], [0]->[2]
[1]->[2]
Which is all we need, and saves the needless duplication if both
loops were to run through everything.

=cut

sub _calculatePod{
    my $self=shift;
    my $locus=shift;

    my $pod=0;
    my $numberOfLoci = scalar(@{$locus});
    $self->logger->debug("Number of pod loci: $numberOfLoci");
    for my $i(0..($numberOfLoci-2)){
        for my $j(($i+1)..($numberOfLoci-1)){
            if(defined $self->_matchedPairs->{$i}->{$j}){
               # $self->logger->debug("Matched pair defined: $i:$j");
                next;
            }
            else{
               # $self->logger->debug("Incrementing pod");
                unless(defined $self->_missingChars->{$locus->[$i]} || defined $self->_missingChars->{$locus->[$j]}){
                    if($locus->[$i] ne $locus->[$j]){
                        $pod++;
                        $self->_matchedPairs->{$i}->{$j}=1;
                    }
                }                
            }
        }
    }
    return $pod;  
}

=head2 _getAllBestFingerprintLoci

Looks at all available loci, and returns all that create the most fingerprints from the data.
We want to choose loci that do not contain "missingCharacters" if possible.
To accommodate this, we push those loci to the back of the array, and unshift any that do not
contain them. In this way, if they exist, the [0] element will contain a "missingCharacter" free
locus.

=cut

sub _getAllBestFingerprintLoci{
    my $self=shift;  

    $self->logger->info("Getting best fingerprints");

    my @chosenLoci;
    my $initialValue = $self->_numberOfFingerprints;
    my $topValue=$initialValue;

    foreach my $locus(sort keys %{$self->_data}){
        #$self->logger->debug("Next locus: $locus");
        if(defined $self->_selectedLoci->{$locus}){
            next;
        }

        my $fingerprintValue = $self->_calculateFingerprint($self->_data->{$locus});

        if($fingerprintValue > $topValue){
            @chosenLoci=($locus);
            $topValue=$fingerprintValue;
            #$self->logger->debug("New best locus: $locus value: $fingerprintValue");
        }
        elsif($fingerprintValue == $topValue && $topValue != $initialValue){
            push @chosenLoci, $locus;
        }    
    }
    $self->logger->info("Best fingerprint value of $topValue");
    $self->logger->info("initialValue: $initialValue topValue: $topValue");
    $self->_numberOfFingerprints($topValue);
    return \@chosenLoci;
}


=head2 _substituteMissingChars

Determines whether or not a locus contains a missing char.
Replaces with a '.'

=cut

sub _substituteMissingChars{
    my $self=shift;
    my $locus=shift;

    foreach my $missingChar(sort keys %{$self->_missingChars}){
        for my $i(0..scalar(@{$locus})-1){           
            if($locus->[$i] eq $missingChar){
                $locus->[$i]='.';
            }
        }
    }
    return $locus;
}

=head2

Returns the number of unique fingerprints, if the current locus is chosen.

=cut

sub _calculateFingerprint{
    my $self = shift;
    my $locus = shift;

    my %tempData;
    foreach my $i(0..(scalar(@{$locus})-1)){
        my $locusValue = $locus->[$i];
        my $value;

        if(defined $self->_missingChars->{$locusValue}){
            $value='.';
            #$self->logger->debug("Missing $value");
        }
        else{
            $value = $locusValue;
        }
        my $key = $self->_currentFingerprints->[$i] . $value;
        $tempData{$key} = 1;
    }
    return $self->_accountForMissingChararcters([keys %tempData]);
}

=head2 _accountForMissingCharacters

We need to ensure that the 'missing characters' are not
seen as unique values, which will give an incorrect, inflated number of
fingerprints.

=cut

sub _accountForMissingChararcters{
    my $self=shift;
    my $prints = shift;

    my $numberOfFingerprints = scalar @{$prints};
    #$self->logger->debug("Accounting for missing chars: number of fingerprints: $numberOfFingerprints");

    for my $i(0..($numberOfFingerprints-2)){
        my $iPrint = $prints->[$i];
        #$self->logger->debug("iprint: $iPrint");
        for my $j(($i+1)..($numberOfFingerprints-1)){
            my $jPrint = $prints->[$j];
            #$self->logger->debug("jPrint: $jPrint");
            if($jPrint =~ m/$iPrint/ || $iPrint =~ m/$jPrint/){
                #$self->logger->debug("Decreasing number of fingerprints from $numberOfFingerprints");
                $numberOfFingerprints--;
                last;
            }
            else{
                #$self->logger->debug("$jPrint not equal to $iPrint");
            }          
        }
    }
    return $numberOfFingerprints;
}


=head2 _storeDataInMemory

Stores the input file as a hash, with the locus name as the key,
and the data as an array reference.

=cut

sub _storeDataInMemory{
    my $self=shift;
    my $file =shift;

    my $inFH=IO::File->new('<' . $file) or die "$!";
    $self->logger->info("Collecting data from $file");

    while(my $line = $inFH->getline){ 
        $line =~ s/\R//g;
        my @la = split('\t',$line);

         if($inFH->input_line_number == 1){
            foreach my $head(@la){
                if($head eq ''){
                    next;
                }
                push @{$self->_dataHeader},$head;
            }
            next;
        }

        my $locus = shift(@la);
        $self->_data->{$locus}=\@la;
    }
    $inFH->close();
}




1;