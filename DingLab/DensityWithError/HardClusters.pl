use strict;
use warnings;

use LWP::Simple;
use FileHandle;
use Data::Dumper qw(Dumper);
use Scalar::Util qw(looks_like_number);
use List::Util qw[min max shuffle];

######################## Command Line Arguments ################################

# quit unless we have the correct number of command-line args
my $num_args = $#ARGV + 1;
if ($num_args != 3) {
    print "\nUsage: HardClusters.pl RD.*.clusters_in_./Results Epsilon MinPts\n\n";
    exit;
}

my $inputFile1 = "./Results/$ARGV[0]";
my $Epsilon = $ARGV[1];
my $MinPts = $ARGV[2];

##########################  Reading from RD.out  #############################

my $this = {};
my @InitialCuts;

open(IN, "<$inputFile1") || die "Can't open $inputFile1: $!";
while (my $line = <IN>) {
	if ( not $line =~ /Cluster/ ) {
		chomp $line;
		my @tabs2 = split(/\t/,$line);
		push @InitialCuts, [$tabs2[0],$tabs2[1],$tabs2[2],$tabs2[7],$tabs2[8],$tabs2[9]];
		# Cluster	Gene/Drug	Mutation/Gene	Epsilon_prime	Avg_density	 Covering_clusters
	}
}

###############################################################################

#print Dumper \@InitialCuts;

for (my $i = 0; $i < scalar @InitialCuts; $i++) {
	$InitialCuts[$i][0] =~ /(\d+)\.(\d+)\.(\d+)/g;
	if ($2 != 0) {
		$this->{"InitialCuts"}->{$1}->{$2} = $InitialCuts[$i][3];
	}
	$this->{"Variants"}->{$InitialCuts[$i][1].":".$InitialCuts[$i][2]}->{"run0"}->{$1}->{$2}->{$3} = $InitialCuts[$i][4];
	$this->{"Memberships"}->{$1}->{$2}->{$3}->{$InitialCuts[$i][1].":".$InitialCuts[$i][2]} = 1;
}

#print Dumper $this;

################################################################################

#print "\nOrderedNodes array=\n";
MainOPTICS($this, "met.pairwise");
#print Dumper $this->{CurrentRDarray};

# Identify super clusters:

for (my $run = 1; $run < 2; $run++) {
	my $scs = 0; # super cluster start
	for (my $i = 1; $i < scalar @{$this->{CurrentRDarray}}; $i++) {
		#print "i=$i\n";
		if ( ${$this->{CurrentRDarray}}[$i][1] == 10 ) {
			#print "RD($i)=10\n";
			$scs = $i;	
		}
		else {
			$this->{"CurrentSuperClusters"}->{$scs} = $i;
		}	
	}
	print "\nCurrent SC=\n";
	print Dumper $this->{CurrentSuperClusters};

	foreach my $CurrentSCstart (keys $this->{CurrentSuperClusters}) {
		#if ($this->{CurrentSuperClusters}->{$CurrentSCstart} - $CurrentSCstart >= $MinPts) {
			for (my $i = $CurrentSCstart; $i <= $this->{CurrentSuperClusters}->{$CurrentSCstart}; $i++) {
				# print "current variant=${$this->{CurrentRDarray}}[$i][0]\n";
				# print "In Old=\n";
				# print Dumper keys %{$this->{Variants}->{${$this->{CurrentRDarray}}[$i][0]}->{run0}};
				my @SCArray = keys %{$this->{Variants}->{${$this->{CurrentRDarray}}[$i][0]}->{run0}};
				my $TotHits;
				if (exists $this->{SuperClusterMatching}->{$CurrentSCstart}->{$SCArray[0]}) {
					$TotHits = $this->{SuperClusterMatching}->{$CurrentSCstart}->{$SCArray[0]};
				}
				else {
					$TotHits = 0;
				}
				$TotHits++;
				$this->{"SuperClusterMatching"}->{$CurrentSCstart}->{$SCArray[0]} = $TotHits;
			}			
		#}
		my @SCmatchArray = sort { $this->{SuperClusterMatching}->{$CurrentSCstart}->{$a} <=> $this->{SuperClusterMatching}->{$CurrentSCstart}->{$b} } keys %{$this->{SuperClusterMatching}->{$CurrentSCstart}};
		my $SCmatch = pop @SCmatchArray;
		print "SC map= $CurrentSCstart\t$SCmatch\n";

	}

}
print "SC matching\n";
print Dumper $this->{"SuperClusterMatching"};



####################################################################
##########################  Functions  #############################
####################################################################

# sub MapSuperClusters {
# 	my $this = @_;

# 	foreach my $CurrentSCstart (keys $this->{CurrentSuperClusters})
# }

sub MainOPTICS {
	my ($this, $Pairwisefile)=@_;

	my %SetOfNodes;

	my $file = "./Test/$Pairwisefile";
	open(IN, "<$file") || die "Can't open $file: $!";

	while (my $line = <IN>) {
		chomp $line;
		my @tabs = split(/\t/,$line);
		my @char19 = split("",$tabs[19]);
		my $dis = $char19[0].$char19[1].$char19[2].$char19[3].$char19[4];
		my $key1 = CombineWords($tabs[0],$tabs[4]);
		my $value1 = CombineWords($tabs[9],$tabs[13]);

		$SetOfNodes{$key1}{distances}{$value1} = $dis;
		$SetOfNodes{$value1}{distances}{$key1} = $dis;
	}
	###### For variants in the same residue and chain 

	foreach my $key ( keys %SetOfNodes ) {
		#print "key= $key\n";
		$key =~ /(\w+)\:\D\.(\D+\d+)\D/g;
		my $keyGene = $1;
		my $keyRes = $2;
		my @hits = grep(/$keyGene\:\D\.$keyRes\D/g, keys %SetOfNodes);
		#print Dumper \@hits;
		foreach my $hit (@hits) {
			if ( $hit ne $key ) {
				$SetOfNodes{$key}{distances}{$hit} = "0";
				$SetOfNodes{$hit}{distances}{$key} = "0";
			}
		}
	}

	foreach my $i (keys %SetOfNodes) {
		$SetOfNodes{$i}{processInfo} = "False";
	}

	#print Dumper \%SetOfNodes;
	print "Number of Objects = ";
	print scalar keys %SetOfNodes;
	print "\n";

	my @SetOfCores;
	my @SetOfEdges;
	foreach my $key ( keys %SetOfNodes ) {
		if ( scalar keys $SetOfNodes{$key}{distances} >= $MinPts ) {
			push @SetOfCores, $key;
		}
		else {
			push @SetOfEdges, $key;
		}
	}
	@SetOfCores = shuffle @SetOfCores;
	@SetOfEdges = shuffle @SetOfEdges;
	my @SetOfCoresThenEdges = ( @SetOfCores, @SetOfEdges );

	###########################################################

	my @OrderedNodes;

	################# Main OPTICS function ####################

	foreach my $p ( @SetOfCoresThenEdges ) {
		#print "first p=$p\n";
		if ($SetOfNodes{$p}{processInfo} =~ "False") {
			########## Expand Cluster Order ###########
			my %neighbors; # is a hash with keys neigbor indices whose values are mutual separations
			my %OrderSeeds; # is a hash to add seeds
			%neighbors = %{GetNeighbors($p,$Epsilon,\%SetOfNodes)};
			$SetOfNodes{$p}{processInfo} = "True"; # set as processed
			my $RD = undef;
			my $CD;
			$CD = GetCoreDistance(\%neighbors,$MinPts);
			# print "p=$p and ";
			# print "CD=$CD\n";
			push @OrderedNodes, [$p,$RD,$CD]; # write to the file 
			if (defined $CD) {
				OrderSeedsUpdate(\%neighbors,$p,$CD, \%OrderSeeds, \%SetOfNodes);
				# print "For p=$p, OrderSeeds= \n";
				# print Dumper \%OrderSeeds;
				my $PrevObj = $p; # used to get the current obj. (To check whether variants are at the same location)
				while (scalar keys %OrderSeeds != 0) {
					my @SeedKeys = sort { $OrderSeeds{$a} <=> $OrderSeeds{$b} } keys %OrderSeeds;
					my @SeedValues = @OrderSeeds{@SeedKeys};
					#my $CurrentObject =  $SeedKeys[0]; # CurrentObject is the object having the least RD in OrderSeeds
					my $CurrentObject = GetCurrentObject(\@SeedValues, \@SeedKeys, $PrevObj);
					$PrevObj = $CurrentObject;
					#print "\n\n current object= $CurrentObject\t neighbors=";
					%neighbors = %{GetNeighbors($CurrentObject,$Epsilon,\%SetOfNodes)};
					#print Dumper \%neighbors;
					#print Dumper $SetOfNodes{$CurrentObject}{distances};
					$SetOfNodes{$CurrentObject}{processInfo} = "True"; # set as processed
					$RD = $SeedValues[0];
					$CD = GetCoreDistance(\%neighbors,$MinPts);
					push @OrderedNodes, [$CurrentObject,$RD,$CD]; # write to the file 
					delete $OrderSeeds{$CurrentObject};
					if (defined $CD) {
						#print "\tCurrent object is a core.(CD=$CD)\n Updated Order seeds list\n\t";
						OrderSeedsUpdate(\%neighbors,$CurrentObject,$CD, \%OrderSeeds, \%SetOfNodes);
						#print Dumper \%OrderSeeds;
					}
				}
			}
			# print "p=$p,(undefined CD) OrderedNodes= \n";
			# print Dumper \@OrderedNodes;
		}
	}

	### Replacing undefined RD by 10
	for (my $i = 0; $i < scalar @OrderedNodes; $i++) {
		if (not defined $OrderedNodes[$i][1]) {
			$OrderedNodes[$i][1] = 10;
		}
	}
	$this->{"CurrentRDarray"} = \@OrderedNodes;
	return $this;
}

sub GetNeighbors {
	my ($Obj, $Epsilon, $Set_ref)=@_;
	my %neighborHash;
	foreach my $i (keys %{$Set_ref->{$Obj}->{distances}}) {
			$neighborHash{$i} = "$Set_ref->{$Obj}->{distances}->{$i}";
	}
	return \%neighborHash;
}

sub GetCoreDistance {
	my ($neighbors_ref, $MinPts)=@_;
	my @keys = sort { $neighbors_ref->{$a} <=> $neighbors_ref->{$b} } keys %{$neighbors_ref}; # sort keys according to distances
	my @vals = @{$neighbors_ref}{@keys};
	my $CoreDist;
	if (scalar keys %{$neighbors_ref} >= $MinPts){
			$CoreDist = $vals[$MinPts-1]; # MinPt^th-distance
		}
	else {
		$CoreDist = undef;
	}
	return $CoreDist;
}

sub OrderSeedsUpdate {
	my ($neighbors_ref, $CenterObject, $CD, $OrderSeeds_ref, $Set_ref) = @_;
	my $c_dist = $CD; 
	my %neighborsHash = % { $neighbors_ref };
	my %OrderSeedsHash = % { $OrderSeeds_ref};
	foreach my $q (keys %{$neighbors_ref}) {
		if (${$Set_ref}{$q}{processInfo} =~ "False") {
			my $new_r_dist = max ($c_dist,${$neighbors_ref}{$q});
			if (exists ${$OrderSeeds_ref}{$q}) {
				if ($new_r_dist < ${$OrderSeeds_ref}{$q}) {
					${$OrderSeeds_ref}{$q}="$new_r_dist";
				}
			}
			else {
					${$OrderSeeds_ref}{$q}="$new_r_dist";
				}
		}
	}
}

sub CombineWords {
	my ($word1,$word2)=@_;
	return $word1.":".$word2;
}

sub GetCurrentObject { # To check whether variants are at the same location
	my ($ValueSet, $KeySet, $PrevObj)=@_;
	my @SmallestKeys;
	for (my $i = 0; $i < scalar @$ValueSet; $i++) {
		if ($ValueSet->[$i] == $ValueSet->[0]) {
			push @SmallestKeys, $KeySet->[$i];
		}
	}
	if (scalar @SmallestKeys > 1) { # more than one variant has the smallest RD
		$PrevObj =~ /(\w+)\:\D\.(\D+\d+)\D/g;
		my $keyGene = $1;
		my $keyRes = $2;
		my @hits = grep(/$keyGene\:\D\.$keyRes\D/g, @SmallestKeys);

		if (scalar @hits > 0) {
			unshift @SmallestKeys, $hits[0];
		}
	}
	return shift @SmallestKeys;
}
