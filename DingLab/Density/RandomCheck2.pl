use strict;
use warnings;

use LWP::Simple;
use FileHandle;
use Data::Dumper qw(Dumper);
use Scalar::Util qw(looks_like_number);
use List::Util qw[min max];

######################## Command Line Arguments ################################

# quit unless we have the correct number of command-line args
# my $num_args = $#ARGV + 1;
# if ($num_args != 2) {
#     print "\nUsage: RandomCheck.pl cluster_file1 cluster_file2\n\n";
#     exit;
# }

for (my $file = 1; $file < 5; $file++) {
	my $secondfile = $file+1;
	my $inputFile1 = "./Results/$file/RD.8.4.met.pairwise.SuperClustersID.clusters";
	my $inputFile2 = "./Results/$secondfile/RD.8.4.met.pairwise.SuperClustersID.clusters";

	########################################

	my @InitialSet1;

	open(IN, "<$inputFile1") || die "Can't open $inputFile1: $!";
	while (my $line = <IN>) {
		chomp $line;
		my @tabs1 = split(/\t/,$line);
		push @InitialSet1, [$tabs1[2],$tabs1[0]];
	}

	close(IN);

	my @InitialSet2;

	open(IN, "<$inputFile2") || die "Can't open $inputFile2: $!";
	while (my $line = <IN>) {
		chomp $line;
		my @tabs2 = split(/\t/,$line);
		push @InitialSet2, [$tabs2[2],$tabs2[0]];
	}
	close(IN);

	my $Variants = {};

	for (my $i = 0; $i < scalar @InitialSet1; $i++) {
		if (exists $Variants->{$InitialSet1[$i][0]}->{"1"}) {
			$Variants->{$InitialSet1[$i][0]}->{"1"} = $Variants->{$InitialSet1[$i][0]}->{"1"}.":".$InitialSet1[$i][1];
		}
		else {
			$Variants->{$InitialSet1[$i][0]}->{"1"} = $InitialSet1[$i][1];
		}
	}
	for (my $i = 0; $i < scalar @InitialSet2; $i++) {
		if (exists $Variants->{$InitialSet2[$i][0]}->{"2"}) {
			$Variants->{$InitialSet2[$i][0]}->{"2"} = $Variants->{$InitialSet2[$i][0]}->{"2"}.":".$InitialSet2[$i][1];
		}
		else {
			$Variants->{$InitialSet2[$i][0]}->{"2"} = $InitialSet2[$i][1];	
		}
	}

	#print Dumper $Variants;

	my $OutFile1 = "./Results/Compare$file-$secondfile";
	open (OUT, ">$OutFile1");

	foreach my $variant (keys $Variants) {
		my @first = split(":", $Variants->{$variant}->{"1"});
		my @second = split(":", $Variants->{$variant}->{"2"});

		# print "$variant\n";
		# print Dumper \@first;
		# print "\n";
		if (scalar @first == scalar @second) {
			for (my $i = 1; $i < scalar @first; $i++) {
				
				$first[$i] =~ /\d+\.(\d+)\.(\d+)/g;
				my $firstsub1 = $1;
				my $firstsub2 = $2;
				$second[$i] =~ /\d+\.(\d+)\.(\d+)/g;
				if ($firstsub1 != $1 || $firstsub2 != $2) {
					print OUT "\n$variant is not in the same $firstsub1.$firstsub2 cluster. It's in $1.$2\n";
					for (my $j = 0; $j < scalar @InitialSet1; $j++) {
						if ($InitialSet1[$j][1] eq $first[$i]) {
							$Variants->{$first[$i]}->{"1"}->{$InitialSet1[$j][0]} = 0;
						}
						if ($InitialSet2[$j][1] eq $second[$i]) {
							$Variants->{$first[$i]}->{"2"}->{$InitialSet2[$j][0]} = 0; 
						}
					}
					my $numberOfobj = scalar keys $Variants->{$first[$i]}->{"1"};
					my $n = 0;
					foreach my $key1 (keys $Variants->{$first[$i]}->{"1"}) {
						foreach my $key2 (keys $Variants->{$first[$i]}->{"2"}) {
							if ($key1 eq $key2) {
								$n++;
								delete $Variants->{$first[$i]}->{"1"}->{$key1};
								delete $Variants->{$first[$i]}->{"2"}->{$key2};
								last;
							}
						}
					}

					if ($numberOfobj != $n) {
						print OUT "Looks fine!\n";
					}
					else {
						print OUT "\t Looks suspicious... $first[$i]--$second[$i]\n";
					}
				}
			}
		}
		else {
			print OUT "\nFor $variant, number of clusters doesn't match\n";
		}
	}	
	print OUT "Variants Hash=\n";
	print OUT Dumper $Variants;
}

