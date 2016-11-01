use strict;
use warnings;

use LWP::Simple;
use FileHandle;
use Data::Dumper qw(Dumper);
use Scalar::Util qw(looks_like_number);
use List::Util qw[min max];

######################## Command Line Arguments ################################

# quit unless we have the correct number of command-line args
my $num_args = $#ARGV + 1;
if ($num_args != 2) {
    print "\nUsage: RandomCheck.pl cluster_file1 cluster_file2\n\n";
    exit;
}

my $inputFile1 = "$ARGV[0]";
my $inputFile2 = "$ARGV[1]";

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

print Dumper $Variants;

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
			if ($firstsub1 == $1 && $firstsub2 == $2) {
				print "for $variant, firstsub1.firstsub2 = $1.$2\n";
			}
			else {
				print "$variant is not in the same $firstsub1.$firstsub2 cluster. It's in $1.$2\n";
			}
			
		}
	}
	else {
		print "For $variant, number of clusters doesn't match\n";
	}
}