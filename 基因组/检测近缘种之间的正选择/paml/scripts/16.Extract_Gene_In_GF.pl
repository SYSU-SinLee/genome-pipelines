use autodie;
use strict;

open POSGF,"<","$ARGV[0]";
open GF,"<","$ARGV[1]";
open ATHGENELIST,">","Ath.genelist.txt";
my $Spe=$ARGV[3];
open SPEGENELIST,">","$Spe.genelist.txt";

my @Positive_Selected_OGs;
while (<POSGF>) {
		s/\s+$//;
		my $OG=(split/\s+/)[0];
		$OG=~s/pal2aln_/OG/;
		$OG=~s/\-gb.*$//;
		push @Positive_Selected_OGs,$OG;
}
close POSGF;

while (<GF>) {
	s/\s+$//;
	my $Line=$_;
	my @eles=split/\s+/;
	my $OG=shift @eles;
	$OG=~s/://;
	if (grep {$OG eq $_} @Positive_Selected_OGs) {
		foreach my $ele (@eles) {
			if ($ele=~m/Ath\|(.*)/) {
				my $Ath_Gene=$1;
				print ATHGENELIST "$Ath_Gene\n";
			} elsif ($ele=~m/$Spe\|(.*)/) {
				my $Spe_Gene=$1;
				print SPEGENELIST "$Spe_Gene\n";
			}
		}
	}
}
close GF;