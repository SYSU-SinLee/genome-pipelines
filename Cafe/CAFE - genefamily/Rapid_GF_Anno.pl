use strict;
use autodie;

die "Usage: perl $0 rapid_fams.txt Orthogroups.txt ath.gene.discription\n" if @ARGV != 3 ;

open RGF,"<","$ARGV[0]";
open ORTHO,"<","$ARGV[1]";
open ANNO,"<","$ARGV[2]";

my %Ath_Genes;
while (<ORTHO>) {
		s/\s+$//;
		my @eles=split/\s+/,$_;
		my $OG=shift @eles;
		$OG=~s/://;
		my @Ath_Genes;
		foreach my $ele (@eles) {
			push @Ath_Genes,$ele if $ele=~/Ath/;
		}
		$Ath_Genes{$OG}=[@Ath_Genes];
		undef @Ath_Genes;
}
close ORTHO;


my %Description;
while (<ANNO>) {
	next if /^Locus/;
	next if /^\s+$/;
	s/\s+$//;
	s/\(source\:Araport11\)//;
	s/\s+protein_coding//;
	my @eles=split/\s+/;
	my $Locus=shift @eles;
	my $Gene=shift @eles;
	my $Description=join(' ',@eles);
	$Description{$Gene}=$Description;
}

while (<RGF>) {
	next if (/^#/ || /^Overall/);
	s/\s+$//;
	my @eles=split/\s+/;
	my $spe=shift @eles;
	$spe=~s/\<//;
	$spe=~s/\>//;
	my @OGs=split/,/,(shift @eles);
	(my $out=$spe)=~s/:/\.rapid\.AthAnno/;
	open OUT,">","$out";
	foreach my $OG (@OGs) {
		$OG=~m/(\+|\-)/;
		my $Status=$1;
		$OG=~s/\[.*\]//;
		my @Ath_Genes=@{$Ath_Genes{$OG}};
		foreach my $Ath_Gene (@Ath_Genes) {
			$Ath_Gene=~s/^Ath\|//;
			print OUT "$OG$Status: $Description{$Ath_Gene}\n";
		}
	}

}
