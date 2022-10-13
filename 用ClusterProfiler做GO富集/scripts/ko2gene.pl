open KO,"<","$ARGV[0]";
open GENE,">","$ARGV[1]";

print GENE "KO\tGENE\n";

while (<KO>) {
	s/\s+$//;
	my @eles=split/\t/;
	my $KO=shift @eles;
	for (1..4) {
		shift @eles;
	}
	foreach my $ele (@eles) {
		print GENE "$KO\t$ele\n";
	}
}