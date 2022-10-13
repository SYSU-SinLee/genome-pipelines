open ANNO,"<","$ARGV[0]";
open GENE,">","$ARGV[1]";

print GENE "GO\tGENE\tCLASS\n";


while (<ANNO>) {
	next if /^Gene/;
	s/\s+$//;
	my ($Gene,$BP,$MF,$CC)=split/\t/;
	@BPOGs=split/\) /,$BP;
	foreach my $BPOG (@BPOGs) {
		$BPOG=~m/(?<GO>GO:\d+)\((?<DESC>.*)/;
		print GENE "$+{GO}\t$Gene\tbiological_process\n";
	}
	@MFOGs=split/\) /,$MF;
	foreach my $MFOG (@MFOGs) {
		$MFOG=~m/(?<GO>GO:\d+)\((?<DESC>.*)/;
		print GENE "$+{GO}\t$Gene\tmolecular_function\n";
	}
	@CCOGs=split/\) /,$CC;
	foreach my $CCOG (@CCOGs) {
		$CCOG=~m/(?<GO>GO:\d+)\((?<DESC>.*)/;
		print GENE "$+{GO}\t$Gene\tcellular_component\n";
	}
}

