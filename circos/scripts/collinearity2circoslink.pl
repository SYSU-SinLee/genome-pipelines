# Usage: perl collinearity2circoslink.pl PREFIX.gff PREFIX.collinearity >PREFIX.inks.txt

use strict;
use autodie;
die "Usage: perl collinearity2circoslink.pl PREFIX.gff PREFIX.collinearity >PREFIX.inks.txt\n" if @ARGV==0;
open GFF,"<","$ARGV[0]";
open CLO,"<","$ARGV[1]";

my (%chr_of_gene,%sta_of_gene,%end_of_gene);
while (<GFF>) {
	chomp;
	my ($chr,$gene,$sta,$end)=split/\s+/;
	$chr_of_gene{$gene}=$chr;
	$sta_of_gene{$gene}=$sta;
	$end_of_gene{$gene}=$end;
}
close GFF;

while (<CLO>) {
	next if /^#/;
	chomp;
	$_=~s/^\s*\d+\-\s*\d+:\s+//;
	my ($gene1,$gene2)=(split/\s+/,$_)[0,1];
	print "$chr_of_gene{$gene1} $sta_of_gene{$gene1} $end_of_gene{$gene1} $chr_of_gene{$gene2} $sta_of_gene{$gene2} $end_of_gene{$gene2}\n"; 
}
close CLO;