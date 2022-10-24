use strict;
use warnings;
use autodie;
die "Usage: perl sca_length2circors_karyotype.pl PREFIX.LENGTH PREFIX SCA_NUM >PREFIX.karyotype\n" if @ARGV!=3;
open LEN,"<","$ARGV[0]";
my $prefix=$ARGV[1];
my $sca_num=$ARGV[2];

while (<LEN>) {
	next if /^#/;
	chomp;
	my ($rank,$sca,$length)=(split/\s+/)[0,1,2];
	die if $rank>$sca_num;
	$length-=1;
	my $color=$rank%18;
	$color=($color==0) ? ("18") : ("$color");
	print "chr - $prefix$rank $prefix$rank 0 $length chr$color\n";
}
close LEN;