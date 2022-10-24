use strict;
use warnings;

die "Usage: perl $0 PREFIX.length PREFIX.genedensity [target_sca_num] [SPECIES PREFIX] PREFIX.genedensity.4circos >prefix.gd.range\n" if @ARGV!=5;

open LENGTH,"<$ARGV[0]";
open GD,"<$ARGV[1]";
my $target_sca_num=$ARGV[2];
my $prefix=$ARGV[3];
open OUT,">$ARGV[4]";
my @target_scas;
my %scaname;
while (<LENGTH>) {
	next if /^#/;
	my ($index,$chr)=(split/\s+/,$_)[0,1];
	push @target_scas,$chr if ($index<=$target_sca_num);
	$scaname{$chr}="$prefix$index";
}
close LENGTH;
my @gds;
my @gd_on_target_sca;
my @lines;
while (<GD>) {
	chomp;
	my ($sca,$sta,$end,$gd)=(split/\s+/)[0..3];
	$end-=1;
	$gd=int $gd;
	push @gds,$gd;
	if ($sca~~@target_scas) {
		$sca=$scaname{$sca};
		my @name=($sca,$sta,$end,$gd);
		my $name=\@name;
		push @lines,$name;
	}
}

@lines=sort {$a->[0] cmp $b->[0] or $a->[1]<=>$b->[1]} @lines;
map {print OUT "@$_\n"} @lines;
@gds=sort {$a<=>$b} @gds;
my $max=pop @gds;
my $min=shift @gds;
print "$max\t$min\n";