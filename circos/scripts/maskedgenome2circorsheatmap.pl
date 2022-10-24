use strict;
use warnings;
use Bio::SeqIO;
die "Usage: perl $0 PREFIX.length PREFIX.fa.masked [target_sca_num] [block_size] [species prefix] PREFIX.repeatcontent >PREFIX.repeat.range\n" if @ARGV!=6;

open LENGTH,"<$ARGV[0]";
my $fastq=$ARGV[1];
my $target_sca_num=$ARGV[2];
my $block_size=$ARGV[3];
my $prefix=$ARGV[4];
open OUT,">$ARGV[5]";
my @target_scas;
my %scaname;
my %length;
while (<LENGTH>) {
	next if /^#/;
	my ($index,$chr,$length)=(split/\s+/,$_)[0,1,2];
	push @target_scas,$chr if ($index<=$target_sca_num);
	$scaname{$chr}="$prefix$index";
	$length{$chr}=$length;
}
close LENGTH;
my @values;
my $maskedgenome=Bio::SeqIO->new(-file=>$fastq,-format=>'fasta');
while (my $obj=$maskedgenome->next_seq()) {
	my $id=$obj->id;
	my $seq=$obj->seq;
	if ($id~~@target_scas) {
		my @seq=split//,$seq;
		my $length=$length{$id};
		my $block_num=int ($length/$block_size);
		foreach my $num (1..$block_num) {
			my %num;

			for (1..$block_size) {
				my $base=shift @seq;
				$num{$base}++;
			}
			my $percent=$num{'N'}/$block_size;
			$percent=int ($percent*=100);
			push @values,$percent;
			my $sta=($num-1)*$block_size;
			my $end=$num*$block_size-1;
			print OUT "$scaname{$id}\t$sta\t$end\t$percent\n";
			undef %num;
		}
		my %num;
		foreach my $base (@seq) {
			$num{$base}++;
		}
		my $sta=$block_num*$block_size;
		my $end=$length{$id}-1;
		my $percent=$num{'N'}/$block_size;
		$percent=int ($percent*=100);
		push @values,$percent;
		print OUT "$scaname{$id}\t$sta\t$end\t$percent\n";
		undef %num;
	}
}
@values=sort {$a<=>$b} @values;
my $max=pop @values;
my $min=shift @values;
print "$max\t$min\n";