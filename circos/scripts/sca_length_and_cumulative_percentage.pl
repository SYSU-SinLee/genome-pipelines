use strict;
use warnings;
use Bio::SeqIO;

die "Usage: perl sca_length_and_cumulative_percentage.pl PREFIX.fa >PREFIX.sca_length\n" if @ARGV==0;
my $fa=Bio::SeqIO->new(-file=>$ARGV[0],-format=>'fasta');

my %length;
my $total_length;
while (my $obj=$fa->next_seq()) {
	my $id=$obj->id;
	my $length=$obj->length;
	$length{$id}=$length;
	$total_length+=$length;
}

print "# rank\tchr\tlength\tcumulative_percentage\n";
my @ids=sort {$length{$b}<=>$length{$a}} keys %length;
my $cumulative_length;
my $sort_index;
foreach my $id (@ids) {
	$sort_index++;
	$cumulative_length+=$length{$id};
	my $cumulative_percentage=100*$cumulative_length/$total_length;
	print "$sort_index\t$id\t$length{$id}\t";
	printf "%.2f%%\n",$cumulative_percentage;
}