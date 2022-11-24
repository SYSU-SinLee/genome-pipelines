use 5.011;

while (<>) {
        if (/^Orthogroup/) {
                chomp;
                my $line=$_;
                my @eles=split/\s+/,$line;
                print "Desc\tFamily ID\t";
                shift @eles;
                pop @eles;
                my $lastone=pop @eles;
                map {print "$_\t"} @eles;
                say $lastone;
        } else {
                chomp;
                my $line=$_;
                my @eles=split/\s+/,$line;
                print "(null)\t";
                pop @eles;
                my $lastone=pop @eles;
                map {print "$_\t"} @eles;
                say $lastone;
        }
}
