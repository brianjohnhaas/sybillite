use strict;

sub parse_scaffold_list {
    my ($org, $gff3_file, $tmpDir, $PURGE_CACHE) = @_;

    my @scaffolds;
    
    my $cached_file = $tmpDir . "/$org." . basename($gff3_file) . ".scaffolds.cached";
    if (-s $cached_file && ! $PURGE_CACHE) {
        open (my $fh, $cached_file);
        while (<$fh>) {
            chomp;
            push (@scaffolds, $_);
        }
        close $fh;
    
    } else {
        my %scaffs;
        open (my $fh, $gff3_file) or die "Error, cannot open $gff3_file";
        while (<$fh>) {
            chomp;
            unless (/\w/) { next; }
            if (/^\#/) {
                next;
            }

            my @x = split(/\t/);
            my $contig = $x[0];
            $scaffs{"$org;$contig"} = 1;
        }
        close $fh;
        @scaffolds = sort keys %scaffs;
        
        # write to cache file
        open (my $ofh, ">$cached_file") or die "Error, cannot write to cache: $cached_file";
        foreach my $scaffold (@scaffolds) {
            print $ofh "$scaffold\n";
        }
        close $ofh;
    }

    return(@scaffolds);
}



## CGI.pm was messing this up, so I'm doing my own
sub scaffold_encode {
    my $str = $_[0];
    $str =~ s/\;/___/;
    return $str;
}

sub scaffold_decode {
    my $str = $_[0];
    $str =~ s/___/\;/;
    return $str;
}





1==1;
