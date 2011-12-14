#!/usr/bin/env perl

=head1 DESCRIPTION

This parses a projects "regions of interest" (ROI) JSON files and returns them
as JSON objects, presumably from an AJAX call.

The only argument is a project name, and all files are found by file name
and path convention after that.  If the example project 'Schizos' is passed,
the script assumes the following to be present (relative to this script's
path):

    ./data/Schizos/regions_of_interest/some_label.json
    ./data/Schizos/regions_of_interest/some_other_label.json

Where 'some_label' above will be the ID of a given ROI.

=cut

use strict;
use CGI;
use JSON;

my $cgi = new CGI();
print $cgi->header( -type => 'application/json' );

my $json = JSON->new->allow_nonref;

my $project = $cgi->param('project') || die "ERROR: project is a required parameter";

my %json_data = ( regions => [] );

my $roi_dir = "./data/$project/regions_of_interest/";

if ( -d $roi_dir ) {
    opendir(my $idh, $roi_dir) || die "failed to read directory $roi_dir: $!";
    
    while ( my $file = readdir($idh) ) {
        next unless $file =~ /\.json$/;
        
        local $/;
        open( my $json_fh, "$roi_dir/$file" ) || die "failed to read input file: $!";
        my $json_text = <$json_fh>;
        push @{$json_data{regions}}, decode_json( $json_text );
    }
}

## JSON output
print $json->encode(\%json_data);

exit(0);
