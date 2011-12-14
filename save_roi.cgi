#!/usr/bin/env perl

=head1 DESCRIPTION

This saves a "region of interest" (ROI) as a JSON file and corresponding
image screenshot.

Files are saved by file name and path convention.  If the example project 
'Schizos' and label "awesome region" are passed, the script would create
the following files (relative to this script's path):

    ./data/Schizos/regions_of_interest/awesome_region.json
    ./data/Schizos/regions_of_interest/awesome_region.png
    ./data/Schizos/regions_of_interest/awesome_region_150h.png

There are at least two image files.  The first is just the label passed but
with non-alphanumeric characters and whitespace replaced with underlines.  The 
second is the same basename but with "_Nh.png" appended to the end, where 'N'
is the height of the image in pixels.  Depending on user settings, different 
versions of the image might be used in different parts of the interface.

The JSON looks like this:

    {
        "id" : "awesome_region",
        "project" : "Schizos",
        "label" : "Awesome region",
        "orgsOrder" : [ "SP2", "SO3", "SJ1" ],
        "orgsSelected" : [ "SP2", "SO3" ],
        "scaffold" : "SP2;7000000090838467",
        "range" : "600000-1700000",
        "flank" : 0,
        "pixelsPerKb" : 1,
        "description" : "Here a long region of synteny exists between the two displayed genomes with little rearrangement."
    }

=cut

use strict;
use CGI;
use JSON;
use File::Copy;
use GD;

my $cgi = new CGI();
print $cgi->header( -type => 'application/json' );

my $json = JSON->new->allow_nonref;
   $json->pretty;  ## he is a handsome fellow

my $project     = $cgi->param('project') || die "ERROR: project is a required parameter";
my $label       = $cgi->param('roi_save_label') || die "ERROR: label is a required parameter";
my $id          = $label;  ## this will be transformed
   $id =~ s/[^a-zA-Z0-9]/_/g;

my @orgsOrder = split(',', $cgi->param('orgsOrder') );
my @orgsSelected = split(',', $cgi->param('orgsSelected') );

my %json_data = ( 
    "id" => $id,
    "project" => $project,
    "label" => $label,
    "orgsOrder" => \@orgsOrder,
    "orgsSelected" => \@orgsSelected,
    "scaffold" => $cgi->param('scaffold'),
    "range" => $cgi->param('range'),
    "flank" => $cgi->param('flank') || 0,
    "pixelsPerKb" => $cgi->param('pixelsPerlKb') || 10, ## make configurable default
    "description" => $cgi->param('roi_save_desc') || ''
);



my $roi_dir = "./data/$project/regions_of_interest/";
my $json_file = "$roi_dir/$id.json";
my $full_png = "$roi_dir/$id.png";
my $small_png = "$roi_dir/${id}_150h.png";

open(my $json_fh, ">$json_file") || die "failed to create JSON output file ($json_file): $!";

## JSON output
print $json_fh $json->encode(\%json_data);

## place the full-res file
copy( $cgi->param('current_image'), $full_png );

## create the 150h image file
my $full_gd = GD::Image->newFromPng($full_png) || die "Can't load PNG $full_png: $!";

my ($full_w, $full_h) = $full_gd->getBounds;
my ($small_w, $small_h) = (  int( $full_w * (150/$full_h)), 150  );
my $small_gd = GD::Image->new($small_w, $small_h);
$small_gd->copyResized($full_gd, 0, 0, 0, 0, $small_w, $small_h, $full_gd->getBounds);

open(my $small_fh, ">$small_png") || die "can't create small PNG file: $!";
binmode $small_fh;

print $small_fh $small_gd->png;

my %status = (
    success => 1,
    id => $id
);

## return our status as JASON too, to enable richer error handling
print $json->encode( \%status );

exit(0);






















