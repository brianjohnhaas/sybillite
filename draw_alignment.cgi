#!/usr/bin/env perl

## test range: SJ1: 7000000036978760 : 140000-160000

=head1 DESCRIPTION

To be called via AJAX only, this script draws an alignment and creates one of the
two options, depending on parameters passed:

    - JSON file defining image map and associated image
    - JSON file containing path to available (temporary) SVG

=cut

use strict;
use warnings;
use CGI;
use Carp;
use CGI::Carp qw(fatalsToBrowser);
use Data::Dumper;
use Cwd;
use File::Basename;
use FindBin;
use HTML::Template;
use JSON;

use lib ("$FindBin::Bin/PerlLib", "$FindBin::Bin/../PerlLib");
use IniReader;
use SybilLite;
use SynView;

my $cgi = new CGI();

my %params = $cgi->Vars;
my $params = \%params;

$|++;
umask(0000);
my $tmpDir = $ENV{WEBSERVER_TMP} || cwd() . "/tmp";
my $pid = $$;

## When a user clicks to scroll left or right, this percentage of that side of the previous
##  range will be retained.  This is meant to give the user some reference point of familiarity
##  so scrolling is more fluid.
my $SCROLL_OVERHANG_PCT = 20;

## When a user clicks to zoom in or out, this is the fold increase or decrease in zoom level, as
#   implemented by the number of pixels per kilobyte in the image.
my $ZOOM_LEVEL_CHANGE = 2;

my $SybilLiteConf = new IniReader("conf/SybilLite.conf");
my $project = $params{project} || die "ERROR: project is a required parameter";
my $project_conf_file = $SybilLiteConf->get_value($project, "Conf_file");
my $conf = new IniReader("conf/$project_conf_file");

my $PURGE_CACHE = $params->{PURGE_CACHE} ? 1 : 0;
my $DEBUG = $params->{DEBUG} || 0;

## this will ultimately be transformed into a JSON structure and printed
my %json_data = ();

my $pixelsPerKb = $params->{pixelsPerKb} || 10;
my $showAllVsAll = $params->{showAllVsAll} || 0;
my $export_as = $params->{export_as} || 'none';

## build the list of synteny pairs
my @synteny_labels = $conf->get_section_attributes("Synteny");
my @syn_files;
foreach my $synteny_label (@synteny_labels) {
    my $syn_filename = $conf->get_value("Synteny", $synteny_label);
    push (@syn_files, $syn_filename);
}
my $syn_files_string = join (",", @syn_files);

## build list of gff3 files
my @gff3_files;
my @organisms = split (/,/, $conf->get_value("Meta", "Organisms"));

foreach my $org (@organisms) {
    $org =~ s/\s+//g;
    my $gff3_file = $conf->get_value("Genes", $org);
    push (@gff3_files, "$org" . "::" . $gff3_file);
}

my $gff3_files_string = join(",", @gff3_files);

my $organism_order = $params{orgsOrder} || die "ERROR: orgsOrder is a required parameter";
   $organism_order =~ s/\s//g;
my $organisms_selected = $params->{orgsSelected} || die "ERROR: orgsSelected is a required parameter";
   $organisms_selected =~ s/\s//g;

my $scaffold = scaffold_decode($params->{scaffold});
my $range = $params->{range};

## if not defined Javascript will report the range as NaN-NaN.  Provide a 
#   default instead.
#if ( $range eq 'NaN-NaN' ) {
#    $range = '1-50000';
#}

my $original_range = $range;

my $image_map_file = $tmpDir . "/imap.$$.txt";
my $image_file;

if ( $export_as eq 'none' ) {
    $image_file = $tmpDir . "/sybilLite.$$.png";
} elsif ( $export_as eq 'svg' ) {
    $image_file = $tmpDir . "/sybilLite.$$.svg";
}

my $map_points = [];

if ( $scaffold ) {
    if ( $range ) {
        if ($range !~ /\b\d+-\d+\b/) {
            die "ERROR: Sorry, cannot parse range: [$range] ";
        }
        my ($lend, $rend) = split (/-/, $range);
        if ($lend > $rend) {
            die "ERROR: coordinate range is out of order: $lend > $rend.   Please reexamine. ";
        }
        
        ## if a scroll direction was defined we need to change the range
        my $scroll_direction = $params->{scroll_direction} || 'none';

        if ( $scroll_direction ne 'none' ) {
            my $range_span = abs( $rend - $lend );
            
            ## if we're at the beginning of the molecule, we have to correct the span (since we say 1-N rather than 0-N)
            $range_span++ if ( $lend == 1 );
            
            my $range_overhang = int($range_span * ( $SCROLL_OVERHANG_PCT/100) );
            my $range_shift = $range_span - $range_overhang;
            
            if ( $scroll_direction eq 'right' ) {
                $rend += $range_shift;
                $lend += $range_shift;

            } elsif ( $scroll_direction eq 'left' ) {
                $rend -= $range_shift;
                $lend -= $range_shift;

            } else {
                die("ERROR: invalid scroll direction passed.  Values should be 'none', 'left', or 'right'");
            }
        }
        
        if (my $flank = $params->{flank}) {
            $lend -= $flank;
            $rend += $flank;
            if ($lend < 1) { $lend = 1; }
        }
        
        ## if we've shifted off the left side of the molecule, correct
        if ( $lend < 1 ) {
            $rend = 1 + ( $rend - $lend );
            $lend = 1;
        }
        
        ## set the range back after any corrections
        $range = "${lend}-$rend";
    }
    
    ## did the user request a zoom?
    my $zoom_direction = $params->{zoom_direction} || 'none';
    
    if ( $zoom_direction ne 'none' ) {
    
        if ( $zoom_direction eq 'in' ) {
            $pixelsPerKb *= $ZOOM_LEVEL_CHANGE;
            
        } elsif ( $zoom_direction eq 'out' ) {
            $pixelsPerKb /= $ZOOM_LEVEL_CHANGE;
            
        } else {
            die("ERROR: invalid zoom direction passed.  Values should be 'none', 'in', or 'out'");
        }
    }

    my %synPlotOptions = (   "synGenePairs" => $syn_files_string,
                             "gff3Files" => $gff3_files_string,
                             "refScaffold" => $scaffold,
                             "refScaffoldRange" => $range,
                             "orgsOrder" => $organisms_selected,
                             "pixelsPerKb" => $pixelsPerKb,
                             "imageFile" => $image_file,
                             "imageMapCoordsFile" => $image_map_file,
                             "showAllVsAll" => $showAllVsAll,
                         );

    #use Data::Dumper;
    #print STDERR Dumper(\%synPlotOptions);

    if ( $export_as eq 'svg' ) {
        $synPlotOptions{format} = 'svg';
    
    }

    &SynView::createSyntenyPlot( %synPlotOptions );

    # build the image map
    open (my $fh, $image_map_file) or die "Error, cannot open image map file: $image_map_file";
    while (<$fh>) {
        chomp;
        my ($x1, $y1, $x2, $y2, $name) = split (/\t/);
        push @$map_points, { label => $name, x1 => $x1, y1 => $y1, x2 => $x2, y2 => $y2 };
    }
    close $fh;

    unlink ($image_map_file);
}

$json_data{image_file} = $image_file;
$json_data{map_points} = $map_points;
## the range is passed back so both client and server layers don't have to
#   do all the shift/zoom calculations
$json_data{range} = $range;
$json_data{pixels_per_kb} = $pixelsPerKb;

if ( $export_as eq 'none' ) {

    ## JSON output
    print $cgi->header( -type => 'application/json' );
    my $json = JSON->new->allow_nonref;
    print $json->encode(\%json_data);

} elsif ( $export_as eq 'svg' ) {
    
    print $cgi->header( -type => 'image/svg+xml' );
    open(my $svg_fh, "<$image_file") || die "failed to open SVG image file: $!";
    while (<$svg_fh>) {
        print $_;
    }
    
}























