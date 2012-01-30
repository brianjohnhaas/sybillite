#!/usr/bin/env perl

## test range: SJ1: 7000000036978760 : 140000-160000

=head1 DESCRIPTION

The script controls the main functional interface of SybilLite.  It presents the header,
footer and any appropriate sub-controls, then decides which content is appropriate for the 
result panel.  

The general options are:

    - display any available 'regions of interest' defined if no aligments parameters exist
      yet (this is the initial state)
    - display a multi-alignment


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

use lib ("$FindBin::Bin/PerlLib", "$FindBin::Bin/../PerlLib");
use IniReader;
use SybilLite;

my $cgi = new CGI();
    
## this gets populated later once we actually know which template we're exporting
my $tmpl;

my $SybilLiteConf = new IniReader("conf/SybilLite.conf");

my $tmpDir = $ENV{WEBSERVER_TMP} || cwd() . "/tmp";

my $project = $cgi->param('project') || die "ERROR: project is a required parameter";
my $project_conf_file = $SybilLiteConf->get_value($project, "Conf_file");
my $conf = new IniReader("conf/$project_conf_file");
my $proj_title = $conf->get_value("Meta", "Title"); 

## this gets toggled if the URL API is used to define an alignment region to view
#   keyed on whether scaffold is manually defined.
my $url_api_used = 0;

print $cgi->header( -type => 'text/html' );
$tmpl = HTML::Template->new( filename => 'templates/sybil_lite.tmpl',
                             die_on_bad_params => 0 );

my $default_scaff_value;

if ( $cgi->param('scaffold') ) {
   $default_scaff_value = $cgi->param('scaffold');
   $default_scaff_value =~ s/___/\;/;
   $url_api_used = 1;

} else {
    $default_scaff_value = $conf->get_value("Meta", "DefaultMolecule") || ""
}

$default_scaff_value = scaffold_decode($default_scaff_value);

## build list of gff3 files
my @organisms = split (/,/, $conf->get_value("Meta", "Organisms"));
my @scaffolds;

foreach my $org (@organisms) {
    $org =~ s/\s+//g;
    my $gff3_file = $conf->get_value("Genes", $org);
    my @scaffold_list = &parse_scaffold_list($org, $gff3_file, $tmpDir, 0);
    push (@scaffolds, @scaffold_list);
}

my $organism_order = get_organism_order();
my $molecule_selections = &create_molecule_selections(\@scaffolds, $default_scaff_value);

## create the organism list and, for now, make them all pre-selected unless passed in the URL
my $organism_selections = get_organism_selections();

## handle the ROI view, which can optionally be specified in the project conf file
my $graphical_roi_view = 1; ## set the default

## first check passed parameters, then check the conf file
if ( $cgi->param('ROI_view') && $cgi->param('ROI_view') eq 'table') {
    $graphical_roi_view = 0;
    
} elsif ( $conf->get_value('Display', 'ROI_view') && 
          $conf->get_value('Display', 'ROI_view') eq 'table' ) {

    $graphical_roi_view = 0;
}

my $range = $cgi->param('range') || $conf->get_value("Meta", "DefaultRange") || "";


$tmpl->param( GRAPHICAL_ROI_VIEW => $graphical_roi_view );    ## set to 0 for 'table' view
$tmpl->param( DEFAULT_FLANK => $cgi->param('flank') || 0 );
$tmpl->param( DEFAULT_SCAFF_RANGE => $range );
$tmpl->param( KEEP_IMAGE_FLAG => $cgi->param('keep_image') ? 1 : 0 );
$tmpl->param( MOLECULE_SELECTIONS => $molecule_selections );
$tmpl->param( ORGANISM_ORDER => $organism_order );
$tmpl->param( ORGANISM_SELECTIONS => $organism_selections );
$tmpl->param( PROJ_TITLE => $proj_title );
$tmpl->param( PROJECT => $project );
$tmpl->param( PIXELS_PER_KB => 10 );  ## this should be made configurable later
$tmpl->param( URL_API_USED => $url_api_used );

print $tmpl->output;


exit(0);


## Gets organism order from several possible data sources.  Priority is:
#   1. URL API
#   2. Config file Display:Organism_order
#   3. Config file Meta:Oragnisms
sub get_organism_order {
    my $org_order;
    
    if ( $cgi->param('orgsOrder') ) {
        $org_order = $cgi->param('orgsOrder');
        
    } elsif ( $conf->get_value("Display", "Organism_order") ) {
        $org_order = $conf->get_value("Display", "Organism_order");
        
    } elsif ( $conf->get_value("Meta", "Organisms") ) {
        $org_order = $conf->get_value("Meta", "Organisms");
    
    } else {
        die "Failed to get an organism order\n";
    }
    
    $org_order =~ s/\s//g;
    
    return $org_order;
}


####
#  currently populates the "scaffold" select box
sub create_molecule_selections {
    my ($scaffolds_aref, $default_scaff_value) = @_;
    my $mols = [];
    
    foreach my $scaffold (@$scaffolds_aref) {
        if ($scaffold eq "$default_scaff_value") {
            push @$mols, { selected => 1, label => $scaffold, value => $scaffold };
        } else {
            push @$mols, { selected => 0, label => $scaffold, value => $scaffold };
        }
    }
    
    return $mols;
}

sub get_organism_selections {
    my $org_selected_list;
    
    ## was this passed in the URL API?
    if ( $cgi->param('orgsSelected') ) {
        $org_selected_list = $cgi->param('orgsSelected');
    
    ## else just use the full list
    } else {
        $org_selected_list = $organism_order;
    }
    
    $org_selected_list =~ s/\s//g;
    
    my $org_selections = &create_org_selections($organism_order, $org_selected_list);
    
    return $org_selections;
}

####
sub create_org_selections {
    my ($chosen_org_order, $orgs_selected) = @_;
    my @organisms = split (/,/, $chosen_org_order);
    my %preselected = ();
    
    for ( split(/,/, $orgs_selected) ) {
        $preselected{$_} = 1;
    }
    
    my $org_selections = [];
    
    foreach my $org (@organisms) {
        if ( exists $preselected{$org} ) {
            push @$org_selections, { label => $org, enabled => 1 };
        } else {
            push @$org_selections, { label => $org, enabled => 0 };
        }
    }
    
    return $org_selections;
}
