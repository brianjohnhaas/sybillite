#!/usr/bin/env perl

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use HTML::Template;
use FindBin;

use lib ("$FindBin::Bin/PerlLib", "$FindBin::Bin/../PerlLib");
use IniReader;

my $q = new CGI;
print $q->header( -type => 'text/html' );

my $tmpl = HTML::Template->new( filename => 'templates/project_choice.tmpl',
                                die_on_bad_params => 0,
                              );

my $SybilLiteConf = new IniReader("conf/SybilLite.conf");
my $project_list = $SybilLiteConf->get_value("Projects", "Project_list");

## each element like:
#   { project => '', organisms => '', description => '' };
my $orgs = [];

foreach my $project (split (/,/, $project_list) ) {
	$project =~ s/^\s+|\s+$//;
	
	my $organisms = $SybilLiteConf->get_value("$project", "Organisms");
	my $description = $SybilLiteConf->get_value("$project", "Description");
	
	push @$orgs, {
        project => $project,
        organisms => $organisms,
        description => $description
    };
	
  }

$tmpl->param( ORGS => $orgs );
$tmpl->param( PROJ_TITLE => 'Synteny browser' );

print $tmpl->output;

exit(0);
