#!/usr/bin/env perl

=head1 DESCRIPTION

This CGI should eventually be removed.  Previously, it held all the functional code for the SybilLite
display but has been replaced by more modular scripts.

This script only remains so that any other software linking to SybilLite will still work as intended.
It only checks the parameters passed and redirects to the appropriate (current) page.

=cut

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);

my $cgi = new CGI();

my $this_url = $cgi->self_url();
my $new_url = '';

## if a project defined, pass on all parameters to sybil_lite
if ( $cgi->param('project') ) {
    $new_url = $this_url;
    $new_url =~ s/SybilLite.cgi/sybil_lite.cgi/;

## if a project isn't displayed redirect to the project choice page
} else {
    $this_url =~ /^(.+)\/SybilLite.cgi/;
    $new_url = "$1/project_choice.cgi";
}

print $cgi->redirect(-uri => $new_url );

exit(0);
