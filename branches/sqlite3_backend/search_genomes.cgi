#!/usr/bin/env perl

=head1 DESCRIPTION

This is used to search the GFF3 of one or more genomes in a given project to look
for matches to the user search input.  It relies on a very de-normalized SQLite3 
schema that exists for each project in a file called 'annotation.db'. 

Results are returned in a JSON object.

=cut

use strict;
use CGI;
use DBI;
use FindBin;
use JSON;

use lib ("$FindBin::Bin/PerlLib");
use IniReader;
use SybilLite;

my $cgi = new CGI();
print $cgi->header( -type => 'application/json' );

my $json = JSON->new->allow_nonref;

my $project = $cgi->param('project') || die "ERROR: project is a required parameter";
   $project =~ s/[^A-Za-z0-9\-]/_/g;   ## a wee bit of sanitization
   
my $term    = $cgi->param('term') || die "ERROR: term is a required parameter";

#print STDERR "SEARCHING FOR: ($term)\n";

my %json_data = ( results => [] );
my @json_data = ();

my $db_file = "./data/$project/annotation.db";

## create the db index if it doesn't exist already
if (! -f $db_file ) {
    die "ERROR: couldn't find expected index file: $db_file";
}

## connect to an existing db
my $dbh = DBI->connect( "dbi:SQLite:$db_file" ) || die "Cannot connect: $DBI::errstr";

my $qry = qq{
    SELECT * FROM gene WHERE name LIKE ? OR gene_id LIKE ?
};
my $dsh = $dbh->prepare($qry);
   $dsh->bind_param( 1, "%$term%" );
   $dsh->bind_param( 2, "%$term%" );
   $dsh->execute(); 

my $limit = 7;
my $found = 0;

while ( my $row = $dsh->fetchrow_hashref ) {
    push @json_data, { label => "$$row{org_abbrev}:$$row{gene_id} - $$row{name}", value => $$row{gene_id} };
    last if ++$found == $limit;
}

$dsh->finish();
$dbh->disconnect();

## JSON output
print $json->encode(\@json_data);

exit(0);













