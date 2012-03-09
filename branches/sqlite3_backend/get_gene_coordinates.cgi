#!/usr/bin/env perl

=head1 DESCRIPTION



=cut

use strict;
use CGI;
use DBI;
use JSON;

my $cgi = new CGI();
print $cgi->header( -type => 'application/json' );

my $json = JSON->new->allow_nonref;
my $project = $cgi->param('project') || die "ERROR: project is a required parameter";
my $gene_id = $cgi->param('gene') || die "ERROR: gene_id is a required parameter";

my %json_data = ( found => 0, fmin => 0, fmax => 0, org_abbrev => '', molecule => '' );

my $db_file = "./data/$project/annotation.db";
my $dbh = DBI->connect( "dbi:SQLite:$db_file" ) || die "Cannot connect: $DBI::errstr";

my $qry = qq{
    SELECT *
      FROM gene
     WHERE gene_id = ?
};

my $dsh = $dbh->prepare($qry);
   $dsh->execute( $gene_id );

while (my $row = $dsh->fetchrow_hashref ) {
    
    $json_data{found} = 1;
    
    $json_data{org_abbrev} = $$row{org_abbrev};
    $json_data{molecule}   = $$row{molecule};
    
    if ( $$row{strand} eq '+' ) {
        $json_data{fmin} = $$row{start};
        $json_data{fmax} = $$row{end};
    } else {
        $json_data{fmin} = $$row{end};
        $json_data{fmax} = $$row{start};
    }
    
    last;
}

$dsh->finish();
$dbh->disconnect();

#use Data::Dumper;
#print STDERR "found gene with data: ", Dumper(\%json_data);

## JSON output
print $json->encode(\%json_data);

exit(0);
