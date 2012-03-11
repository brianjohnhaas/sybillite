#!/usr/bin/perl -w

=head1 DESCRIPTION

This script creates the SQLite3 file for a given project from its GFF.

Initialization of the search database assumes that the 9th column of each GFF3
file has ID and Name attributes.

=head1 INPUT

You pass a project name (such as 'Schizos') as the first parameter on the command
line and the SybilLite.conf file will use this to do everything else.  Therefore,
you must have set up that project in the Sybillite.conf first before running this
indexing script.

To allow incomplete paths to be defined in the conf file, this script must be run
from within the conf directory.

=head1 OUTPUT

Within the project directory a single SQLite3 database file called 'annotation.db'
will be created.

=cut

use strict;
use DBI;
use FindBin;

use lib ("$FindBin::Bin/../PerlLib");
use IniReader;
use SybilLite;
use Data::Dumper;

my $project = shift || die "ERROR: pass a project as the first parameter";
   $project =~ s/[^A-Za-z0-9\-]/_/g;   ## a wee bit of sanitization

## load configurations
my $SybilLiteConfFile = "$FindBin::Bin/../conf/SybilLite.conf";
unless (-e $SybilLiteConfFile) {
    die "Error, cannot find SybilLite main configuration file: $SybilLiteConfFile";
}
my $SybilLiteConf = new IniReader($SybilLiteConfFile);
my $project_conf_file = $SybilLiteConf->get_value($project, "Conf_file");
unless ($project_conf_file) {
    die Dumper($SybilLiteConf) . "Error, cannot find project configuration file for project $project"; 
}

unless (-e $project_conf_file) {
    die "Error, cannot locate project conf file for project: $project_conf_file";
}

my $conf = new IniReader("$FindBin::Bin/../conf/$project_conf_file");


## key = org abbreviation, value = GFF3 path
my $orgs = get_organisms_from_conf( $conf );

my $data_dir = "$FindBin::Bin/../data/$project/";
unless (-d $data_dir) {
    mkdir $data_dir or die "Error, cannot mkdir $data_dir";
    chmod(0775, $data_dir);
}

## prep regions of interest folder
my $roi_dir = "$data_dir/regions_of_interest";
unless (-d $roi_dir) {
    mkdir ($roi_dir) or die "Error, cannot mkdir $roi_dir";
    chmod (0777, $roi_dir); # webserver needs to be able to write to it.
}


my $db_file = "$data_dir/annotation.db";

## create the db index if it doesn't exist already
my $next_gene_id = 1;

if ( -f $db_file ) {
    
    print STDERR "ERROR: the db file ($db_file) already exists.  Cowardly refusing to " .
                 "overwrite it.  Remove or rename it and try again.\n";
    exit(1);
    
} else {
    
    print STDERR "INFO: creating database file: ($db_file)\n";
    
    my $dbh = DBI->connect( "dbi:SQLite:$db_file" ) || die "Cannot connect: $DBI::errstr";
    
    ## without this statement you could go make a nice dinner and enjoy it with your
    #   family in the time it would take to populate a medium project.  Notes why:
    #       http://www.sqlite.org/faq.html#q19
    $dbh->do("PRAGMA synchronous=OFF");
    
    initialize_database( $dbh, $db_file, $orgs );
    $dbh->disconnect();
}



exit(0);



sub get_organisms_from_conf {
    my $config = shift;

    my @organisms = split (/,/, $conf->get_value("Meta", "Organisms"));
    
    my $data = {};
    
    print STDERR "Got organisms: @organisms\n";

    foreach my $org (@organisms) {
        $org =~ s/\s+//g;
        
        $$data{$org} = $conf->get_value("Genes", $org);
    }
    
    return $data;
}


sub initialize_database {
    my ($dbh, $file, $orgs) = @_;
    
    my $create_gene_ddl = qq{
        CREATE TABLE gene (
            id          INTEGER,
            gene_id     TEXT,
            alias       TEXT,
            org_abbrev  TEXT,
            molecule    TEXT,
            start       TEXT,
            end         TEXT,
            strand      INTEGER,
            name        TEXT
        )
    };
    $dbh->do( $create_gene_ddl );
    
    my $insert_gene_dml = qq{
        INSERT INTO gene
        VALUES ( ?, ?, ?, ?, ?, ?, ?, ?, ? )
    };
    my $insert_gene_dsh = $dbh->prepare($insert_gene_dml);
    
    ## load each of the GFF3 files
    for my $abbreviation ( keys %$orgs ) {
        print STDERR "INFO: calling load_gff3_file( $insert_gene_dsh, $abbreviation, '$FindBin::Bin/../$$orgs{$abbreviation}' );\n";
        load_gff3_file( $insert_gene_dsh, $abbreviation, "$FindBin::Bin/../$$orgs{$abbreviation}" );
    }
    
    $insert_gene_dsh->finish();
    
    ## add the indexes
    $dbh->do("CREATE UNIQUE INDEX idx_gene_id ON gene (gene_id)");
    $dbh->do("CREATE INDEX idx_molecule ON gene (molecule)");
    $dbh->do("CREATE INDEX idx_name ON gene (name)");
}

sub load_gff3_file {
    my ($dsh, $org_abbrev, $file) = @_;
    
    open(my $ifh, $file) || die "can't read GFF3 file ($file): $!";
    
    while ( my $line = <$ifh> ) {
        chomp $line;
        my @cols = split("\t", $line);
        
        next unless scalar @cols >= 9;
        next unless $cols[2] eq 'gene';
        
        my @insert_data = (
            $next_gene_id++, ## id
            '', ## gene_id
            '', ## alias
            $org_abbrev,
            $cols[0], ## molecule
            $cols[3], ## start
            $cols[4], ## end
            $cols[6], ## strand
            '', ## name
        );
        
        my @atts = split(';', $cols[8]);
        for my $att ( @atts ) {
            my ($k, $v) = split('=', $att);
            
            if ( $k eq 'ID' ) {
                $insert_data[1] = $v;
                
            } elsif ( $k eq 'Name' ) {
                $insert_data[8] = $v;
            
            } elsif ( $k eq 'Alias' ) {
                $insert_data[2] = $v;
            }
        }
        
        $dsh->execute( @insert_data );
    }
}










