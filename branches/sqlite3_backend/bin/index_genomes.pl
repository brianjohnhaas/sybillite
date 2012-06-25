#!/usr/bin/env perl 

=head1 DESCRIPTION

This script creates the SQLite3 file for a given project from its GFF.

Initialization of the search database assumes that the 9th column of each GFF3
file has ID and Name attributes.

=head1 INPUT

You pass a project name (such as 'Schizos') as the first parameter on the command
line and the SybilLite.conf file will use this to do everything else.  Therefore,
you must have set up that project in the Sybillite.conf first before running this
indexing script.

Instead of a single project name, if you pass the string 'ALL' each of the projects
in the SybilLite.conf file will be indexed.  Errors processing any individual project
will be reported but won't be fatal, so they won't stop processing of the other
projects.

=head1 OUTPUT

Within the project directory a single SQLite3 database file called 'annotation.db'
will be created.  If it already exists it will be removed and replaced.

=cut

use strict;
use warnings;

use DBI;
use FindBin;
use File::Copy;

use lib ("$FindBin::Bin/../PerlLib");
use IniReader;
use SybilLite;
use Data::Dumper;

my $project = shift || die "ERROR: pass a project as the first parameter, or ALL to index all genomes";

## load configurations
my $SybilLiteConfFile = "$FindBin::Bin/../conf/SybilLite.conf";
unless (-e $SybilLiteConfFile) {
    die "Error, cannot find SybilLite main configuration file: $SybilLiteConfFile";
}

my $SybilLiteConf = new IniReader($SybilLiteConfFile);
my @projects = ();

if ( $project eq 'ALL' ) {
    my $project_list = $SybilLiteConf->get_value('Projects', 'Project_list');
    @projects = split(',', $project_list);
    
} else {
    push @projects, $project;
}

for my $p ( @projects ) {
    process_project( $SybilLiteConf, $p );
}

exit(0);


sub process_project {
    my ($sybil_conf, $proj) = @_;
    
    $proj =~ s/\s//g;
    $proj =~ s/[^A-Za-z0-9\-]/_/g;   ## a wee bit of sanitization

    print STDERR "INFO: processing project: ($proj)\n";

    my $project_conf_file = "$FindBin::Bin/../conf/" . $sybil_conf->get_value($proj, "Conf_file");

    unless ($project_conf_file) { 
        print STDERR "ERROR: Can't find project configuration file for project $proj\n";
        return 0;
    }

    unless (-e $project_conf_file) {
        print STDERR "ERROR: Can't locate project conf file for project: $project_conf_file\n";
        return 0;
    }

    my $conf = new IniReader($project_conf_file);

    ## key = org abbreviation, value = GFF3 path
    my $orgs = get_organisms_from_conf( $conf );

    my $data_dir = "$FindBin::Bin/../data/$proj/";
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
    my $db_file_tmp = "$db_file.tmp";  ## written here first, then moved over
    
    my $coord_files = get_coordfiles_from_conf( $conf );

    print STDERR "INFO:\tcreating database file: ($db_file)\n";
    
    if ( -e $db_file_tmp ) {
        unlink($db_file_tmp);
    }

    my $dbh = DBI->connect( "dbi:SQLite:$db_file_tmp" ) || die "Cannot connect: $DBI::errstr";

    ## without this statement you could go make a nice dinner and enjoy it with your
    #   family in the time it would take to populate a medium project.  Notes why:
    #       http://www.sqlite.org/faq.html#q19
    $dbh->do("PRAGMA synchronous=OFF");

    initialize_database( $dbh, $orgs, $coord_files );
    $dbh->disconnect();

    move( $db_file_tmp, $db_file );
    
    return 1;
}

sub get_coordfiles_from_conf {
    my $config = shift;
    
    my @synteny_labels = $config->get_section_attributes("Synteny");
    my @syn_files;
    
    foreach my $synteny_label (@synteny_labels) {
        my $syn_filename = $config->get_value("Synteny", $synteny_label);
        push (@syn_files, $syn_filename);
    }
    
    return \@syn_files;
}

sub get_organisms_from_conf {
    my $conf = shift;

    my @organisms = split (/,/, $conf->get_value("Meta", "Organisms"));
    
    my $data = {};
    
    print STDERR "INFO:\tGot organisms: @organisms\n";

    foreach my $org (@organisms) {
        $org =~ s/\s+//g;
        
        $$data{$org} = $conf->get_value("Genes", $org);
    }
    
    return $data;
}

## schema changes here should be reflected in documentation updates
#   in doc/SCHEMA
sub initialize_database {
    my ($dbh, $orgs, $coordfiles) = @_;
    
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
        print STDERR "INFO:\tloading GFF3 file: $$orgs{$abbreviation}\n";
        load_gff3_file( $insert_gene_dsh, $abbreviation, "$FindBin::Bin/../$$orgs{$abbreviation}" );
    }
    
    $insert_gene_dsh->finish();
    
    ## add the indexes
    $dbh->do("CREATE UNIQUE INDEX idx_gene_id ON gene (gene_id)");
    $dbh->do("CREATE INDEX idx_molecule ON gene (molecule)");
    $dbh->do("CREATE INDEX idx_name ON gene (name)");
    
    my $create_aligncoords_ddl = qq{
        CREATE TABLE aligncoords (
            id          INTEGER,
            qry_org     TEXT,
            qry_mol     TEXT,
            qry_acc     TEXT,
            qry_fmin    INTEGER,
            qry_fmax    INTEGER,
            qry_strand  INTEGER,
            ref_org     TEXT,
            ref_mol     TEXT,
            ref_acc     TEXT,
            ref_fmin    INTEGER,
            ref_fmax    INTEGER,
            ref_strand  INTEGER,
            e_value     TEXT
        )
    };
    $dbh->do( $create_aligncoords_ddl );
    
    my $insert_alncoord_dml = qq{
        INSERT INTO aligncoords
        VALUES ( ?,?,?,?,?,?,?,?,?,?,?,?,?,? );
    };
    my $insert_alncoord_dsh = $dbh->prepare($insert_alncoord_dml);
    
    for my $file ( @$coordfiles ) {
        print STDERR "loading alignment coordinates from file: $file\n";
        load_aligncoord_file( $insert_alncoord_dsh, "$FindBin::Bin/../$file" );
    }
    
    $insert_alncoord_dsh->finish();
    
    ## add the indexes
    $dbh->do("CREATE UNIQUE INDEX idx_aligncoords_id ON aligncoords (id)");
    $dbh->do("CREATE INDEX idx_qry_org ON aligncoords (qry_org)");
    $dbh->do("CREATE INDEX idx_qry_mol ON aligncoords (qry_mol)");
}

sub load_aligncoord_file {
    my ($dsh, $file) = @_;
    
    open(my $ifh, $file) || die "can't read aligncoord file ($file): $!";
    
    my $next_aligncoords_id = 1;
    
    while ( my $line = <$ifh> ) {
        chomp $line;
        my @cols = split("\t", $line);
        
        next unless scalar @cols >= 14;
        
        my ($qry_fmin, $qry_fmax, $qry_strand) = tigr_to_interbase($cols[3], $cols[4]);
        my ($ref_fmin, $ref_fmax, $ref_strand) = tigr_to_interbase($cols[10], $cols[11]);
        
        my @insert_data = (
            $next_aligncoords_id++,     # id
            $cols[0],                   # qry_org
            $cols[1],                   # qry_mol
            $cols[2],                   # qry_acc
            $qry_fmin,                  
            $qry_fmax,
            $qry_strand,
            $cols[7],                   # ref_org
            $cols[8],                   # ref_mol
            $cols[9],                   # ref_acc
            $ref_fmin,
            $ref_fmax,
            $ref_strand,
            $cols[13],                  # e_value
        );
        
        $dsh->execute( @insert_data );
        $next_aligncoords_id++;
    }
}

sub tigr_to_interbase {
    my ($end5, $end3) = @_;
    my ($fmin, $fmax, $strand);
    
    if ( $end5 <= $end3 ) {
        $fmin = $end5 - 1;
        $fmax = $end3;
        $strand = 1;
        
    } else {
        $fmin = $end3 - 1;
        $fmax = $end5;
        $strand = -1;
    }
    
    return ($fmin, $fmax, $strand);
}

sub load_gff3_file {
    my ($dsh, $org_abbrev, $file) = @_;
    
    ## create the db index if it doesn't exist already
    my $next_gene_id = 1;
    
    if ( ! -e $file || ! -r $file ) {
        print STDERR "ERROR:\tCan't read and skipping GFF3 file: ($file): $!";
        return;
    }
    
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










