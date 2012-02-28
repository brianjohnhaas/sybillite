#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use URI::Escape;
use FindBin;
use lib ("/util/lib/perl5/site_perl/5.8.8/", "$FindBin::Bin/PerlLib", "$FindBin::Bin/../PerlLib");
use DPchain;

## Modified from J. Crabtree's testMultiPanel.pl script, bhaas Wed Sep 24 14:28:13 EDT 2008

my $usage = <<_EOUSAGE_;

#############################################################################
#
#  --synGenePairs  : list of files containing dagchainer results (ie. "fileA,fileB,...")
#
#  --gff3Files      : list of gff3 files containing the protein-coding gene coordinates (ie. "OrgA::fileA.gff3,OrgB::fileB.gff3")
#
#  --refScaffold    : accession of the reference scaffold molecule.
#  --refScaffoldRange  : range to view on the reference scaffold (ie. "1-10000") 
#
#  --orgsOrder      : order of organisms, including the reference
#
#  --pixelsPerKb    :  number of pixels per kilobase
#    
#  --miscFeaturesGff3    : misc features to include in the image.  Provide list of gff3 files to consider.
#
#  --imageFile           : file to save image; otherwise, output to stdout
#  --imageMapCoords      : file to which the image mapped coordinates will be written
#
#  --hideMatches         : pairs of organisms for which matches should be hidden.  Format:   orgA~orgB,orgC~orgD ...  hides between A-B and separately between C-D
#
##########################################################################################################


_EOUSAGE_

	;

my ($synGenePairs, $gff3Files, $refScaffold, $refScaffoldRange, $orgsOrder, 
	$miscFeaturesGff3, $imageMapCoordsFile, $imageFilename, $hideMatchesString);
my $pixelsPerKb = 10;
my $help;

&GetOptions("synGenePairs=s" => \$synGenePairs,
			"gff3Files=s" => \$gff3Files,
			"refScaffold=s" => \$refScaffold,
			"pixelsPerKb=f" => \$pixelsPerKb,
			"refScaffoldRange=s" => \$refScaffoldRange,
			"orgsOrder=s" => \$orgsOrder,
			"miscFeaturesGff3=s" => \$miscFeaturesGff3,
			"imageMapCoords=s" => \$imageMapCoordsFile,
			"imageFile=s" => \$imageFilename,
			"help|h" => \$help,
			"hideMatches=s" => \$hideMatchesString,
			);

if ($help) { 
	die $usage;
}

unless ($synGenePairs && $gff3Files && $refScaffold) {
	die $usage;
}

my %HIDE_MATCHES;
if ($hideMatchesString) {
  my @org_pairs = split (/,/, $hideMatchesString);
  foreach my $org_pair (@org_pairs) {
	my ($orgA, $orgB) = split (/\~/, $org_pair);
	$orgA =~ s/\s+//g;
	$orgB =~ s/\s+//g;
	$HIDE_MATCHES{$orgA}->{$orgB} = 1;
	$HIDE_MATCHES{$orgB}->{$orgA} = 1;
  }
}


# prepare for graphics drawing.
use Bio::Graphics::Panel;
use Bio::SeqFeature::Generic;
use FindBin;
use lib ("$FindBin::Bin");
use MultiPanel;


main: {

	
	## parse the input files
	my %syn_gene_pairs = &parse_syn_gene_pairs($synGenePairs, $refScaffold);
	 
	my %scaffold_to_gene_structs = &parse_gene_coords($gff3Files);
	
	
	# transpose the scaffold to gene data
	my %gene_acc_to_gene_struct = &transpose_scaffold_genes_list_to_gene_acc_lookup(\%scaffold_to_gene_structs);
	
	# get scaffold lengths (cheap; just taking the right-most coordinate of any parsed feature on the scaffold for now.
	my %scaffold_lengths = &parse_scaffold_lengths(\%scaffold_to_gene_structs);	
	
	
	## assign reference scaffold range coordinates.
	my ($refScaffLend, $refScaffRend) = &assign_reference_scaffold_range_coordinates($refScaffold, $refScaffoldRange, \%scaffold_lengths);
		
	## pull out the syntenic scaffolds and genes found syntenic to the reference scaffold range of interest.
	my %syntenic_scaffolds_to_genes = &extract_scaffolds_syntenic_to_ref_scaffold($refScaffold, [$refScaffLend, $refScaffRend], \%scaffold_to_gene_structs, \%gene_acc_to_gene_struct, \%syn_gene_pairs);
	

	## Extract syn region ranges:
	my ($ref_org, $trash) = split (/;/, $refScaffold);
	my %syn_org_to_scaffolds = ( $ref_org => [$refScaffold] ); # init to ref info.
	my %syn_contig_to_seq_range  = ( $refScaffold => { lend => $refScaffLend, rend => $refScaffRend, ref_lend => $refScaffLend, ref_rend => $refScaffRend } ); # init reference info.
	&assign_syntenic_contig_seq_ranges(\%syntenic_scaffolds_to_genes, \%syn_org_to_scaffolds, \%syn_gene_pairs, \%syn_contig_to_seq_range, \%gene_acc_to_gene_struct, [$refScaffLend, $refScaffRend]);	
	
	
	## organize the order of organism tiers for display
	
	my @top_orgs;
	my @bottom_orgs;
	&organize_order_of_organism_tiers(\%syn_org_to_scaffolds, $ref_org, $refScaffold, $orgsOrder, \@top_orgs, \@bottom_orgs);

	## get other misc features if they are provided:
	my %scaffold_to_misc_features;
	if ($miscFeaturesGff3) {
		%scaffold_to_misc_features = &parse_misc_features_in_range($miscFeaturesGff3, \%syn_contig_to_seq_range);
	}
			
	
	#####################################################################
	## draw them out.

	my $yposn = 50;
	my $PANEL_VGAP = 30;
	my $pad_left = 50;
	my $pad_right = 50;
	my $pad_bottom = 50;
	
	my %gene_acc_to_feats;
	my @panels;

	my %scaff_orients;


	my $get_glyph_type_sref = sub { my ($feature,$option_name,$part_no,$total_parts,$glyph) = @_;
									return($feature->{__glyph_type});
						 
	};
	


	foreach my $org_tier (\@top_orgs, [$ref_org], \@bottom_orgs) { 
		
		my @orgs = @$org_tier;
		foreach my $org (@orgs) {
			unless ($org) { next; } ## should look into this.  Shouldn't have undefined values here.
			print STDERR "-processing $org\n";
			
			unless (exists $syn_org_to_scaffolds{$org}) { next; }
			
			my @scaffolds = @{$syn_org_to_scaffolds{$org}};
			
			my @scaff_tiers;
			if (scalar @scaffolds == 1 && $scaffolds[0] eq $ref_org) {
				## no tiering necessary
				@scaff_tiers = ([$ref_org]);
			}
			else {
				## tiers, plus sets the draw_lend and draw_rend atts of %syn_contig_to_seq_range components
				@scaff_tiers = &tier_scaffolds_compute_positions(\@scaffolds, \%syn_contig_to_seq_range, $refScaffold, { pad_left => $pad_left, pad_right => $pad_right});
			}
			
			foreach my $scaff_tier (@scaff_tiers) {
				
				my $prev_panel;
				foreach my $scaffold (@$scaff_tier) {
					
					my $range_coords_href = $syn_contig_to_seq_range{$scaffold};
					my ($rangeLend, $rangeRend) = ($range_coords_href->{lend}, $range_coords_href->{rend});
					my $length = $rangeRend - $rangeLend + 1;
					
					my $syn_orient = ($scaffold eq $refScaffold) ? '+' : &estimate_synteny_orientation($scaffold, [$rangeLend, $rangeRend], 
																									   $refScaffold, [$range_coords_href->{ref_lend}, $range_coords_href->{ref_rend}],
																									   \%scaffold_to_gene_structs,
																									   \%syn_gene_pairs,
																									   \%gene_acc_to_gene_struct);
					
					$scaff_orients{$scaffold} = $syn_orient; # store for later so we know if the gene matches are same vs. inverted
					
					my ($draw_lend, $draw_rend);

					if ($scaffold eq $refScaffold) {
						($draw_lend, $draw_rend) = (1, int ($length * $pixelsPerKb/1000 + 0.5));
					}
					else {
						($draw_lend, $draw_rend) = ($range_coords_href->{draw_lend}, $range_coords_href->{draw_rend});
					}
					
					my $draw_length = $draw_rend + $pad_left + $pad_right;
					

					## Drawing the scaffold ticker
					print STDERR "-drawing $scaffold, range $rangeLend-$rangeRend length: $length, draw length: $draw_length, draw_lend: $draw_lend to draw_rend: $draw_rend\n";
										
					my $reverse_flag = ($syn_orient eq '-') ? 1 : 0;
					my $panel = Bio::Graphics::Panel->new(-length => $length, -offset => $rangeLend, -width => $draw_length, 
														  -pad_top => $yposn, -pad_left => $draw_lend + $pad_left, 
														  -pad_right => $pad_right, -pad_bottom => $pad_bottom, -flip => $reverse_flag,
														  -key_style => "between");
					$prev_panel = $panel;
					
					my $scale = Bio::SeqFeature::Generic->new(-primary_tag=> "$scaffold-scale", -start=>$rangeLend, -end => $rangeRend);
					$panel->add_track($scale, -glyph => 'anchored_arrow', -tick=> 2, -key => $scaffold);
					
					#my $scaff_glyph = Bio::SeqFeature::Generic->new(-primary_tag=> $scaffold, -start=>$rangeLend, -end => $rangeRend);
					#$panel->add_track([$scaff_glyph], -glyph => 'generic');
					
					#####################################
					## Draw each of the gene features
					my @feats;
					foreach my $gene (@{$scaffold_to_gene_structs{$scaffold}}) {
						my ($lend, $rend, $name) = ($gene->{lend}, $gene->{rend}, $gene->{name});
						unless ($lend < $rangeRend && $rend > $rangeLend) { next; } # no overlap
						
						my $gene_orient = $gene->{orient};
						my $strand = ($gene_orient eq '+') ? 1:-1; #($syn_orient eq $gene_orient) ? 1 : -1;
						
						my $acc = $gene->{acc};
						my $feat = Bio::Graphics::Feature->new(-primary_tag => $acc, -start => $lend, -end => $rend, -strand => $strand, -name => $name, -type => 'gene');
						$gene_acc_to_feats{$acc} = $feat;
						push (@feats, $feat);
						#$panel->add_track($feat, -glyph => 'generic', -label => 1, -description => 1, -key => $acc);
					}
					$panel->add_track([@feats], -glyph => "transcript", -label => 0);
				  
					#####################################
					## Draw each of the misc features:
					
					## -each misc features file gets a separate track.
					
					if (my $tracks_href = $scaffold_to_misc_features{$scaffold}) {
						
						foreach my $track (keys %$tracks_href) {

							my @misc_feats;
							
							my $track_color = "black";
							my $track_bump_flag = 0;
							
							
							foreach my $feature (@{$tracks_href->{$track}}) {
								my ($lend, $rend, $name) = ($feature->{lend}, $feature->{rend}, $feature->{name});
								my $feature_orient = $feature->{orient};
								
								# print STDERR "-adding misc feature to $scaffold : $lend-$rend $name\n";
								
								my $glyph_type = "generic";
								my $strand = 0;
								
								if ($name =~ /inkler|CRN/) { 
									$track_color = "green";
									$track_bump_flag = 1;
									$glyph_type = "transcript";
									$strand = ($feature_orient eq '+') ? 1: -1; #($syn_orient eq $feature->{orient}) ? 1:-1;
								}
								elsif ($name =~ /rxlr/i) {
									$track_color = 'red';
									$track_bump_flag = 1;
									$glyph_type = "transcript";
									$strand = ($feature_orient eq '+') ? 1 : -1; #($syn_orient eq $feature->{orient}) ? 1:-1;
								}
																
								my $feat = Bio::Graphics::Feature->new(-start => $lend, -end => $rend, -strand => $strand, -name => $name, -type => "misc_feature");
								$feat->{__glyph_type} = $glyph_type;
								push (@misc_feats, $feat);
								
							}
							
							$panel->add_track([@misc_feats], -glyph => $get_glyph_type_sref, -label => 0, -bump => $track_bump_flag, -bgcolor => $track_color);
						}
						
					}
					push (@panels, $panel);
					#$yposn = $panel->height() + $PANEL_VGAP;
				} # end scaffold
				
				## stack panels between tiers and organisms
				$yposn = $prev_panel->height() + $PANEL_VGAP;
			
			} # end scaff tier
		} # end org
	} # end org_tier
	
	
	###########################################
	## describe the matches between genes:
	
	my @matches;
	#goto skip_matches;
	my %seen;
	my %orient_swap = ( '+' => '-', '-' => '+');
	foreach my $gene_acc (keys %syn_gene_pairs) {
		my $syn_gene_accs_href = $syn_gene_pairs{$gene_acc};
		foreach my $syn_acc (keys %$syn_gene_accs_href) {
			my $pair = join ("_", sort ($gene_acc, $syn_acc));
			if ($seen{$pair}) { next; }
			$seen{$pair} = 1;
			
			my $featA = $gene_acc_to_feats{$gene_acc};
			my $featB = $gene_acc_to_feats{$syn_acc};
			
			unless ($featA && $featB) { next; } # not within selected range.
			
			my $geneA = $gene_acc_to_gene_struct{$gene_acc};
			my $geneB = $gene_acc_to_gene_struct{$syn_acc};
			
			my $scaff_A = $geneA->{scaffold};
			my $scaff_B = $geneB->{scaffold};
			
			my ($orgA, $molA) = split (/;/, $scaff_A);
			my ($orgB, $molB) = split (/;/, $scaff_B);

			if ($HIDE_MATCHES{$orgA}->{$orgB}) { next; }
			

			my $geneA_orient = $geneA->{orient};
			my $geneB_orient = $geneB->{orient};
			my $scaffA_orient = $scaff_orients{$scaff_A};
			my $scaffB_orient = $scaff_orients{$scaff_B};
			
			if ($scaffA_orient eq '-') { 
				$geneA_orient = $orient_swap{$geneA_orient};
			}
			if ($scaffB_orient eq '-') {
				$geneB_orient = $orient_swap{$geneB_orient};
			}
			
			my ($reverse_flag, $color) = ($geneA_orient eq $geneB_orient) ? (0, 'pink') : (1, 'blue');
			
			my $match = { 'feat1' => $featA, 
						  'feat2' => $featB,
						  'bg' => $color,
						  'fg' => 'black',
						  'reverse' => $reverse_flag,
					  };
			push (@matches, $match);
		}
	}
	
    skip_matches:
	
	my $multi_panel = new MultiPanel([@panels], [@matches]);
	
	my $gd = $multi_panel->gd();
	if ($imageFilename) {
		open (my $ofh, ">$imageFilename") or die "Error, cannot write image to $imageFilename";
		print $ofh $gd->png();
		close $ofh;
	}
	else {
		print $gd->png();
	}

	
	if ($imageMapCoordsFile) {
		## Write the image map info:
		open (my $ofh, ">$imageMapCoordsFile") or die "Error, cannot write image map to $imageMapCoordsFile";
		foreach my $panel (@panels) {
			my @boxes = $panel->boxes();
			foreach my $box (@boxes) {
				my ($feature, $x1, $y1, $x2, $y2, $track) = @$box;
				my $name = $feature->{name};
				$name = uri_unescape($name);
				print $ofh "$x1\t$y1\t$x2\t$y2\t$name\n" if $name;
			}
		}
		close $ofh;
	}
	
	
	exit(0);
}


####
sub parse_syn_gene_pairs {
	my ($fileListString, $refScaffold) = @_;
	
	my %synPairs;
	## store A-> B and B-> A links.
	
	foreach my $file (split (/,/, $fileListString)) {
		$file =~ s/\s//g;
		
		open (my $fh, $file) or die "Error, cannot open file $file";
		while (<$fh>) {
			chomp;
			if (/^\#/) { next; }
			unless (/\w/) { next; }
			my @x = split (/\t/);
			#if (lc($x[0]) eq lc($refScaffold) || lc($x[4]) eq lc($refScaffold)) {
				my ($accA, $accB) = ($x[1], $x[5]);
				$synPairs{$accA}->{$accB} = 1;
				$synPairs{$accB}->{$accA} = 1;
			#}
		}
		close $fh;
	}

	return(%synPairs);
}


####
sub parse_gene_coords {
	my ($fileListString) = @_;

	my %scaffold_to_gene_structs;

	foreach my $org_file_info (split (/,/, $fileListString)) {
	  my ($org, $file) = split (/::/, $org_file_info);
	  $file =~ s/\s+//g;

		open (my $fh, $file) or die "Error, cannot find file $file";
		while (<$fh>) {
			chomp;
			if (/^\#/) { next; }
			unless (/\w/) { next; }
			
			my @x = split (/\t/);
			if ($x[2] eq 'gene') {
				my ($scaffold, $lend, $rend, $orient, $gene_info) = ($x[0], $x[3], $x[4], $x[6], $x[8]);
				
				$scaffold = "$org;$scaffold";
				
				$gene_info =~ /ID=([^\;\s]+)/ or die "Error, no gene ID for $_";
				my $gene_id = $1;
				if ($gene_info =~ /Alias=([^;\s]+)/) {
					$gene_id = $1;
				}
				my $name = $gene_id;
				if ($gene_info =~ /Name=([^;]+)/) {
					$name = $1;
					if ($gene_id ne $name) {
						$name = "$gene_id $name";
					}
				}
				
				push (@{$scaffold_to_gene_structs{$scaffold}}, { acc => $gene_id,
																 name => $name,
																 lend => $lend,
																 rend => $rend,
																 orient => $orient,
																 scaffold => $scaffold,
															 } );
			}
		}
		close $fh;
	}

	return (%scaffold_to_gene_structs);
}


####
sub parse_scaffold_lengths {
	my ($scaffold_to_gene_structs_href) = @_;

	## for now, take the right most gene coordinate as an approximate scaffold length

	my %scaffold_lengths;

	foreach my $scaffold (keys %$scaffold_to_gene_structs_href) {
		my @gene_structs = @{$scaffold_to_gene_structs_href->{$scaffold}};
		@gene_structs = sort {$a->{rend}<=>$b->{rend}} @gene_structs;
		my $max_rend = $gene_structs[$#gene_structs]->{rend};
		$scaffold_lengths{$scaffold} = $max_rend;
	}


	return(%scaffold_lengths);
}


####
sub extract_scaffolds_syntenic_to_ref_scaffold {
	my ($refScaffold, $refScaffCoords_aref, $scaffold_to_gene_structs_href, $gene_acc_to_gene_struct_href, $syn_gene_pairs_href) = @_;
	
	my ($ref_lend, $ref_rend) = @$refScaffCoords_aref;
	
	my @ref_genes = @{$scaffold_to_gene_structs_href->{$refScaffold}};
	
	my %other_scaffolds;
	foreach my $ref_gene (@ref_genes) {
		my $ref_gene_acc = $ref_gene->{acc};

		unless ($ref_gene->{lend} < $ref_rend && $ref_gene->{rend} > $ref_lend) { 
			# no overlap to reference scaffold range of interest.
			next; 
		}

		
		if (my $syn_href = $syn_gene_pairs_href->{$ref_gene_acc}) {
			my @syn_gene_accs = keys %$syn_href;
			foreach my $syn_gene_acc (@syn_gene_accs) {
				my $syn_gene = $gene_acc_to_gene_struct_href->{$syn_gene_acc};
				my $scaffold = $syn_gene->{scaffold};
				push (@{$other_scaffolds{$scaffold}}, $syn_gene);
			}
		}
	}
	return (%other_scaffolds);
}


####
sub compute_reference_synteny_range {
	my ($syn_genes_aref, $synGenePairs_href, $gene_acc_to_gene_struct_href, $refScaffold, $refScaffCoords_aref) = @_;
	
	my ($refScaffLend, $refScaffRend) = @$refScaffCoords_aref;
	
	my @ref_coords;
	
	foreach my $gene (@$syn_genes_aref) {
		my $acc = $gene->{acc};
	   
		my $syn_pairs_href = $synGenePairs_href->{$acc};
		if (ref $syn_pairs_href) {
			foreach my $other_acc (keys %$syn_pairs_href) {
				my $gene_struct = $gene_acc_to_gene_struct_href->{$other_acc};
				if ($gene_struct->{scaffold} eq $refScaffold
					&& $gene_struct->{lend} > $refScaffLend && $gene_struct->{rend} < $refScaffRend) {
					push (@ref_coords, $gene_struct->{lend}, $gene_struct->{rend});
				}
			}
		}
	}

	@ref_coords = sort {$a<=>$b} @ref_coords;
	my $lend = shift @ref_coords;
	my $rend = pop @ref_coords;

	return ($lend, $rend);
}


####
sub tier_scaffolds_compute_positions {
	my ($scaffolds_aref, $syn_contig_to_seq_range_href, $refScaffold, $padding_href) = @_;

	my $ref_scaffold_coords_info_href = $syn_contig_to_seq_range_href->{$refScaffold};
	my $ref_scaff_lend = $ref_scaffold_coords_info_href->{lend};

	my $min_draw_left = 1;

	my @eles;
	foreach my $scaffold (@$scaffolds_aref) {
		my $coords_info_href = $syn_contig_to_seq_range_href->{$scaffold};
		$coords_info_href->{scaffold} = $scaffold; # should have stored earlier

		my $lend = $coords_info_href->{lend};
		my $rend = $coords_info_href->{rend};
		my $length = $rend - $lend + 1;
		my $draw_length = int ($length * $pixelsPerKb / 1000 + 0.5);
		
		my $ref_lend = $coords_info_href->{ref_lend};
		my $ref_rend = $coords_info_href->{ref_rend};
		
		my $ref_midpt = int( ($ref_lend + $ref_rend)/2 - $ref_scaff_lend);  # remember, only examining a range within the reference scaffold!
		## what is the ref midpt in pixels?
		my $ref_midpt_pixels = int ($ref_midpt * $pixelsPerKb/1000 + 0.5);
		
		my $draw_left = int ($ref_midpt_pixels - $draw_length/2);
		#if ($draw_left < 1) { $draw_left = 1; }
		my $draw_right = $draw_left + $draw_length - 1;
		
		if ($draw_left < $min_draw_left) { $min_draw_left = $draw_left; }

		## store draw coords
		$coords_info_href->{draw_lend} = $draw_left;
		$coords_info_href->{draw_rend} = $draw_right;
		
		push (@eles, $coords_info_href);
	}
	
	## adjust for coordinates off scale
	if ($min_draw_left <= 0) {
	  my $adj = abs($min_draw_left) + 1;
	  foreach my $ele (@eles) {
		$ele->{draw_lend} += $adj;
		$ele->{draw_rend} += $adj;
	  }
	}
	

	## compute drawing coordinates based on midpoint match to reference scaffold
	
	## use DP to find the longest non-overlapping series of entries:
	if (scalar @$scaffolds_aref == 1) { 
		## nothing to do here, since there's only one scaffold.
		return ($scaffolds_aref);
	}
	
	my $get_base_score_sref = sub { my ($ele) = @_;
									return ($ele->{draw_rend} - $ele->{draw_lend} + 1);
								};
	my $are_chainable_sref = sub { my ($eleA, $eleB) = @_;
								   if ($eleB->{draw_lend} - $eleA->{draw_rend} >= $padding_href->{pad_left}) { return (1); } else { return(0); } 
							   };

	my %Eles;
	foreach my $ele (@eles) {
		$Eles{ $ele->{scaffold} } = $ele;
	}

	my @scaff_tiers;
	
	while (%Eles) {
		my @eles = sort {$a->{draw_lend}<=>$b->{draw_lend}} values %Eles;
		
		my @tier;
		my @chain = &DPchain::find_highest_scoring_chain(\@eles, $get_base_score_sref, $are_chainable_sref);
		foreach my $ele (@chain) {
			my $scaffold = $ele->{scaffold};
			push (@tier, $scaffold);
			delete $Eles{$scaffold};
		}
		push (@scaff_tiers, [@tier]);
	}
		
	return (@scaff_tiers);
}
	

####
sub estimate_synteny_orientation {
	my ($syn_scaffold, $syn_coords_aref, 
		$ref_scaffold, $ref_coords_aref, 
		$scaffold_to_gene_structs_href, 
		$syn_gene_pairs_href,
		$gene_acc_to_gene_struct_href) = @_;
	
	my ($syn_lend, $syn_rend) = @$syn_coords_aref;
	my ($ref_lend, $ref_rend) = @$ref_coords_aref;
	
	my %orient_counts;
	
	my @orient_orders;

	foreach my $syn_gene (sort {$a->{lend}<=>$b->{lend}} @{$scaffold_to_gene_structs_href->{$syn_scaffold}}) {
		my $syn_acc = $syn_gene->{acc};
		my $syn_gene_lend = $syn_gene->{lend};
		my $syn_gene_rend = $syn_gene->{rend};
		unless ($syn_gene_lend >= $syn_lend && $syn_gene_rend <= $syn_rend) {
			next; # out of range
		}

		
		if (my $other_genes_href = $syn_gene_pairs_href->{$syn_acc}) {
			
			foreach my $other_gene_acc (keys %$other_genes_href) {
				my $other_gene = $gene_acc_to_gene_struct_href->{$other_gene_acc};

				if ($other_gene->{scaffold} eq $ref_scaffold &&
					$other_gene->{lend} >= $ref_lend && $other_gene->{rend} <= $ref_rend) {
					# within range
					
					my $orient_compare = ($syn_gene->{orient} eq $other_gene->{orient}) ? '+' : '-';
					$orient_counts{$orient_compare}++;
					push (@orient_orders, $orient_compare);
				}
			}
		}
	}

	## see if first and last gene are in the same orient.  If so, assume it.
	if (0 && $orient_orders[0] eq $orient_orders[$#orient_orders]) { ## taking majority vote for now, exclusively.  should make this a future option.
		return ($orient_orders[0]);
	}
	else {
		## take a majority vote across the region:

		my @orients = sort {$orient_counts{$a}<=>$orient_counts{$b}} keys %orient_counts;
		my $most_supported_orient = pop @orients;
		
		return ($most_supported_orient);
	}
}
					


####
sub transpose_scaffold_genes_list_to_gene_acc_lookup {
	my ($scaffold_to_gene_structs_href) = @_;
	
	my %gene_acc_to_gene_struct;

	foreach my $gene_list_aref (values %$scaffold_to_gene_structs_href) {
		foreach my $gene_struct (@$gene_list_aref) {
			my $acc = $gene_struct->{acc};
			$gene_acc_to_gene_struct{$acc} = $gene_struct;
			
		}
	}

	return (%gene_acc_to_gene_struct);
}

####
sub assign_syntenic_contig_seq_ranges {
	my ($syntenic_scaffolds_to_genes_href, $syn_org_to_scaffolds_href, $syn_gene_pairs_href, $syn_contig_to_seq_range_href, $gene_acc_to_gene_struct_href, $refScaffCoords_aref) = @_;
	
	# store the syntenic scaffold info
	foreach my $syn_scaffold (keys %$syntenic_scaffolds_to_genes_href) {
		my @syn_genes = sort {$a->{rend}<=>$b->{rend}} @{$syntenic_scaffolds_to_genes_href->{$syn_scaffold}};
		my $syn_lend = $syn_genes[0]->{lend};
		my $syn_rend = $syn_genes[$#syn_genes]->{rend};
		
		my ($ref_lend, $ref_rend) = &compute_reference_synteny_range(\@syn_genes, $syn_gene_pairs_href, $gene_acc_to_gene_struct_href, $refScaffold, $refScaffCoords_aref);
		
		$syn_contig_to_seq_range_href->{$syn_scaffold} = { lend => $syn_lend, 
														   rend => $syn_rend,
														   ref_lend => $ref_lend,
														   ref_rend => $ref_rend,
													   };
		
		my ($org, $scaff_name) = split (/;/, $syn_scaffold);
		push (@{$syn_org_to_scaffolds_href->{$org}}, $syn_scaffold);
		
	}
	
	
	return;

}


####
sub organize_order_of_organism_tiers {
	my ($syn_org_to_scaffolds_href, $ref_org, $refScaffold, $orgsOrder, $top_orgs_aref, $bottom_orgs_aref) = @_;
	
	
	if ($orgsOrder) {
		$orgsOrder =~ s/\s+//g; 
		my @orgs = split (/,/, $orgsOrder);
		my $found_ref_org = 0;
		foreach my $org (@orgs) {
		  if ($org eq $ref_org) {
			$found_ref_org = 1;
		  }
		  elsif (! $found_ref_org) {
			push (@$top_orgs_aref, $org);
		  }
		  else {
			push (@$bottom_orgs_aref, $org);
		  }
		}
	  }
	else {
	  ## just take one and make it the top org:
	  my @syn_orgs = grep { $_ !~ /\b$ref_org\b/ } keys %$syn_org_to_scaffolds_href;
	  my $top = shift @syn_orgs;
	  @$top_orgs_aref = ($top);
	  @$bottom_orgs_aref = @syn_orgs;
	}
	
	return;
}

####
sub assign_reference_scaffold_range_coordinates {
	my ($refScaffold, $refScaffoldRange, $scaffold_lengths_href) = @_;
	
	my ($refScaffLend, $refScaffRend);
	
	if ($refScaffoldRange) {
		$refScaffoldRange =~ s/\s//g;
		($refScaffLend, $refScaffRend) = split (/-/, $refScaffoldRange);
		unless ($refScaffLend =~ /^\d+$/ && $refScaffRend =~ /^\d+$/ 
				&& $refScaffLend < $refScaffRend) {
			die "Error, cannot properly parse reference scaffold coordinates: $refScaffoldRange";
		}
	}
	else {
		($refScaffLend, $refScaffRend) = (1, $scaffold_lengths_href->{$refScaffold});
		unless ($refScaffRend) {
			die "Error, no genes parsed from $refScaffold reference scaffold!";
		}
	}
	
	return($refScaffLend, $refScaffRend);
}


####
sub parse_misc_features_in_range {
	my ($miscFeaturesGff3, $syn_contig_to_seq_range_href) = @_;

	my %scaffold_to_misc_features;
	
	foreach my $file (split (/,/, $miscFeaturesGff3)) {
		$file =~ s/\s//g;

		open (my $fh, $file) or die "Error, cannot open file $file";
		while (<$fh>) {
			my $line = $_;
			chomp;
			if (/^\#/) { next; }
			unless (/\w/) { next; }
			my @x = split (/\t/);
			my ($scaffold, $lend, $rend, $orient, $name_info) = ($x[0], $x[3], $x[4], $x[6], $x[8]);
			
			if (my $range_info_href = $syn_contig_to_seq_range_href->{$scaffold}) {
				
				my ($range_lend, $range_rend) = ($range_info_href->{lend}, $range_info_href->{rend});
				if ($lend < $range_rend && $rend > $range_lend) {
					## within range:
					
					push (@{$scaffold_to_misc_features{$scaffold}->{$file}}, { name => $name_info,
																			   lend => $lend,
																			   rend => $rend,
																			   orient => $orient,
																			   scaffold => $scaffold,
																		   } );
					
					# print STDERR "Misc feat in range: $line";
				}
			}
		}
		close $fh;
	}

	
	return (%scaffold_to_misc_features);
}



																	

