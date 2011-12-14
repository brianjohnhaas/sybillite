#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use URI::Escape;
use FindBin;
use lib ("$FindBin::Bin/PerlLib");
use IniReader;

# prepare for graphics drawing.
use Bio::Graphics::Panel;
use Bio::SeqFeature::Generic;
use MultiPanel;


main: {

	my $usage = "usage: $0 synView.conf\n\n";
	
	my $conf_file = $ARGV[0] or die $usage;
	
	my $conf = new IniReader($conf_file);

	my $annotations_gff3 = $conf->get_value("main", "annotations") or die "Error, no annotations file specified in conf file";
	
	my %scaffold_to_gene_structs = &parse_gene_coords($annotations_gff3);
	
	
	# transpose the scaffold to gene data
	my %gene_acc_to_gene_struct = &transpose_scaffold_genes_list_to_gene_acc_lookup(\%scaffold_to_gene_structs);


	my $match_pairs_file = $conf->get_value("main", "connected_gene_pairs") or die "Error, no connected_gene_pairs file specified in conf file";
	my %match_pairs = &parse_match_pairs($match_pairs_file);
	
	# get scaffold lengths (cheap; just taking the right-most coordinate of any parsed feature on the scaffold for now.
	my %scaffold_lengths = &parse_scaffold_lengths(\%scaffold_to_gene_structs);	
	
	
	#####################################################################
	## draw them out.

	
	my $PANEL_VGAP = $conf->get_value("main", "pad_top") || 0;
	my $pad_left = $conf->get_value("main", "pad_left") || 0;
	my $pad_right = $conf->get_value("main", "pad_right") || 0;
	my $pad_bottom = $conf->get_value("main", "pad_bottom") || 0;
	my $pixelsPerKb = $conf->get_value("main", "pixelsPerKb") or die "Error, must specify pixelsPerKb in conf ";
	
	my $yposn = $PANEL_VGAP;
	
	my @panels;

	my @region_tokens = grep {/^Region_/} $conf->get_section_headings();
	
	my %region_ranks;
	{ # assign region ranks
		foreach my $region (@region_tokens) {
			my $rank = $conf->get_value("$region", "rank");
			$region_ranks{$region} = $rank;
		}
	}


	my %scaff_orients;
	my %gene_acc_to_feats;
	
	my $prev_panel;

	my $tier_order = 0;

	foreach my $region (sort {$region_ranks{$a}<=>$region_ranks{$b}} @region_tokens) {
		
		$tier_order++;

		my $molecule = $conf->get_value($region, "molecule");
		my $range = $conf->get_value($region, "range");
		
		my ($lend, $rend) = split (/-/, $range);
		
		my @gene_structs = &get_genes_in_range($lend, $rend, $scaffold_to_gene_structs{$molecule});
		
		my $length = $rend - $lend + 1;
				
		$scaff_orients{$molecule} = '+'; #$syn_orient; # store for later so we know if the gene matches are same vs. inverted
		
		my ($draw_lend, $draw_rend) =  (1, int ($length * $pixelsPerKb/1000 + 0.5));
			
		my $draw_length = $draw_rend + $pad_left + $pad_right;
		
		
		## Drawing the scaffold ticker
		print STDERR "-drawing $molecule, range $lend-$rend length: $length, draw length: $draw_length, draw_lend: $draw_lend to draw_rend: $draw_rend\n";
		
		my $reverse_flag = 0; #($syn_orient eq '-') ? 1 : 0;
		my $panel = Bio::Graphics::Panel->new(-length => $length, -offset => $lend, -width => $draw_length, 
											  -pad_top => $yposn, -pad_left => $draw_lend + $pad_left, 
											  -pad_right => $pad_right, -pad_bottom => $pad_bottom, -flip => $reverse_flag,
											  -key_style => "between");
		

		push (@panels, $panel);
		
		$prev_panel = $panel;
		
		my $scale = Bio::SeqFeature::Generic->new(-primary_tag=> "$molecule-scale", -start=>$lend, -end => $rend);
		$panel->add_track($scale, -glyph => 'anchored_arrow', -tick=> 2, -key => $molecule);
		
		#my $scaff_glyph = Bio::SeqFeature::Generic->new(-primary_tag=> $scaffold, -start=>$rangeLend, -end => $rangeRend);
		#$panel->add_track([$scaff_glyph], -glyph => 'generic');
		
		#####################################
		## Draw each of the gene features
		my @feats;
		foreach my $gene (@gene_structs) {
			my ($lend, $rend) = ($gene->{lend}, $gene->{rend});
			
			my $name = $gene->{alias_id};
			
			my $gene_orient = $gene->{orient};
			my $strand = ($gene_orient eq '+') ? 1:-1; #($syn_orient eq $gene_orient) ? 1 : -1;
			
			my $acc = $gene->{acc};
			my $feat = Bio::Graphics::Feature->new(-primary_tag => $acc, -start => $lend, -end => $rend, -strand => $strand, -name => $name, -type => 'gene');
			
			$feat->{__tier_order} = $tier_order; # bundling a hidden attribute

			$gene_acc_to_feats{$acc} = $feat;
			push (@feats, $feat);
			#$panel->add_track($feat, -glyph => 'generic', -label => 1, -description => 1, -key => $acc);
		}
		$panel->add_track([@feats], -glyph => "transcript", -label => 1);
		
						
		## stack panels between tiers and organisms
		$yposn = $prev_panel->height() + $PANEL_VGAP;
		
	} # end scaff tier

	
	
	###########################################
	## describe the matches between genes:
	
	my @matches;
	#goto skip_matches;

	my %seen;
	my %orient_swap = ( '+' => '-', '-' => '+');
	
	my $match_draw_style = $conf->get_value("main", "matchDrawing");

	foreach my $gene_acc (keys %gene_acc_to_feats) {
		
		if (exists ($match_pairs{$gene_acc}) ) {
			my @matching_accs = keys (%{$match_pairs{$gene_acc}});
			foreach my $matching_acc (@matching_accs) {
				
				if ($matching_acc eq $gene_acc) { next; } # no self matches
				
				unless (exists ($gene_acc_to_feats{$matching_acc})) { next; }
								
				my $pair = join ("_", sort ($gene_acc, $matching_acc));
				if ($seen{$pair}) { next; }
				$seen{$pair} = 1;
				
				my $featA = $gene_acc_to_feats{$gene_acc};
				my $featB = $gene_acc_to_feats{$matching_acc};
				
				unless ($featA && $featB) { die "Error, no feats stored for both $gene_acc and $matching_acc matching pairs"; } # not within selected range.
				
				my $tier_diff = abs ($featA->{__tier_order} - $featB->{__tier_order});
				if ($tier_diff == 0) { next; } # no show matches on same tier level.
				
				if ($match_draw_style eq "top_down" && $tier_diff != 1) { next; }
				

				my $geneA = $gene_acc_to_gene_struct{$gene_acc};
				my $geneB = $gene_acc_to_gene_struct{$matching_acc};
				
				my $scaff_A = $geneA->{scaffold};
				my $scaff_B = $geneB->{scaffold};
				
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
	}
	
    skip_matches:
	
	my $multi_panel = new MultiPanel([@panels], [@matches]);
	
	my $gd = $multi_panel->gd();
	print $gd->png();
	
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

	foreach my $file (split (/,/, $fileListString)) {
		
		open (my $fh, $file) or die "Error, cannot find file $file";
		while (<$fh>) {
			chomp;
			if (/^\#/) { next; }
			unless (/\w/) { next; }
			
			my @x = split (/\t/);
			if ($x[2] eq 'gene') {
				my ($scaffold, $lend, $rend, $orient, $gene_info) = ($x[0], $x[3], $x[4], $x[6], $x[8]);
								
				$gene_info =~ /ID=([^\;\s]+)/ or die "Error, no gene ID for $_";
				my $gene_id = $1;
				my $alias_id;
				if ($gene_info =~ /Alias=([^;\s]+)/) {
					$alias_id = $1;
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
																 alias_id => $alias_id,
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
			
			if (my $alias = $gene_struct->{alias_id}) {
				$gene_acc_to_gene_struct{$alias} = $gene_struct;
			}
			

		}
	}

	return (%gene_acc_to_gene_struct);
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
			my ($gene_id, $alias, $name);
			
			
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


####
sub parse_match_pairs {
	my ($file) = @_;

	my %pairs;

	open (my $fh, $file) or die $!;
	while (<$fh>) {
		chomp;
		my ($accA, $accB) = split (/\s+/);
		
		$pairs{$accA}->{$accB} = 1;
		$pairs{$accB}->{$accA} = 1;
	}
	close $fh;

	return(%pairs);
}

####
sub get_genes_in_range {
	my ($lend, $rend, $gene_list_aref) = @_;

	my @genes;
	
	if (ref $gene_list_aref) {
		foreach my $gene (@$gene_list_aref) {
			if ($gene->{lend} >= $lend && $gene->{rend} <= $rend) {
				push (@genes, $gene);
			}
		}

	}

	return(@genes);
}
