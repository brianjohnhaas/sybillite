=head1 NAME

Sybil::Graphics::MultiPanel

=head1 SYNOPSIS

Draws the contents of several Bio::Graphics::Panels into a single image, 
optionally highlighting links/alignments between features in different 
Panels.

=head1 DESCRIPTION

NOTE: this package relies on a small change to the source code for 
Bio::Graphics::Panel, which allows the same GD::Image to be passed 
to the draw() method of multiple Bio::Graphics::Panel objects:

  # JC: changed colorAllocate to colorResolve to allow same $gd to be reused
  my $idx = $gd->colorResolve(@{$COLORS{$name}});

(line 382 of Bio/Graphics/Panel.pm in bioperl 1.302)

Note that Perl will generate a warning for the redefinition of the gd() 
method when run with -w.

=cut

use Bio::Graphics::Panel;
package Bio::Graphics::Panel;

# The following method overrides Bio::Graphics::Panel::gd, replacing calls 
# to colorAllocate with equivalent calls to colorResolve.  This allows the 
# same GD::Image object to be safely reused across several calls to gd()
# without disrupting the color palette.
#
# This method was copied verbatim from the source for Bio::Graphics::Panel
# in Bioperl 1.4.  The Bio::Graphics package was written by Lincoln Stein
# and is copyrighted by Cold Spring Harbor Laboratory.  As per the Artistic
# License, our changes are annotated below, with the keyword "SYBIL":
#
sub gd {
  my $self        = shift;
  my $existing_gd = shift;

  local $^W = 0;  # can't track down the uninitialized variable warning

  return $self->{gd} if $self->{gd};

  $self->setup_fonts;

  unless ($existing_gd) {
    # Encapsulating this eval within a BEGIN block
    # adds nothing since $image_class is undefined at compile time.
    # gd exported functions should all use ();
    # BEGIN {
    my $image_class = $self->image_class;
    eval "use $image_class; 1" or die $@;
    # }
  }

  my $width  = $self->width;
  my $height = $self->height;

  my $pkg = $self->image_package;
  my $gd  = $existing_gd || $pkg->new($width,$height,
				      ($self->{truecolor} && $pkg->can('isTrueColor') ? 1 : ())
				     );

  my %translation_table;

  # SYBIL: edits begin here
  my @colorNames = $self->color_names();
  for my $name ('white','black',@colorNames) {
      # calling color_name_to_rgb and color_names (above) to avoid accessing %COLOR, 
      # a lexically-scoped global in the original package
      my @rgb = $self->color_name_to_rgb($name);
      # changed colorAllocate to colorResolve to allow same $gd to be reused
    my $idx = $gd->colorResolve(@rgb);
    $translation_table{$name} = $idx;
  }
  # SYBIL: edits end here

  $self->{translations} = \%translation_table;
  $self->{gd}           = $gd;
  if ($self->bgcolor) {
    $gd->fill(0,0,$self->bgcolor);
  } elsif (eval {$gd->isTrueColor}) {
    $gd->fill(0,0,$translation_table{'white'});
  }

  my $pl = $self->pad_left;
  my $pt = $self->pad_top;
  my $offset = $pt;
  my $keyheight   = $self->{key_font}->height;
  my $bottom_key  = $self->{key_style} eq 'bottom';
  my $between_key = $self->{key_style} eq 'between';
  my $left_key    = $self->{key_style} eq 'left';
  my $right_key   = $self->{key_style} eq 'right';
  my $empty_track_style = $self->empty_track_style;
  my $spacing = $self->spacing;

  # we draw in two steps, once for background of tracks, and once for
  # the contents.  This allows the grid to sit on top of the track background.
  for my $track (@{$self->{tracks}}) {
    my $draw_between = $between_key && $track->option('key');
    next if !$track->parts && ($empty_track_style eq 'suppress'
			   or  $empty_track_style eq 'key' && $bottom_key);
    $gd->filledRectangle($pl,
			 $offset,
			 $width-$self->pad_right,
			 $offset+$track->layout_height
			 + ($between_key ? $self->{key_font}->height : 0),
			 $track->tkcolor)
      if defined $track->tkcolor;
    $offset += $keyheight if $draw_between;
    $offset += $track->layout_height + $spacing;
  }

  $self->draw_grid($gd)  if $self->{grid};

  $offset = $pt;
  for my $track (@{$self->{tracks}}) {
    my $draw_between = $between_key && $track->option('key');
    my $has_parts = $track->parts;
    next if !$has_parts && ($empty_track_style eq 'suppress'
			or  $empty_track_style eq 'key' && $bottom_key);

    if ($draw_between) {
      $offset += $self->draw_between_key($gd,$track,$offset);
    }

    elsif ($self->{key_style} =~ /^(left|right)$/) {
      $self->draw_side_key($gd,$track,$offset,$self->{key_style});
    }

    $self->draw_empty($gd,$offset,$empty_track_style)
      if !$has_parts && $empty_track_style=~/^(line|dashed)$/;

    $track->draw($gd,0,$offset,0,1);
    $self->track_position($track,$offset);
    $offset += $track->layout_height + $spacing;
  }


  $self->draw_bottom_key($gd,$pl,$offset) if $self->{key_style} eq 'bottom';
  return $self->{gd} = $gd;
}

# We now return you to your regularly-scheduled Perl package:

package MultiPanel;

# ------------------------------------------------------------------
# Globals
# ------------------------------------------------------------------

my $Y_OFFSET = 2;

# ------------------------------------------------------------------
# Constructor
# ------------------------------------------------------------------

# $panels   - listref of Bio::Graphics::Panel
# $matches  - listref of hashrefs with the following keys:
#
sub new {
    my($invocant, $panels, $matches) = @_;
    my $class = ref($invocant) || $invocant;
    my $self = {};
    bless($self, $class);

    $self->{'panels'} = $panels;
    $self->{'matches'} = $matches;

    return $self;
}

# ------------------------------------------------------------------
# Public methods
# ------------------------------------------------------------------

sub gd {
    my($self) = @_;
   
    # create a GD::Image large enough to hold all of the component Panels
    my($width, $height) = $self->_getBounds();
    my $gdImageClass = $self->_getGdImageClass();
    my $gdImagePackage = $gdImageClass . '::Image';
    my $gdImg = $gdImagePackage->new($width, $height);

    # use a dummy panel to set the color palette of the GD::Image
    my $dummyPanel = Bio::Graphics::Panel->new(-length => 100, -offset => 0, -width => 10 );
    $dummyPanel->gd($gdImg);

    # draw matches/links between objects in the different panels
    $self->_drawMatches($gdImg);

    # draw panels
    map { $_->gd($gdImg); } @{$self->{'panels'}};
    
    return $gdImg;
}

# ------------------------------------------------------------------
# "Private" methods
# ------------------------------------------------------------------

sub _getBounds {
    my($self) = @_;
    my $panels = $self->{'panels'};
    my $width = undef;
    my $height = undef;

    foreach my $panel (@$panels) {
	my $w = $panel->width();
	my $h = $panel->height();
	$width = $w if (!defined($width) || ($w > $width));
	$height = $h if (!defined($height) || ($h > $height));
    }

    return ($width, $height);
}

sub _getGdImageClass {
    my($self) = @_;
    my $panels = $self->{'panels'};
    my $class = undef;

    # check that all the panels are using the same GD class
    foreach my $panel (@$panels) {
	my $ic = $panel->image_class();
	if (defined($class)) {
	    die "all Panels must have the same image_class()" unless ($class eq $ic);
	} else {
	    $class = $ic;
	}
    }
    return $class;
}

sub _getBoundingBoxes {
    my($self, $feature) = @_;
    my $feat2Boxes = $self->_getFeatToBoundingBoxHash();
    return $feat2Boxes->{$feature};
}

sub _resetFeatToBoundingBoxHash {
    my($self) = @_;
    $self->{'feat_2_boxes'} = undef;
}

sub _getFeatToBoundingBoxHash {
    my($self) = @_;

    if (!defined( $self->{'feat_2_boxes'})) {
	# cache mapping from $feature to bounding boxes
	my $feat2Boxes = $self->{'feat_2_boxes'} = {};
	my $panels = $self->{'panels'};
	my $np = scalar(@$panels);

	for (my $p = 0;$p < $np;++$p) {
	    my $panel = $panels->[$p];

	    # boxes() returns a listref of [ $feature, $x1, $y1, $x2, $y2, $track ]
	    my $boxes = $panel->boxes();
	    
	    foreach my $box (@$boxes) {
		my $list = $feat2Boxes->{$box->[0]};
		$list = $feat2Boxes->{$box->[0]} = [] if (!defined($list));
		push(@$list, {'panel' => $p, 'box' => $box});
	    }
	}
    }

    return $self->{'feat_2_boxes'};
}

sub _gdColorIndex {
    my($self, $gdImg, $color) = @_;
    my $panel = $self->{'panels'}->[0];
    my @rgb;

    # NOTE: we cannot call translate_color at this point because the Panel
    # will not have initialized its translation table.

    # a. hex color name
    if ($color =~ /^\#/) {
	my @hexColors = ($color =~ /^\#?(\S\S)(\S\S)(\S\S)$/);
	@rgb = map { hex($_) } @hexColors;
    } 

    # b. assume it's a bioperl color name
    else {
	@rgb = $panel->color_name_to_rgb($color);
    }
    return $gdImg->colorResolve(@rgb);
}

sub _drawMatches {
    my($self, $gdImg) = @_;
    my $imgClass = $self->_getGdImageClass();
    my $gdPolygonPackage = $imgClass . '::Polygon';
    my $panels = $self->{'panels'};

    # sort matches by priority (default = 0; higher numbers are drawn first)
    my $matches = $self->{'matches'};
    foreach my $match (@$matches) {
	$match->{'priority'} = 0 if (!defined($match->{'priority'}));
    }
    my @sortedMatches = sort { $b->{'priority'} <=> $a->{'priority'} } @$matches;
    
    foreach my $match (@sortedMatches) {
	my($feat1,$feat2,$x1,$y1,$x2,$y2,$bg,$fg,$reverse) = 
	    map {$match->{$_}} ('feat1','feat2','x1','y1','x2','y2','bg','fg','reverse' );

	my $bgColor = $self->_gdColorIndex($gdImg, $bg);
	my $fgColor = $self->_gdColorIndex($gdImg, $fg);

	my $bbs1 = $self->_getBoundingBoxes($feat1);
	my $bbs2 = $self->_getBoundingBoxes($feat2);

	# since the same feature may appear in multiple Panels, matches are drawn between
	# all possible pairs of features
	#
	foreach my $bb1 (@$bbs1) {
	    foreach my $bb2 (@$bbs2) {
		my $f1upper = ($bb1->{'box'}->[4] < $bb2->{'box'}->[2]);
		my $upper = $f1upper ? $bb1 : $bb2;
		my $lower = $f1upper ? $bb2 : $bb1;
		my $uc = $f1upper ? [$x1,$y1] : [$x2,$y2];
		my $lc = $f1upper ? [$x2,$y2] : [$x1,$y1];

		my $b1 = $upper->{'box'};
		my($f1,$ux1,$uy1,$ux2,$uy2,$t1) = @$b1;
		my $p1 = $upper->{'panel'};

		my $b2 = $lower->{'box'};
		my($f2,$lx1,$ly1,$lx2,$ly2,$t2) = @$b2;
		my $p2 = $lower->{'panel'};

		# adjust coordinates according to $uc, $lc (if defined)
		($ux1) = $panels->[$p1]->location2pixel($uc->[0]) if (defined($uc->[0]));
		($ux2) = $panels->[$p1]->location2pixel($uc->[1]) if (defined($uc->[1]));

		($lx1) = $panels->[$p2]->location2pixel($lc->[0]) if (defined($lc->[0]));
		($lx2) = $panels->[$p2]->location2pixel($lc->[1]) if (defined($lc->[1]));

		# draw match
		my $poly = $gdPolygonPackage->new();
		$poly->addPt($ux1, $uy2+$Y_OFFSET);
		$poly->addPt($ux2, $uy2+$Y_OFFSET);

		if ($reverse) {
		    $poly->addPt($lx1, $ly1-$Y_OFFSET);
		    $poly->addPt($lx2, $ly1-$Y_OFFSET);
		} 
		else {
		    $poly->addPt($lx2, $ly1-$Y_OFFSET);
		    $poly->addPt($lx1, $ly1-$Y_OFFSET);
		}

		$gdImg->filledPolygon($poly, $bgColor);
		$gdImg->polygon($poly, $fgColor);
	    }
	}
    }
}

1;
