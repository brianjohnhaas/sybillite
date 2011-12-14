#!/usr/local/bin/perl -w

# ---------------------------------------------------------------
# Script that demonstrates how to use MultiPanel.pm, a 
# generic module for rendering "matches" between features 
# in distinct Bio::Graphics::Panels.
#
# Run the script like so:
#  ./testMultiPanel.pl > test.png
#
# If you get a compile-time error somewhere in the gd() method
# at the beginning of MultiPanel.pm, it's probably because 
# you're not using the same version of Bioperl (1.4, I think)
# from which I lifted the method.  If this happens, either
# upgrade/downgrade to 1.4, or replace the gd() method in 
# MultiPanel.pm with one that matches your Bioperl version.
# Look for the keyword "SYBIL" to see which part of the 
# method needs to be altered to work with MultiPanel.pm
#
# Created: Thu Mar 10 23:19:27 EST 2005
#
# Jonathan Crabtree
# ---------------------------------------------------------------

# bioperl modules
use Bio::Graphics::Panel;
use Bio::SeqFeature::Generic;

# module used to draw matches between bioperl features in different panels
use MultiPanel;

# STEP 1: Create as many Bio::Graphics::Panels as your application requires,
# using pad_top, pad_left, pad_right, pad_bottom, width, and height to 
# position them in a non-overlapping fashion.  In the following code, for
# example, two panels are created, with their pad_top values set to ensure
# that the second (or at least its content) appears below the first.
#
my $PANEL_VGAP = 100;
my $yposn = 50;

# first panel
#
my $panel1 = Bio::Graphics::Panel->new(-length => 1000, -offset => -100, -width => 600, -pad_top => $yposn, -pad_left => 50);
my $scale1 = Bio::SeqFeature::Generic->new(-primary_tag => 'misc', -start => -100, -end => 900);
my $feat1 = Bio::SeqFeature::Generic->new(-primary_tag => 'misc', -start => -75, -end => 100);
my $feat2 = Bio::SeqFeature::Generic->new(-primary_tag => 'misc', -start => 200, -end => 350);
$panel1->add_track($scale1, -glyph => 'anchored_arrow', -tick => 2, );
$panel1->add_track([$feat1, $feat2], -glyph => 'generic', );

$yposn += $panel1->height() - $yposn + $PANEL_VGAP;

# second panel
#
my $panel2 = Bio::Graphics::Panel->new(-length => 1000, -offset => 2000, -width => 600, -pad_top => $yposn, -pad_right => 50, -pad_bottom => 100);
my $scale2 = Bio::SeqFeature::Generic->new(-primary_tag => 'misc', -start => 0, -end => 3000);
my $feat10 = Bio::SeqFeature::Generic->new(-primary_tag => 'misc', -start => 2500, -end => 2750);
my $feat11 = Bio::SeqFeature::Generic->new(-primary_tag => 'misc', -start => 2200, -end => 2400);
$panel2->add_track($scale2, -glyph => 'anchored_arrow', -tick => 2, );
$panel2->add_track([$feat10, $feat11], -glyph => 'generic', );

# STEP 2: Construct an array of "matches" that describe which pairs of features 
# correspond to one another.  Typically these pairs of features will be
# in different panels.  A "match" is simply a hashref with the following keys,
# some of which are optional:
#
# feat1       =>  a Bio::SeqFeatureI object that was added to one of the panels
# feat2       =>  the Bio::SeqFeatureI object that matches/is linked with feat1
# x1,y1(o)    =>  a subsequence of feat1, expressed in the parent coordinate system 
# x2,y2(o)    =>  a subsequence of feat1, expressed in the parent coordinate system
# bg          =>  fill/background color for the match
#                 either a bioperl symbolic color (e.g., 'red','pink','blue') or a
#                 hex-encoded color string of the form '#ff00f0'
# fg          =>  outline/foreground color for the match
# reverse     =>  whether the match should be drawn as a reverse match (see example image)
# priority(o) =>  matches with higher priority are drawn first  (i.e., they will appear *below* the others)
#
# (o) denotes optional parameters
#
# Note also that MultiPanel.pm does very little error and bounds checking.  x1,y1,x2,
# and y2 may, for example, extend outside of the range occupied by feat1 and feat2.

# a match between $feat2 (panel 1) and $feat10 (panel 2)
# no coordinates (x1,y1,x2,y2) are given, so *all* of $feat2 matches all of $feat10
#
my $m1 = { 'feat1' => $feat2, 'feat2' => $feat10, 'bg' => 'pink', 'fg' => 'black', 'reverse' => 0 };

# a reverse match between a portion of $feat1 (panel 1) and a portion of $feat11 (panel 2)
#
my $m2 = { 'feat1' => $feat1, 'feat2' => $feat11, 'bg' => 'lightblue', 'fg' => '#0000ff', 'reverse' => 1,
	   'x1'=> 25, 'y1' => 100, 'x2' => 2200, 'y2' => 2300 };

# listref of matches
#
my $matches = [ $m1, $m2 ];

# STEP 3: Use the listref of Bio::Graphics::Panels and the listref of matches to
# create a new MultiPanel.
#
my $multiPanel = new MultiPanel([$panel1, $panel2], $matches);

# STEP 4: Render the image.
#
my $gd = $multiPanel->gd();
print $gd->png();

# Note that there are some additional restrictions on the component panels passed
# to the MultiPanel.  First, each of the panels must use the same image_class
# (which is "GD" by default; in the above example the image_class is not specified
# so the condition is met.)  Second, setting a background or tkcolor for any of
# the tracks in the component Panels will overwrite the matches drawn by the 
# MultiPanel.  If the image needs a colored background, you'll need to modify the
# source code for MultiPanel.pm.
