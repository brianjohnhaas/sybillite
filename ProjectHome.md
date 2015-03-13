# SybilLite #

SybilLite is a web interface used to display comparative genomics data.  It can be used to display output from tools such as [DagChainer](http://sourceforge.net/projects/dagchainer/).

The inputs expected are GFF3 files for each of the genomes in your comparison as well as a simple tab-delimited file describing the alignments generated between them.  An indexing utility provided reduces these to a SQLite3 database on disk, which is read by the interface.

Here is an example of the front page of the demo project:

<div>
<blockquote><img src='http://sybillite.googlecode.com/svn/trunk/docs/screenshots/project_view.png' />
</div></blockquote>


Users can then manually enter an organism, molecule, and range to use as the reference for display, or choose one of the user-pre-defined "regions of interest".  This brings up the primary comparative display, for example:

<div>
<blockquote><img src='http://sybillite.googlecode.com/svn/trunk/docs/screenshots/example_inversion.png' />
</div>