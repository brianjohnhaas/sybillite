# Introduction #

SybilLite is a lightweight web-based comparative genomics viewer.  It uses GFF3 as input from tools such as [DAGchainer](http://sourceforge.net/projects/dagchainer/), which identifies chains of gene pairs sharing conserved order between genomic regions, by identifying paths through a directed acyclic graph (DAG).


# Installation #

The [Installation Guide](InstallationGuide.md) provides both a quick setup and detailed tutorial, depending on your preference.  The detailed setup lists the required software as well as example commands and configuration file changes needed to get everything running.


# Input data #

SybilLite uses a SQLite database to hold all the comparative data for the display, but comes with scripts to generate this from GFF3.  The [Input Data Guide](InputDataGuide.md) has more details about this process, including details about what exactly is expected in the GFF3 and how to convert it to a SQLite file.


# Using the interface #

SybilLite comes with an example comparative dataset from a small collection of [Schistosoma genomes](http://en.wikipedia.org/wiki/Schistosoma).  The [Interface Walkthrough](InterfaceWalkthrough.md) guides new users through the features of the interface using this example dataset.


# Getting help #

There are a few different resources available to you if you need help with SybilLite.  If you have a question, you can write the [User Support Group](http://groups.google.com/group/sybillite-users).  For bugs or feature requests, please create a ticket in the [Issue tracking system](http://code.google.com/p/sybillite/issues/list).