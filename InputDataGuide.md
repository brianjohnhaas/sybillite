# Overview #

SybilLite relies on a combination of file formats organized following a specific directory and file name convention.  No relational databases are currently used.

The SybilLite interface is driven by the concept of a 'project', and each project has its own directory and set of files to support it.  Adding a new project to the system is as easy as creating a new folder and populating it correctly.  How this is done is described below, but the basic starting point is a GFF3 file.

For those who prefer to dive right in feel free to use the example Schizos project directory and files as your template.  Below you'll find detailed descriptions of the things found there.


## Directory layout convention ##

In the root of the install you'll find a 'data' directory.  Under here there is one directory per project, which is usually one group of organisms for which you have comparative data.  For our example project, you'll find the following:

```
    data/Schizos/
    data/Schizos/annotation.db
    data/Schizos/data.pairs.aligncoords
    data/Schizos/regions_of_interest/
    data/Schizos/SJ1_CALLGENES_Final_2.gff3
    data/Schizos/SO3_CALLGENES_FINAL_1.gff3
    data/Schizos/SP2_Sanger_072008_1.gff3
```


## Expected GFF3 format ##


## Converting to SQLite3 ##


## SQLite3 schema ##


## Alignment data ##

Information on the generation of the data.pairs.aligncoords file from DAGChainer goes here.

## Regions of Interest graphics ##
