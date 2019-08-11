# Assembly Dereplicator

__Ryan R. Wick<sup>1</sup> and Kathryn E. Holt<sup>1,2</sup>__
<br>
<sub>1. Department of Infectious Diseases, Central Clinical School, Monash University, Melbourne, Victoria 3004, Australia<br>2. London School of Hygiene & Tropical Medicine, London WC1E 7HT, UK</sub>

[![Build Status](https://travis-ci.com/rrwick/Assembly-dereplicator.svg?branch=master)](https://travis-ci.com/rrwick/Assembly-dereplicator) [![License GPL v3](https://img.shields.io/badge/license-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3365572.svg)](https://doi.org/10.5281/zenodo.3365572)



## Table of contents

* [Introduction](#introduction)
* [Requirements](#requirements)
* [Installation](#installation)
* [Method](#method)
   * [Clustering](#clustering)
   * [Batches](#batches)
* [Quick usage](#quick-usage)
* [Full usage](#full-usage)
* [License](#license)



## Introduction

This repo contains a stand-alone Python script ([`dereplicator.py`](dereplicator.py)) to solve a problem I occasionally run into: dereplicating a group of bacterial genome assemblies. Dereplication means removing assemblies for which there are sufficiently close relatives (defined by a distance threshold), resulting in a smaller set where the assemblies are more unique.

<p align="center"><img src="images/trees.png" alt="Trees" width="70%"></p>

As an example, imagine you have 10000 genome assemblies for a particular taxon and want to do some analysis on them, maybe building a pan genome. You know there is redundancy in this set because some of the genomes come from outbreaks and are nearly identical to each other. So instead of doing the analysis on all 10000 assemblies, you can dereplicate them to a smaller set (i.e. remove near-identical redundant genomes) so your analysis will be faster.



## Requirements

You'll need Python 3.5 or later to run `dereplicator.py`. It's a stand-alone script and has no Python package dependencies (i.e. it only uses the standard library).

There is one external dependency, [Mash](https://github.com/marbl/Mash), which you'll need installed and callable on your command line. If you can run `mash sketch -h` and `mash dist -h` without getting an error, you should be good to go!



## Installation

Since Assembly Dereplicator is just a single script, no installation is required. You can simply clone it and run it:
```
git clone https://github.com/rrwick/Assembly-dereplicator
Assembly-dereplicator/dereplicator.py --help
```

If you plan on using it often, you can copy it to someplace in your PATH variable for easier access:
```
git clone https://github.com/rrwick/Assembly-dereplicator
cp Assembly-dereplicator/dereplicator.py ~/.local/bin
dereplicator.py --help
```

If you'd like to double-check that everything works as intended, you can run this repo's [automated tests](test).



## Method

### Clustering

The basic idea of this script is to cluster assemblies and choose a single representative assembly per cluster. It first uses Mash to get all-vs-all distances between your assemblies. It then builds a graph where the assemblies are vertices and any two assemblies below the distance threshold are connected by an edge. Each [connected component](https://en.wikipedia.org/wiki/Component_(graph_theory)) of the graph is a cluster, and the assembly in that cluster with the largest N50 is chosen as the representative (so complete assemblies are preferred over draft assemblies). This graph-based approach is equivalent to [single-linkage clustering](https://en.wikipedia.org/wiki/Single-linkage_clustering).


### Batches

Despite Mash being very efficient, an all-vs-all comparison can be cumbersome when you have a large number of assemblies (e.g. ten thousand or more). This script therefore clusters in smaller randomly-selected batches (500 at a time by default but configurable with `--batch_size`). This is repeated iteratively until an iteration fails to remove any assemblies. When this happens, the script assumes that most replication has been cleared out and conducts a final dereplication round with all remaining assemblies.

I added this batch-based dereplication mainly to help with RAM usage on large datasets, especially if the assembly collection contains a lot of redundancy. But it can also help a bit with a common criticism of single-linkage clustering: long thin clusters. By clustering subsets at a time, it increases the chance that large clusters will be broken apart. Note that this means changing the batch size may change the number of dereplicated assemblies.



## Quick usage

Assembly Dereplicator is easy to run. Just give it an input directory containing your assemblies in FASTA format, an output directory and a distance threshold:
```
dereplicator.py --threshold 0.01 assembly_dir output_dir
```



## Full usage

```
usage: dereplicator.py [--threshold THRESHOLD] [--sketch_size SKETCH_SIZE]
                       [--batch_size BATCH_SIZE] [--threads THREADS] [-h] [--version]
                       in_dir out_dir

Assembly dereplicator

Positional arguments:
  in_dir                     Directory containing all assemblies
  out_dir                    Directory where dereplicated assemblies will be copied

Settings:
  --threshold THRESHOLD      Mash distance clustering threshold (default: 0.005)
  --sketch_size SKETCH_SIZE  Mash assembly sketch size (default: 10000)
  --batch_size BATCH_SIZE    Dereplication iterations will occur on random batches of
                             this many assemblies (default: 500)
  --threads THREADS          Number of CPU threads for Mash (default: 12)

Other:
  -h, --help                 Show this help message and exit
  --version                  Show program's version number and exit
```

Positional arguments:
* `in_dir`: directory containing input assemblies, will be searched recursively. The input assemblies must be in FASTA format and the script will find them in the input directory based on extension. Any of the following are acceptable: `.fasta`, `.fasta.gz`, `.fna`, `.fna.gz`, `.fa` and `.fa.gz`.
* `out_dir`: directory where dereplicated assemblies (i.e. the representative assembly for each cluster) will be copied. If this directory doesn't exist, it will be created. Copied assemblies will have the same filename they had in the input directory.

Settings:
* `--threshold` is the Mash distance below which two assemblies will cluster together and is the most important parameter. Mash distance roughly corresponds to one minus average nucleotide identity, so a threshold of 0.005 (the default) means that assemblies closer than \~99.5% identity will cluster together. Setting it to a small value (e.g. 0.001) will result in more dereplicated assemblies – i.e. only very close relatives will cluster. Setting it to a large value (e.g. 0.02) will result in fewer dereplicated assemblies, perhaps just one per species.
* `--sketch_size` controls the [Mash sketch size](https://mash.readthedocs.io/en/latest/sketches.html#sketch-size). A smaller value (e.g. 1000) will make the process faster but slightly less accurate.
* `--batch_size` controls how many assemblies will be clustered at once (see [Batches](#batches)). You probably don't need to adjust this, but it might be worth playing with if you are trying to minimise RAM usage. If you have plenty of RAM and just want to cluster in a single iteration, set this to a very high number (more than your number of assemblies).
* `--threads` controls how many threads Mash uses for its sketching and distances. Mash scales well in parallel, so use lots of threads if you have them!



## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
