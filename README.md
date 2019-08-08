# Assembly dereplicator

## Introduction

This repo contains a stand-alone Python script to solve a common problem I ran into: dereplicating a group of bacterial genome assemblies.

As an example, imagine that you had 500 genome assemblies for a particular bacterial species and wanted to do some analysis on them, perhaps building a phylogeny. You knew that there was redundancy in this set, i.e. some of your genomes were very similar to each other. So instead of doing the analysis on all 500 assemblies, it might make more sense to reduce your set of assemblies to a smaller set where everything is sufficiently unique.

<p align="center"><img src="images/trees.png" alt="Trees" width="70%"></p>



## Requirements

You'll need Python 3.5 or later to run the `dereplicator.py` script. It's stand-alone and has no package dependencies (i.e. it only uses the Python standard library).

There is one external dependency, [Mash](https://github.com/marbl/Mash), which you'll need installed and callable on your command line. If you can run `mash sketch -h` and `mash dist -h` without getting an error, you should be good to go!



## Method

The basic idea of this script is to cluster assemblies and choose a single representative assembly per cluster. It first uses Mash to get all-vs-all distances between your assemblies. It then builds a graph where the assemblies are vertices and any two assemblies below the distance threshold are connected by an edge. Each [connected component](https://en.wikipedia.org/wiki/Component_(graph_theory)) of the graph is a cluster, and the assembly in that cluster with the largest N50 is chosen as the representative (the use of N50 makes this script prefer completed assemblies over draft assemblies). This graph-based approach is equivalent to a [single-linkage clustering](https://en.wikipedia.org/wiki/Single-linkage_clustering) approach.

However, despite Mash being very fast, the all-vs-all comparison can be cumbersome when you have a large number (e.g. tens of thousands) of assemblies. So this script does clustering in smaller randomly-selected batches: 500 at a time by default (configurable with `--batch_size`). This batch-based dereplication is repeated iteratively until an iteration fails to remove any assemblies. When this happens, most replication has been cleared out, and so the script conducts a final dereplication round with everything that remains.



## Usage

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

The input assemblies must be in fasta format and the script will find them in the input directory based on extension. Any of the following are acceptable: `.fasta`, `.fasta.gz`, `.fna`, `.fna.gz`, `.fa` or `.fa.gz`. Dereplicated assemblies (i.e. the representative assembly for each cluster) will be copied to the output directory with their same filename.

Using a smaller Mash sketch size (e.g. 1000) will make the process faster but slightly less accurate. Mash scales well in parallel, so use lots of threads if you have them!



## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
