# Tn-seq

Analyze a Tn-seq experiment, in which a bacterial transposon insertion mutant library is subjected to some selective pressure, and abundance of mutants relative to a control (or relative to an expected distribution) is measured as an indicator of mutant fitness.

## Requirements

This package is deployed on [Docker](https://hub.docker.com/r/khturner/tn-seq/), so the easiest and most reproducible way to run the software is with `docker`, which can be easily installed on a [variety of platforms](https://docs.docker.com/engine/installation/)[.

## Use

### `process_reads.py`: Process sequence reads to detect and quantify inserts

Once you have verified your docker installation, the first step is to process a fastq file containing your transposon-adjacent reads to generate a list of sites and their abundance in the library. This is done as follows:

`docker run --rm -v $PWD:/data khturner/tn-seq:latest process_reads.py [OPTIONS]`

The [Docker documentation](https://docs.docker.com/engine/reference/run/) has more information on what that command does. The `-v $PWD:/data` bit is important for our purposes, as it mounts your current directory inside the docker container at `/data`, so assuming you are in the directory where your reads and your reference genome are located, this will expose them to the analysis scripts. See the [Docker documentation on volumes](https://docs.docker.com/engine/tutorials/dockervolumes/) for more information. `process_reads.py` is run like:

```
TODO: input here
```

Here is what a typical run of this step looks like, using data from [Gislason *et al*.](http://www.theuselessweb.com/):

```
TODO: input here
```

### `obligate_essentiality_analysis.py`: Generate a list of essential genes

The goal of this script is to test how different the observed distribution of transposon insertion mutant abundances is from randomly sampled, "fitness-neutral" pseudodata. Genes that contain significantly fewer transposon-adjacent reads than expected by chance are likely essential. See [Turner *et al.*](http://www.pnas.org/content/112/13/4110.abstract) for more details on this method, specifically the [SI Materials and Methods](http://www.pnas.org/lookup/suppl/doi:10.1073/pnas.1419677112/-/DCSupplemental/pnas.201419677SI.pdf?targetid=nameddest=STXT). `obligate_essentiality_analysis.py` is run like:

```
TODO: input here
```

Some important options here include:
* `xxx` - something something

Here is what a typical run of this step looks like, using data from [Gislason *et al*.](http://www.theuselessweb.com/):

```
TODO: input here
```


