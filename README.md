# haplotype wide association study (hwas)


## Example

In R
```
library(hwas)

dir.create("analysis_dir")
setwd("analysis_dir")

hwas::init(< phenotype of interest >,
           < geno root directory that contains all genotypes >,
           < path to sample file relative to geno root >,
           < TRUE / FALSE whether the sample file >,
           < TRUE / FALSE whether the positions are inclusive >)
```

The initialization configures the following directory structure
and configuration file.

```
|- analysis_dir/
   |--- config.yaml
   |--- .scripts/
   |    |--- compute_grm.R
   |    |--- compute_herit.R
   |    |--- process_pos.R
   |    |--- process_traits.R
   |    |--- unique_samples.R
   |
   |--- preprocess_data/
   |    |--- *traits.tsv*
   |    |--- *samples*
   |
   |--- postprocess_data/
   |    |--- **covariates.tsv**
   |    |--- **phenotype.tsv**
   |    |--- **samples**
   |    |--- pos/
   |         |--- **chr01.tsv**
   |         |--- **chr02.tsv**
   |         |--- ...
   |         |--- **chr20.tsv**
   |
   |--- results/
        |--- **heritability**
        |--- grms/
        |    |--- **chr01.RData**
        |    |--- **chr02.RData**
        |    |--- ...
        |    |--- **chr20.RData**
        |
        |--- lod/
        |    |--- **chr01.bed**
        |    |--- ...
        |    |--- **chr20.bed**
        |
        |--- blup/
             |--- **chr01.tsv**
             |--- ...
             |--- **chr20.tsv**
```

*NOTE ON DIR STRUCTURE*
- *traits.tsv* a file with the sample id, covariates, and trait
      measurements.  The user is responsible for making this file,
      and it is considered raw data.
- **covariates.tsv**, and any other name between ** are produced
      by pipeline.  


## Installation and Compiling

As the package depends on the systems htslib
you need to tell R where to find the header and
library files.  To do this set the following 
environment variables

```
export HTSLIB_LIBS=-L<PATH_TO_LIB>
export HTSLIB_CFLAGS=-isystem<PATH_TO_HEADER_DIR>
```

Note here the use of `isystem` instead of `-I` to 
specify the path of the header files.  The reason 
for this is that we want the header files in R 
packages to be discovered before packages locally
installed on our system.


## Features outstanding

* initialization functionality
* heritability script
* lod script
* blup script
* test against qtl2


## AI Disclaimer

The AI Claude 4.7 Opus by Anthropic was used to review 
code, architectural recommendations / discussions, and
in very few cases contributed code.  Any code contributed
by Claude will be made known in the code comments or
in the git logs.


## Copyright notice

Portions from the linear mixed model fitting code are:

* Copyright (C) 2020 Karl Browman
* Copyright (C) 1995, 1996 Robert Gentleman and Ross Ihaka,
* Copyright (C) 1998-2014 The R Core Team

