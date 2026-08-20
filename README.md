# haplotype wide association study (hwas)

This package contains tools for running an HWAS interactively
(e.g. interactive R console or Jupyter notebook) or semi-
automatically by a programmatic pipeline.  The data are read
from the vcf, vcf.gz, or bcf file format with htslib [1, 2] 
and the statistics by code from Karl Broman's QTL2 [3, 4].


**Contents**

[HWAS tools](#hwas-tooling)
[Pipeline Instructions](#pipeline-instructions)
[Compiling and Installation](#compiling-and-installation)
[Testing C++ Code](#testing-cpp)
[Features Outstanding](#features-outstanding)
[AI Disclosure](#ai-disclosure)
[Copyright](#copyright)
[References](#references)

## HWAS tools

Please checkout the package vignettes for examples using the
`hwas` package tools.


## Compiling and Installation

**HTSLIB** 

This package depends on the systems htslib you need to tell R 
where to find the header and library files.  To do this set 
the following environment variables

```
export HTSLIB_LIBS=-L<PATH_TO_LIB>
export HTSLIB_CFLAGS=-isystem<PATH_TO_HEADER_DIR>
```

Note, we use `-isystem` instead of `-I` to specify the path 
to header files, because we want the header files in R 
packages to be discovered before packages locally
installed on our system.

**R LIBRARY TREE**

When building the package be mindful that vignettes are also built
from source.  When R builds the package vignette it doesn't read a
users `Rprofile` file.  Consequently, the `hwas` package dependency
needs to be installed either in the default R library path or that
the environment variable `R_LIBS` contains the path to the relevant
library directory.

## Testing C++ code

The `hwas` package uses C/C++ libraries and is connected to R by
Rcpp.  The result a combination of Rcpp independent and dependent
C++ code.  We've created GoogleTest based unit tests for the Rcpp
independent C++ code.  These tests can be run using the `Makefile`
in the root of the source package.  Simply,

```
make tests
```

and the code will be built and tested.


## Features outstanding

* initialization functionality
* heritability script
* lod script
* blup script
* test against qtl2


## AI Disclaimer

The AI, Claude 4.7- Opus by Anthropic was used to review 
code, architectural recommendations / discussions, and
in very few cases contributed code.  Any code contributed
by Claude will be made known in the code comments or
in the git logs.


## Copyright Notice

Portions from the linear mixed model fitting code are:

* Copyright (C) 2020 Karl Browman
* Copyright (C) 1995, 1996 Robert Gentleman and Ross Ihaka,
* Copyright (C) 1998-2014 The R Core Team

## References

[1] HTSLIB citation
[2] https://github.com/samtools/htslib
[3] QTL2 citation
[4] https://github.com/rqtl/qtl2



