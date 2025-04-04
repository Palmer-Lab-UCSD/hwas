# :construction: Under construction :building_construction:

# Haplotype Wide Association Study (HWAS) pipeline


Here is a collection of scripts that are used for haplotype 
wide association studies (HWAS) in the Palmer Lab.  Scripts
are not meant to be completely generalizable, but instead for
specific use by the Palmer Lab given the Palmer Lab infrastructure.

## Dependencies

* R >= 4.1.1
    - qtl2: package for association testing
* python >= 3.10
    - psycopg 
* hgrm: Lab tool for computing kinship matrices
* bcftools >= 1.21
* tabix (htslib) >= 1.21
* bgzip (htslib) >= 1.21


## Interactive use of `hwas`

`hwas` consists of several steps, each of which are implemented
as a subcommand.  To start, use the subcommand `init`

```
hwas init schema_name phenotype
```

which will create directories:

```
- schema_name/
    |
    |- phenotype/
        |
        |- 001/
            |
            |- config
```

The `hwas` tools use the written `config` file to specify or change
parameters.  The config file is divided into sections, in each section
there are option and value pairs.  To see what sections exists simply

```
hwas config
```

which will list the sections.  To see the option and value pairs in
a section

```
hwas config -s section_name
```

or to see the value of an individual section and option pair

```
hwas config -s section_name -o option_name
```

and lastly, to change the value of a section and option pair

```
hwas config -s section_name -o option_name -v value
```

A subset of requisite parameters may be specified in the analysis step,
for example let's query data from a postgres data base with the 
Palmer Lab organization,

```
hwas query --dbname my_db_name \
    --host hostname \
    --port port_number \
    --db_user me
```

if the database requires a password you'll need to set the password
in the enivronmental variable `PALMER_DB_PW` prior to querying the data.
The query subcommand will produce two files, let's assume that the current
working directory is `path_to_data/schema_name/phenotype/001`,

```
./
|- config
|- covariates.csv
|- phenotype.csv
```

The data queried from the database consists of all relevant samples, 
in which any one sample may or may not be genotyped.
To find the samples with covariates, phenotype, and genotypes 
available run

```
hwas intersect path_to_vcf/genotypes.vcf
```

This command will overwrite the 
`covariates.csv` and `phenotype.csv` files with the data of 
samples specified in the `samples` file.  

```
./
|- config
|- covariates.csv
|- phenotype.csv
|- samples
```

The `samples` file
is critical, as it is repeatedly used in subsequent steps,
so please do not modify or delete.

Next we can compute the genetic relationship matrix `G` of each
sample $i$ and $j$ using the expected haplotype counts $X_{ih}$ 
and defined as

$$
G_{ij} = \sum_{m=1}^{M_\text{loci}} \sum_{h = 1}^{N_\text{founders}} X_{imh} X_{jmh}
$$

with $h$ and $m$ being the indices over $N_\text{founders}$ ancestral 
founder haplotypes and $M_\text{loci}$ genetic loci.  To compute
this hgrm for `chr12` we would use the `hgrm` submodule as follows

```
hwas hgrm chr12
```

This command makes the direcotry `hgrm`, if it doesn't already exist, and
writes $G$ to file `chr12.mat`.  This is a comma delimited text file with
$N_\text{samples}$ rows and $N_\text{samples}$ columns.  The order of samples
is that in the `samples` file.  Now suppose we next compute the `hgrm` of
`chr1`

```
hwas hgrm chr1
```

would result in the directory structure

```
./
|- config
|- covariates.csv
|- phenotype.csv
|- samples
|- hgrm/
    |- chr1.mat
    |- chr12.mat
```


## Pipeline use of `hwas` on hpc system with `SLURM`



# Acknowledgment


This code has been reviewed by Claude, the AI assistant from Anthropic. 
The code was designed and implemented by Robert Vogel, code recommendations
that were provided by Claude were adapted and implemented by Robert Vogel.

