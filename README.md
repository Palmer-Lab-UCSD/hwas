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
a sectino

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
in an enivronmental variable.

```
export password_env_var=my_password_to_db

hwas query --dbname my_db_name \
    --host hostname \
    --port port_number \
    --db_user me \
    --db_pw_env password_env_var
```

This command will produce two files, let's assume that the current
working directory is `path_to_data/schema_name/phenotype/001`,

```
./
|- config
|- covariates.csv
|- phenotype.csv
```

The data queried from the database consists of all relevant samples, 
which may or may not be identical to that which we have genotypes.
To find the samples with covariates, phenotype, and genotypes 
available run

```
hwas intersect path_to_vcf/genotypes.vcf
```

This will create a file `samples` and will subset the samples
in the data tables `covariates.csv` and `phenotype.csv`

```
./
|- config
|- covariates.csv
|- phenotype.csv
|- samples
```




## Pipeline use of `hwas` on hpc system with `SLURM`



# Acknowledgment


This code has been reviewed by Claude, the AI assistant from Anthropic. 
The code was designed and implemented by Robert Vogel, code recommendations
that were provided by Claude were adapted and implemented by Robert Vogel.

