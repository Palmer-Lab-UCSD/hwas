"""Data types and tables specific to the Palmer Lab

By: Robert Vogel
Affiliation: Palmer Lab at UCSD


The Palmer Lab data base is organized, roughly, as follows.
Experimental data produced by a single lab, award, or other
logical unit is organized in a schema.  Each schema contains,
at a minimum, two tables 'descriptions' and 'gwas_phenotypes'.
The 'descriptions' table enumerates all measurements in the 
'gwas_phenotypes' table, and documents which are to be used
as covariates vs. experimental phenotypes.  Moreover, for each
experimental phenotype the set of measured covariates to be
used are documented.

Example

Schema name:
    p50_investigator

Table names:
    descriptions,
    gwas_phenotypes,
    etc.
"""


from dataclasses import dataclass


PHENOTYPE_TABLENAME = "gwas_phenotypes"
METADATA_TABLENAME = "descriptions"
SAMPLE_COLNAME = "rfid"
COVARIATE_DELIMITER = ','

@dataclass
class Metadata:
    measure: str
    description: str
    trait_covariate: str
    covariates: str


@dataclass
class PhenotypeRecord:
    rfid:str
    measurement:str
