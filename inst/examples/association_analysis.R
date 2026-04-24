
library(yaml)
library(hwas)

config <- yaml::yaml.load_file("config.yaml")

phenotypes <- read.table(config$phenotype$output, 
                         sep='\t',
                         header=TRUE)

covariates <- read.table(config$covariates$output)

grm



