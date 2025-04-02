# Determine the set of valid samples
# 
# Typical usage example:
#

library(argparse)



main <- function(covariates_file, phenotype_file, vcf_samples_file, samp_id, me) {

    if (!is.null(me)) {
        source(file.path(dirname(me),"trait_io.R"))
    }

    genotyped_samples <- read.table(vcf_samples_file, header=FALSE, col.names=samp_id)

    print(sprintf("N samples vcf: %d", dim(genotyped_samples)[1]))

    covariates <- read_data(covariates_file)
    phenotype <- read_data(phenotype_file)

    print(sprintf("N samples covariate table: %d", dim(covariates$data)[1]))
    print(sprintf("N samples phenotype table: %d", dim(phenotype$data)[1]))

    out <- merge(covariates$data, phenotype$data, 
                 by=samp_id, sort=TRUE, all=FALSE)
    out <- merge(out, genotyped_samples, by=samp_id, sort=TRUE, all=FALSE)

    # read.table converts missing data as follows:
    #   empty cell  ->  ''
    #   NA          ->  NA 
    #   NULL        ->  'NULL'

    is_not_valid_matrix <- ((out == "") 
                            + is.na(out) 
                            + t(apply(out, 1, is.nan)))

    invalid_sample_counts <- colSums(is_not_valid_matrix)

    s = "Invalid sample counts by column:"
    for (cname in names(invalid_sample_counts))
        s <- sprintf("%s    %s=%d", s, cname, invalid_sample_counts[cname])
    print(s)

    is_valid_row <- rowSums(is_not_valid_matrix) == 0

    print(sprintf("N samples output %d", sum(is_valid_row)))

    sample_ids <- out[is_valid_row,samp_id]

    rownames(covariates$data) <- covariates$data[,samp_id]
    rownames(phenotype$data) <- phenotype$data[,samp_id]

    covariates$data <- covariates$data[sample_ids,]
    phenotype$data <- phenotype$data[sample_ids,]

    write_data(covariates_file, covariates)
    write_data(phenotype_file, phenotype)

    # write sample ids of those in the intersection of genotypes,
    # covariates, and phenotype values
    write.table(covariates$data[,samp_id],
                file.path(dirname(covariates_file), "samples"),
                col.names=FALSE,
                row.names=FALSE,
                quote=FALSE)
}



parser <- argparse::argument_parser(
        argparse::argument_def(
            ref = "--covariate",
            type = "character",
            help = "Covariate file"),
        argparse::argument_def(
            ref = "--phenotype",
            type = "character",
            help = "Phenotype file"),
        argparse::argument_def(
            ref = "--vcf",
            type = "character",
            help = "Variant call file (vcf)"),
        argparse::argument_def(
            ref = "--id",
            type = "character",
            help = "Column name that represents sample identifiers.")
)



if (!interactive())
{
    out <- parser(commandArgs(), return_program_name=TRUE)
    args <- out$args
    main(args$covariate, args$phenotype, args$vcf, args$id, out$program)
}
