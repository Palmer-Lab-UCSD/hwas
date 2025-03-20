library(argparse)



main <- function(covariates_file, phenotype_file, vcf_file, samp_id, me) {

    if (!is.null(me)) {
        source(file.path(dirname(me),"trait_io.R"))
    }

    covariates <- read_data(covariates_file)$data
    phenotype <- read_data(phenotype_file)$data

    out <- merge(covariates, phenotype, by=samp_id)

    is_valid_row <- apply(out, 1,
                          function (x) (all(x != "")
                            && !is.element(NaN, x)
                            && !is.element(NA, x)))

    print(out[is_valid_row,samp_id])
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
