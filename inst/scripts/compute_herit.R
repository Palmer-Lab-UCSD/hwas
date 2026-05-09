library(yaml)
library(hwas)

CONFIG <- "config.yaml"
DATADIR <- "postprocess_data"
GRMDIR <- "grms"
CHROM <- "chr12"
FORMAT <- "HD"


main <- function(config) {

    load(file.path(GRMDIR, paste0(CHROM, ".RData")))

    pheno <- read.table(config$harmonized$phenotype,
                        header=TRUE,
                        row.names="sample_id")
    covars <- read.table(config$harmonized$covariates,
                         header=TRUE,
                         row.names="sample_id")

    out <- hwas::est_herit(pheno, s / 4, 
                            addcovar=covars, 
                            reml=TRUE, 
                            cores=8)
    print(out)
}


if (!interactive()) {

    if (!file.exists(CONFIG))
        stop("Samples can't be found")

    main(yaml::yaml.load_file(CONFIG))
}
