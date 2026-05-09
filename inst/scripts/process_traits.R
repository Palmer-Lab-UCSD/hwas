# Phenotype data formatting for hwas
#
# Palmer Lab at UCSD
# 

library(yaml)

char2binary <- function(tablevals) {
    u <- sort(unique(tablevals))
    if (length(u) != 2)
        return(NULL)

    tdata <- rep(0, length(tablevals))
    tdata[tablevals == u[1]] <- 1

    return(structure(tdata,
           factors=u,
           class="char2binary"))
}


zscores <- function(tablevals) {
    m <- mean(tablevals)
    s <- sd(tablevals)
    return((tablevals - m) / s)
}


quantile_norm_gauss <- function(tablevals) {
    u <- rank(tablevals, ties.method = "average")
    u <- u / (length(u) + 1)
    return(qnorm(u, mean = 0, sd = 1))
}


normalize <- function(tablevals, tcfg) {
    if (tcfg$type == "character")
        tablevals <- char2binary(tablevals);

    if (tcfg$normalize == "zscore")
        return(zscores(tablevals))

    if (tcfg$normalize == "quantile")
        return(quantile_norm_gauss(tablevals))

    return(NULL)
}

samps_without_missing_data <- function(d, samp_col_id) {
    is_valid <- apply(d, 1, 
                      function(x) { return(!any(is.na(x))) })
    samps <- as.vector(d[, samp_col_id])
    return(samps[is_valid])
}


get_covariate_names <- function(covariates) {
    cov_names <- vector(mode="character", 
                        length=length(covariates))

    for (i in seq(length(covariates)))
        cov_names[i] <- covariates[[i]]$name

    return(cov_names)
}

# @title the effective formatting script
#
#
main <- function(config) {
    trait_table <- read.table(config$traits$filename, 
                              header=TRUE,
                              sep="\t",
                              na.string=config$traits$na_string)

    pheno_table <- trait_table[,c(config$traits$sample_id, 
                                  config$traits$phenotype$name)]

    psamps <- samps_without_missing_data(pheno_table,
                                         config$traits$sample_id)

    if (length(psamps) == 0)
        stop("No phenotyped samples detected.  Please make sure
              that phenotype data columns and cfg.yaml agree.")

    if (anyDuplicated(psamps) != 0)
        stop("Phenotype data have duplicate sample entries, please,
              resolve before proceeding.")

    # cov is short for covariates
    cov_names <- get_covariate_names(config$traits$covariates)
    cov_table <- trait_table[,c(config$traits$sample_id, cov_names)]

    csamps <- samps_without_missing_data(cov_table,
                                         config$traits$sample_id)


    if (length(csamps) == 0)
        stop("No phenotyped samples detected.  Please make sure
             that phenotype data columns and cfg.yaml agree.")

    if (anyDuplicated(csamps) != 0)
        stop("Covariate data have duplicate sample entries, please,
             resolve before proceeding.")

    # validate genotype sample id file
    gsamps <- readLines(config$genotypes$samples)
    if (length(gsamps) == 0)
        stop("No genotyped samples detected.  Please make sure
             that genotype sample file is not empty.")

    if (anyDuplicated(gsamps) != 0)
        stop("Genotype sample file includes duplicate sample entries,
             please, resolve before proceeding.")


    # Find samples that are in the intersection of those genotyped
    # and phenotyped
    sample_ids <- intersect(gsamps, psamps)
    sample_ids <- intersect(sample_ids, csamps)

    if (anyDuplicated(sample_ids) != 0)
        stop("How did I get here?") 

    
    writeLines(sample_ids, config$harmonized$samples)

    row.names(cov_table) <- cov_table[, config$traits$sample_id]
    cov_table <- cov_table[sample_ids,]

    for (i in seq(length(config$traits$covariates))) {
        cov_name <- config$traits$covariates[[i]]$name
        cov_table[, cov_name] <- normalize(cov_table[,cov_name], 
                                       config$traits$covariates[[i]])
    }

    write.table(cov_table,
                file=config$harmonized$covariates,
                sep="\t",
                quote=FALSE,
                row.names=FALSE)

    row.names(pheno_table) <- pheno_table[,config$traits$sample_id]
    pheno_table <- pheno_table[sample_ids,]

    pheno_name <- config$traits$phenotype$name
    pheno_table[,pheno_name] <- normalize(pheno_table[,pheno_name], 
                                          config$traits$phenotype)

    write.table(pheno_table,
                file=config$harmonized$phenotype,
                sep="\t",
                quote=FALSE,
                row.names=FALSE)
}


# Parse configuration file and run main when called as a command
# line utility
if (!interactive()) {
    if (!file.exists("config.yaml"))
        stop("Can't find the required config, config.yaml")

    config = yaml::yaml.load_file("config.yaml")
    main(config)
}
