
library(parallel)
library(yaml)

all_list <- function(in_list) {
    for (l in in_list)
        if (!l) return(FALSE)

    return(TRUE)
}

get_all_pos <- function(tmp_dir, seq_dir, cfg) {
    geno_fname <- file.path(seq_dir, cfg$dir, cfg$bcf)

    tmpname <- file.path(tmp_dir, sprintf("%s_all", cfg$name))

    if (!file.exists(geno_fname))
        return(NULL)

    if (file.exists(tmpname))
        return(NULL)

    system2("bcftools", c("query", 
                          "-o", tmpname,
                          "-f", "'%CHROM:%POS\n'",
                          geno_fname), 
            stdout = NULL)

    return(readLines(tmpname))
}

get_exclude_set <- function(seq_dir, cfg) {
    fname <- file.path(seq_dir, cfg$dir, cfg$pos_exclude)
    if (!file.exists(fname))
        return(NULL)

    return(readLines(fname))
}


get_chrom_posits <- function(chr, tmp_dir, seq_dir, out_dir) {
    cat(sprintf("%s\n", chr$name))

    if (!dir.exists(tmp_dir)) return(FALSE)
    if (!dir.exists(seq_dir)) return(FALSE)
    if (!dir.exists(out_dir)) return(FALSE)

    exclude_pos <- get_exclude_set(seq_dir, chr)
    if (is.null(exclude_pos)) return(FALSE)

    all_pos <- get_all_pos(tmp_dir, seq_dir, chr)
    if (is.null(all_pos)) return(FALSE)

    include_pos <- setdiff(all_pos, exclude_pos)
    tmp <- strsplit(include_pos, ":")
    
    output <- vector(mode="character", length=length(tmp))
    for (i in seq(length(output)))
        output[i] <- sprintf("%s\t%s", tmp[[i]][1], tmp[[i]][2])
    

    out_fname <- file.path(out_dir,
                         sprintf("%s.tsv", chr$name))

    writeLines(output, out_fname)
                   
    cat(sprintf("Completed %s\n", chr$name))

    return(TRUE)
}


main <- function(cfg) {

    tmp_outdir <- tempdir()

    pos_dir <- file.path(cfg$dest$dir, cfg$dest$pos_dir)
    if (!dir.exists(pos_dir))
        dir.create(pos_dir, mode = "0750")

    out <- mclapply(cfg$genotypes$chrom, 
                    get_chrom_posits,
                    tmp_dir = tmp_outdir,
                    seq_dir = cfg$genotypes$dir,
                    out_dir = pos_dir,
                    mc.cores = 8)

    if (!all_list(out))
        cat("ERROR")
}


if (!interactive()) {
    if (!file.exists("config.yaml"))
        stop("Can't find config")

    main(yaml::yaml.load_file("config.yaml"))
}
