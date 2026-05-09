
library(parallel)
library(yaml)
library(hwas)

all_list <- function(in_list) {
    for (l in in_list)
        if (!l) return(FALSE)

    return(TRUE)
}

init_env_vars <- function() {
    nodes <- Sys.getenv("SLURM_NODELIST")
    ncores <- Sys.getenv("SLURM_CPUS_PER_TASK")
    if (nodes == "" || ncores == "") {
        hpc <- FALSE
        ncores <- as.integer(ncores)
        nodes <- ifelse(Sys.getenv("HOSTNAME") == "",
                        "unknown",
                        Sys.getenv("HOSTNAME"))
        ncores <- ifelse(detectCores() >= 1, detectCores(), 1)
    } else {
        hpc <- TRUE
        ncores <- as.integer(ncores)
        ncores <- ifelse(ncores >= 1, ncores, 1)
    }
        
    tmp_dir <- Sys.getenv("TMPDIR")
    if (tmp_dir == "" || !dir.exists(tmp_dir)) {
        if (!dir.create("tmp", mode="0750")) 
            stop("attempted to make tmp dir failed")

        tmp_dir <- file.path(getwd(), "tmp")
    }

    slurm_info <- list(nodes = nodes,
                       account = Sys.getenv("SLURM_JOB_ACCOUNT"),
                       qos = Sys.getenv("SLURM_JOB_QOS"),
                       partition = Sys.getenv("SLURM_JOB_PARTITION"),
                       user = Sys.getenv("SLURM_JOB_USER"),
                       cpus_per_task = Sys.getenv("SLURM_CPUS_PER_TASK"),
                       mem_per_cpu = Sys.getenv("SLURM_MEM_PER_CPU"),
                       ntasks = Sys.getenv("SLURM_NTASKS"))


    return(structure(list(tmp_dir = tmp_dir,
                          ncores = ncores,
                          slurm_info = slurm_info,
                          sys_info = Sys.info()),
            class="sys"))
}


transfer <- function(bcfname, 
                     sample_filename,
                     pos_filename,
                     output_bcf,
                     ncores) {

    if (!file.exists(bcfname))
        return(FALSE)
    if (!file.exists(sample_filename))
        return(FALSE)
    if (!file.exists(pos_filename))
        return(FALSE)
    if (file.exists(output_bcf))
        return(TRUE)
    if (ncores <= 0 || ncores > detectCores())
        ncores <- 1

    success <- system2("bcftools", c("view",
                                     "--regions-file", pos_filename,
                                     "--samples-file", sample_filename,
                                     "--output", output_bcf, 
                                     "--output-type", "u",
                                     "--threads", ncores,
                                     bcfname),
                        stdout = NULL)
    return(success == 0)
}

calculate_grm <- function(chrom, trait_cfg, geno_cfg, env_vars, bcf_cores) {
    ################################################
    # PREPARE VARIABLES
    trait_dir <- file.path(trait_cfg$dir, 
                           trait_cfg$phenotype$name, 
                           trait_cfg$phenotype$version)
    if (!dir.exists(trait_dir))
        return(FALSE)

    sample_filename <- file.path(trait_dir, 
                                 trait_cfg$harmonized$dir,
                                 trait_cfg$harmonized$sample_file)
    if (!file.exists(sample_filename))
        return(FALSE)

    pos_filename <- file.path(geno_cfg$dest$dir,
                              geno_cfg$dest$pos_dir,
                              sprintf("%s.tsv", chrom$name))
    if (!file.exists(pos_filename))
        return(FALSE)

    bcfname <- file.path(geno_cfg$genotypes$dir,
                         chrom$dir,
                         chrom$bcf)
    if (!file.exists(bcfname))
        return(FALSE)

    ################################################
    # SUBSET BCF AND WRITE TO TMP DIR, COMPUTE NODE

    status <- TRUE
    output_bcf <- file.path(env_vars$tmp_dir,
                            sprintf("%s.bcf", chrom$name))
    if (!file.exists(output_bcf)) {
        status <- transfer(bcfname, 
                           sample_filename,
                           pos_filename,
                           output_bcf,
                           bcf_cores)
    }
    if (!status)
        return(FALSE)

    ################################################
    # COMPUTE GRM
    bid <- hwas::bopen(output_bcf, "r")
    if (is.null(bid))
        return(FALSE)

    grmatrix <- hwas::calc_grm(bid, trait_cfg$grm$genotype_format)

    rownames(grmatrix) <- hwas::sample_names(bid)
    colnames(grmatrix) <- hwas::sample_names(bid)

    grmc <- structure(grmatrix,
                   bcf              = bcfname,
                   pos_filename     = pos_filename, 
                   sample_filename  = sample_filename,
                   chrom            = chrom$name,
                   env_vars         = env_vars,
                   class            = "grm")

    grmdir <- file.path(trait_dir, trait_cfg$grm$dir)
    save(grmc, file=file.path(grmdir, 
                              paste0(chrom$name, ".RData")))

    return(TRUE)
}


main <- function(trait_cfg, geno_cfg) {

    env_vars <- init_env_vars()

    # TODO: This is totall ad hoc and I need to return to
    total_cores <- env_vars$ncores
    bcf_cores <- floor(total_cores / 4)
    if (bcf_cores == 0)
        bcf_cores = 1
    chrm_cores <- 4 + total_cores %% 4

    # On multiprocessing, I can divide process by the number
    # of course specified.  So if ncores = 8 and mc.cores is
    # 8, then bcftools for each of the 8 chroms will run on 
    # 1 core.  If mc.cores = 4, then bcftools will run on
    # 2 cores.
    out <- mclapply(geno_cfg$genotypes$chrom,
                    calculate_grm,
                    trait_cfg = trait_cfg,
                    geno_cfg  = geno_cfg,
                    env_vars  = env_vars,
                    bcf_cores = bcf_cores,
                    mc.cores  = chrm_cores)#env_vars$ncores)
    if (!all_list(out))
        cat("ERROR\n")
}


if (!interactive()) {
    if(!file.exists("config.yaml"))
        stop("Config is not found")

    trait_cfg <- yaml::yaml.load_file("config.yaml")

    geno_cfg_filename <- file.path(trait_cfg$genotypes$dir,
                                   trait_cfg$genotypes$config)

    if (!file.exists(geno_cfg_filename))
        stop("Genotypes can't be found")

    geno_cfg <- yaml::yaml.load_file(geno_cfg_filename)

    main(trait_cfg, geno_cfg)
}
