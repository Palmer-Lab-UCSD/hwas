
.MAGIC_CONST <- "__UNIT_TEST_RESULTS__"
.PRINT_BARRIER <- "===============================================\n"

glob.dir <- function(pattern, dir_name) {
    dnames <- list.dirs(dir_name, recursive = FALSE)
    if (length(dnames) == 0)
        return(character(0))

    idxes <- grep(pattern, dnames)
    return(dnames[idxes])
}

main <- function() {
    dname <- glob.dir("*.Rcheck", ".")
    if (length(dname) != 1)
        cat("No unique *.Rcheck directory found\n")
    
    test_dir <- file.path(dname, "tests")
    result_files <- list.files(test_dir,
                               pattern = "*.Rout")
    regex_str <- sprintf("^%s$", .MAGIC_CONST)
    for (filename in result_files) {
        fname <- file.path(test_dir, filename)
        cat(.PRINT_BARRIER)
        cat(sprintf("TEST FILE: %s\n", fname))
        fid <- file(fname, "r")

        start_results <- FALSE
        for (tline in readLines(fid)) {
            if (start_results)
                cat(sprintf("%s\n", tline))
            if (!start_results && regexpr(regex_str, tline, perl=TRUE) > 0)
                start_results <- TRUE
        }

        close(fid)

    }
}


if (!interactive())
    main()
