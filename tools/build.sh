#/bin/bash
#
#

for opt in $@; do
    if [ $opt = "-h" ] || [ $opt = "--help" ]; then
        echo "Build hwas package"
        echo ""
        echo "Options"
        echo "  --only-compile  compile code but do not make manual, run"
        echo "                  tests, or build vignettes."

        exit 0
    fi
done


only_compile=1
if [ $# -ge 1 ] && [ "$1" = "--only-compile" ]; then
    only_compile=0
fi


echo "=========================================================="
echo "COMPILE ATTRIBUTES"
echo "=========================================================="

R --no-echo -e "library(Rcpp); Rcpp::compileAttributes()"

echo "=========================================================="
echo "BUILD"
echo "=========================================================="

if [ $only_compile -eq 1 ]; then
    R CMD build .
else
    R CMD build --no-build-vignettes --no-manual .
fi

echo "=========================================================="
echo "CHECK"
echo "=========================================================="

hwas_tarball="$(find -E . -iregex './hwas_[0-9]+.[0-9]+.[a-zA-Z0-9+]+.tar.gz' -type f)"
if [ ! -f "$hwas_tarball" ]; then
    echo "Error: No $hwas_tarball built"
    exit 1
fi


if [ $only_compile -eq 1 ]; then
    R CMD check "$hwas_tarball"
    success=$?
else
    R CMD check --no-tests \
        --no-examples \
        --no-manual \
        --no-vignettes "$hwas_tarball"
fi

if [ $only_compile -eq 1 ] && [ $success -eq 0 ]; then
    echo "=========================================================="
    echo "UNIT TESTS"
    echo "=========================================================="
    
    Rscript "inst/utils/print_unittest_results.R"
fi
