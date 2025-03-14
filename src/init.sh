#!/usr/bin/bash
#
# Write parameters to be used by pipeline to file.


author="Robert Vogel"
author_contact="rmvogel<at>ucsd<dot>edu"

maintainer="Robert Vogel"
maintainer_contact="rmvogel<at>ucsd<dot>edu"



# magic numbers
num_random=15


function err {
    echo "Runtime Error:"
    echo "$1"
    exit 1
}


function is_valid_option {
    [ ${1#--} != $1 ] && \
        err "$1 is an invalide argument value"
}



function print_help {
    echo "Configure HWAS run"
    echo ""
    echo "./init.sh --account <account name> --qos <qos name>"
    echo ""
    echo "This tool configures and initializes scripts and directories"
    echo "necessary to run a haplotype wide association study (HWAS)."
    echo "All entities will be written to the current working directory."
    echo "For questions or reporting bugs submit an issue on GitHub"
    echo ""
    echo "#TODO PUT IN GITHUB URL"
    echo ""
    echo " or reach out to the maintainer:"
    echo ""
    echo "$maintainer $maintainer_contact"
    echo ""
    echo ""
    echo "DEPENDS ON"
    echo ""
    echo "* bcftools"
    echo "* hgrm"
    echo "* R"
    echo "* R package: qtl2"
    echo ""
    echo ""
    echo "REQUIRED ARGUMENTS"
    echo ""
    echo "--account <identifier>"
    echo "  Account SLURM parameter to be used for all jobs."
    echo "--qos"
    echo "  Quality of service SLURM parameter to be used for all jobs"
    echo "--phenotype"
    echo "  file with phenotypes"
    echo "--vcf_dir"
    echo "  path to directory containing vcf_dir file per chromosome"
    echo ""
    echo ""
    echo "OPTIONAL ARGUMENTS"
    echo ""
    echo "--hgrm"
    echo "  path to directory of haplotype specific genetic relationship"
    echo "  (hgrm) matrices.  Each matrix is the hgrm computed by leave-"
    echo "  one-chromosome-out (loco) method.  For example, chr12_loco.hgrm"
    echo "  is the hgrm using all genetic markers except those of chromosome"
    echo "  12. If none is specified, then the loco hgrm matrices are"
    echo "  computed"

}


# ===========================================================================
# parse command line input

hgrm="__null__"
prev="__null__"
seed="$(tr -dc '0-9' < /dev/urandom | head -c $num_random)"
tmp_path="/tmp"
for curr in $@; do

    if [ $curr = "--help" ]; then
        print_help
        exit 0
    fi


    if [ $prev = "--account" ]; then
        is_valid_option $curr
        account=$curr
    elif [ $prev = "--seed" ]; then 
        is_valid_option $curr
        seed=$curr
    elif [ $prev = "--qos" ]; then
        is_valid_option $curr
        qos=$curr
    elif [ $prev = "--vcf_dir" ]; then
        is_valid_option $curr
        vcf_dir=$curr
    elif [ $prev = "--proj" ]; then
        is_valid_option $curr
        proj=$curr
    elif [ $prev = "--phenotype" ]; then
        is_valid_option $curr
        phenotype=$curr
    elif [ $prev = "--hgrm" ]; then
        is_valid_option $curr
        hgrm=$curr
    elif [ $prev = "--tmp_path" ]; then
        is_valid_option $curr
        tmp_path=$curr
    fi

    prev=$curr
done


# ===========================================================================
# Validate inputs

# if [ ! -d "$vcf_dir" ]; then
#     err "Haplotype vcf_dir path not found"
# elif [ $hgrm != "__null__" ] && [ ! -d $hgrm ]; then
#     err "internal error, please report bug."
# elif [ $qos = "" ]; then
#     err "--qos is a required parameter and was not specified"
# elif [ $account = "" ]; then
#    err "--acount is a required parameter and was not specified"
# fi


# ===========================================================================
# make config file and required directories

if [ -d results ]; then
    err "results directory already exists."
elif [ -d logs ]; then
    err "logs directory already exists."
elif [ $hgrm = "__null__" ] && [ -d hgrm ]; then
    # TODO: I don't understand this logic
    err "hgrm already exists."
fi


if [ ! -d "$tmp_path" ]; then
    err "$tmp_path directory does not exist or are permission constrained"
fi


# # Set up temporary directory
# tmp_dir="${tmp_dir%/}" \
#     "/$(tr -dc 'a-zA-Z' < /dev/urandom | head -c $tmp_dir_char_count)"
# 
# # count the number of attempts and give up after 10
# declare -i i=0
# while [ -d $tmp_dir ]; do
# 
#     tmp_dir="${tmp_dir%/}" \
#         "/$(tr -dc 'a-zA-Z' < /dev/urandom | head -c $tmp_dir_char_count)"
# 
#     i=i+1
# 
#     if [ i -eq 10 ]; then
#         err "All previous 10 direcotries already exist."
#     fi
# done



mkdir -m 750 results 
mkdir -m 750 logs
mkdir -m 750 hgrm
# mkdir -m 750 "$tmp_dir"


echo "# Do not edit file by hand" > config_hwas
printf "account\t%s\n" "$account" >> config_hwas
printf "qos\t%s\n" "$qos" >> config_hwas
printf "vcf_dir\t%s\n" "$vcf_dir" >> config_hwas
printf "seed\t%s\n" "$seed" >> config_hwas
printf "project\t%s\n" "$proj" >> config_hwas
printf "phenotype\t%s\n" "$phenotype" >> config_hwas



# ===========================================================================
# The samples, rats, we have phenotypes are a subset of those we have
# genotypes.  We therefore need to subset our genotypes stored in the vcf 
# format to the samples we have phenotype measurements.
#


# bcftools query -l "$vcf" > "${tmp_dir}/genotyped_sample_names.txt"

