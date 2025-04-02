"""hwas pipeline

This module is the entry point for the haplotype wide study association
analysis (hwas) pipeline developed by the Palmer Lab at UCSD.  Components
of the pipeline are represented as subcommands.  See subcommand module
for more module specific information.  

The pipeline uses a combination of custom tools and those commonly used
by the statistical genetics community.  It is not fully generalizable to
any system as it makes some assumptions that are specific to the Palmer
Lab data organization.


"""

import sys
import argparse
import os
import logging


from . import _constants


logger = logging.getLogger(__name__)




def _parse_args(args):

    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers(title = "subcommands",
                                      description = ("Each subcommand runs scripts"
                                                     " for that corresponding step"
                                                     " in the HWAS pipeline."),
                                      dest = "subcommand")
    
    
    
    
    # -----------------------------------------------------------------------------
    # Initialize pipeline
    
    parser_init = subparser.add_parser("init",
                                        help = "Initialize current directory for pipeline.")
    
    # required 
    parser_init.add_argument("schema",
                            type = str,
                            help = ("Name of schema / project for which phenotype"
                                     " are stored."))
    parser_init.add_argument("phenotype",
                            type = str,
                            help = ("Name of the phenotype/measurement to compute"
                                    " genetic associations."))
    
    # -----------------------------------------------------------------------------
    # pipeline: generate scripts from templates and specified configuration
    
    # options
    parser_pipeline = subparser.add_parser("pipeline",
                                help = "Generate scripts from templates.")

    parser_init.add_argument("--partition",
                             type = str,
                             default = os.environ.get(_constants.ENV_SLURM_PARTITION),
                             help = "SLURM partition to schedule job.")
    parser_init.add_argument("--account",
                             type = str,
                             default = os.environ.get(_constants.ENV_SLURM_ACCOUNT),
                             help = "SLURM account to bill the job.")
    parser_init.add_argument("--qos",
                             type = str,
                             default = os.environ.get(_constants.ENV_SLURM_QOS),
                             help = "SLURM quality-of-service")
    parser_init.add_argument("--bin",
                             type = str,
                             default = os.environ.get(_constants.ENV_BIN),
                             help = "Path to where program binary files are found.")
    parser_init.add_argument("--vcf",
                             type = str,
                             default = None,
                             help = "Path to the VCF file")


    # -----------------------------------------------------------------------------
    # Query database
    
    parser_query = subparser.add_parser("query",
                        help = ("Query phenotype and covariates from a"
                                " specified database."))

    parser_query.add_argument("--dbname",
            type=str,
            default=os.environ.get(_constants.ENV_DB_NAME),
            help="Name of data base to query from.")
    parser_query.add_argument("--host",
            type=str,
            default=os.environ.get(_constants.ENV_DB_HOST),
            help="Hostname of database")
    parser_query.add_argument("--port",
            type=str,
            default=os.environ.get(_constants.ENV_DB_PORT),
            help="Port in which to connect to db.")
    parser_query.add_argument("--db_user",
            type=str,
            default=os.environ.get(_constants.ENV_DB_USERNAME),
            help="User name for logging into database.")
    parser_query.add_argument("--db_pw_env",
            type = str,
            default = _constants.ENV_DB_PW,
            help = ("Name of environment variable that stores the"
                   " database password."))
    parser_query.add_argument("--phenotype_file",
            type = str,
            default = _constants.FILENAME_PHENOTYPE,
            help = ("Filename to write phenotype data queried from database"
                   f" , default {_constants.FILENAME_PHENOTYPE}."))
    parser_query.add_argument("--covariates_file",
            type = str,
            default = _constants.FILENAME_COVARIATES,
            help = ("Filename to write covariate data queried from database"
                   f" , default {_constants.FILENAME_PHENOTYPE}."))


    # -----------------------------------------------------------------------------
    # Find the set of samples that do not have missing data 
    # at the intersection of the covariates, phenotype and vcf samples

    parser_intersect = subparser.add_parser("intersect",
                        help = ("Find the set of samples that intersect"
                                " those with genotype data with non"
                                "-missing covariates and phenotypes."))
    
    parser_intersect.add_argument("vcf",
            type=str,
            help = "The path and name of the variant call file (vcf).")


    # -----------------------------------------------------------------------------
    # set and get configuration values

    parser_config = subparser.add_parser("config",
                        help = "Get and set configuration parameters")

    parser_config.add_argument("-s", "--section",
            type = str,
            default = None,
            help = ("Specify section for options to print, if None print"
                    " enumerate sections."))

    parser_config.add_argument("-o", "--option",
            type = str,
            default = None,
            help = ("Specify the option to set or print.  If None,"
                    " enumerate all options in section."))
    
    parser_config.add_argument("-v", "--value",
            type = str,
            default = None,
            help = "Set section.option to value")


    # -----------------------------------------------------------------------------
    # Compute the hgrm: a single chromsome, all chromosomes, or LOCO

    parser_hgrm = subparser.add_parser("hgrm",
                        help = ("Compute the GRM from expected haplotype"
                                " counts and write to file."))

    
    parser_hgrm.add_argument("--chrm",
            type = str,
            default = None,
            help = "Chromosome to compute hgrm")
    parser_hgrm.add_argument("--hgrm",
            type = str,
            default = None,
            help = "Path to `hgrm` binary if not on path.")
    parser_hgrm.add_argument("--tempdir",
            type = str,
            default = None,
            help = "Temporary directory for intermediate vcf and matrix file")
    parser_hgrm.add_argument("--vcf",
            type = str,
            default = None,
            help = ("Filename of vcf to be analyzed for which genotypes are"
                    " extracted."))

    return parser.parse_args(args)
    
    
    
def main(input_args=None):


    if input_args is None:
        input_args = sys.argv[1:]

    args = _parse_args(input_args)

    if args.subcommand == "config":
        from . import _config
        print(args)
        _config.interface(**args.__dict__)
        sys.exit(0)

    elif args.subcommand == "init":
        from . import _init
        _init.interface(args.schema, args.phenotype) 
        sys.exit(0)

    elif args.subcommand == "query":
        from . import _query
        args.cmd = ' '.join(input_args)
        _query.interface(**args.__dict__)
        sys.exit(0)
    
    elif args.subcommand == "intersect":
        from . import _intersect
        _intersect.interface(args.vcf)
        sys.exit(0)

    elif args.subcommand == "hgrm":
        from . import _hgrm
        _hgrm.interface(**args.__dict__)
        sys.exit(0)

    # TODO Think carefully about when I want to activate the logger
    logging.basicConfig(format = ("%(asctime)s\t%(levelname)s"
                                  "\t%(filename)s"
                                  "\t%(message)s"),
                        level = logging.INFO)

    logger.info(' '.join(sys.argv))


    if args.subcommand == "pipeline":
        from . import _pipeline
        _pipeline.interface(**args.__dict__)

    sys.exit(0)



