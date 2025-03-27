"""hwas pipeline

This module is the entry point for the haplotype wide study association
analysis (hwas) pipeline developed by the Palmer Lab at UCSD.  Components
of the pipeline are represented as subcommands.  See subcommand module
for more module specific information.  

The pipeline uses a combination of custom tools and those commonly used
by the statistical genetics community.  It is not fully generalizable to
any system as it makes some assumptions that are specific to the Palmer
Lab data organization.


Typical usage example:
    # Let us assume the `$` denote a `bash` prompt

    $ hwas init --dbname my_database \
                --bin path_to_binaries \
                --
                data_schema \
                phenotype_measurement_name
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
    
    # options
    parser_init.add_argument("--config",
                             type = str,
                             default = None,
                             help = ("If a configuration file exists, use it for default"
                                     " parameter values."))
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
    # Compute the hgrm: a single chromsome, all chromosomes, or LOCO
    
    
    
    

    return parser.parse_args(args)
    
    
    
def main(input_args=None):

    logging.basicConfig(format = ("%(asctime)s\t%(levelname)s"
                                  "\t%(filename)s"
                                  "\t%(message)s"),
                        level = logging.INFO)

    logger.info(' '.join(sys.argv))

    if input_args is None:
        input_args = sys.argv[1:]

    args = _parse_args(input_args)
    
    if args.subcommand == "init":
        from . import _init
        _init.run(args.schema,
                  args.phenotype,
                  args.partition,
                  args.account,
                  args.qos,
                  args.config, 
                  args.bin)
    
    elif args.subcommand == "query":
        from . import _query
        _query.run(args.dbname,
                   args.host,
                   args.port,
                   args.db_user,
                   args.db_pw_env,
                   cmd = ' '.join(input_args))
    
    elif args.subcommand == "intersect":
        from . import _intersect
        _intersect.run(args.vcf)

    elif args.subcommand == "hgrm":
        raise NotImplementedError


