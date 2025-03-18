"""Entry point to HWAS pipeline
"""
import sys
import argparse
import os






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
    parser_init.add_argument("--account",
                             type = str,
                             default = None,
                             help = "SLURM account to bill the job.")
    parser_init.add_argument("--qos",
                             type = str,
                             default = None,
                             help = "SLURM quality-of-service")
    parser_init.add_argument("--db_pw_env_var",
                             type = str,
                             default = None,
                             help = ("Environment variable that stores the"
                                    " database password."))


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
                                        help = ("Query phenotype and covariates"
                                                " from a specified database."))

    parser_query.add_argument("--dbname",
            type=str,
            default=None,
            help="Name of data base to query from.")
    parser_query.add_argument("--host",
            type=str,
            default=None,
            help="Hostname of database")
    parser_query.add_argument("--port",
            type=str,
            default=None,
            help="Port in which to connect to db.")
    parser_query.add_argument("--user",
            type=str,
            default=None,
            help="User name for logging into database.")

    parser_query.add_argument("--schema",
            default=None,
            type=str,
            help="Name of schema to query from.")
    parser_query.add_argument("--phenotype",
            type=str,
            default=None,
            help="Name of phenotype to query.")


    # -----------------------------------------------------------------------------
    # Compute the hgrm: a single chromsome, all chromosomes, or LOCO
    
    
    
    
    
    
    

    return parser.parse_args(args)
    
    
    
def main(input_args=None):

    if input_args is None:
        input_args = sys.argv[1:]

    args = _parse_args(input_args)
    
    if args.subcommand == "init":
        from . import init
        init.run(args.config, args.account, args.qos, args.db_pw_env_var,
                 args.schema, args.phenotype)
    
    elif args.subcommand == "query":
        from . import query
        query.run(args.dbname, args.host, args.port, args.user,
                  args.schema, args.phenotype, cmd = ' '.join(input_args))
    
    elif args.subcommand == "hgrm":
        raise NotImplementedError


