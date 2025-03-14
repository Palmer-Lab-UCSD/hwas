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


    # required 
    parser_init.add_argument("schema_name",
                             type = str,
                             help = ("Name of schema / project for which phenotype"
                                     " are stored."))
    parser_init.add_argument("phenotype",
                             type = str,
                             help = ("Name of the phenotype/measurement to compute"
                                     " genetic associations."))
    
    
    # -----------------------------------------------------------------------------
    # Compute the hgrm: a single chromsome, all chromosomes, or LOCO
    
    
    
    
    
    # -----------------------------------------------------------------------------
    # Query database
    
    

    return parser.parse_args(args)
    
    
    
def main(args=None):

    if args is None:
        args = sys.argv[1:]

    args = _parse_args(args)
    
    if args.subcommand == "init":
        from . import _init
        _init.run(args.config, args.account, args.qos, args.schema_name, args.phenotype)
    
    elif args.subcommand == "query":
        raise NotImplementedError
        # _query.run()
    
    elif args.subcommand == "hgrm":
        raise NotImplementedError


