""" Find sample set intersection

By: Robert Vogel
Affiliation: Palmer Lab at UCSD


"""
import importlib.resources

import subprocess


from . import _settings
from . import _constants



def run(vcf: str) -> None:

    args = _settings.SetIntersectParameters(vcf)

    rscript_path = importlib.resources.files('hwas.R').joinpath('trait_intersect.R')

    subprocess.run(["Rscript", rscript_path,
                   "--covariate", args.covariates_file,
                   "--phenotype", args.phenotype_file,
                   "--vcf", args.vcf,
                   "--id", _constants.SAMPLE_COLNAME],
                   check = True)


