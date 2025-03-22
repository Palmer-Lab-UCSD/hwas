"""Find sample set intersectionWARE.

By: Robert Vogel
Affiliation: Palmer Lab at UCSD



DESCRIPTION

Find the intersection of samples with genotypes, rfid, valid covariate
measurements, and valid phenotype measurements.  Meta data are printed
to logger.
"""
import os
import importlib.resources
import logging
import tempfile
import subprocess


from . import _settings
from . import _constants

logger = logging.getLogger(__name__)

def run(vcf: str) -> None:

    args = _settings.SetIntersectParameters(vcf)

    with tempfile.NamedTemporaryFile(delete_on_close=False) as ftmp:

        out = subprocess.run([args.bcftools,"query","--list-samples",vcf],
                             check = True,
                             stdout=ftmp,
                             stderr=subprocess.PIPE)
        ftmp.close()

        rscript_path = importlib.resources.files('hwas.R').joinpath('trait_intersect.R')

        if not os.path.isfile(rscript_path):
            raise FileNotFoundError(rscript_path)

        out = subprocess.run([args.rscript, rscript_path,
                              "--covariate", args.covariates_file,
                              "--phenotype", args.phenotype_file,
                              "--vcf", ftmp.name,
                              "--id", _constants.SAMPLE_COLNAME],
                              check = True,
                              stdout=subprocess.PIPE,
                              stderr=subprocess.PIPE)

        for w in out.stdout.decode('utf-8').strip().split('\n'):
            logging.info(w)


