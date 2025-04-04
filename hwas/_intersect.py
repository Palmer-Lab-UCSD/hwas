"""Find sample set intersectio

Find the intersection of samples with genotypes, rfid, valid covariate
measurements, and valid phenotype measurements.  Meta data are printed
to logger.


"""
import os
import importlib.resources
import logging
import tempfile
import subprocess


from . import _config
from . import _constants

logger = logging.getLogger(__name__)

def interface(vcf: str) -> None:

    args = _config.get_config_section(_constants.FILENAME_CONFIG,
                                      "intersect")

    if not os.path.isabs(vcf):
        vcf = os.path.abspath(vcf)
    args.vcf = vcf                                          # type: ignore[attr-defined]

    _config.update_config_section(args)

    with tempfile.NamedTemporaryFile(delete_on_close=False) as ftmp:

        out = subprocess.run(["bcftools","query",
                              "--list-samples", args.vcf],  # type: ignore[attr-defined]
                             check = True,
                             text = True,
                             stdout=ftmp,
                             stderr=subprocess.PIPE)

        if out.returncode != 0:
            raise subprocess.CalledProcessError(returncode = out.returncode,
                                                cmd = out.args,
                                                output = out.stdout,
                                                stderr = out.stderr)

        ftmp.close()

        rscript_path = importlib.resources.files('hwas.R').joinpath('trait_intersect.R')

        if not os.path.isfile(str(rscript_path)):
            raise FileNotFoundError(rscript_path)

        out = subprocess.run(["Rscript", str(rscript_path),
                              "--covariate", args.covariates_file,  # type: ignore[attr-defined]
                              "--phenotype", args.phenotype_file,   # type: ignore[attr-defined]
                              "--vcf", ftmp.name,
                              "--id", _constants.SAMPLE_COLNAME,
                              "--sample_filename", _constants.FILENAME_SAMPLES],
                              text = True,
                              capture_output = True)

        if out.returncode != 0:
            raise subprocess.CalledProcessError(returncode = out.returncode,
                                                cmd = out.args,
                                                output = out.stdout,
                                                stderr = out.stderr)

        for w in out.stdout.strip().split('\n'):
            logging.info(w)


