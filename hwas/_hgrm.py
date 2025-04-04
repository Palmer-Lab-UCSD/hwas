"""Compute and write to file the haplotype grm.



"""

import os
import subprocess
import tempfile


from . import _constants
from . import _config


def interface(**kwargs) -> None:
    pars = _config.get_config_section(_constants.FILENAME_CONFIG,
                                      "hgrm")

    pars.update(**kwargs)
    chrm = kwargs["chrm"]
    matrix_filename = f"{chrm}.{_constants.MATRIX_SUFFIX}"

    if not os.path.isdir(hgrm_path := os.path.join(pars.path, pars.hgrm_dir)):
        os.mkdir(hgrm_path, mode = _constants.DEFAULT_DIRMODE)

    if os.path.isfile(os.path.join(hgrm_path, matrix_filename)):
        raise FileExistsError(f"Matrix file {matrix_filename} alread exists.")


    if pars.temp_dir is None:
        temp_dir = tempfile.TemporaryDirectory()
        pars.temp_dir = temp_dir.name

    if not pars.is_specification_complete():
        print(pars)
        raise ValueError("A required value wasn't specified, see above list.")



    # I require that the vcf be indexed for retrieving chromosome data,
    # if the vcf isn't bgzipped, then bgzip,

    if pars.vcf.endswith(".vcf"):
        out_vcf = os.path.join(pars.temp_dir, f"{os.path.basename(pars.vcf)}.gz")
        out = subprocess.run(["bgzip", "-o", out_vcf, pars.vcf],
                            check = True,
                            capture_output = True,
                            text = True)

        if out.returncode != 0:
            raise subprocess.CalledProcessError(returncode = out.returncode,
                                                cmd = out.args,
                                                output = out.stdout,
                                                stderr = out.stderr)

        pars.vcf = out_vcf

    # if no tabix file exists, then make

    if pars.vcf.endswith(".vcf.gz") and not os.path.isfile(f"{pars.vcf}.tbi"):

        out = subprocess.run(["tabix", pars.vcf],
                            check = True,
                            capture_output = True,
                            text = True)

        if out.returncode != 0:
            raise subprocess.CalledProcessError(returncode = out.returncode,
                                                cmd = out.args,
                                                output = out.stdout,
                                                stderr = out.stderr)


    chrm_vcf = os.path.join(pars.temp_dir, f"{chrm}.vcf")

    out = subprocess.run(["bcftools", "view",
                        "--samples-file", pars.samples_filename,
                        "-r", chrm,
                        "-o", chrm_vcf,
                        pars.vcf],
                        capture_output = True,
                        text = True)

    if out.returncode != 0:
        raise subprocess.CalledProcessError(returncode = out.returncode,
                                            cmd = out.args,
                                            output = out.stdout,
                                            stderr = out.stderr)




    with open(os.path.join(hgrm_path, f"{chrm}.mat"), "w") as fid:

        subprocess.run(["hgrm", chrm_vcf],
                    check = True,
                    stdout = fid)
