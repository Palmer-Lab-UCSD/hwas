
import os
import datetime
from collections.abc import Iterable
import psycopg as pg


from . import _constants


def make_output_metadata(dbname: str,
                         schema: str,
                         phenotype: str,
                         cmd: str) -> str:

    date = datetime.datetime.now(datetime.UTC)
    
    return (f"{_constants.DEFAULT_META_PREFIX}date"
              f"={date.year}-{date.month:02}-{date.day:02}"
              f"-{date.hour:02}:{date.minute:02}:{date.second:02}\n"
        f"{_constants.DEFAULT_META_PREFIX}timezone={date.tzinfo}\n"
        f"{_constants.DEFAULT_META_PREFIX}user={os.environ['USER']}\n"
        f"{_constants.DEFAULT_META_PREFIX}database={dbname}\n"
        f"{_constants.DEFAULT_META_PREFIX}schema={schema}\n"
        f"{_constants.DEFAULT_META_PREFIX}phenotype={phenotype}\n"
        f"{_constants.DEFAULT_META_PREFIX}pipeline_version={_constants.VERSION}\n"
        f"{_constants.DEFAULT_META_PREFIX}input_command={cmd}\n")


def write_to_file(filename: str,
                  cur: pg.Cursor,
                  colnames: Iterable[str],
                  meta_data: str | None = None,
                  delimiter: str = _constants.OUTPUT_DELIMITER) -> None:

    if os.path.exists(filename):
        raise FileExistsError(f"The file {filename} already exists")


    with open(filename, 'w') as fid:

        if meta_data is not None:
            fid.write(meta_data)

        fid.write(_constants.DEFAULT_HEADER_PREFIX
                    + _constants.OUTPUT_DELIMITER.join(colnames)
                    + "\n")

        for record in cur:
            s = ""
            for cname in colnames:
                s += f"{record[cname]}{_constants.OUTPUT_DELIMITER}"
                
            s = s.removesuffix(_constants.OUTPUT_DELIMITER)
            fid.write(f"{s}\n")
