import re
from unittest import mock,TestCase,main

from hwas import _db



class TestMetaData(TestCase):
    def setUp(self):
        self.phenotype = "test_pheno"
        self.schema = "test_schema"
        self.dbname = "test_dbname"
        self.cmd = "this is the command used at the command line"

    def test_output(self):
        output = _db.make_output_metadata(self.dbname,
                                          self.schema,
                                          self.phenotype,
                                          self.cmd)

        # The end of expected `output` string is a new line character
        # and consequently, the last element of the for loop would
        # not match the regex.  To avoid this loop over elements 
        # 0:n-2
        output = output.split('\n')

        for tline in output[:-1]:
            self.assertIsNotNone(re.match("^##\\w+\\s*=\\s*[-:\\w\\s.]+$", tline))


if __name__ == "__main__":
    main()
