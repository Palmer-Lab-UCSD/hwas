"""

By: Robert Vogel
Affiliation: Palmer Lab at UCSD
"""
from unittest import mock, TestCase, main
import collections

from hwas import _db
from hwas import _constants



class TestCovariate(TestCase):
    def _metadata(self, measure, description, trait_covariate, covariates):
        return dict(measure = measure,
                    description = description,
                    trait_covariate = trait_covariate,
                    covariates = covariates)

    def setUp(self):

        self.cur = mock.MagicMock()
        self.out = mock.MagicMock()
        self.out.rowcount = 1
        self.out.fetchone.return_value = self._metadata("measurement",
                                                "This is the description.",
                                                "trait",
                                                "batch,date,meta,testing")
        self.cur.execute.return_value = self.out

    @mock.patch('hwas._db.is_covariate')
    def test_valid(self, mock_is_covariate):
        mock_is_covariate.return_value = True
        _db.get_covariate_names(self.cur, "test", "random")


    @mock.patch('hwas._db.is_covariate')
    def test_not_covariate(self, mock_is_covariate):
        mock_is_covariate.return_value = False

        with self.assertRaises(ValueError):
            _db.get_covariate_names(self.cur, "test", "random")


    @mock.patch('hwas._db.is_covariate')
    def test_measure_is_covariate(self, mock_is_covariate):
        mock_is_covariate = True
        self.out.fetchone.return_value = self._metadata("measurement",
                                                "This is the description.",
                                                "covariate_trait",
                                                "batch,date,meta,testing")
        with self.assertRaises(ValueError):
            _db.get_covariate_names(self.cur, "test", "random")

    @mock.patch('hwas._db.is_covariate')
    def test_measure_not_unique(self, mock_is_covariate):
        mock_is_covariate = True
        self.out.rowcount = 2

        with self.assertRaises(ValueError):
            _db.get_covariate_names(self.cur, "test", "random")



if __name__ == "__main__":
    main()


