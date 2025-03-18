

import os
from unittest import mock, TestCase, main
from dataclasses import dataclass

from hwas import _db




@dataclass
class MockOut:
    rowcount: int


class TestIsSchemaUnique(TestCase):
    def _execute(self, sql_str: str, parameters: tuple[str]) -> mock.MagicMock:

        if "unique_proj" in parameters:
            return MockOut(1)
        elif "not_unique" in parameters:
            return MockOut(3)

        return MockOut(0)

    def setUp(self):
        self.cur = mock.MagicMock()
        self.cur.execute = self._execute

    def test_schema_unique(self):
        self.assertTrue(_db.is_schema_unique(self.cur, "unique_proj"))

    def test_schema_duplicate(self):
        self.assertFalse(_db.is_schema_unique(self.cur, "not_unique"))

    def test_schema_absent(self):
        self.assertFalse(_db.is_schema_unique(self.cur, "absent"))

    def test_schema_number(self):
        self.assertFalse(_db.is_schema_unique(self.cur, 5))

    def test_schema_list(self):
        self.assertFalse(_db.is_schema_unique(self.cur, ["cat", "dog"]))



class TestIsTableUnique(TestCase):
    def _execute(self, sql_str: str, parameters: tuple[str]) -> mock.MagicMock:

        if "unique_proj" in parameters:
            return MockOut(1)
        elif "not_unique" in parameters:
            return MockOut(3)

        return MockOut(0)

    def setUp(self):
        self.cur = mock.MagicMock()
        self.cur.execute = self._execute

    def test_schema_table_unique(self):
        self.assertTrue(_db.is_table_unique(self.cur, "schema", "unique_proj"))

    def test_schema_table_not_unique(self):
        self.assertFalse(_db.is_table_unique(self.cur, "schema", "not_unique"))

    def test_schema_table_absent(self):
        self.assertFalse(_db.is_table_unique(self.cur, "schema", "absent"))



class TestIsCovariate(TestCase):
    def _execute(self, sql_str: str, parameters: tuple[str]) -> mock.MagicMock:
        if _db.COVARIATE_TYPE_TOKEN not in parameters:
            return MockOut(0)

        if "not_unique" in sql_str.as_string():
            return MockOut(3)
        elif "unique_measurement" in sql_str.as_string():
            return MockOut(1)

        return MockOut(0)

    def setUp(self):
        self.cur = mock.MagicMock()
        self.cur.execute = self._execute

    def test_covariate_unique(self):
        self.assertTrue(_db.is_covariate(self.cur, "schema", "unique_measurement"))

    def test_covariate_not_unique(self):
        self.assertFalse(_db.is_covariate(self.cur, "schema", "not_unique"))

    def test_covariate_absent(self):
        self.assertFalse(_db.is_covariate(self.cur, "schema", "absent"))






if __name__ == "__main__":
    main()
