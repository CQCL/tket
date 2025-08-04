import unittest
import pytket


class TestExample(unittest.TestCase):
    def test_foo(self):
        self.assertEqual("foo".upper(), "FOO")
