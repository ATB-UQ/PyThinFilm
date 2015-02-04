import Deposition
import unittest


class Test(unittest.TestCase):

    def testBasicRun(self):
        Deposition.runDeposition("runConfig.yml")


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()