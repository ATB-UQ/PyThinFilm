import os
import shutil
import unittest
from pathlib import Path

import yaml
from pkg_resources import resource_filename

from PyThinFilm.pytf import run_deposition
from PyThinFilm.common import PACKAGE_MAME


class TestVacuumDeposition(unittest.TestCase):

    def _quicktest_setup(self):
        test_dir = Path(resource_filename(PACKAGE_MAME, "test"))
        run_config_file = test_dir / "quicktest.yml"
        with open(run_config_file) as fh:
            run_config = yaml.safe_load(fh)
        if os.path.exists(run_config["work_directory"]):
            shutil.rmtree(run_config["work_directory"])
        return run_config

    def _test_single_core(self):
        run_config = self._quicktest_setup()
        run_deposition(run_config, "quick_test", 1, debug=True)

    def test_multi_core(self):
        run_config = self._quicktest_setup()
        run_deposition(run_config, "quick_test", 8, debug=True)


if __name__ == "__main__":
    unittest.main()

