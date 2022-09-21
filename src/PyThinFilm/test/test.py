import os
import shutil
import unittest
from pathlib import Path

import yaml
from pkg_resources import resource_filename

from PyThinFilm.pytf import main
from PyThinFilm.common import PACKAGE_MAME


class TestVacuumDeposition(unittest.TestCase):

    def _test_setup(self, test_config):
        test_dir = Path(resource_filename(PACKAGE_MAME, "test"))
        run_config_file = test_dir / test_config
        with open(run_config_file) as fh:
            run_config = yaml.safe_load(fh)
        if os.path.exists(run_config["work_directory"]):
            shutil.rmtree(run_config["work_directory"])
        return run_config

    def test_quick_single_core(self):
        run_config = self._test_setup("quick_test.yml")
        main(run_config, 1)

    def test_multicore(self):
        run_config = self._test_setup("multicore_test.yml")
        main(run_config, 8)

    def test_fullerene(self):
        run_config = self._test_setup("fullerene_test.yml")
        main(run_config, 8)

    def test_solvent(self):
        run_config = self._test_setup("solvent_test.yml")
        main(run_config, 8)


if __name__ == "__main__":
    unittest.main()

