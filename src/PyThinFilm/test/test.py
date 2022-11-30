import os
import shutil
import unittest
from pathlib import Path
from psutil import cpu_count

import yaml
from pkg_resources import resource_filename

from PyThinFilm.pytf import main
from PyThinFilm.common import PACKAGE_MAME

N_CORES = cpu_count(logical=False)


def setup_test(test_config):
    test_dir = Path(resource_filename(PACKAGE_MAME, "test"))
    run_config_file = test_dir / test_config
    with open(run_config_file) as fh:
        run_config = yaml.safe_load(fh)
    if os.path.exists(run_config["work_directory"]):
        shutil.rmtree(run_config["work_directory"])
    return run_config


class TestPyThinFilm(unittest.TestCase):

    def test_quick_single_core(self):
        run_config = setup_test("quick_test.yml")
        main(run_config)

    def test_multicore(self):
        run_config = setup_test("multicore_test.yml")
        run_config["n_cores"] = N_CORES
        main(run_config)

    def test_fullerene(self):
        run_config = setup_test("fullerene_test.yml")
        run_config["n_cores"] = N_CORES
        main(run_config)

    def test_solvent(self):
        run_config = setup_test("solvent_test.yml")
        run_config["n_cores"] = N_CORES
        main(run_config)

    def test_annealing(self):
        run_config = setup_test("annealing_test.yml")
        run_config["n_cores"] = N_CORES
        main(run_config)

    def test_equilibration(self):
        run_config = setup_test("equilibration_test.yml")
        run_config["n_cores"] = N_CORES
        main(run_config)

    def test_solvent_accel(self):
        run_config = setup_test("solvent_accel_test.yml")
        run_config["n_cores"] = N_CORES
        main(run_config)

    def test_solvent_self_insert(self):
        run_config = setup_test("solvent_self_insert_test.yml")
        run_config["n_cores"] = N_CORES
        main(run_config)


if __name__ == "__main__":
    unittest.main()
