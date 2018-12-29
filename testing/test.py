import sys
import os

TESTING_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_PARENT_DIR = os.path.join(TESTING_DIR, "../../")
sys.path.append(PROJECT_PARENT_DIR)

from vacuum_deposition.deposition import runDeposition

TEST_CONFIG = sys.argv[1]
TEST_CORES = sys.argv[2]
GROMACS_BIN_PATH = "/usr/local/gromacs-4.6/bin"
GROMACS_LIB_PATH = "/usr/local/gromacs-4.6/lib64"
GROMACS_C_INCLUDE_PATH = "/usr/local/gromacs-4.6/include"

def setup_env_variables():
    def append_env(variable, value):
        if variable in os.environ:
            os.environ[variable] += ":{0}".format(value)
        else:
            os.environ[variable] = value
    append_env("PATH", GROMACS_BIN_PATH)
    append_env("LD_LIBRARY_PATH", GROMACS_LIB_PATH)
    append_env("C_INCLUDE_PATH", GROMACS_C_INCLUDE_PATH)

if __name__ == "__main__":
    setup_env_variables()
    runDeposition(TEST_CONFIG, TEST_CORES, debug=True, continuation=True,)
