from pathlib import Path

import pkg_resources

PACKAGE_MAME = "PyThinFilm"
TEMPLATE_DIR = Path(pkg_resources.resource_filename(PACKAGE_MAME, "templates"))
RESOURCES_DIR = Path(pkg_resources.resource_filename(PACKAGE_MAME, "resources"))
DEFAULT_SETTING = Path(pkg_resources.resource_filename(PACKAGE_MAME, "default_settings.yml"))

# grompp produces at least 1 dubious warning message related to the GROMOS force field
GPP_TEMPLATE = "{GMX_EXEC} {grompp} -maxwarn 1 -f {MDP_FILE} -c {initial} -r {restraints} -p {top} -o {tpr} "
GROMPP = "grompp"
TPBCONV = "tpbconv"
MDRUN = "mdrun"
MDRUN_TEMPLATE = "{GMX_EXEC} {mdrun} -s {tpr} -x {xtc} -e {edr} -g {log} -cpo {cpo} -c {final}"
MDRUN_TEMPLATE_GPU = "{GMX_EXEC} {mdrun} -ntomp 1 " \
                     "-s {tpr} -x {xtc} -e {edr} -g {log} -cpo {cpo} -c {final}"
K_B = 0.00831451  # kJ / (mol K)

ROOT_DIRS = ["topology", "trajectory", "log", "tpr", "energy", "restraints", "checkpoint", "final-coordinates",
             "control", "input-coordinates", "stdout", "stderr", "deposition-log"]

VACUUM_DEPOSITION = "vacuum_deposition"
SOLVENT_EVAPORATION = "solvent_evaporation"
THERMAL_ANNEALING = "thermal_annealing"
EQUILIBRATION = "equilibration"
SIMULATION_TYPES = [VACUUM_DEPOSITION, SOLVENT_EVAPORATION, THERMAL_ANNEALING, EQUILIBRATION]