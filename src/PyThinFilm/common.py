from pathlib import Path

import pkg_resources

PACKAGE_MAME = "PyThinFilm"
TEMPLATE_DIR = Path(pkg_resources.resource_filename(PACKAGE_MAME, "templates"))
RESOURCES_DIR = Path(pkg_resources.resource_filename(PACKAGE_MAME, "resources"))
DEFAULT_SETTING = Path(pkg_resources.resource_filename(PACKAGE_MAME, "default_settings.yml"))

K_B = 0.00831451  # kJ / (mol K)

ROOT_DIRS = ["topology", "trajectory", "log", "tpr", "energy", "restraints", "checkpoint", "final-coordinates",
             "control", "input-coordinates", "stdout", "stderr", "deposition-log"]

VACUUM_DEPOSITION = "vacuum_deposition"
SOLVENT_EVAPORATION = "solvent_evaporation"
THERMAL_ANNEALING = "thermal_annealing"
EQUILIBRATION = "equilibration"
SIMULATION_TYPES = [VACUUM_DEPOSITION, SOLVENT_EVAPORATION, THERMAL_ANNEALING, EQUILIBRATION]