import logging
import click

from PyThinFilm.deposition import Deposition
from PyThinFilm.common import VACUUM_DEPOSITION, THERMAL_ANNEALING


def main(config, n_cores, debug=False):

    deposition = Deposition(config, n_cores, debug)
    # First check to see if any simulation cycles are required.
    if deposition.run_ID > deposition.last_run_ID:
        logging.info(f"{deposition.run_ID - 1}/{deposition.last_run_ID} depositions completed, "
                     f"there are no more depositions to run")
        return

    while deposition.run_ID <= deposition.last_run_ID:
        # Some housekeeping for the new cycle.
        deposition.setup_logging()
        logging.info(f"Running {deposition.type} simulation #{deposition.run_ID}")
        logging.debug(f"Settings: \n{deposition.run_config_summary()}")
        deposition.resize_box()

        # Remove molecules from the gas phase; in some cases molecules can remain in the gas phase
        # during vacuum deposition.
        deposition.remove_gas_phase_residues()

        if deposition.type == VACUUM_DEPOSITION:
            # Get the next molecule and insert it into the deposition model with random position,
            # orientation and velocity.
            deposition.insert_residues()
        elif deposition.type == THERMAL_ANNEALING:
            # Set temperature based on temperature_list
            deposition.set_temperature_thermal_annealing()

        # Create run directory and run setup gromacs simulation files.
        deposition.init_gromacs_simulation()

        # Perform MD simulation
        deposition.run_gromacs_simulation()

        # Cycle complete, increment run ID.
        deposition.run_ID += 1
        logging.debug(f"Run cycle completed, ID incremented: {deposition.run_ID}")

    # The specified number of MD simulation cycles has been reached
    if deposition.run_ID > deposition.last_run_ID:
        logging.info("Finished {0} cycles".format(deposition.last_run_ID))


@click.command(short_help="Run a PyThinFilm simulation protocol.")
@click.argument('config', nargs=1, type=click.Path(exists=True))
@click.option('-n', '--n_cores', default=1, help="Number of cores. Using n_cores > 1 requires MPI.")
@click.option('-d', '--debug', is_flag=True, default=False, help="Print debugging information.")
def cli(config, n_cores, debug):
    main(config, n_cores, debug)


if __name__ == "__main__":
    cli()
