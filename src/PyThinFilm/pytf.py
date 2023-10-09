import logging
import click

from PyThinFilm.deposition import Deposition
from PyThinFilm.common import VACUUM_DEPOSITION, SOLVENT_EVAPORATION, THERMAL_ANNEALING, EQUILIBRATION


def main(config, debug=False):

    deposition = Deposition(config, debug)
    # First check to see if any simulation cycles are required.
    if deposition.run_ID > deposition.last_run_ID:
        msg = f"{deposition.run_ID - 1}/{deposition.last_run_ID} cycles completed, " \
              f"there are no more cycles to run."
        logging.error(msg)
        raise Exception(msg)

    while deposition.run_ID <= deposition.last_run_ID and not deposition.batch_complete:
        # Some housekeeping for the new cycle.
        deposition.init_cycle()

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
        elif deposition.type == EQUILIBRATION:
            deposition.equilibration()
        elif deposition.type == SOLVENT_EVAPORATION:
            # Delete solvent molecules from lower section of the skin if enabled
            deposition.solvent_delete()
            # Insert extra slab of solution below the skin if enabled and
            # conditions are met.
            if not deposition.should_abort:
                deposition.insert_soln_layer()
            # Either of the above may set the should_abort flag, in which case
            # mdrun should be aborted.
            if deposition.should_abort:
                logging.info(f"Aborting mdrun on ID: {deposition.run_ID}")
                break

        # Perform MD simulation
        deposition.run_gromacs_simulation()

        # Cycle complete, increment run ID.
        deposition.run_ID += 1
        logging.debug(f"Run cycle completed, ID incremented: {deposition.run_ID}")
        deposition.update_batch_count()

    # The specified number of MD simulation cycles has been reached
    if deposition.run_ID > deposition.last_run_ID:
        logging.info("Finished {0} cycles".format(deposition.last_run_ID))


@click.command(short_help="Run a PyThinFilm simulation protocol.")
@click.argument('config', nargs=1, type=click.Path(exists=True))
@click.option('-d', '--debug', is_flag=True, default=False, help="Print debugging information.")
@click.help_option("--help", "-h")
def cli(config, debug):
    main(config, debug)


if __name__ == "__main__":
    cli()
