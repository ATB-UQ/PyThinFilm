import logging
import click

from PyThinFilm.deposition import Deposition


def main(config, debug=False):

    deposition = Deposition(config, debug)
    # First check to see if any simulation cycles are required.
    if deposition.run_ID > deposition.last_run_ID:
        msg = f"{deposition.run_ID - 1}/{deposition.last_run_ID} cycles completed, " \
              f"there are no more cycles to run."
        logging.error(msg)
        raise Exception(msg)

    while deposition.cycle():
        pass

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
