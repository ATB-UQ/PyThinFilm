import argparse
import logging

from PyThinFilm.deposition import Deposition


def run_deposition(run_config, name, max_cores, debug=False):
    deposition = Deposition(name,
                            run_config,
                            max_cores,
                            debug,
                            )

    if deposition.run_ID > deposition.last_run_ID:
        logging.info("{}/{} depositions completed, there are no more depositions to run".format(deposition.run_ID, deposition.last_run_ID))
        return

    while deposition.run_ID <= deposition.last_run_ID:

        deposition.setup_logging()

        initial_layer_height = deposition.max_layer_height(
            density_fraction_cutoff = deposition.run_config["density_fraction_cutoff"],
        )

        for i in range(deposition.remove_top_molecule):
            residue_id = deposition.top_molecule(deposition.solvent_name)
            logging.info("removing top molecule {0}".format(residue_id))
            deposition.remove_residue_with_id(residue_id)

        deposition.resize_box()

        if deposition.highest_z() > deposition.run_config["escape_tolerance"]:
            # Iterate over the residues and remove the ones that have left the layer
            layer_height = deposition.max_layer_height(density_fraction_cutoff=deposition.run_config['density_fraction_cutoff'])
            logging.info("Layer height {0}".format(layer_height))
            leaving = [ residue for residue in deposition.model.residues[1:] \
                       if deposition.has_residue_left_layer(
                               residue.id,
                               minimum_layer_height=initial_layer_height,
                               layer_height=layer_height,
                               density_fraction_cutoff = deposition.run_config["density_fraction_cutoff"],
                       )
                       ]
            deposition.remove_residues(leaving)
        logging.info("[DEPOSITION] Running deposition with parameters: {parameters_dict}".format(rid=deposition.run_ID, parameters_dict=deposition.runParameters()))
        # Get the next molecule and insert it into the deposition model with random position, orientation and velocity
        deposition.insert_residues()

        # create run directory and run setup make file
        deposition.runSetup()

        # deposition.zero_substrate_velocity()

        # Write updated model to run directory
        deposition.write_init_configuration()
        logging.debug("creating restraints file")
        deposition.write_restraints_file() #NOTE! This modifies substrate positions in deposition.model.

        actualMixture = ",".join([" {0}:{1}".format(r["res_name"], r["count"]) for r in deposition.mixture.values()])
        # Do first Run
        logging.info("    Current mixture is: {0}".format(actualMixture))
        deposition.runSystem()

        # Increment run ID
        deposition.run_ID += 1
        logging.debug("Run ID incremented: '{0}'".format(deposition.run_ID))

    if deposition.run_ID > deposition.last_run_ID:
        logging.info("Finished deposition of {0} molecules".format(deposition.last_run_ID))


def parse_command_line():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-n', '--name', help="name of simulation.")
    parser.add_argument('--debug', dest='debug', action='store_true')
    parser.add_argument('--max-cores', dest='max_cores', default = 1, type=int,
            help='{int} Provide the maximum number of cores to use with mpi.')
    args = parser.parse_args()

    run_deposition(args.input,
                   args.name,
                   args.max_cores,
                   debug=args.debug,
                   )


if __name__ == "__main__":
    parse_command_line()
