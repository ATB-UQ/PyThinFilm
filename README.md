# Launching a Simulation

## Run config file
The run config files use a markup language called [YAML](http://www.yaml.org).
It is a simple but powerful key-values format.
If you just modify the one provided, you won't need to learn all the details.

The variable in there have been given very explicit names.

## Running multiple simulations
Running multiple simulations is as simple as calling the script multiple times, as long as the
different runConfig.yml files have different `work_directory` (otherwise you are heading towards trouble).

## Things that can go wrong

### Missing python libraries
The following python libraries are necessary:

* pmx (modified version)
* jinja2
