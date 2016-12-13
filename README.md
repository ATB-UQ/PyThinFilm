# Launching a Simulation

## Installation on Raijin

    # First of all, add your Raijin ssh key located at ~/.ssh/id_rsa.pub to your Gitlab's account ssh keys.
    # If you don't have one yet, create ssh keys
    [[ ! -d ~/ssh/id_rsa.pub ]] && cd ~ ; ssh-keygen; cd -

    # Then
    git clone ssh://git@scmb-gitlab.biosci.uq.edu.au:2023/ATB-Dependencies/pmx.git
    cd pmx && make install-all-raijin
    cd ..
    git clone ssh://git@scmb-gitlab.biosci.uq.edu.au:2023/MD-protocols/vacuum_deposition.git
    cd vacuum_deposition
    module load gromacs/4.0.7
    make test

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
