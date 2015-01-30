SHELL=/bin/bash

%/topo.top: 
	n=$(shell echo $@ | egrep -o '[0-9]+') erb templates/topo.top > $@

%/run.pbs: 
	n=$(shell echo $@ | egrep -o '[0-9]+') erb templates/run.pbs > $@
	chmod +x $@

%/run.mdp:
	cp templates/depo-50ps-run.mdp $@

%/gph80OK.itp:
	cp templates/gph80OK.itp $@

%/cbp-massmod.itp:
	cp templates/cbp-massmod.itp $@

%/cbp300-init.pdb:
	cp templates/cbp300-init.pdb $@

%:
	[[ ! -d $@ ]] && mkdir $@
	make $@/topo.top $@/run.pbs $@/run.mdp $@/cbp-massmod.itp $@/gph80OK.itp $@/cbp300-init.pdb


