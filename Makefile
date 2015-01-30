SHELL=/bin/bash

%/topo.top: 
	n=$(shell echo $@ | egrep -o '[0-9]+') erb templates/topo.top > $@

%/run.pbs: 
	n=$(shell echo $@ | egrep -o '[0-9]+') erb templates/run.pbs > $@

%/run.mdp:
	cp templates/depo-50ps-run.mdp $@

%:
	mkdir $@
	make $@/topo.top $@/run.pbs $@/run.mdp


