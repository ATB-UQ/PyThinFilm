SHELL=/bin/bash

%/topo.top: 
	n=$(shell echo $@ | egrep -o '[0-9]+') erb $(template_dir)/topo.top > $@

%/run.pbs: 
	n=$(shell echo $@ | egrep -o '[0-9]+') erb $(template_dir)/run.pbs > $@
	chmod +x $@

%/run.mdp:
	n=$(shell echo $@ | egrep -o '[0-9]+') erb $(template_dir)/depo-50ps-run.mdp > $@
	
%/gph80OK.itp:
	cp $(template_dir)/gph80OK.itp $@

%/cbp-massmod.itp:
	cp $(template_dir)/cbp-massmod.itp $@

%:
	[[ ! -d $@ ]] && mkdir $@
	make $@/topo.top $@/run.pbs $@/run.mdp $@/cbp-massmod.itp $@/gph80OK.itp


