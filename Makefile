SHELL=/bin/bash

default :
	python MovieGenerator.py -i slowDepConfig.yml -b 100:100 -d --fast --keep_png
	scp /mddata/uqmstroe/phaseTransitionData/slowRun/100/md.mp4 uqbcaron@scmb-momar01d.md.smms.uq.edu.au:~/Downloads
