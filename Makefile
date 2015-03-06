SHELL=/bin/bash

deposition-test:
	python Deposition.py -i quickTestRunConfig.yml --debug

movie1:
	python batchedMovieGenerator.py -i slowDepConfig.yml -b 1:10
	python MovieGenerator.py -i slowDepConfig.yml -b 1:10 -c
	scp /mddata/uqmstroe/phaseTransitionData/slowRun/md_tot.mp4 uqbcaron@scmb-momar01d.md.smms.uq.edu.au:~/Downloads/movie1.mp4

movie2:
	python batchedMovieGenerator.py -i slowDepConfig.yml -b 101:130
	python MovieGenerator.py -i slowDepConfig.yml -b 101:130 -c
	scp /mddata/uqmstroe/phaseTransitionData/slowRun/md_tot.mp4 uqbcaron@scmb-momar01d.md.smms.uq.edu.au:~/Downloads/movie2.mp4

movie3:
	python MovieGenerator.py -i slowDepConfig.yml -b 1001:1001 -d -kp
	python MovieGenerator.py -i slowDepConfig.yml -b 1001:1001 -c
	scp /mddata/uqmstroe/phaseTransitionData/slowRun/md_tot.mp4 uqbcaron@scmb-momar01d.md.smms.uq.edu.au:~/Downloads/movie3.mp4
