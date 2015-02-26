SHELL=/bin/bash

default :
	python MovieGenerator.py -i slowDepConfig.yml -b 100:100 -d --fast
