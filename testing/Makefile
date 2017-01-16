
TESTCASES:=quicktest test IPA restarting CLI MPI MPI2 exactcomposition

.PHONY : all $(TESTCASES)

all : $(TESTCASES)

quicktest :
	python test.py quicktest.yml

test :
	python test.py test.yml

IPA :
	python test.py IPA.yml

restarting :
	python test.py restarting.yml 
	python test.py restarting.yml 
	python test.py restarting.yml 

MPI : 
	python test.py MPI.yml 

MPI2 : 
	python ../deposition.py -i MPI2.yml \
	    --remove-mol-bouncing \
	    --remove-mol-leaving \
	    --continuation \
	    --mpi True --max-cores 6

CLI :
	python ../deposition.py -i CLI.yml \
	    --remove-mol-bouncing \
	    --remove-mol-leaving \
	    --continuation \
	    --debug \
	    --start 0

exactcomposition :
	python test.py exactcomposition.yml

clean :
	rm -rf $(foreach X,$(TESTCASES),$X_working)

