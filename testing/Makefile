
TESTCASES:=quicktest IPA CLI MPI solution

.PHONY : all $(TESTCASES)

all : $(TESTCASES)

quicktest :
	python test.py quicktest.yml

IPA :
	python test.py IPA.yml

restarting :
	python test.py restarting.yml 
	python test.py restarting.yml 
	python test.py restarting.yml || echo done

MPI : 
	python test.py MPI.yml 

CLI :
	python ../deposition.py -i CLI.yml \
	    --remove-mol-bouncing \
	    --remove-mol-leaving \
	    --continuation \
	    --debug \
	    --start 0

solution:
	mkdir -p solution_working/1
	cp configurations/solution.gro solution_working/1/end.gro
	python test.py solution.yml 


clean :
	rm -rf $(foreach X,$(TESTCASES),$X_working)

