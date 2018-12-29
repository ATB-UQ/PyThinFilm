
TESTCASES:=quicktest IPA CLI MPI solution

.PHONY : all $(TESTCASES)

all : $(TESTCASES)

quicktest :
	python test.py quicktest.yml 1

IPA :
	python test.py IPA.yml 1

restarting :
	python test.py restarting.yml 
	python test.py restarting.yml 
	python test.py restarting.yml || echo done

MPI : 
	python test.py MPI.yml 2

CLI :
	python ../deposition.py -i CLI.yml \
	    --max-cores 1 \
	    --debug \
	    --start 0

solution:
	mkdir -p solution_working/1
	cp configurations/solution.gro solution_working/1/end.gro
	python test.py solution.yml 2


clean :
	rm -rf $(foreach X,$(TESTCASES),$X_working)

