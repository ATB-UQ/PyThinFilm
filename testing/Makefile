
TESTCASES:=quicktest test IPA restarting CLI

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

CLI :
	python ../deposition.py -i CLI.yml \
	    --remove-mol-bouncing \
	    --remove-mol-leaving \
	    --continuation \
	    --debug \
	    --start 0


clean :
	rm -rf $(foreach X,$(TESTCASES),$X_working)

