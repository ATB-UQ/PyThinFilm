
TESTCASES:=quicktest test IPA restarting

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

clean :
	rm -rf $(foreach X,$(TESTCASES),$X_working)

