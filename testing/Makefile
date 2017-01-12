
TESTCASES:=quicktest test IPA

.PHONY : all $(TESTCASES)

all : $(TESTCASES)

quicktest :
	python test.py quicktest.yml

test :
	python test.py test.yml

IPA :
	python test.py IPA.yml

clean :
	rm -rf $(foreach X,$(TESTCASES),$X_working)

