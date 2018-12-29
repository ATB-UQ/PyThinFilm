SHELL=/bin/bash

test:
	$(MAKE) -C testing

clean_tests:
	$(MAKE) -C testing clean

