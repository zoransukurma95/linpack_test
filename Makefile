all:
	@if [ ! -d vers/ ] ; then mkdir vers/ ; cp src/* vers/ ; fi ;
	@rsync -a src/ vers/ ;
	$(MAKE) -C vers/ 
	@if [ ! -d bin/ ] ; then mkdir bin/ ; fi ;
	cp vers/linpack_test bin/linpack_test

.PHONY: clean

clean: 
	rm -rf vers/
