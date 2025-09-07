.PHONY: all default clean

default:
	cd examples && $(MAKE)
	@echo Please find your binaries in the example/ directory.

all:
	cd examples && $(MAKE) all
	@echo "Please find your binaries in the example/ and"\
		"smoluchowski/ directory."

doxydocs:
	doxygen docs/doxygen/smolsolver.conf

cleandocs:
	rm -rf docs/doxygen/html docs/doxygen/man

clean:
	cd examples && $(MAKE) clean
	cd smoluchowski && $(MAKE) clean
	$(MAKE) cleandocs
