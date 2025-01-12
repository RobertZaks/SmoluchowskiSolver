default:
	cd examples && $(MAKE)
	@echo Please find your binaries in the example/ directory.

all:
	cd examples && $(MAKE) all
	@echo Please find your binaries in the example/ and smoluchowski/ directory.

clean:
	cd examples && $(MAKE) clean
	cd smoluchowski && $(MAKE) clean
