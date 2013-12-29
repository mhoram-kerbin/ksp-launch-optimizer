all: launch simple documentation.pdf

include ../Makefile_linux.inc


c: clean
	rm -f documentation.pdf

KD = $(USERHOME)/ksp-physics-documentation/src/
KAPRI = -I$(KD) $(KD)Cartesian.o $(KD)Engine.o $(KD)Planet.o $(KD)Vector.o

LAUNCH = launch   $(SNOPT_WRAPPER)

LAUNCH_O = $(LAUNCH:%=$(EXAMPLESDIR)/%.o)

SIMPLE = simple   $(SNOPT_WRAPPER)

SIMPLE_O = $(SIMPLE:%=$(EXAMPLESDIR)/%.o)

launch.o: launch.cxx launch.hh setup.hh Makefile
	$(CXX) -c $(CXXFLAGS) $< -o $@


launch: $(LAUNCH_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(KAPRI) $(LDFLAGS)

simple.o: simple.cxx simple.hh setup.hh Makefile
	$(CXX) -c $(CXXFLAGS) $< -o $@

simple: $(SIMPLE_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(KAPRI) $(LDFLAGS)

documentation.pdf: documentation.tex
	latexmk -pdf documentation.tex
	latexmk -c documentation.tex
