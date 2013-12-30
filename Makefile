all: alpha documentation.pdf

include ../Makefile_linux.inc


c: clean
	rm -f *.pdf *.dat *.txt gnuplot.scp *.out launch simple mesh_statistics_launch.tex *.tap

KD = $(USERHOME)/ksp-physics-documentation/src/
KAPRI = -I$(KD) $(KD)Cartesian.o $(KD)Engine.o $(KD)Planet.o $(KD)Vector.o

LAUNCH = launch   $(SNOPT_WRAPPER)

LAUNCH_O = $(LAUNCH:%=$(EXAMPLESDIR)/%.o)

ALPHA = alpha   $(SNOPT_WRAPPER)

ALPHA_O = $(ALPHA:%=$(EXAMPLESDIR)/%.o)

SIMPLE = simple   $(SNOPT_WRAPPER)

SIMPLE_O = $(SIMPLE:%=$(EXAMPLESDIR)/%.o)

alpha.o: alpha.cxx alpha.hh setup.hh Makefile
	$(CXX) -c $(CXXFLAGS) $< -o $@

alpha: $(ALPHA_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)

launch.o: launch.cxx launch.hh setup.hh Makefile
	$(CXX) -c $(CXXFLAGS) $< -o $@

launch: $(LAUNCH_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(KAPRI) $(LDFLAGS)

simple.o: simple.cxx simple.hh setup.hh Makefile
	$(CXX) -c $(CXXFLAGS) $< -o $@

simple: $(SIMPLE_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)

documentation.pdf: documentation.tex
	latexmk -pdf documentation.tex
	latexmk -c documentation.tex
