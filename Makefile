regular: ksp-launch-optimizer.pdf alpha

all: alpha launch ksp-launch-optimizer.pdf

include ../Makefile_linux.inc


clean:
	rm -f *.o *.pdf *.dat *.txt gnuplot.scp *.out mesh_statistics_*.tex *.tap *.png *.dat launch alpha

KD = $(USERHOME)/ksp-physics-documentation/src/
KAPRI = -I$(KD) $(KD)Cartesian.o $(KD)Engine.o $(KD)Planet.o $(KD)Vector.o

LAUNCH = launch   $(SNOPT_WRAPPER)

LAUNCH_O = $(LAUNCH:%=$(EXAMPLESDIR)/%.o)

ALPHA = alpha   $(SNOPT_WRAPPER)

ALPHA_O = $(ALPHA:%=$(EXAMPLESDIR)/%.o)

alpha.o: alpha.cxx alpha.hh setup.hh Makefile
	$(CXX) -c $(CXXFLAGS) $< -o $@

alpha: $(ALPHA_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(LDFLAGS)

launch.o: launch.cxx launch.hh setup.hh Makefile
	$(CXX) -c $(CXXFLAGS) $< -o $@

launch: $(LAUNCH_O) $(PSOPT_LIBS) $(DMATRIX_LIBS) $(SPARSE_LIBS)
	$(CXX) $(CXXFLAGS) $^ -o $@ -L$(LIBDIR) $(ALL_LIBRARIES) $(KAPRI) $(LDFLAGS)

ksp-launch-optimizer.pdf: ksp-launch-optimizer.tex
	latexmk -pdf ksp-launch-optimizer.tex
	latexmk -c ksp-launch-optimizer.tex
