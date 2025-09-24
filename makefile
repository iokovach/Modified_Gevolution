# programming environment
COMPILER     := ~/bin/mpic++
INCLUDE      := -I/u/iok2/fftw-3.3.10 -I/projects/illinois/eng/physics/sheltonj/iostuff/hdf5-1.14.6/lib -I/projects/illinois/eng/physics/sheltonj/iostuff/gev_dependencies/class/external/heating -I/projects/illinois/eng/physics/sheltonj/iostuff/gev_dependencies/class/include -I/projects/illinois/eng/physics/sheltonj/iostuff/gev_dependencies/class/external/RecfastCLASS -L/projects/illinois/eng/physics/sheltonj/iostuff/gev_dependencies/class -I/projects/illinois/eng/physics/sheltonj/iostuff/gev_dependencies/class/external/HyRec2020 -L/projects/illinois/eng/physics/sheltonj/iostuff/gev_dependencies/lib -I/projects/illinois/eng/physics/sheltonj/iostuff/gev_dependencies/include -I/projects/illinois/eng/physics/sheltonj/iostuff/gev_dependencies/LATfield2 -L/projects/illinois/eng/physics/sheltonj/iostuff/hdf5-1.14.6/lib -I/projects/illinois/eng/physics/sheltonj/iostuff/hdf5-1.14.6/include -I/projects/illinois/eng/physics/sheltonj/iostuff/gev_dependencies/fftw/api 
LIB          := -lfftw3 -lm -lhdf5 -lgsl -lgslcblas -lclass
HPXCXXLIB    := -lhealpix_cxx -lcfitsio

# target and source
EXEC         := gevolution
SOURCE       := main.cpp
HEADERS      := $(wildcard *.hpp)

# mandatory compiler settings (LATfield2)
DLATFIELD2   := -DFFT3D -DHDF5

# optional compiler settings (LATfield2)
#DLATFIELD2   += -DH5_HAVE_PARALLEL
#DLATFIELD2   += -DEXTERNAL_IO # enables I/O server (use with care)
#DLATFIELD2   += -DSINGLE      # switches to single precision, use LIB -lfftw3f

# optional compiler settings (gevolution)
DGEVOLUTION  := -DPHINONLINEAR
DGEVOLUTION  += -DBENCHMARK
DGEVOLUTION  += -DEXACT_OUTPUT_REDSHIFTS
DGEVOLUTION  += -DORIGINALMETRIC
DGEVOLUTION  += -DVELOCITY      # enables velocity field utilities
#DGEVOLUTION  += -DCOLORTERMINAL
#DGEVOLUTION  += -DCHECK_B
#DGEVOLUTION  += -DHAVE_CLASS    # requires LIB -lclass
#DGEVOLUTION  += -DHAVE_HEALPIX  # requires LIB -lchealpix

# further compiler options
OPT          := -O3 -std=c++11

$(EXEC): $(SOURCE) $(HEADERS) makefile
	$(COMPILER) $< -o $@ $(OPT) $(DLATFIELD2) $(DGEVOLUTION) $(INCLUDE) $(LIB)
	
lccat: lccat.cpp
	$(COMPILER) $< -o $@ $(OPT) $(DGEVOLUTION) $(INCLUDE)
	
lcmap: lcmap.cpp
	$(COMPILER) $< -o $@ $(OPT) -fopenmp $(DGEVOLUTION) $(INCLUDE) $(LIB) $(HPXCXXLIB)

clean:
	-rm -f $(EXEC) lccat lcmap

