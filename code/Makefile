FC = gfortran
FCFLAGS = -O3 -fopenmp -g -ffree-line-length-none -fbounds-check
HEALPIX = /Users/ciaranohare/Work/Zaragoza/Healpix_3.40# if using healpix this will need to be linked
FCFLAGS += -I$(HEALPIX)/include -L$(HEALPIX)/lib -lhealpix -lcfitsio # and this is how you link it

OBJFILES = params.o neldermead.o util.o LabFuncs.o WIMPFuncs.o NeutrinoFuncs.o like.o
PROGRAMS = runCYGNUS_All runCYGNUS_Ideal runCYGNUS_thresholds runMollweide runCYGNUS_Exposure runCYGNUS_timeline runSkymaps runNuFloor


all: clean $(OBJFILES) $(PROGRAMS)

%.o: %.f95
	$(FC) $(FCFLAGS) -c $<

neldermead.o: neldermead.f
		$(FC) $(FCFLAGS) -c neldermead.f

$(PROGRAMS): $(OBJFILES)
	$(FC) $(FCFLAGS) $@.f95  -o $@ $^

clean:
	rm -f *.o *.mod
