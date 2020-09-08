FC = gfortran
FCFLAGS = -O3 -fopenmp -g -ffree-line-length-none -fbounds-check
OBJFILES = params.o neldermead.o util.o LabFuncs.o WIMPFuncs.o NeutrinoFuncs.o like.o
PROGRAMS = runCYGNUS_All runCYGNUS_Final runCYGNUS_thresholds runCYGNUS_Exposure runCYGNUS_timeline runCYGNUS_Final2 runNuFloor


all: clean $(OBJFILES) $(PROGRAMS)

%.o: %.f95
	$(FC) $(FCFLAGS) -c $<

neldermead.o: neldermead.f
		$(FC) $(FCFLAGS) -c neldermead.f

$(PROGRAMS): $(OBJFILES)
	$(FC) $(FCFLAGS) $@.f95  -o $@ $^

clean:
	rm -f *.o *.mod
