FC=gfortran
#FFLAGS=-g
# edrop
#---------
#
OBJ = film.o bacoli.o bacoli-aux.o driver-gridmesh.o d1mach_i1mach.o
#
FFLAGS = -O3 -g
#
PDEtest: $(OBJ)
	$(FC) $(FFLAGS) -o PDEtest $(OBJ)	
film.o: film.f
	$(FC) $(FFLAGS) -c film.f
bacoli.o: bacoli.f
	$(FC) $(FFLAGS) -c bacoli.f
bacoli-aux.o: bacoli-aux.f
	$(FC) $(FFLAGS) -c bacoli-aux.f
driver-gridmesh.o: driver-gridmesh.f
	$(FC) $(FFLAGS) -c driver-gridmesh.f
d1mach_i1mach.o: d1mach_i1mach.f
	$(FC) $(FFLAGS) -c d1mach_i1mach.f