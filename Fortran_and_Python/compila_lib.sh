gfortran -c -fPIC -O3 src/precision.f90 -o src/precision.o
gfortran -c -fPIC -O3 src/nrtype.f90 -o src/nrtype.o
gfortran -c -fPIC -O3 src/nrutil.f90 -o src/nrutil.o
gfortran -c -fPIC -O3 src/cg_aux_mod.f90 -o src/cg_aux_mod.o
gfortran -c -fPIC -O3 src/conjgrad_mod.f90 -o src/conjgrad_mod.o
gfortran -c -fPIC -O3 src/straight_ray_module.f90 -o src/straight_ray_module.o
gfortran -c -fPIC -O3 src/isotropic_ray_module.f90 -o src/isotropic_ray_module.o
gfortran -c -fPIC -O3 src/refracted_ray_module.f90 -o src/refracted_ray_module.o
gfortran -c -fPIC -O3 src/count_file_lines_mod.f90 -o src/count_file_lines_mod.o
gfortran -c -fPIC -O3 src/ray_tracing_mod.f90 -o src/ray_tracing_mod.o

gfortran -O3 -shared -fPIC \
	 src/precision.o \
	 src/nrtype.o \
	 src/nrutil.o \
	 src/cg_aux_mod.o \
	 src/conjgrad_mod.o \
	 src/straight_ray_module.o \
	 src/isotropic_ray_module.o \
	 src/refracted_ray_module.o \
	 src/count_file_lines_mod.o \
	 src/ray_tracing_mod.o \
	 ray_tracing_wrap_mod.f90 -o libray_tracing.so

rm *.mod
