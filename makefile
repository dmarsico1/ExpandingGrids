.PHONY: clean

smr_test: pointer_type.o grid.o netcdf_data.o therm.o init.o interp.o advect.o diffuse.o Boundary.o matrix.o poisson.o average.o rk_time_step.o step.o smr_test.o
		gfortran -o smr_test pointer_type.o grid.o netcdf_data.o therm.o init.o interp.o advect.o diffuse.o Boundary.o matrix.o poisson.o average.o \
		rk_time_step.o step.o smr_test.o\
		-I/usr/local/opt/netcdf/nclude\
		-L/usr/local/opt/netcdf/lib -lnetcdff -lnetcdf

pointer_type.o: pointer_type.f95
		gfortran -c pointer_type.f95

grid.o: grid.f95
		gfortran -c grid.f95


netcdf_data.o: netcdf_data.f95
	gfortran -c netcdf_data.f95 -I /usr/local/opt/netcdf/include -L/usr/local/opt/netcdf/lib -lnetcdff -lnetcdf

therm.o: therm.f95
			gfortran -c therm.f95

init.o: init.f95
		gfortran -c init.f95

interp.o: interp.f95
			gfortran -c interp.f95

advect.o: advect.f95
			gfortran -c advect.f95

diffuse.o: diffuse.f95
			gfortran -c diffuse.f95

matrix.o: matrix.f95
			gfortran -c matrix.f95

poisson.o: poisson.f95
			gfortran -c poisson.f95

Boundary.o: Boundary.f95
			gfortran -c Boundary.f95

average.o: average.f95
			gfortran -c average.f95

rk_time_step.o: rk_time_step.f95
				gfortran -c rk_time_step.f95

step.o: step.f95
		gfortran -c step.f95

smr_test.o: smr_test.f95
		gfortran -c smr_test.f95

clean:
		rm -f pointer_type.o grid.o netcdf_data.o matrix.o init.o interp.o advect.o diffuse.o\
				Boundary.o rk_time_step.o step.o smr_test.o poisson.o therm.o average.o\
				pointer_type.mod grid.mod netcdf_data.mod matrix.mod init.mod interp.mod advect.mod diffuse.mod\
				Boundary.mod rk_time_step.mod step.mod smr_test.mod poisson.mod therm.mod average.mod smr_test
