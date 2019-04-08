aaaa:
	gfortran module/array_mod.f90 energyband.f90 -llapack
	./a.out
	rm ./a.out
	gnuplot -persist "plot/plot.gp"
	rm array.mod constants.mod hamiltonian.mod
test:
	gfortran array_mod.f90 energyband.f90 -llapack
	./a.out
	rm ./a.out
