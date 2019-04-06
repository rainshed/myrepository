aaaa:
	gfortran array_mod.f90 energyband.f90 -llapack
	./a.out
	rm ./a.out
	gnuplot -persist "plot/plot.gp"
	rm data/energyband.dat array.mod constants.mod hamiltonian
test:
	gfortran array_mod.f90 energyband.f90 -llapack
	./a.out
	rm ./a.out
