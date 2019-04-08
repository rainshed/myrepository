test:
	gfortran module/array_mod.f90 antiferr.f90 -llapack
	./a.out
	rm ./a.out
	rm constants.mod array.mod calc.mod
	
aaa:
	gfortran module/array_mod.f90 antiferr.f90 -llapack
	./a.out
	gnuplot --persist plot/antiferr1.gp
	gnuplot --persist plot/antiferr2.gp
	rm ./a.out
	rm constants.mod array.mod calc.mod
