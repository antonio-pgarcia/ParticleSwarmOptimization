!********************************************************************
!** @Program: cna2009
!**
!** Antonio P. Garcia
!** antonio.garcia2@gmail.com
!** Calculo Numerico Avanzado
!** Facultad de Informatica - Universidad Politecnica de Madrid
!** 2009
!**
!** 
!** 
!********************************************************************
program cna2009
	use ParticleSwarmOptimization
	use testcases
	implicit none
		
	print*,'---------- Calculo Numerico Avanzado/2009 ----------'
	print*
	
	call pswo_run( nr_func, 500, 2 )
	!call pswo_run( rosenbrock, 100, 2 )
	!call pswo_run( rosenbrock_n, 5000, 4 )

end program



