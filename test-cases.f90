!********************************************************************
!** @Module: testcases
!**
!** Antonio P. Garcia
!** antonio.garcia2@gmail.com
!** Calculo Numerico Avanzado
!** Facultad de Informatica - Universidad Politecnica de Madrid
!** 2009
!**
!** Funciones para la validacion del algoritmo. 
!** 
!********************************************************************
module testcases
	implicit none
	
	contains 
	
	!********************************************************************
	!** @subroutine: nr_func
	!**    
	!** Adaptado de Numerical Recipes
	!** Stoer/Bulirsh
	!******************************************************************** 
	subroutine nr_func(  vector, i_dimension, fvalue ) 
		double precision, intent(in out), dimension(:)	::	vector 
		integer, intent(in out) 									:: i_dimension
		double precision, intent(out) 							:: fvalue
		integer :: i
		
 		fvalue= 0.0
 		fvalue= 10.*( vector(2)**2 * (3. - vector(1)) - vector(1)**2  * (3. + vector(1)))**2 + & 
 					& (2. + vector(1))**2 / (1. + (2. + vector(1))**2)
	end subroutine nr_func 	
	
	!********************************************************************
	!** @subroutine: nr_dfunc
	!**    
	!** Adaptado de Numerical Recipes
	!** Stoer/Bulirsh
	!******************************************************************** 
	subroutine nr_dfunc( vector, df )
		double precision, intent(in out), dimension(:)	::	vector 
		double precision, intent(out), dimension(:)	::	df
		integer :: i_dimension
		
		i_dimension= size( vector )

		df(1)= 20. * (vector(2)**2 * (3. -vector(1)) - vector(1)**2 * (3. + vector(1))) * & 
						& (- vector(2)**2 - 6. * vector(1) - 3. * vector(1)**2) + & 
						& 2. * (2. + vector(1))  / (1. + (2. + vector(1)) ** 2) - 2. * &
						& (2. + vector(1))**3 / (1. + (2. + vector(1))**2)**2
		df(2)= 40. * (vector(2)**2 * (3. -vector(1)) - vector(1)**2 * (3. +vector(1))) * vector(2) * (3. - vector(1))
	end subroutine nr_dfunc	
	
	!********************************************************************
	!** @subroutine: rosenbrock
	!**    
	!** Funcion de rosenbrock de 2 variables
	!**
	!******************************************************************** 
	subroutine rosenbrock(  vector, i_dimension, fvalue ) 
		double precision, intent(in out), dimension(:)	::	vector 
		integer, intent(in out) 									:: i_dimension
		double precision, intent(out) 							:: fvalue
		integer :: i
		
 		fvalue= 0.0
  		fvalue= (100.D00*(vector(2) - vector(1)**2)**2 + (vector(1)-1.D00)**2)
	end subroutine rosenbrock
	
	!********************************************************************
	!** @subroutine: rosenbrock_n
	!**    
	!** Funcion de rosenbrock de n variables
	!**
	!******************************************************************** 
	subroutine rosenbrock_n( vector, i_dimension, fvalue ) 
		double precision, intent(in out), dimension(:)	::	vector 
		integer, intent(in out) 									:: i_dimension
		double precision, intent(out) 							:: fvalue
		integer :: i
			
 		fvalue= 0.0
 		do i= 1, i_dimension - 1 
 			fvalue= fvalue + (100.D00*(vector(i+1)-vector(i)**2)**2 + (vector(i)-1.D00)**2)
 		end do
	end subroutine rosenbrock_n
	
end module testcases
