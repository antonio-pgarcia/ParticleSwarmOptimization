!********************************************************************
!** @Module: ParticleSwarmOptimization
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
module ParticleSwarmOptimization
	use testcases
    	implicit NONE
   
	integer, parameter :: c_iteration_max= 10000	
    	integer, parameter :: c_neighbors_max= 	40
	integer, parameter :: c_pop_max= 10000
    	integer, parameter :: c_selected_max= 15
	integer, parameter :: c_var_max= 100

	double precision, dimension(c_pop_max) 			:: a_functionvalues
	double precision, dimension(c_pop_max, c_var_max) 	:: a_particles  
	double precision, dimension(c_pop_max, c_var_max) 	:: a_speeds  
	double precision, dimension(c_pop_max, c_var_max) 	:: a_pbest
	double precision, dimension(c_var_max) 			:: a_gbest
	double precision 													:: f_best
	integer 																:: p_bestindex
		
 	integer :: i_pop_size
    	integer :: i_var_size
		
    	double precision, parameter :: C1				= .5D00
    	double precision, parameter :: C2				= .6D00
    	double precision, parameter :: INERTIA_W	= .5D00
    	double precision, parameter :: EPSILON		= 1.D-08 
    	double precision, parameter :: STEADY		= 10
   
    	double precision, dimension(c_var_max) 	:: m_previous
    	double precision 			:: ff_best
    	integer 				:: i_steady_count

contains
       	
    	
    	
        !********************************************************************
	!** @Function: random
	!**    
	!** Wrapper generador de numeros pseudo aleatorios
	!** Modificar 
	!** 
	!********************************************************************
        function random( )
        	double precision :: random
             	call random_number( random )
	end function random
        
	!********************************************************************
	!** @Function: pswo_run
	!**    
	!** 
	!** 
	!** 
	!********************************************************************
	subroutine pswo_run( stub_function, i_pop, i_var )
		interface stub_function
			subroutine stub_function(  vector, i_dimension, fvalue ) 
				double precision, intent(in out), dimension(:)	::	vector 
				integer, intent(in out) 									:: i_dimension
				double precision, intent(out) 							:: fvalue
			end subroutine stub_function
		end interface stub_function
		integer, intent(in) :: i_pop
        	integer, intent(in) :: i_var
		integer :: i	

    		i_pop_size= i_pop
        	i_var_size= i_var
    	
    		call pswo_init( )	
    		
     		do i= 1, c_iteration_max
  			call pswo_evaluate( stub_function )
     			call findbestneighbor()
     			call pswo_updatespeeds()
     			call pswo_updatepositions()
     			if ( pswo_coverged( a_gbest ) ) then
     				exit 
     			end if	
		end do			
		call pswo_results(stub_function, i)

	end subroutine pswo_run
		
        !********************************************************************
	!** @Function: pswo_init
	!**    
	!** Inicializacion de la poblacion de particulas
	!**
	!********************************************************************
    	subroutine pswo_init( )
            integer :: i
            integer :: j
            
        	call random_seed
			do i= 1, i_pop_size
            			do j= 1, i_var_size
					a_particles(i,j)= ( random() - 0.5D00) * 10
					a_pbest(i,j)= 0.0
					a_gbest(j)= 0
                		end do
				a_functionvalues(i)= 1.0D30
			end do
			
			do i= 1, i_pop_size
            			do j= 1, i_var_size
                			a_speeds(i, j)= ( random() - 0.5D00)
                		end do
				a_functionvalues(i)= 1.0D30
			end do
			
			ff_best= 1.0D30 
			f_best= 1.0D30
			p_bestindex= 0
			m_previous= 0
    	 		i_steady_count= 0
			
				
		end subroutine pswo_init

		!********************************************************************
		!** @Function: pswo_evaluate
		!**    
		!** 
		!**
		!********************************************************************
		subroutine pswo_evaluate( stub_function )
			interface stub_function
				subroutine stub_function(  vector, i_dimension, fvalue ) 
					double precision, intent(in out), dimension(:)	::	vector 
					integer, intent(in out) 									:: i_dimension
					double precision, intent(out) 							:: fvalue
				end subroutine stub_function
			end interface stub_function
			
			integer 												:: i
			integer 												:: j
			double precision, dimension(c_var_max)	:: a_particle 
         	double precision, dimension(c_var_max)	:: a_speed
         	double precision  								:: f_value
         	
			do i= 1, i_pop_size	 
				do j= 1, i_var_size  				
						a_particle(j)= a_particles( i, j )
						a_speed(j)= a_speeds( i, j )
				end do
				call stub_function ( a_particle, i_var_size, f_value )
				
				if( f_value < a_functionvalues( i ) ) then 
					a_functionvalues( i )= f_value
					do j= 1, i_var_size
						a_gbest(j)= a_particle(j) 
						a_pbest(i, j)= a_particle(j) 
					end do
				end if 
			end do	
			
		end subroutine pswo_evaluate
        
		!********************************************************************
		!** @Function: findbestneighbor
		!**    
		!** 
		!**
		!******************************************************************** 
		subroutine findbestneighbor()
			integer :: i
			integer :: k
			integer :: j
			integer :: i_index
			
			do i= 1, i_pop_size
				do k=1, c_neighbors_max
					i_index= int( random() * i_pop_size)+1
					if( f_best > a_functionvalues(i_index) ) then
						f_best= a_functionvalues(i_index)
						p_bestindex= i_index
						do j= 1, i_var_size
							a_gbest(j)= a_particles(i_index,j)
						end do	
					end if
				end do
			end do				
         end subroutine findbestneighbor
         
        !********************************************************************
		!** @Function: updatespeeds
		!**    
		!** 
		!**
		!******************************************************************** 
		subroutine pswo_updatespeeds()
			integer :: i
			integer :: j
			do i= 1, i_pop_size
				do j= 1, i_var_size
 					a_speeds( i, j )= INERTIA_W * a_speeds( i, j ) + C1 * random() * ( a_pbest(i, j)  - a_particles( i, j ) )  
 					if( p_bestindex > 0 ) then 
 						a_speeds( i, j )= a_speeds( i, j ) + C2 * &
 						& random() * ( a_pbest(p_bestindex, j) - a_particles( i, j) )
 					end if
 				end do
 			end do
		end subroutine pswo_updatespeeds
    
		              
		!********************************************************************
		!** @Function: updatepositions
		!**    
		!** 
		!**
		!******************************************************************** 
		subroutine pswo_updatepositions()
			integer :: i
			integer :: j
			do i= 1, i_pop_size
				do j= 1, i_var_size
 					a_particles( i, j )= a_particles( i, j )  +  a_speeds( i, j ) * ( 1.D00 + random() )
 				end do
 			end do
 			
		end subroutine pswo_updatepositions
		
		!********************************************************************
		!** @Function: pswo_coverged
		!**    
		!** 
		!**
		!********************************************************************
		function pswo_coverged( m_gbest )
			double precision, dimension(c_var_max) 	:: m_gbest
			logical 												:: pswo_coverged
			integer 												:: i
			integer 												:: i_var_count
			
			i_var_count= 0
			do i= 1, i_var_size
				if( dabs( m_gbest(i)  - m_previous(i) ) < EPSILON ) then
					i_var_count= i_var_count + 1
				end if	
				m_previous(i)= m_gbest(i)
			end do

			if( i_var_count == i_var_size ) then 
				i_steady_count= i_steady_count + 1
			end if	
			
			if ( i_steady_count > STEADY ) then 
				pswo_coverged= .true. 
			else 	
				pswo_coverged= .false. 
			end if	
		end function pswo_coverged
		
		!********************************************************************
		!** @Function: pswo_results
		!**    
		!** 
		!** 
		!** 
		!********************************************************************
		subroutine pswo_results( stub_function, iterations )
			interface stub_function
				subroutine stub_function(  vector, i_dimension, fvalue ) 
					double precision, intent(in out), dimension(:)	::	vector 
					integer, intent(in out) 									:: i_dimension
					double precision, intent(out) 							:: fvalue
				end subroutine stub_function
			end interface stub_function
			integer :: iterations
			double precision :: x
			integer :: i
			call stub_function ( a_gbest, i_var_size, x )			
			
			print*,'----------------------------------------------------------------------'
			print*,'-- RESULTADOS FINALES'
			print*,'--' 
			print*,'-- Iterations= ' , iterations
			print*,'-- Particle = {', (a_gbest(i), i= 1, i_var_size), '}'
			print*,'-- Function = ', x
			print*,'----------------------------------------------------------------------'
		end subroutine pswo_results
		
end module

