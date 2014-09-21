program quasinewton
	integer ndim
	real gtol
	parameter(ndim= 2, gtol= 1.0E-4)
	common /stats/ nfunc, ndfunc
	integer iter, ndfunc, nfunc
	real fret, p(ndim)
	external func, dfunc
	write(*,'(/1x,a)') 'True minimun is at ( -2.0,+-0.89442719)'
	nfunc= 0
	ndfunc= 0
	!p(1)= .1
	p(1)=  -1
	p(2)= 4.2
	write(*,'(/1x,a,2(f7.4,a))') 'Starting vector: (', p(1),',',p(2),')'
	call dfpmin(p,ndim,gtol,iter,fret,func,dfunc)
	write(*,'(1x,a,i3)') ' Iterations:',iter
	write(*,'(1x,a,i3)') ' Func. evals:',nfunc
	write(*,'(1x,a,i3)') ' Deriv. evals:',ndfunc
	write(*,'(1x,a,2(f9.6,a))') 'Solution vector: (', p(1),',',p(2),')'
	write(*,'(1x,a,e14.6)') ' Func. value at solution',fret
end 



real function func(x)
	integer ndfunc, nfunc
	common /stats/ nfunc, ndfunc
	real x(*)
	nfunc=nfunc+1
	func= 10.*( x(2)**2 * (3. -x(1)) - x(1)**2 * (3. + x(1)))**2 + (2. + x(1))**2 / (1. + (2. + x(1))**2)
end

subroutine dfunc(x, df)
	integer nmax
	parameter (nmax= 50)
	integer ndfunc,nfunc
	common /stats/ nfunc, ndfunc
	real x(*), df(nmax)
	ndfunc= ndfunc+1
	
	df(1)= 20. * (x(2)**2 * (3. -x(1)) - x(1)**2 * (3. + x(1))) * (-x(2)**2 - 6. * x(1) - 3. * x(1)**2) + &
	&2. * (2. + x(1))  / (1. + (2. + x(1)) ** 2) - 2. * (2. + x(1))**3 / (1. + (2. + x(1))**2)**2
	df(2)= 40. * (x(2)**2 * (3. -x(1)) - x(1)**2 * (3. +x(1))) * x(2) * (3. - x(1))
	return
end	
	
	
		
!----------------------------------------
! suroutine: dfpmin     
!----------------------------------------
subroutine dfpmin(p,n,gtol,iter,fret,func,dfunc)

      integer iter
      integer n 
      integer nmax
      integer itmax
      
      real fret,gtol,p(n),func,EPS,STPMX,TOLX
      
      parameter (nmax=50, itmax=200, stpmx=100., eps=3.e-8,TOLX=4.*EPS)
      EXTERNAL dfunc,func
!    USES dfunc,func,lnsrch
      
      LOGICAL check
      REAL den,fac,fad,fae,fp,stpmax,sum,sumdg,sumxi,temp,dg(NMAX),hdg(NMAX),hessin(NMAX,NMAX),pnew(NMAX)
      
      
      integer i
      integer j
      integer its
      
      real m_sum
      real m_hessian(nmax, nmax)
      real m_stpmax
      real g(nmax)
      real xi(nmax)
      
      real test
      
!----------------------------------------
      m_sum= 0.0
!----------------------------------------
      
	fp= func(p)
	call dfunc(p,g)

!----------------------------------------
! Hessian initialization
!----------------------------------------      
	do i=1, n
    	do j=1, n
        	m_hessian(i, j)= 0.0
		end do
        m_hessian(i, i)= 1.0
        xi(i)= -g(i)
        m_sum= m_sum + p(i)**2
	end do

	m_stpmax= stpmx * max( sqrt(m_sum) , float(n) )	
	
	
	do its= 1, itmax
		iter= its
		
		call lnsrch(n, p, fp, g, xi, pnew,fret, m_stpmax, check, func)
		
		do i= 1, n
			xi(i)= pnew(i) - p(i)
			p(i)= pnew(i)
		end do

		test= 0.0
		
		do i= 1, n
			temp= abs( xi(i) ) / max( abs(p(i)) , 1.0 )
			if( temp > test ) then
				test= temp
			end if	   
		end do
			
		if( test < tolx ) then
			return
		end if			
		
		do i= 1, n
			dg(i) = g(i)
		end do
		
		call dfunc(p, g)
		test= 0.0;
		den= max(fret, 1.0)
		
		do i= 1, n
			temp= abs( g(i) ) * max( abs(p(i)) , 1.0) / den
			if (temp > test) then
					test= temp;
			end if					
		end do
		
		if ( test < gtol ) then
			return
		end if		
	
		do i= 1, n
			dg(i)= g(i) - dg(i)
		end do		
		
		do i= 1, n
			hdg(i)= 0.0
			
			do j= 1, n
				hdg(i) = hdg(i) + m_hessian(i,j) * dg(j)
			end do
			
		end do
		fac= 0.0
        fae= 0.0
        sumdg= 0.0
        sumxi= 0.0			

		do i= 1, n
			fac= fac + dg(i) * xi(i)
          	fae= fae + dg(i) * hdg(i)
          	sumdg= sumdg + dg(i)**2
          	sumxi= sumxi + xi(i)**2 
		end do
		
		if ( fac**2 > (eps * sumdg * sumxi) ) then
			fac= 1.0 / fac
          	fad= 1.0 / fae
          	
          	do i= 1, n
          		dg(i)= fac * xi(i) - fad * hdg(i)
          	end do
		end if          	
        
		do i= 1, n
			do j= 1, n
				m_hessian(i,j) = m_hessian(i,j) + fac * xi(i) * xi(j) - fad * hdg(i) * hdg(j) + fae * dg(i) * dg(j)
			end do
		end do

		do i= 1, n
			xi(i)= 0.0	
			do j= 1, n				    
				xi(i)= xi(i) - m_hessian(i,j) * g(j)
			end do
		end do
				
	end do

 	return
end 


!----------------------------------------
! suroutine: lnsrch
!----------------------------------------      
subroutine lnsrch(n,xold,fold,g,p,x,f,stpmax,check,func)
      INTEGER n
      LOGICAL check
      REAL f,fold,stpmax,g(n),p(n),x(n),xold(n),func,ALF,TOLX
      PARAMETER (ALF=1.e-4,TOLX=1.e-7)
      EXTERNAL func
!   USES func
      INTEGER i
      REAL a,alam,alam2,alamin,b,disc,f2,rhs1,rhs2,slope,sum,temp,test,tmplam
      check=.false.
      sum=0.
      
	do i=1, n
        sum= sum + p(i) * p(i)
	end do
	
    sum= sqrt(sum)
    
	if(sum > stpmax) then
        do i=1,n
          p(i)= p(i) * stpmax / sum
    	end do
	end if
      
    slope=0.
	do i=1, n
        slope = slope + g(i) * p(i)
	end do	
	
	
   	if( slope >= 0. ) then 
   		print*, 'roundoff problem in lnsrch'
   	end if
   			
      test=0.
      do 14 i=1,n
        temp=abs(p(i))/max(abs(xold(i)),1.)
        if(temp.gt.test)test=temp
14    continue
      alamin=TOLX/test
      alam=1.
1     continue
        do 15 i=1,n
          x(i)=xold(i)+alam*p(i)
15      continue
        f=func(x)
        if(alam.lt.alamin)then
          do 16 i=1,n
            x(i)=xold(i)
16        continue
          check=.true.
          return
        else if(f.le.fold+ALF*alam*slope)then
          return
        else
          if(alam.eq.1.)then
            tmplam=-slope/(2.*(f-fold-slope))
          else
            rhs1=f-fold-alam*slope
            rhs2=f2-fold-alam2*slope
            a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
            b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
            if(a.eq.0.)then
              tmplam=-slope/(2.*b)
            else
              disc=b*b-3.*a*slope
              if(disc.lt.0.)then
                tmplam=.5*alam
              else if(b.le.0.)then
                tmplam=(-b+sqrt(disc))/(3.*a)
              else
                tmplam=-slope/(b+sqrt(disc))
              endif
            endif
            if(tmplam.gt..5*alam)tmplam=.5*alam
          endif
        endif
        alam2=alam
        f2=f
        alam=max(tmplam,.1*alam)
      goto 1
      END

