
!constans set
!----------------------------------------------------------------------
module constants
	implicit none
!	integer,parameter :: dp = selected_real_kind(8)
	real(8),parameter :: a=1       !Lattice constant
	real(8),parameter :: Sz=1  !anverage z component of spin
	real(8),parameter :: S=1   !spin of system
	real(8),parameter :: hbar=1    !Plank constant
	real(8),parameter :: N=1       !Number of atoms
	real(8),parameter :: delta=0.01   !self energy
	real(8),parameter :: z=4       !coordination number
	real(8),parameter :: E0=0      !energy of ground state
	real(8),parameter :: J=1       !exchange constant 
	real(8),parameter :: beta = 1  !temperature 
	real(8),parameter :: pi=3.1415926        
	real(8),parameter :: k = 2d0*pi/a
	real(8),parameter :: Sz1 = 0.5,Sz2=-0.5 
!	real(8) :: Sz1,Sz2
	!-----------------------------------------------------------------------
end module
!-------------------------------------------------------------------





!some calculation
!---------------------------------------------------------------------------------
module calc
	use constants
contains
!given k and omega,calculate the spectrum
!--------------------------------------------------------------------------
	!real function spe(k1,omega)
	!	real(dp) :: k1
	!	real(dp) :: omega
	!	real(dp) :: Ek
	!	Ek = J*S*z*sqrt(1-(4*cos(k1*a)/z)**2) !this formula is dependent on crystal structure
	!	spe = 2d0*hbar*Sz*delta/((hbar*omega-Ek)**2+delta**2)
	!end function spe
!---------------------------------------------------------------------------
	function green_fun(kx,ky,omega) result(green) 
		!-------------------------------------------------------------
		! for dgeev
		integer,parameter :: order = 2
		complex(8) :: P(order,order),Eigenvalue(order),VL(order,order),&
			&Eigenvector(order,order), Work(68)
		integer :: Lwork = 68,info
		real(8) :: Rwork(2*order)

		!-------------------------------------------------------------
		real(8) :: Jk,J0,kx,ky
		complex(8) :: omega,green(2,2),c_delta
		complex(8) :: det !detmination of eigenvector
		integer :: m,n,l,k
		complex(8) :: F(2,2)
		F(1,1:2) = (/2*Sz1,0d0/)
		F(2,:) = (/0d0,2*Sz2/)
		Jk = J*(2*cos(kx*a)+2*cos(ky*a))
		J0 = z*J
		P(1,1) = complex(-J0*Sz2,0)
		P(1,2) = complex(Jk*Sz1,0)
		P(2,1) = complex(Jk*Sz2,0)
		P(2,2) = complex(-J0*Sz1,0)
		call zgeev('V','V',order,P,order,Eigenvalue,VL,order,Eigenvector,order,WORK,LWORK,RWORK,INFO)
		if (info .ne. 0) then
			write(*,*) 'dgeev error'
		end if
		p = Eigenvector!P is Eigenvector matrix
		!compute the inverse of Eigenvectors
		!compute green function
		do m = 1,2
			do n = 1,2
				green(m,n)=0
				do l=1,2
					do k=1,2
					c_delta = (0,delta)
					green(m,n) =green(m,n)+(P(m,l)*VL(l,k))*F(k,n)/(omega-Eigenvalue(l)+c_delta)	
					end do
				end do
				green(m,n) = P(m,n) + green(m,n)
			end do
		end do
		



	end function	



end module calc
!------------------------------------------------------------------------------------------------

!main program
!------------------------------------------------------------------------------------
program main
	use constants
	use array
	use calc
	implicit none
!	call calc_spin(-0.5d0,0.5d0,Sz1,Sz2)
!	call iteration(0.5d0,-0.5d0,Sz1,Sz2)
!	write(*,*) Sz1,Sz2
	!------------------------------------------------------------------------------------------------
	real(dp),parameter :: k_inter=0.01,o_inter = 0.01 !k_inter is interval of k, o_inter is interval of omega
	real(dp),parameter :: o_max = 4d0, o_min=0d0  !o_max is max value of omega at omega axial
	integer,parameter ::  k_num = int(k/k_inter)+1    !number of k point
	integer,parameter ::  o_num = int((o_max-o_min)/o_inter)+1  !number of omega point
	real(dp) :: k1(k_num),omega(o_num)
	real(dp) :: spec1(k_num,o_num),spec2(k_num,o_num)
	integer :: knum, onum
	!---------------------------------------------------------------------------------------------------
	!for spec
	complex(8) :: green(2,2),c_omega(o_num)


	!------------------------------------------------------------------------------------------------------
	k1 = vec(-k,k,k_inter)
	omega = vec(o_min,o_max,o_inter)
	open(unit=10,file='data/antiferr1.dat')
	open(unit=11,file='data/antiferr2.dat')
	do knum = 1,size(k1)
		do onum = 1,size(omega)
			c_omega(onum) = complex(omega(onum),0)
			green = -green_fun(k1(knum),0d0,c_omega(onum))
			spec1(knum,onum) = aimag(green(1,1))
			spec2(knum,onum) = aimag(green(2,2))
			write(10,*) k1(knum),omega(onum),spec1(knum,onum)
			write(11,*) k1(knum),omega(onum),spec2(knum,onum)
		end do
	end do

!	use calc
!	complex(8) :: test(2,2) 
!	test = green_fun(1d0,1d0,(1d0,0d0))
!	write(*,*) test

end             
!--------------------------------------------------------------------------------------



subroutine iteration(Sz1_init,Sz2_init,Sz1_new,Sz2_new)
	use array
	use constants
	implicit none

	!-------------------------------------------------------------
	! for dgeev
	integer,parameter :: order = 2
	complex(8) :: P(order,order),Eigenvalue(order),VL(order,order),&
		&VR(order,order), Work(68),inv(order,order)
	integer :: Lwork = 68,info
	real(8) :: Rwork(2*order)

	!-------------------------------------------------------------
	real(8),intent(in) :: Sz1_init,Sz2_init
	real(8),intent(out) :: Sz1_new,Sz2_new
	real(8) :: Sz1_old,Sz2_old

	real(8) :: Jk,J0,kx,ky
	complex(8) :: omega,green(2,2),c_delta,Phi(2)!Phi is a summation which is used to calculate Sz

	!k_num is number of point of k, knum is used to cycle
	integer :: m,l,h,knum 
	integer,parameter :: k_num=500
	real(8) :: k1(k_num) 
	Sz1_old=0
	Sz2_old=0
	Sz1_new=Sz1_init
	Sz2_new=Sz2_init
	k1 = vecn(-k/2d0,k/2d0,k_num)
	
	!compute the Sz1 and Sz2 by iteration
	do while(abs(Sz1_old-Sz1_new) > 0.00001 .and. abs(Sz2_old-Sz2_new) > 0.00001)

		Sz1_old = Sz1_new
		Sz2_old = Sz2_new

		
		Phi = 0
		do l = 1,2
			do knum = 1,k_num
				do m = 1,2
				!write(*,*) k_num
				Jk = J*(2*cos(k1(knum)*a))  !There only calculate one dimension case
				J0 = z*J
				!write(*,*) k1(k_num) 
				!write(*,*) cos(k1(knum)*a) 
				!write(*,*) Jk
				P(1,1) = complex(-J0*Sz2_old,0)
				P(1,2) = complex(Jk*Sz1_old,0)
				P(2,1) = complex(Jk*Sz2_old,0)
				P(2,2) = complex(-J0*Sz1_old,0)
				!write(*,*) P

				call zgeev('V','V',order,P,order,Eigenvalue,VL,order,VR,order,WORK,LWORK,RWORK,INFO)
				!VR is eigenvector matrix,VL is the inverse of VR
				
				if (info .ne. 0) then
					write(*,*) 'dgeev error'
				end if
				!write(*,*) VR
				!write(*,*) VR(1,1),VR(2,2)
				!write(*,*) Jk*Sz1/(sqrt((Jk*Sz1_old)**2+(sqrt(J0**2-Jk**2)*Sz1_old-J0*Sz1_old)**2))
				inv = inverse(VR,2)

				Phi(l) = Phi(l)+2*VR(l,m)*inv(m,l)/(exp(beta*Eigenvalue(m))-1)	
			!	write(*,*)VR(l,m)*inv(m,l)/(exp(beta*Eigenvalue(m))-1)
			!	write(*,*) exp(beta*Eigenvalue(m))-1
				!write(*,*) (exp(beta*Eigenvalue(m))-1)
				end do
				
		!	write(*,*) 'eigenvalue',Eigenvalue
		!	write(*,*) sqrt(J0**2-Jk**2)*Sz1_new
			end do
			Phi(l) = Phi(l)/k_num
	!		write(*,*) Phi(l)
		end do


		stop !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



		write(*,*) Phi(1)+Phi(2)
		Sz1_new = 1/(2*(2*Phi(1)+1))
		Sz2_new = 1/(2*(2*Phi(2)+1))

		write(*,*) Sz1_new,Sz2_new
		
	end do
end subroutine	
!debug	
!-----------------------------------------------------------------------------
subroutine calc_spin(Sz1_init,Sz2_init,Sz1_new,Sz2_new)
	use array
	use constants
	implicit none
	
	!-----input and output-------
	real(8),intent(in) :: Sz1_init,Sz2_init
	real(8),intent(out) :: Sz1_new,Sz2_new
	!----------------------------

	!--------local var-----------
	integer,parameter :: order = 2
	complex(8) :: P(2,2),Phi(order)
	real(8) :: J0 = z*J, Jk, Sz1_old, Sz2_old
	integer,parameter :: k_num = 300
	integer :: knum,l,m
	real(8) :: k1(k_num) 

	!-----------------------------------------

	!--------for lapack-----------------------
	integer, parameter :: Lwork = 100
	real(8) :: Rwork(2*order)
	complex(8) ::Eigenvalue(order)
	complex(8) :: VL(order,order),VR(order,order),Work(Lwork),inv(order,order)
	integer :: info
	!--------------------------------------------------------
	Sz1_new = Sz1_init
	Sz2_new = Sz2_init
	Sz1_old = 0
	Sz2_old = 0

	do  while(abs(Sz1_old-Sz1_new)>0.0001 .or. abs(Sz2_old-Sz2_new)>0.0001)

		open(unit=20,file='aaa')
	Sz1_old = Sz1_new
	Sz2_old = Sz2_new

	k1 = vecn(-k/2d0,k/2d0,k_num)
	Phi = 0
	do l=1,order
		do m=1,order
			do knum=1,k_num
			Phi = 0	
				Jk = J*(2*cos(k1(knum)*a)) !This line must be modified	when go to general
				P(1,1) = -J0*Sz2_old
				P(1,2) = Jk*Sz1_old
				P(2,1) = Jk*Sz2_old
				P(2,2) = -J0*Sz1_old
			!	write(*,*) P
			
				call zgeev('V','V',order,P,order,Eigenvalue,VL,order,VR,order,WORK,LWORK,RWORK,INFO)
				if (info /= 0 ) then 
					write(*,*) 'zgeev error',info
					stop
				end if
				write(*,*) Eigenvalue
				write(20,*) k1(knum),real(eigenvalue(1)),real(eigenvalue(2)) 
		!		write(*,*) sqrt(J0**2-Jk**2)*Sz1
!				write(*,*) Eigenvalue
			!	if (abs(Eigenvalue(1)+Eigenvalue(2))>0.1) then
		!			stop
			!	endif
			!	stop
				inv = inverse(VR,2)
		!		write(*,*) matmul(inv,VR)
				Phi(l)=Phi(l) + 2*VR(l,m)*inv(m,l)/(exp(beta*Eigenvalue(m))-1) 
	!			write(*,*) VR(l,m)*inv(m,l)/(exp(beta*Eigenvalue(m))-1)
			!	write(*,*)'Phi',Phi
			end do
		end do
		Phi(l) = Phi(l)/k_num
!	write(*,*) 'phi',Phi
!	write(*,*) Sz1_new,Sz2_new
		
	end do
!	write(*,*) Phi
	Sz1_new = 1/(2*(2*Phi(1)+1))
	Sz2_new = 1/(2*(2*Phi(2)+1))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	stop!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!	write(*,*) Sz1_new,Sz2_new
	end do

end subroutine
!------------------------------------------------------------------------------

subroutine spin(Sz1_init,Sz2_init,Sz1_new,Sz2_new)
	use array
	use constants
	implicit none

	!----input and output-------
	real(8),intent(in) :: Sz1_init,Sz2_init
	real(8),intent(out) :: Sz1_new,Sz2_new

	!----local vars----------
	integer,parameter :: k_num=300
	real(8) :: k1(k_num),Phi(2)
	integer :: knum,l,m

	real(8) :: Sz1_old,Sz2_old,Jk,J0
	J0 = z*j
	Sz1_new = Sz1_init
	Sz2_new = Sz2_init
	

!	do while(abs(Sz1_old-Sz1_new)>0.001 .or. abs(Sz2_old-Sz2_new)>0,001)
	do  while(abs(Sz1_old-Sz1_new)>0.0001 .or. abs(Sz2_old-Sz2_new)>0.0001)

	Sz1_old = Sz1_new
	Sz2_old = Sz2_new

	k1 = vecn(-k/2d0,k/2d0,k_num)
	Phi = 0
	do l = 1,2
		do knum = 1,k_num
			Jk=2*J*cos(k1(knum)*a)
			Phi(l) = Phi(l) + 1/sqrt(1-Jk**2/J0**2)/tanh(beta*sqrt(J0**2-Jk**2)/2)
		end do
	end do
	Sz1_new = 1/(2*Phi(1))
	Sz2_new = 1/(2*Phi(2))
!	write(*,*) Phi

	end do

	write(*,*) Sz1_new,Sz2_new
end subroutine spin
