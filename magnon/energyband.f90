module constants
	implicit none
	real(8),parameter :: J1=1d0
	real(8),parameter :: J2=0.134d0*J1
	real(8),parameter :: pi=acos(-1d0)
	complex(8),parameter :: ima = complex(0,1d0)
end module

module Hamitonian
	use constants
	implicit none
	contains
	
	!-----Hamitonian of nest neighbor-----------------
	function Hn(ka,kb,kc) result(Ham_n)
		implicit none
		real(8) :: ka,kb,kc
		complex(8) :: Ham_n(12,12) 
		integer :: i,j
		Ham_n = 0
!		Ham_n(1,9) = J1*exp(-ima*(ka+kb))
!		Ham_n(1,10) = J1
!		Ham_n(1,11) = J1*exp(ima*kc)
!		Ham_n(1,12) = J1
!		Ham_n(2,9) = J1
!		Ham_n(2,10) = J1
!		Ham_n(2,11) = J1
!		Ham_n(2,12) = J1*exp(ima*ka) 
!		Ham_n(3,7) = J1*exp(-ima*(ka+kb))
!		Ham_n(3,8) = J1
!		Ham_n(3,11) = J1
!		Ham_n(3,12) = J1
!		Ham_n(4,7) = J1
!		Ham_n(4,8) = J1
!		Ham_n(4,11) = J1*exp(ima*kb)
!		Ham_n(4,12) = J1*exp(-ima*kc)
!		Ham_n(5,7) = J1*(ima*kc)
!		Ham_n(5,8) = J1
!		Ham_n(5,9) = J1
!		Ham_n(5,10) = J1*(ima*kb)
!		Ham_n(6,7) = J1
!		Ham_n(6,8) = J1*(ima*ka)
!		Ham_n(6,9) = J1
!		Ham_n(6,10) = J1*(-ima*kc)
!		
!		do i = 7,12
!			do j = 1,6
!			Ham_n(i,j) = conjg(Ham_n(j,i))
!			end do
!		end do
		Ham_n(1,9) = J1*complex(cos(-(ka+kb)),sin(-(ka+kb)))!J1*exp(-ima*(ka+kb))
		Ham_n(1,10) = J1
		Ham_n(1,11) = J1*complex(cos(kc),sin(kc))!exp(ima*kc)
		Ham_n(1,12) = J1
		Ham_n(2,9) = J1
		Ham_n(2,10) = J1
		Ham_n(2,11) = J1
		Ham_n(2,12) = J1*complex(cos(ka),sin(ka))!exp(ima*ka)
		Ham_n(3,11) = J1
		Ham_n(3,12) = J1
		Ham_n(4,11) = J1*complex(cos(kb),sin(kb))!exp(ima*kb)
		Ham_n(4,12) = J1*complex(cos(-kc),sin(-kc))!exp(-ima*kc)

		do i = 1,5
			do j = i+7,12
				Ham_n(j-6,i+6) = Ham_n(i,j)
			end do
		end do

		do i = 1,6
			do j = 7,12
				Ham_n(j,i) = conjg(Ham_n(i,j))
			end do
		end do
		
	end function
	!--------------------------------------------------------------------------


	!-----Hamitonian of next nearest neighbor---------------------------------------------
	function Hnn(ka,kb,kc) result(Ham_nn)
		implicit none
		real(8) :: ka,kb,kc
		complex(8) :: Ham_nn(12,12)
		integer :: i,j
		Ham_nn = 0
		Ham_nn(1,3) = J2*complex(cos(kc),sin(kc))!exp(ima*kc)
		Ham_nn(1,4) = J2*complex(cos(kc),sin(kc))!exp(ima*kc)
		Ham_nn(1,5) = J2*complex(cos(-1d0*(ka+kb)),sin(-1d0*(ka+kb)))!exp(-ima*(ka+kb))
		Ham_nn(1,6) = J2*complex(cos(kc),sin(kc))!exp(ima*kc)
		Ham_nn(2,3) = J2*complex(cos(ka),sin(ka))!exp(ima*ka)
		Ham_nn(2,4) = J2*complex(cos(-kb),sin(-kb))!exp(-ima*kb)
		Ham_nn(2,5) = J2*complex(cos(-kb),sin(-kb))!exp(-ima*kb)
		Ham_nn(2,6) = J2
		Ham_nn(3,5) = J2*complex(cos(-(ka+kb+kc)),sin(-(ka+kb+kc)))!exp(-ima*(ka+kb+kc))
		Ham_nn(3,6) = J2*complex(cos(ka),sin(ka))!exp(-ima*ka)
		Ham_nn(4,5) = J2
		Ham_nn(4,6) = J2

		do i = 1,5
			do j = i+1,6,1
				Ham_nn(j,i) = conjg(Ham_nn(i,j))
				Ham_nn(i+6,j+6) = conjg(Ham_nn(i,j))
				Ham_nn(j+6,i+6) = conjg(Ham_nn(i+6,j+6))
			end do
		end do

!		do i = 1,6
!			do j = 1,6
!				Ham_nn(i+6,j+6) = conjg(Ham_nn(i,j))
!			end do
!		end do
!
!		!----Hnn is hermite---------------------
!		do i = 1,12
!			do j = i+1,12
!				Ham_nn(j,i) = conjg(Ham_nn(i,j))
!			end do
!		end do


	end function
	!-------------------------------------------------------

	!----get the whole Hamitonian-------------------
	function Ham(ka,kb,kc) result(H)
		implicit none
		real(8) :: ka,kb,kc
		complex(8) :: H(12,12),A(12,12),B(12,12)
		integer :: i,j
		write(*,*) ka,kb,kc!!!!!!!!!!!!!!!!!!!!
		H = Hn(ka,kb,kc) + Hnn(ka,kb,kc)
		A = Hn(0d0,0d0,0d0)
		B = Hnn(0d0,0d0,0d0)
		do i = 1,12
			do j = 1,12
				H(i,i) = H(i,i) + A(i,j) - B(i,j)
			end do
		end do

	end function 
	!-------------------------------------------------------


	!----generate high symmetry line-----------------------
	subroutine high_sym_line(k,knum)
		use array
		implicit none
		integer,intent(in) :: knum
		real(8),intent(out) :: k(knum,3)
		
		!---TH--------
		k(1:knum/3,1) = vecn(0d0,0.5d0*pi*2,knum/3)
		k(1:knum/3,2) = vecn(0d0,0.5d0*pi*2,knum/3)
		k(1:knum/3,3) = vecn(0d0,-0.5d0*pi*2,knum/3)
		!-------------

		!---HP-------
		k(knum/3+1:2*knum/3,1) = vecn(0.5d0*pi*2,0.25d0*pi*2,knum/3)
		k(knum/3+1:2*knum/3,2) = vecn(0.5d0*pi*2,0.25d0*pi*2,knum/3)
		k(knum/3+1:2*knum/3,3) = vecn(-0.5d0*pi*2,0.25d0*pi*2,knum/3)
		!------------

		!---PT-------
		k(2*knum/3+1:knum,1) = vecn(0.25d0*pi*2,0d0,knum/3)
		k(2*knum/3+1:knum,2) = vecn(0.25d0*pi*2,0d0,knum/3)
		k(2*knum/3+1:knum,3) = vecn(0.25d0*pi*2,0d0,knum/3)
		!------------
	end subroutine
	!--------------------------------------------------------

	!-----diagonalization----------------------
	function diag(H) result(eigen_E)
		implicit none
		complex(8) :: eigen_E(12)
		complex(8) :: H(12,12)

		!----for zheev-----------
		integer,parameter :: N = 12, Lwork =24 
		complex(8) :: Work(Lwork),VL(N,N),VR(N,N)
		real(8) :: Rwork(3*N-2)
		integer :: info
		!------------------------

		call ZGEEV('N','N',N,H,N,eigen_E,VL,N,VR,N,Work,Lwork,Rwork,info)

		if (info /= 0) then
			print*,'zgeev error'
		end if

	end function
	!------------------------------------------

	subroutine sort(array,length)
		integer,intent(in) :: length
		real(8) :: array(length),b
		integer :: m,j
		logical :: exchange

		exchange = .true.
		
		do m = length,1,-1
			do j = 1,m-1
				if (array(j)>array(j+1)) then
					b = array(j)
					array(j) = array(j+1)
					array(j+1) = b
					exchange = .false.
				end if 
			end do
			if (exchange) then
				exit
			end if
		end do

	end subroutine

end module


program main
	use Hamitonian
	implicit none
	integer,parameter :: knum=60
	real(8) :: E(12)
	complex(8) :: motion(12,12),test(12,12)
	real(8) :: k(knum,3)
	integer :: i,m,j
	call high_sym_line(k,knum)

	open(unit=10,file='data/energyband.dat',status='unknown')
	motion = 0
	do i = 1,6
		motion(i,i) = 1
	end do
	do i = 7,12
		motion(i,i) = -1
	end do


!	do i = 1,knum
!		E = dble(diag(matmul(motion,Ham(k(i,1),k(i,2),k(i,3)))))
!		write(*,*) E
!		call sort(E,12)
!		write(*,*) E
!		stop !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!		write(10,*) i,E(12),E(11),E(10),E(9),E(8),E(7)!,E(6),E(5),E(4),E(3),E(2),E(1)
!	end do 
	
!	write(*,*) diag(motion)
!	E = diag(Ham(0.5d0,0.5d0,-0.5d0))
!	write(*,*) E
	test = Ham(2*pi/3,pi/2,pi/5)
	open(unit=11,file='test.txt')
	do i = 1,12
		do j = 1,12
			write(11,*) test(i,j)
		end do 
	end do
end 
