
!constans set
!----------------------------------------------------------------------
module constants
	implicit none
	integer,parameter :: dp = selected_real_kind(8)
	real(dp),parameter :: a=1       !Lattice constant
	real(dp),parameter :: Sz=1  !anverage z component of spin
	real(dp),parameter :: S=1   !spin of system
	real(dp),parameter :: hbar=1    !Plank constant
	real(dp),parameter :: N=1       !Number of atoms
	real(dp),parameter :: delta=1   !self energy
	real(dp),parameter :: z=4       !coordination number
	real(dp),parameter :: E0=0      !energy of ground state
	real(dp),parameter :: J=1       !exchange constant 
	real(dp),parameter :: pi=3.1415926        
	real(dp),parameter :: k = 2*pi/a
	!-----------------------------------------------------------------------
end module
!-------------------------------------------------------------------






program main
	use constants
		real(8) :: Jk,J0,kx=1,ky=1,Sz1=1,Sz2=-1
		real(8) :: P(2,2)
		integer,parameter :: oder =2
		integer :: LWORK=68
		integer :: INFO
		real(8) :: Eigenvalue_real(2),Eigenvalue_image(oder),VL(2,2),Eigenvecter(2,2),inverse(2,2),WORK
		real(8) :: det !detmination of eigenvector
		Jk = 2*cos(kx*a)+2*cos(ky*a)
		J0 = z*J
		P(1,1) = 1!-J0*Sz1
		P(1,2) = 0!Jk*Sz1
		P(2,1) = 0!Jk*Sz2
		P(2,2) = 1!-J0*Sz1
		call dgeev('V','V',oder,P,oder,Eigenvalue_real,Eigenvalue_image,VL,oder,Eignvector,oder,WORK,LWORK,INFO)
		! to be continue
		write(*,*)info
		!write(*,*)work	
		write(*,*)'P',P
		write(*,*)Eigenvalue_real		
		!write(*,*)Eigenvalue_image
		write(*,*)Eigenvecter
		write(*,*)VL

end program
