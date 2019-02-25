
!constans set
!----------------------------------------------------------------------
module constants
	implicit none
	integer,parameter :: dp = selected_real_kind(8)
	real(dp),parameter :: a=1       !Lattice constant
	real(dp),parameter :: Sz=0.5d0  !anverage z component of spin
	real(dp),parameter :: S=0.5d0   !spin of system
	real(dp),parameter :: hbar=1    !Plank constant
	real(dp),parameter :: N=1       !Number of atoms
	real(dp),parameter :: delta=1   !self energy
	real(dp),parameter :: z=6       !coordination number
	real(dp),parameter :: E0=0      !energy of ground state
	real(dp),parameter :: J=1       !exchange constant 
	real(dp),parameter :: pi=3.1415926        
	real(dp) :: k = 2*pi/a
end module
!-------------------------------------------------------------------




!operation with array
!-----------------------------------------------------------------------------------------------
module array
	!real,external :: vec
	integer,parameter :: dp = selected_real_kind(8)
	contains

!create a array with given initial value, end value and interval
!---------------------------------------------------------------------
	function vec(init,ed,inter) result(arr)
		implicit none
		integer :: i
		real(dp) :: init,ed,inter
		integer :: num
		real(dp) :: arr(int((ed-init)/inter+1))
		num = (ed-init)/inter+1
		do i=1,num
			arr(i) = init+(i-1)*inter
		end do
	end function
!-------------------------------------------------------------------------

end module array
!---------------------------------------------------------------------------------------------



module calc
	use constants
contains
!given k and omega,calculate the spectrum
!--------------------------------------------------------------------------
	real function spe(k1,omega)
		real(dp) :: k1
		real(dp) :: omega
		real(dp) :: Ek
		Ek = J*S*z*sqrt(1-(4*cos(k1*a)/z)**2) !this formula is different with different crystal structure
		spe = 2d0*hbar*Sz*delta/((hbar*omega-Ek)**2+delta**2)
	end function spe
!---------------------------------------------------------------------------


end module calc
!------------------------------------------------------------------------------------------------
program main
	use constants
	use array
	use calc
	implicit none
	
end
