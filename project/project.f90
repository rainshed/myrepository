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

	real(dp) :: k
	k = 2*pi/a
end module
!--------------------------------------------------------------------
module array
	real,external :: vec
	integer,parameter :: dp = selected_real_kind(8)
	contains
	function vec(init,ed,inter)
		implicit none
		real(dp) :: init,ed,inter
		real(dp) :: num = (ed-init)/inter+1
		real(dp) :: 





end module array

program main
	implicit none
	use constants

end
