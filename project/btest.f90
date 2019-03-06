
!constans set
!----------------------------------------------------------------------
module constants
	implicit none
	integer,parameter :: dp = selected_real_kind(8)
	real(dp),parameter :: a=1       !Lattice constant
	real(dp),parameter :: Sz=1  !anverage z component of spin
	real(dp),parameter :: S=1   !spin of system
	real(dp),parameter :: hbar=1    !Plank constant
!	real(dp),parameter :: N=1       !Number of atoms
	real(dp),parameter :: delta=1   !self energy
	real(dp),parameter :: z=4       !coordination number
	real(dp),parameter :: E0=0      !energy of ground state
	real(dp),parameter :: J=1       !exchange constant 
	real(dp),parameter :: pi=3.1415926        
	real(dp),parameter :: Sz1=1,Sz2=-1
	!-----------------------------------------------------------------------
end module
!-------------------------------------------------------------------




!operation with array
!-----------------------------------------------------------------------------------------------
module array
	!real,external :: vec
	!integer,parameter :: dp = selected_real_kind(8)
	use constants
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


program main
	use constants
!---------------------------------------------------------------------------
		real(8) :: Jk,J0,kx=1d0,ky=1d0,RWOEK(4)
		integer,parameter :: order = 2
		complex(8) :: omega,P(order,order),green(2,2)
		integer :: LWORK=100,WORK(2),info
		complex(8) :: Eigenvalue,VL,Eigenvecter(order,order)
		complex(8) :: det !detmination of eigenvector
		complex(8) :: inverse(order,order) !inverse of Eigenvecter
		integer :: m,n,l,k
		complex(8) :: F(2,2)
		F(1,1:2) = (/2*Sz1,0d0/)
		F(2,:) = (/2*Sz2,0d0/)
		Jk = 2*cos(kx*a)+2*cos(ky*a)
		J0 = z*J
		P(1,1) =1! complex(-J0*Sz1,0)
		P(1,2) =2! complex(Jk*Sz1,0)
		P(2,1) =2! complex(Jk*Sz2,0)
		P(2,2) =1! complex(-J0*Sz1,0)
		call dgeev('N','V',order,P,order,Eigenvalue,VL,order,Eigenvector,order,WORK,LWORK,RWORK,INFO)
		if (info .ne. 0) then
			write(*,*) 'dgeev error'
		end if
		p = Eigenvector!P is Eigenvector matrix
		write(*,*) P
		!compute the inverse of Eigenvectors
	!	det = P(1,1)*P(2,2)-P(1,2)*P(2,1)
	!	inverse(1,1) = P(2,2)/det
	!	inverse(2,1) = -P(1,2)/det
	!	inverse(1,2) = -P(2,1)/det
	!	inverse(2,2) = P(1,1)/det
	!	!compute green function
	!	do m = 1,2
	!		do n = 1,2
	!			green(m,n)=0
	!			do l=1,2
	!				do k=1,2
	!				green(m,n) =green(m,n)+(P(m,l)*inverse(l,k))*F(l,k)/(omega-Eigenvalue(l))	
	!				end do
	!			end do
	!			green(m,n) = P(m,n) + green(m,n)
	!		end do
	!	end do
end
!------------------------------------------------------------------------------------------------
