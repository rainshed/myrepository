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

!create a array with given initial value,end value and length of array
!-------------------------------------------------------------------------
	function vecn(init,ed,num) result(arr)
		implicit none
		integer :: i,num
		real(dp) :: init,ed,inter,arr(num)
		inter = (ed-init)/(real(num)-1d0)
		do i = 1,num
		  arr(i) = init + (i-1)*inter
		end do
	end function vecn
!--------------------------------------------------------------------------

!compute the inverse of complex matrix A
!---------------------------------------------------------------------------
	function inverse(A,N) result(inv)
		implicit none
		integer :: N
		complex(8) :: A(N,N),work(N),inv(N,N)	
		integer :: info1,info2,IPIV(N),Lwork =1000 
		inv = A
		call zgetrf(N,N,inv,N,IPIV,info1)
		if (info1 /= 0) then 
			write(*,*) 'zgetrf goes wrong'
		end if
		call zgetri(N,inv,N,IPIV,work,lwork,info2)
		if (info2 /= 0) then
			write(*,*) 'zgetri goes wrong'
		end if
	end function inverse
!----------------------------------------------------------------------------

end module array
!---------------------------------------------------------------------------------------------
