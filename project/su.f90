program main
	implicit none
	real(8) :: b(2),RWOEK(4)
	complex(8) :: A(2,2),WOEK(4)
	integer:: i,j,ok,LWORK
	A(1,:) = [(0d0,0d0),(0d0,1d0)]
	A(2,:) = [(0d0,-1d0),(0d0,0d0)]
	do i=1,2
	write(*,*) (A(i,j),j=1,2)
	end do
	LWOEK =66
	call zheev('V','L',2,A,2,b,WORK,LWORK,RWORK,ok)
	do i=1,2
	 write(*,*) (A(i,j),j=1,2)
	end do
end 
