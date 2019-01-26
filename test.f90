program main
  implicit none 
  integer :: g,beta
  parameter(g=20,beta=1)
  real,external :: random01 
  integer  i
  integer  j
  real(kind=8), external :: ranw
  integer,parameter :: rows=20
  integer,parameter :: cols=20
  integer :: Lattice(rows,cols)
  call RanGen_Ini 
  do i=1,rows
   do j=1,cols
    Lattice(i,j) =1
   end do 
  end do 
  write(*,*) Lattice
end program

function random01( )
  call RanGen_Ini 
  real a
  call random_number(a)
  if (a<0.5) then
    random01 = 1
  return
  else
    random01 = 2
  return
end if
end function
  
!function sample2(vector)
!  implicit none
!  vector = vector(random_number())
!end function
!------------------------------------------------------------------
  subroutine RanGen_Ini
    !
    ! Purpose
    ! =======
    !   Initialize the Random Number Generator. 
    !
    integer :: clock

    ! ... local vars ... 
    integer :: i, n
    integer, allocatable :: seed(:)
    
    call SYSTEM_CLOCK(clock)
    call RANDOM_SEED(size = n)
    ALLOCATE(SEED(n))
    do i = 1, n
       SEED(i) = clock + 37*(i-1)
    end do
    call RANDOM_SEED(PUT = SEED)
    DEALLOCATE(SEED)

  end subroutine RanGen_Ini

  !------------------------------------------------------------------
  real(kind=8) function ranw()
    !
    ! Purpose
    ! =======
    !   random number generator, use the f90 intrinsic routin, random numbers
    !   distuributes in (0, 1)

    call random_number(ranw)
  end function ranw

