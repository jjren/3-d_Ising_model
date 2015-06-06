module array_mod
	use kinds_mod
	implicit none

	integer(kind=i1),allocatable :: szinit(:,:,:),sz(:,:,:)

	contains
!================================================
!================================================

subroutine allocatearray
	use para_mod

	implicit none
	integer :: error
	
	if(initmode=='s') then
		allocate(szinit(nx,ny,nz),stat=error)
		if(error/=0) stop
	end if

	allocate(sz(nx,ny,nz),stat=error)
	if(error/=0) stop
	
	return

end subroutine allocatearray

!================================================
!================================================

end module array_mod
