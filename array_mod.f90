module array_mod
	use kinds_mod
	implicit none

	integer(kind=i1),allocatable :: szinit(:,:,:),sz(:,:,:)
	integer(kind=i4),allocatable :: sum_ordpara(:),sum_ordpara2(:)
	real(kind=r8),allocatable :: sum_energy(:),sum_energy2(:)
	integer(kind=i4),allocatable :: laststep_ordpara(:)

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
	
	! average order parameter and energy on every process
	allocate(sum_ordpara(0:nsteps),stat=error)
	if(error/=0) stop
	sum_ordpara=0
	allocate(sum_ordpara2(0:nsteps),stat=error)
	if(error/=0) stop
	sum_ordpara2=0
	allocate(sum_energy(0:nsteps),stat=error)
	if(error/=0) stop
	sum_energy=0.0D0
	allocate(sum_energy2(0:nsteps),stat=error)
	if(error/=0) stop
	sum_energy2=0.0D0
	allocate(laststep_ordpara(-nx*ny*nz:nx*ny*nz),stat=error)
	laststep_ordpara=0
	if(error/=0) stop
return

end subroutine allocatearray

!================================================
!================================================

end module array_mod
