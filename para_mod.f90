module para_mod

	use kinds_mod
	implicit none
	
	integer :: nx,ny,nz                 ! nx is the x dimension nsite
	integer :: periodx,periody,periodz  ! if have period boundary conditon in the xyz dimension(=1)
	integer :: ntrajs,nsteps            ! number of trajectorys
	character(len=1) :: initmode        ! the initial condition mode
	real(kind=r8) :: couplingJ,temperature  

	integer :: localntrajs              ! every process ntrajs


end module para_mod
