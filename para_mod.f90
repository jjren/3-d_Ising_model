include 'mkl_vsl.f90'
module para_mod

      USE MKL_VSL_TYPE
      USE MKL_VSL
	use kinds_mod
	implicit none
	
	integer :: nx,ny,nz                 ! nx is the x dimension nsite
	integer :: periodx,periody,periodz  ! if have period boundary conditon in the xyz dimension(=1)
	integer :: ntrajs,nsteps            ! number of trajectorys
	character(len=1) :: initmode        ! the initial condition mode
	real(kind=r8) :: couplingJ,temperature  
	integer :: seed                     ! random_seed
	integer :: localntrajs              ! every process ntrajs
	integer :: intratrajavera           ! intra traj average Nconfigurations
	
	! vsl para
	integer :: errcode,method
	type(VSL_STREAM_STATE) :: stream

end module para_mod
