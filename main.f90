Program Main
! this is a Monte Carlo Program treating 3-d Ising model
	use Ising
	use mpi

	implicit none

	real(kind=r8) :: starttime,endtime
	integer :: itraj
	integer :: ierr
	
	! initiate the MPI enviroment
	call init_communicate
	starttime=MPI_WTIME()
	
	! initiate the random seed
	! parallel random
	call sleep(myid)
	call random_seed()

	! Read Parameters
	call ReadInput
	
	! the same initial configuration
	if(initmode=='s') then
		call InitCond_global
	end if

	! every process do
	! no communicate
	do itraj=1,localntrajs,1
		call MC_starter
	end do

	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	endtime=MPI_WTIME()
	call master_print_message(endtime-starttime,"RUMTIME:")
	call exit_Ising(0,"Program Ising end successfully")

end program 

