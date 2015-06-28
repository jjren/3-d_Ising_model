Program Main

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!% this is a Monte Carlo Program treating 3-d Ising model%
!% jiajunren0522@126.com                                 %
!% 2015/6/27                                             %
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

	use Ising
	use mpi
      USE MKL_VSL_TYPE
      USE MKL_VSL
	use blas95
	use f95_precision

	implicit none

	real(kind=r8) :: starttime,endtime
	integer :: itraj,i
	real(kind=r8),allocatable :: &
		sum0_energy(:)   ,  &         ! 0 process store sum of every step energy
		sum0_rordpara(:) ,  &         ! 0 process store sum of every step |order parameter|
		sum0_energy2(:)  ,  &         ! 0 process store sum of every step energy^2
		sum0_rordpara2(:)             ! 0 process store sum of every step order parameter^2
	integer(kind=i4),allocatable :: &
		sum0_ordpara(:)  ,  &
		sum0_ordpara2(:) ,  &
		laststep_ordpara0(:)          ! 0 process store last step order paramter 
	real(kind=r8) :: resM,resE,tmp,resM2,resE2    ! result |M| |E| |M|^2 |E|^2
	integer :: ierr
	
	! initiate the MPI enviroment
	call init_communicate
	starttime=MPI_WTIME()
	
	! Read Parameters
	call ReadInput

	! initiate the random parallel seed
	method= VSL_BRNG_MT2203  ! can change this method MKL lib
	errcode = vslnewstream(stream,method+myid,seed)

	! the same initial configuration
	if(initmode=='s') then
		call InitCond_global
	end if

	! the I/O file
	call fileopen

	! every process do no communicate
	do itraj=1,localntrajs,1
		call MC_starter
	end do
	
	call fileclose

	! 0 process output
	if(myid==0) then
		allocate(sum0_ordpara(0:nsteps))
		allocate(sum0_rordpara(0:nsteps))
		allocate(sum0_ordpara2(0:nsteps))
		allocate(sum0_rordpara2(0:nsteps))
		allocate(sum0_energy(0:nsteps))
		allocate(sum0_energy2(0:nsteps))
		allocate(laststep_ordpara0(-nx*ny*nz:nx*ny*nz))   ! in each |M| fold how many trajectories
	end if

	call MPI_REDUCE(sum_ordpara(0),sum0_ordpara(0),nsteps+1,MPI_INTEGER4,MPI_SUM,0,MPI_COMM_WORLD,ierr)
	call MPI_REDUCE(sum_energy(0),sum0_energy(0),nsteps+1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
	call MPI_REDUCE(sum_ordpara2(0),sum0_ordpara2(0),nsteps+1,MPI_INTEGER4,MPI_SUM,0,MPI_COMM_WORLD,ierr)
	call MPI_REDUCE(sum_energy2(0),sum0_energy2(0),nsteps+1,MPI_REAL8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
	call MPI_REDUCE(laststep_ordpara(-nx*ny*nz),laststep_ordpara0(-nx*ny*nz),2*nx*ny*nz+1,MPI_INTEGER4,MPI_SUM,0,MPI_COMM_WORLD,ierr)
	
	if(myid==0) then
		tmp=DBLE(nx*ny*nz)
		open(unit=11,file="stdout.out",status="replace")
		do i=0,nsteps,1
			! the inter traj average
			sum0_rordpara(i)=DBLE(sum0_ordpara(i))/DBLE(ntrajs)/tmp       ! ordpara/nsite
			sum0_energy(i)=sum0_energy(i)/DBLE(ntrajs)/tmp
			sum0_rordpara2(i)=DBLE(sum0_ordpara2(i))/DBLE(ntrajs)/tmp/tmp ! ordpara^2/nsite^2
			sum0_energy2(i)=sum0_energy2(i)/DBLE(ntrajs)/tmp/tmp    
			if(mod(i,1000)==0) then  ! output every 1000 steps
				write(11,*) sum0_rordpara(i),sum0_energy(i)
			end if
		end do
		
		! the intra traj average
		resM=sum(sum0_rordpara(nsteps-intratrajavera+1:nsteps))
		resM=resM/DBLE(intratrajavera)

		resE=sum(sum0_energy(nsteps-intratrajavera+1:nsteps))
		resE=resE/DBLE(intratrajavera)
		!
		resM2=asum(sum0_rordpara2(nsteps-intratrajavera+1:nsteps))
		resM2=resM2/DBLE(intratrajavera)

		resE2=asum(sum0_energy2(nsteps-intratrajavera+1:nsteps))
		resE2=resE2/DBLE(intratrajavera)

		write(11,*) "===================="
		do i=-nx*ny*nz,nx*ny*nz,1
			write(11,*) i,laststep_ordpara0(i)
		end do
		write(11,*) "seed=",seed
		write(11,*) "resM and resE============"
		write(11,*) resM,resE
		write(11,*) "resM2 and resE2============"
		write(11,*) resM2,resE2
		close(11)
	end if

	errcode = vsldeletestream(stream)
	
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	endtime=MPI_WTIME()
	call master_print_message(endtime-starttime,"RUMTIME:")
	call exit_Ising(0,"Program Ising end successfully")

end program 

