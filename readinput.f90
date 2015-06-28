subroutine ReadInput
	
	use communicate
	use mpi
	use para_mod
	use array_mod

	implicit none

	integer :: position1,packsize
	character(len=1),allocatable :: packbuf(:)
	integer :: ierr  ! MPI flag

	
	call master_print_message("Enter ReadInput Subroutine")

	packsize=1000
	allocate(packbuf(packsize))

	if(myid==0) then
		open(unit=10,file="inp",status="old")
		read(10,*) nx,ny,nz
		read(10,*) periodx,periody,periodz
		read(10,*) ntrajs
		read(10,*) nsteps
		read(10,*) initmode
		read(10,*) couplingJ
		read(10,*) temperature
		read(10,*) seed
		read(10,*) intratrajavera
		close(10)

		position1=0
		call MPI_PACK(nx,1,MPI_integer,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(ny,1,MPI_integer,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(nz,1,MPI_integer,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(periodx,1,MPI_integer,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(periody,1,MPI_integer,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(periodz,1,MPI_integer,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(ntrajs,1,MPI_integer,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(nsteps,1,MPI_integer,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(initmode,1,MPI_character,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(couplingJ,1,MPI_real8,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(temperature,1,MPI_real8,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(seed,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
		call MPI_PACK(intratrajavera,1,MPI_integer4,packbuf,packsize,position1,MPI_COMM_WORLD,ierr)
	end if

	call MPI_BCAST(packbuf,packsize,MPI_PACKED,0,MPI_COMM_WORLD,ierr)

	if(myid/=0) then
		position1=0
		call MPI_UNPACK(packbuf,packsize,position1,nx,1,MPI_integer,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,ny,1,MPI_integer,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,nz,1,MPI_integer,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,periodx,1,MPI_integer,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,periody,1,MPI_integer,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,periodz,1,MPI_integer,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,ntrajs,1,MPI_integer,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,nsteps,1,MPI_integer,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,initmode,1,MPI_character,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,couplingJ,1,MPI_real8,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,temperature,1,MPI_real8,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,seed,1,MPI_integer4,MPI_COMM_WORLD,ierr)
		call MPI_UNPACK(packbuf,packsize,position1,intratrajavera,1,MPI_integer4,MPI_COMM_WORLD,ierr)
	end if

	! set the local ntrajs
	if(myid<MOD(ntrajs,nprocs)) then
		localntrajs=ntrajs/nprocs+1
	else
		localntrajs=ntrajs/nprocs
	end if

	if(myid==0) then
		write(*,*) "============================="
		write(*,*) "Input Infomation"
		write(*,*) "nx,ny,nz=",nx,ny,nz
		write(*,*) "periodic boundary",periodx,periody,periodz
		write(*,*) "initial condition mode   ",initmode
		write(*,*) "number of trajectorys",ntrajs
		write(*,*) "number of steps",nsteps
		write(*,*) "coupling J=",couplingJ
		write(*,*) "temperature=",temperature
		write(*,*) "seed=",seed
		write(*,*) "intratrajavera=",intratrajavera
		write(*,*) "============================="
	end if
	
	! allocate work array
	call allocatearray

	deallocate(packbuf)

return

end subroutine ReadInput


