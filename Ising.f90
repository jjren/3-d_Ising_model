Module Ising
! the core module to do Monte Carlo
	use para_mod
	use kinds_mod
	use array_mod
	use communicate
	use exit_mod
	implicit none
	
	integer,parameter :: nflips=1  ! how many site flip every step
	integer :: flip(3,nflips)      ! flip site index xyz
	real(kind=r8) :: energy        ! every step energy

	contains

!================================================
!================================================

subroutine MC_starter
	
	implicit none
	
	integer :: istep
	real(kind=r8) :: deltaE     ! the energy diffrence
	
	call master_print_message("enter MC_starter")
	
	! set the initial condition
	call InitCond_local
	
	! calculate the energy
	call CalcEnergy(sz)
	call StandardIO(0)

	do istep=1,nsteps,1
		
		! Random flip
		call RandFlip
		
		! calculate the energy
		call CalcFlipEnergy(deltaE)
		
		! Metropolis algorithm if accept the flip
		call Metropolis(deltaE)

		! output
		call StandardIO(istep)
		
	end do

return

end subroutine MC_starter

!================================================
!================================================

subroutine RandFlip
! random flip some site

	implicit none

	integer :: i

	real(kind=r8) :: deltax,deltay,deltaz,random
	integer :: index1
	
	deltax=1.0D0/nx
	deltay=1.0D0/ny
	deltaz=1.0D0/nz

	do i=1,nflips,1
		call random_number(random)
		index1=FLOOR(random/deltax)+1
		if(index1>nx) then
			write(*,*) "========="
			write(*,*) "index1>nx"
			write(*,*) "========="
			stop
		end if
		flip(1,i)=index1

		call random_number(random)
		index1=FLOOR(random/deltay)+1
		if(index1>ny) then
			write(*,*) "========="
			write(*,*) "index1>ny"
			write(*,*) "========="
			stop
		end if
		flip(2,i)=index1

		call random_number(random)
		index1=FLOOR(random/deltaz)+1
		if(index1>nz) then
			write(*,*) "========="
			write(*,*) "index1>nz"
			write(*,*) "========="
			stop
		end if
		flip(3,i)=index1
	end do
	return

end subroutine RandFlip

!================================================
!================================================

subroutine CalcFlipEnergy(deltaE)

! calculate the flip energy site by site
	implicit none

	real(kind=r8) :: deltaE
	integer :: i

	deltaE=0.0D0
	do i=1,nflips,1
		
		! x dim
		if(flip(1,i)<nx) then
			deltaE=deltaE+2.0D0*couplingJ*DBLE(sz(flip(1,i),flip(2,i),flip(3,i))*(-1)*sz(flip(1,i)+1,flip(2,i),flip(3,i)))
		end if
		if(flip(1,i)>1) then
			deltaE=deltaE+2.0D0*couplingJ*DBLE(sz(flip(1,i),flip(2,i),flip(3,i))*(-1)*sz(flip(1,i)-1,flip(2,i),flip(3,i)))
		end if
		if((flip(1,i)==1 .or. flip(1,i)==nx) .and. periodx==1) then
			deltaE=deltaE+2.0D0*couplingJ*DBLE(sz(nx,flip(2,i),flip(3,i))*(-1)*sz(1,flip(2,i),flip(3,i)))
		end if
		
		! y dim
		if(flip(2,i)<ny) then
			deltaE=deltaE+2.0D0*couplingJ*DBLE(sz(flip(1,i),flip(2,i),flip(3,i))*(-1)*sz(flip(1,i),flip(2,i)+1,flip(3,i)))
		end if
		if(flip(2,i)>1) then
			deltaE=deltaE+2.0D0*couplingJ*DBLE(sz(flip(1,i),flip(2,i),flip(3,i))*(-1)*sz(flip(1,i),flip(2,i)-1,flip(3,i)))
		end if
		if((flip(2,i)==1 .or. flip(2,i)==ny) .and. periody==1) then
			deltaE=deltaE+2.0D0*couplingJ*DBLE(sz(flip(1,i),ny,flip(3,i))*(-1)*sz(flip(1,i),1,flip(3,i)))
		end if

		! z dim
		if(flip(3,i)<nz) then
			deltaE=deltaE+2.0D0*couplingJ*DBLE(sz(flip(1,i),flip(2,i),flip(3,i))*(-1)*sz(flip(1,i),flip(2,i),flip(3,i)+1))
		end if
		if(flip(3,i)>1) then
			deltaE=deltaE+2.0D0*couplingJ*DBLE(sz(flip(1,i),flip(2,i),flip(3,i))*(-1)*sz(flip(1,i),flip(2,i),flip(3,i)-1))
		end if
		if((flip(3,i)==1 .or. flip(3,i)==nz) .and. periodz==1) then
			deltaE=deltaE+2.0D0*couplingJ*DBLE(sz(flip(1,i),flip(2,i),nz)*(-1)*sz(flip(1,i),flip(2,i),1))
		end if
		
		! flip this site
		sz(flip(1,i),flip(2,i),flip(3,i))=sz(flip(1,i),flip(2,i),flip(3,i))*-1
	
	end do
return

end subroutine CalcFlipEnergy

!================================================
!================================================

subroutine Metropolis(deltaE)
! MC metropolis algorithm 
! update the energy

	implicit none
	
	real(kind=r8) :: deltaE
	
	! local 
	real(kind=r8) :: expnegTE,random
	logical :: ifaccept
	integer :: i

	ifaccept=.false.
	if(deltaE<0.0D0) then
		ifaccept=.true.
	else
		expnegTE=exp(-1.0D0/temperature*deltaE)  ! exp^(-beta*deltaE)
		call random_number(random)
		if(random<expnegTE) then
			ifaccept=.true.
		end if
	end if

	if(ifaccept==.true.) then
		energy=energy+deltaE
	else
		! recover the last step sz
		do i=1,nflips,1
			sz(flip(1,i),flip(2,i),flip(3,i))=sz(flip(1,i),flip(2,i),flip(3,i))*(-1)
		end do
	end if
return

end subroutine Metropolis

!================================================
!================================================

subroutine CalcEnergy(szlocal)

! nearest neighbour coupling
	implicit none
	
	integer :: ix,iy,iz
	integer(kind=i1) :: szlocal(nx,ny,nz)

	energy=0.0D0

	do iz=1,nz,1
	do iy=1,ny,1
	do ix=1,nx,1
		if(iz<nz) then
			energy=energy+couplingJ*DBLE(szlocal(ix,iy,iz)*szlocal(ix,iy,iz+1))
		end if
		if(iy<ny) then
			energy=energy+couplingJ*DBLE(szlocal(ix,iy,iz)*szlocal(ix,iy+1,iz))
		end if
		if(ix<nx) then
			energy=energy+couplingJ*DBLE(szlocal(ix,iy,iz)*szlocal(ix+1,iy,iz))
		end if

		if(periodz==1 .and. iz==nz) then
			energy=energy+couplingJ*DBLE(szlocal(ix,iy,nz)*szlocal(ix,iy,1))
		end if
		if(periody==1 .and. iy==ny) then
			energy=energy+couplingJ*DBLE(szlocal(ix,ny,iz)*szlocal(ix,1,iz))
		end if
		if(periodx==1 .and. ix==nx) then
			energy=energy+couplingJ*DBLE(szlocal(nx,iy,iz)*szlocal(1,iy,iz))
		end if
	end do
	end do
	end do
return

end subroutine CalcEnergy

!================================================
!================================================

subroutine RandomCond(szlocal,mx,my,mz)
	
	implicit none

	! local
	integer :: mx,my,mz
	integer(kind=i1) szlocal(mx,my,mz)
	integer :: ix,iy,iz
	real(kind=r8) :: random

	do iz=1,mz,1
	do iy=1,my,1
	do ix=1,mx,1
		call random_number(random)
		if(random>0.5) then
			szlocal(ix,iy,iz)=1
		else
			szlocal(ix,iy,iz)=-1
		end if
	end do
	end do
	end do

return
end subroutine RandomCond

!================================================
!================================================

subroutine InitCond_global

! set the initial condition in the InitMode same
	
	use mpi
	implicit none
	integer :: ierr

	if(myid==0) then
		call master_print_message("Enter InitCond subroutine")
		call RandomCond(szinit,nx,ny,nz)
	end if

	call MPI_BCAST(szinit,nx*ny*nz,MPI_integer1,0,MPI_COMM_WORLD,ierr)

return
end subroutine InitCond_global

!================================================
!================================================

subroutine InitCond_local

	implicit none

	if(initmode=='s') then
		sz=szinit
	else if(initmode=='r') then
		call RandomCond(sz,nx,ny,nz)
	else
		call exit_Ising(sigAbort,"InitMode Wrong!")
	end if
return

end subroutine InitCond_local

!================================================
!================================================

subroutine StandardIO(istep)

	implicit none

	integer :: istep
	! local
	integer :: gap=100
	character(len=50) :: filename
	integer :: ordpara
	integer :: i,j,k
	logical :: alive
	
	if(mod(istep,gap)==0) then
		
		write(filename,'(i3.3,a10)') myid,'energy.tmp'
		inquire(file=trim(filename),exist=alive)
		if(alive) then
			open(unit=101,file=trim(filename),status="old",position="append")
		else
			open(unit=101,file=trim(filename),status="replace")
		end if
		write(101,*) istep,energy
		close(101)

		write(filename,'(i3.3,a10)') myid,'config.tmp'
		inquire(file=trim(filename),exist=alive)
		if(alive) then
			open(unit=102,file=trim(filename),status="old",position="append")
		else
			open(unit=102,file=trim(filename),status="replace")
		end if
		
		write(102,*) "=============",istep,"=============="
		ordpara=0
		do i=1,nz,1
		do j=1,ny,1
		do k=1,nx,1
			ordpara=ordpara+sz(k,j,i)
			write(102,*) sz(k,j,i)
		end do
		end do
		end do
		close(102)
		
		write(filename,'(i3.3,a11)') myid,'ordpara.tmp'
		inquire(file=trim(filename),exist=alive)
		if(alive) then
			open(unit=103,file=trim(filename),status="old",position="append")
		else
			open(unit=103,file=trim(filename),status="replace")
		end if
		write(103,*) istep,ordpara
		close(103)
	end if

	return
end subroutine StandardIO

!================================================
!================================================
end module Ising
