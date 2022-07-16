

! Ahsan Zeb, 20 oct 2018
! HTC model with N molecules, m exciton-photon excitations, M+1 vibrational levels.

! We will use permutation symmetry:
! blocks with respect to number of up spins:
!	pth block: p up spins.
!					p up sites are permutaiton symmetric; N-p down sites are perm symm. 
!					m1=Min(m,N)+1 total blocks; i.e., p=0,1,2,3,...,m1
! further details in doc.


	program chi2
	use modmain
	use maps, only: getmap, writemap
	use bases, only: PermSymBasis, writebasis, writebasisf
	use hamiltonian, only: MakeHhtc, HamParts
	use diag, only: diagonalise
	use mpi
	use correlation, only: tcorr, rwallnodesx
	use hamiltonian1, only: MakeH1,Ham1Parts,glueHblocks,glueHblocksD
	use correlationd, only: tcorrd
	use dmat, only: cdms, rwallnodes, cdms1, dmphot
	
	implicit none

	double precision:: tottime, dtau, ed,ew
	logical :: newm, goforcdms
	! mpi
	integer:: ierr, nj
	integer(kind=MPI_OFFSET_KIND),
     .           dimension(:), allocatable :: disp 

	! 
	integer :: i, ijob, ij1, njl,k, j,ll,il,idark,nnzcoef
	integer :: ix,jux,ind 
!======================================================================
	call MPI_INIT(ierr)
	!find out MY process ID, and how many processes were started.
	call MPI_COMM_RANK (MPI_COMM_WORLD, node, ierr)
	call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)

	if(node==0) 	call timestamp(node) ! time stamp

	!write(*,*) "Node = ",node," num_procs, ntp = ",num_procs, ntp 
	! readinput file
	call input(node)

	! just write basis and map?
	if(writebmap) then
		call MPI_FINALIZE(ierr) ! done only at master node
		if(node==0) then
			call getmap(n,mv)
			call PermSymBasis(n,mv)
			! write output files
			call writemap()
			call writebasis()
		endif
		stop "htc: basis and map written to files... done!"
	endif


	! things need to be done only once	
	! allocate space for basis, and calc maps, etc.
	!---------- called once for a given n,mv ---------
	call getmap(n,mv) !,ntot) ! ntot output
	call PermSymBasis(n,mv)

	call setparam()
	
	call setjobs(node)

	!total jobs on this node
	nj = jobs(node)%njobs

	!allocate eig
	allocate(eig(nj))

c	m1max = 0
c	do i=1,nj
c		ijob = jobs(node)%i1+i-1
c		m1max = max(m1max,param(ijob)%m);
c	enddo
c	m1max=min(nsym,m1max);

	ij1 = 0;
	! until its useful to increase m1max, set it to the required m1.	
	newm = .true. ! new m ?
	!+++++++++++++++++++++START IJONB LOOP ++++++++++++++++++++++++++
	do i=1,nj,1

		!write(*,*)'jobs(node)%i1 , i2= ',jobs(node)%i1,jobs(node)%i2
		if(node==0)write(*,'(a,i5,a,i5)')'Node 0: job ',i,' out of ',nj
		ijob = jobs(node)%i1+i-1
		
		m = param(ijob)%m;

		! set max m1 possible; basis, map, data can be reused for smaller m.
		m1max=min(nsym,m);
		!write(6,*)'main: nsym, n, ndummy, Nex = ',nsym, n, ndummy, m
		!write(6,*)'main: m1max=',m1max

	! task: 
	! Response functions: 100 series
	! 101 abs, 102 emission
	! 103 hopping up, 104 hopping down
	! 105 density up, 106 density down
	! Matrix Elements: 200 series
	! 201 abs, 202 emission
	! 203 hopping without vibration, default N in range [1:2N:1]
	! Density Matrices: 300 series
	! 301 conditional vibrational dms like objects in symmetric space, mode=1
	! 302 conditional vibrational dms like objects in full space needed for hopping, mode=2
	! 310 photon reduced dm in the condensate state

	! call internal routines to do the job
	select case(task)
	case(101) ! light absorption by the condensate
	 call absorption(i) !(i,ijob, n, nsym,m,m1max, mv, newm,chi)
	case(102) ! PL/light emission from the condensate
	 call emission(i)
	case(103,104)
	 call hoppingresponse(task)
	case(105,106)
	 call densityresponse(task)
	case(201)
	 call absmatelem(i)
	case(202)
	 call plmatelem(i)
	case(203)
	 call hopmatelem(i)
	case(301, 302)
	 call getcdms(i, ijob, mode, newm, ij1)
	case(303)
	 !call dmphot()
	 	stop "main: task 303 dmphot() not available yet.... DO IT!"
	case default
		stop "main: task not recongnised....!"
	end select
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	if(1==1)  then
	! temporary, to write evals. run in serial
	write(6,'(a5x100000f25.15)') "EVALS: ",
     .  	eig(i)%eval(1:min(20,eig(i)%n2))
	write(*,*)'eig(i)%n2 = ',eig(i)%n2
	endif

	!newm = setnewm()
	if(i < nj) then
		if (param(ijob+1)%m .ne. m ) newm = .true.
	endif

	end do ! i jobs
	!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	call MPI_FINALIZE(ierr)

	! combine all cdm files written by diff nodes???

	! read files	written by all nodes and write a single file.
	if(node==0) then
		!if(writewfs) then
		!	call iowfcombine(nact)
		!endif
		call writeout(task)
		write(*,*)"chi2: everything done.... " 
		call timestamp(node)
	endif

!======================================================================
	contains
!======================================================================
	subroutine printnode(node)
	implicit none
	integer, intent(in) :: node
	if(node==0) then
		write(*,'("main: this is the master node ",I5)')node
	else
		write(*,'("main: this is the slave node ",I5)') node
	endif
	return
	end 	subroutine printnode
!=============================================	
	subroutine CheckError(ierr, node, oper)
	implicit none
	integer, intent(in):: ierr, node
	character(len=*):: oper
	if(ierr .ne. 0) then
		write(*,'("main: mpi operation ",a," at node ",I3," ier = ",I3)')
     .          oper, node, ierr
		stop
	endif
	
	! comment
	!write(*,'("main: mpi operation ",a," at node ",I3," ier = ",I3)')
  !   .          oper, node, ierr
	return
	end 	subroutine CheckError
!=============================================

!=============================================
	subroutine setparam()
	! sets param(:).
	! every ndoes sets its own copy, and uses its part of the job,
	! node=node will use param(jobs(node)%i1:jobs(node)%i2)
	implicit none
	integer :: i,i1,i2,i3,i4,i5

	! total number of jobs
	njobs = nm * nwr * ndel * nlam * nwv;
	allocate(param(njobs))

	!if(node==0) write(*,*) 'njobs = ',njobs

	i = 0;
	do i1=1,nm,1
		do i2=1,nwr,1
			do i3=1,ndel,1
				do i4=1,nlam,1
					do i5=1,nwv,1
						i = i + 1;
						param(i)%m = ms(i1)
						param(i)%wr = wrs(i2)
						param(i)%del = dels(i3)
						param(i)%lam = lams(i4)
						param(i)%wv = wvs(i5)
						param(i)%lamd = lamd
					end do
				end do
			end do
		end do
	end do

	return
	end 	subroutine setparam
!=============================================

!=============================================
	subroutine setjobs(node)
	implicit none
	integer, intent(in) :: node
	integer :: i,j,k

	allocate(jobs(0:num_procs-1))
	
	jobs(:)%njobs = 0; ! set on all processes,
									! independently, and used only for the relevant node
	jobs(:)%i1 = 0; 	jobs(:)%i2 = 0;
	
	! assign jobs one by one
	k =0; ! do i=1,large can be replaced by an open do loop, or a do while loop.
	do i=1, njobs ! i=1, a large dummy number(something equal or bigger than njobs)
		do j=0,num_procs-1
			k = k + 1;
			jobs(j)%njobs = jobs(j)%njobs + 1
			if(k == njobs) exit
		enddo
		if(k == njobs) exit
	end do

	!if(node==0) write(*,*) 'jobs(:)%njobs = ',jobs(:)%njobs

	! assign interval range, only jobs(node) will actually be used later by this node.
	jobs(0)%i1= 1;
	jobs(0)%i2 = jobs(0)%i1 + jobs(0)%njobs - 1; 
	do j=1,num_procs-1,1
		jobs(j)%i1 = jobs(j-1)%i2+1
		jobs(j)%i2 = jobs(j)%i1 + jobs(j)%njobs - 1
	enddo

	return
	end 	subroutine setjobs

!=============================================
	subroutine checkifgoforcdms(nj,i,goforcdms)
	implicit none
	integer, intent(in) :: nj,i
	logical, intent(out) :: goforcdms
	integer :: valueRSS, k

	goforcdms = .false.;
	! memory limit read from input file
	
	! estimate the size of memory allocated
	!call system_mem_usage(valueRSS)
	do k=1,nj
		valueRSS = valueRSS + size(eig(k)%evec)*8; ! 8 bytes for a double precision
	end do
	!write(*,*)'size of eig, memlim = ',valueRSS, memlim
	!if (valueRSS == -1 ) stop "main: valueRSS not fouond!"
	if ((i.eq.nj) .or. (valueRSS .ge. memlim)) then
		goforcdms = .true.
	endif

	return
	end 	subroutine checkifgoforcdms
!=============================================
	subroutine system_mem_usage(valueRSS)
	implicit none
	!use ifport !if on intel compiler
	integer, intent(out) :: valueRSS
	character(len=200):: filename=' '
	character(len=80) :: line
	character(len=8)  :: pid_char=' '
	integer :: pid
	logical :: ifxst	

	valueRSS=-1    ! return negative number if not found

	!--- get process ID
	pid=getpid()
	write(pid_char,'(I8)') pid
	filename='/proc/'//trim(adjustl(pid_char))//'/status'

	!--- read system file
	inquire (file=filename,exist=ifxst)
	if (.not.ifxst) then
		write (*,*) 'system file does not exist'
		return
	endif

	open(unit=100, file=filename, action='read')
	do
		read (100,'(a)',end=120) line
		if (line(1:6).eq.'VmRSS:') then
			read (line(7:),*) valueRSS
			exit
		endif
	enddo
120	continue
	close(100)
	return
	end subroutine system_mem_usage
!=============================================
!=============================================================
!=============================================
! write lowest nev eigenvectors, and eigvals etc
	subroutine iowfnode(i,ij) ! i,ij: local and global job indeces
	implicit none
	integer, intent(in) :: i,ij
	character :: rank*30, fname*100
	
	write(rank,'(i6.6)') node
	fname = 'wf'//'-'//trim(rank)
	! write unformatted file
	if(i==1) then
		open(1,file=trim(fname), form="unformatted", action="write")
	else
		open(1,file=trim(fname), form="unformatted", action="write",
     .                                      position="append")
	endif
	write(1)param(ij)%m,param(ij)%wr,
     .      param(ij)%del,param(ij)%lam, param(ij)%wv
	!...................... basis%pntr ............
	write(1) basis%pntr(0:m1max+1)
	!...................... Vx and its coordinates ............
	write(1) Hg%nnz ! size
	!write(*,*) Hg%coo1
	write(1) Hg%coo1 ! coo of Vx (the same as that of Hg)
	write(1) Hg%coo2 
	write(1) Vx
	!...............................................
	write(1) eig(i)%ntot, eig(i)%n2
	write(1) eig(i)%eval 
	write(1) eig(i)%evec	
	close(1)
	return
	end 	subroutine iowfnode
	!=============================================================






	!=============================================================
	subroutine iowfcombine()
	implicit none
	integer :: i,ij,k,j, ntot, nev, nsize
	character :: rank*30, fname*100
	!local
	integer :: nnz
	integer, dimension(:), allocatable :: coo,pntr
	double precision, dimension(:), allocatable :: val,vec,vx
	double precision, dimension(4) :: rvals
	character :: dir*30
	logical :: ex

		if(n<10) then
			write(dir,'(I1)') n
		elseif(n<100) then
			write(dir,'(I2)') n
		elseif(n<1000) then
			write(dir,'(I3)') n
		elseif(n<10000) then
			write(dir,'(I4)') n
		else
			write(*,*)'ERROR: n too large!!!'
			stop
		endif

	dir = 'n-'//trim(dir);
	inquire (file=trim(dir), exist=ex)
	if(.not. ex) call system('mkdir '//trim(dir))

	! write output file with data from all nodes
	fname = trim(dir)//'/wfs.dat'
	open(10,file=trim(fname), form="unformatted", 
     .                          action="write",position='append')

	nsize = 0;
	allocate(val(1))
	allocate(vec(1))
	allocate(coo(1))
	allocate(vx(1))

	do i=0,min(njobs,num_procs)-1
		write(rank,'(i6.6)') i
		fname = 'wf-'//trim(rank)
		! read unformatted file
		open(1,file=trim(fname), form="unformatted", action="read")
		do ij=1,jobs(i)%njobs
			read(1) j, rvals
			write(10) j, rvals

			allocate(pntr(0:min(j,n)+1))
			read(1) pntr
			write(10) pntr
			deallocate(pntr)
			
			read(1) nnz
			write(10) nnz

			if(size(coo) /= nnz) then
				deallocate(coo, vx)
				allocate(coo(nnz))
				allocate(vx(nnz))
			endif

			read(1) coo
			write(10) coo
			read(1) coo
			write(10) coo
			read(1) vx
			write(10) vx

			read(1) ntot, nev
			write(10) ntot, nev

			if( nsize /= ntot*nev) then
				nsize=ntot*nev
				deallocate(val,vec)
				allocate(val(nev))
				allocate(vec(nsize))		
			endif
			read(1) val 
			read(1) vec	
			!write(*,*)'ntot = ', ntot
			!write(*,'(a)') 'htc: main: eigenvalues:'
			!write(*,'(100000f4.1)') val
			!write(*,'(a)') 'htc: main: eigenvectors:'
			!write(*,'(100000f4.1)') vec	

			! write to final output file
			write(10) val 
			write(10) vec	
		end do
		close(1, status='delete')
	end do
	close(10)

	return
	end 	subroutine iowfcombine
	!=============================================================	
	subroutine checktrace(l,a)
	implicit none
	integer, intent(in) :: l
	double precision, dimension(l,l), intent(in) :: a
	double precision :: tr
	integer :: i
	
	tr =0.0d0
	do i=1,l
		tr = tr + a(i,i)
	end do

	write(*,*) 'tr = ',tr
	
	return
	end 	subroutine checktrace
	!=============================================================	
	! hopping matrix elements, without vibration case, i.e., mv=0. 
	subroutine matelem(i,n,m)
	implicit none
	integer, intent(in) :: i,n,m
	double precision, dimension(min(nev,ntotg),2) :: amp
	integer :: nev1, j,k
	double precision :: x1,x2
	
	nev1 = min(nev,ntotg);

	do j=1,nev1
		amp(j,1) = sum(eig(i)%evec(1:ntotdn,1) * 
     .                   eig(i)%evec(1:ntotdn,j))**2
		amp(j,2) = sum(eig(i)%evec(ntotdn+1:ntotg,1) *
     .                   eig(i)%evec(ntotdn+1:ntotg,j))**2
	end do

	open(12,file='matelem.dat',action='write',position='append')
		! output file data order: LP dn, LP up, all requested excited up=dn
		write(12,'(2i10,2x,10000f25.15)') n, m,
     .             amp(1,1), amp(:,2) ! for non-diagonal transitions, both chan give the same amp, just diff sign. 
	close(12)

	open(12,file='evals.dat',action='write',position='append')
		write(12,'(10000f25.15)') eig(i)%eval(1:nev1)
	close(12)



	if(1==0) then
	open(13,file='evecs.dat',action='write',position='append')
	do j=1,nev1
		!write(13,'(i5,3x,f15.10)')j,eig(i)%eval(j)
		write(13,'(a,3x,10000f8.3)')'dn:',
     .  eig(i)%evec(1:ntotdn,1)*eig(i)%evec(1:ntotdn,j)
		!write(13,'(a,3x,10000f8.3)')'up:',eig(i)%evec(ntotdn+1:ntotg,j)
	end do
	write(13,*)'--------------------------'
	do j=1,nev1
		!write(13,'(i5,3x,f15.10)')j,eig(i)%eval(j)
		!write(13,'(a,3x,10000f8.3)')'dn:',eig(i)%evec(1:ntotdn,j)
		write(13,'(a,3x,10000f8.3)')'up:',
     .  eig(i)%evec(ntotdn+1:ntotg,1)*eig(i)%evec(ntotdn+1:ntotg,j)
	end do
	write(13,*)'==========================='

	close(13)



	! angles?

	do j=1,nev1
		x1 = dsqrt(dabs(sum(eig(i)%evec(1:ntotdn,1)*
     .                  eig(i)%evec(1:ntotdn,1) )));
		x2 = dsqrt(dabs(sum(eig(i)%evec(1:ntotdn,j)*
     .                  eig(i)%evec(1:ntotdn,j) )));

		if(dabs(x1*x2) .gt. 1.0d-7) then
		amp(j,1) = sum(eig(i)%evec(1:ntotdn,1)*
     .             eig(i)%evec(1:ntotdn,j) )/(x1*x2)
		else
		amp(j,1) = 0.0d0
		endif

		x1 = dsqrt(sum(eig(i)%evec(ntotdn+1:ntotg,1)*
     .                  eig(i)%evec(ntotdn+1:ntotg,1) ));
		x2 = dsqrt(sum(eig(i)%evec(ntotdn+1:ntotg,j)*
     .                  eig(i)%evec(ntotdn+1:ntotg,j) ));

		if(dabs(x1*x2) .gt. 1.0d-7) then
		amp(j,2) = sum(eig(i)%evec(ntotdn+1:ntotg,1)*
     .             eig(i)%evec(ntotdn+1:ntotg,j) )/(x1*x2)
		else
		amp(j,2) = 0.0d0
		endif
    
	end do

	x1 = 3.141592653589793d0; 
	open(13,file='theta.dat',action='write',position='append')
	write(13,'(a,3x,1000f10.3)')'dn:', (dacos(amp(j,1))/x1, j=1,nev1)
	write(13,'(a,3x,1000f10.3)')'up:', (dacos(amp(j,2))/x1, j=1,nev1)
	write(13,*)'==========================='
	close(13)

	open(13,file='cos-theta.dat',action='write',position='append')
	write(13,'(a,3x,1000f10.3)')'dn:', (amp(j,1), j=1,nev1)
	write(13,'(a,3x,1000f10.3)')'up:', (amp(j,2), j=1,nev1)
	write(13,*)'==========================='
	close(13)


	open(13,file='thetaDUP.dat',action='write',position='append')
	write(13,'(10f10.5)') ((dacos(dabs(amp(j,k)))*180/x1, k=1,2),
     .       j=2,nev1)
	close(13)

	endif ! 1==0

	return
	end 	subroutine matelem
	!=============================================================	







	!.....................................................
	! optical absorption from a polariton condensate
	!.....................................................
	subroutine absorption(i) !(i,ijob,n, nsym,m,m1max, mv, chi, newm)
	implicit none
	integer, intent(in) :: i
	ddiagOK = .true.
c	if(m==0) then
c	 call setm0variables(i,n,mv,newm) ! set init state manually.
c	 call seteig0(i)
c	else
	 ! calc H at N_ex=m to get its LP_0 eigenstate
	 if(newm) then
	  call HamParts(nsym,m,mv)
		! make up block for the additional site
		!if(mode>1) call Ham1Parts(nsym,m,mv)
		!newm = .false. because we are going to calc m+1 next anyway.
	 endif
	 call MakeHhtc(nsym, ijob, mode) ! ijob===> sets wr,delta,lambda,wv values
	 call diagonalise(i)
c	endif
	! save LP_0 for time evolution
	call seteig0(i)

	m1max = min(m+1, nsym); 
	! calc H at N_ex=m+1 for time evolution
	call HamParts(nsym,m+1,mv)
	ddiagOK = .false.
	call MakeHhtc(nsym, ijob, mode)

	!write(6,*)'main: tcorr(..,ntot), ntot = ',ntot

	! calc \xi(0) = a^+ LP_0 state for time evolution
	! time evolve using H at N_ex=m+1
	! ...??? ntotg = ntot; ! set to avoid
	! Hf%ntot is current hilbert space dimension.
	call tcorr(dt,w1,w2,nt,nw,i,101, Hf%ntot, m) ! Hf = Hhtc in mode 1
	return
	end subroutine absorption
	!.....................................................


	!.....................................................
	! optical emission from a polariton condensate
	!.....................................................
	subroutine emission(i) !(i,ijob,n, nsym,m,m1max, mv, chi, newm)
	implicit none
	integer, intent(in) :: i
	ddiagOK = .true.
	if(m==0) then
	 write(6,'(a)') "Error: N_ex > 0  for Emission task 102..."
	 stop
	endif
	! calc H at N_ex=m to get its LP_0 eigenstate
	if(newm) then
		call HamParts(nsym,m,mv)
		! make up block for the additional site
		!if(mode>1) call Ham1Parts(nsym,m,mv)
		!newm = .false. because we are going to calc m+1 next anyway.
	endif
	call MakeHhtc(nsym, ijob, mode) ! ijob===> sets wr,delta,lambda,wv values
	call diagonalise(i)

	! save LP_0 for time evolution
	call seteig0(i)

	ddiagOK = .false.
	m1max = min(m-1, nsym);
	! calc H at N_ex=m-1 for time evolution
	call HamParts(nsym,m-1,mv)
	call MakeHhtc(nsym, ijob, mode)

	! calc \xi(0) = a^ LP_0 state for time evolution
	! time evolve using H at N_ex=m-1
	call tcorr(dt,w1,w2,nt,nw,i, 102, Hf%ntot, m)
	return
	end subroutine emission
	!.....................................................
	! density-density response function of the polariton condensate
	!.....................................................
	subroutine densityResponse(task)
	implicit none
	integer, intent(in) :: task
	ddiagOK = .false.;
	! parts of hamiltonian
	if(newm) then
		call HamParts(nsym,m,mv)
		! make up block for the additional site
		call Ham1Parts(nsym,m,mv)
		newm = .false.
	endif
	call MakeHhtc(nsym, ijob, mode) ! ijob===> sets wr,delta,lambda,wv values
	call MakeH1(nsym, ijob) ! uses Hhtc 
	! create full hamiltonin
	call glueHblocksD(nsym, ijob, mv) ! uses H1 and Hhtc
	! diagonalise
	call diagonalise(i)
	write(*,*)'main: density-density response calculaton...'
	call tcorr(dt,w1,w2,nt,nw,i,task, Hf%ntot, m)
	return
	end 	subroutine densityResponse
	!.....................................................
	! density-density response function of the polariton condensate
	!.....................................................
	subroutine hoppingResponse(task)
	implicit none
	integer, intent(in) :: task
	ddiagOK = .false.
	! parts of hamiltonian
	if(newm) then
		call HamParts(nsym,m,mv)
		! make up block for the additional site
		call Ham1Parts(nsym,m,mv)
		newm = .false.
	endif
	call MakeHhtc(nsym, ijob, mode) ! ijob===> sets wr,delta,lambda,wv values
	call MakeH1(nsym, ijob) ! uses Hhtc 
	! create full hamiltonin
	call glueHblocksD(nsym, ijob, mv) ! uses H1 and Hhtc
	! diagonalise
	call diagonalise(i)

	! Hopping response calculaton...
	write(*,*)'main: Hopping response calculaton...'
	! time evolve etc using correlationd module...
	call tcorrd(dt,w1,w2,nt,nw,i, task)

	return
	end 	subroutine hoppingResponse
	!.....................................................
	! writes output: node=0 combines all output files
	!.....................................................
	subroutine writeout(task)
	implicit none
	integer, intent(in) :: task
	select case(task)
	case(101) ! light absorption by the condensate
	 call rwallnodesx('temp-t','absorption-t',nt)
	 call rwallnodesx('temp-w','absorption-w',nw)
	case(102) ! PL/light emission from the condensate
	 call rwallnodesx('temp-t','emission-t',nt)
	 call rwallnodesx('temp-w','emission-w',nw)
	case(103,104)
	 if(task == 103) then
		call rwallnodesx('temp-t','hopping-up-t',nt)
		call rwallnodesx('temp-w','hopping-up-w',nw)
	 elseif(task==104) then
		call rwallnodesx('temp-t','hopping-dn-t',nt)
		call rwallnodesx('temp-w','hopping-dn-w',nw)
	 endif
	case(105,106)
	 if(task == 105) then
		call rwallnodesx('temp-t','density-up-t',nt)
		call rwallnodesx('temp-w','density-up-w',nw)
	 elseif(task == 106) then
		call rwallnodesx('temp-t','density-dn-t',nt)
		call rwallnodesx('temp-w','density-dn-w',nw)
	 endif
	case(301,302)
	 call rwallnodes('dmup',nact)
	 call rwallnodes('dmdn',nact)
	end select
	return
	end subroutine writeout
	!.....................................................
	! response functions calculations	
	!.....................................................
	subroutine response(task)
	implicit none
	integer, intent(in) :: task
	select case(task)
	case(101) ! light absorption by the condensate
	 call absorption(i) !(i,ijob, n, nsym,m,m1max, mv, newm,chi)
	case(102) ! PL/light emission from the condensate
	 call emission(i)
	case(103,104)
	 call hoppingresponse(task)
	case(105,106)
	 call densityresponse(task)
	end select
	return
	end subroutine response

	!.....................................................
	! matrix elements for optical absorption from a polariton condensate
	!.....................................................
	subroutine absmatelem(i) !(i,ijob,n, nsym,m,m1max, mv, chi, newm)
	implicit none
	integer, intent(in) :: i
	integer :: n1,n2, nev1, i1,i2,p
	double precision, dimension(nev) :: val
	double precision :: fac
	
	! calc H at N_ex=m to get its LP_0 eigenstate
	if(newm) then
		call HamParts(nsym,m,mv)
		! make up block for the additional site
		!if(mode>1) call Ham1Parts(nsym,m,mv)
		!newm = .false. because we are going to calc m+1 next anyway.
	endif
	call MakeHhtc(nsym, ijob, mode) ! ijob===> sets wr,delta,lambda,wv values
	call diagonalise(i)
	! only LP0 state needed, so can set nev=1 above....? and reset back again??

	call seteig0(i)

	!write(6,*) '===========>> '
	m1max = min(m+1, nsym); ! used in various ham routines...
	! calc H at N_ex=m+1 and final states
	call HamParts(nsym,m+1,mv)
	call MakeHhtc(nsym, ijob, mode)
	call diagonalise(i)

	! calc diploe matrix elements and write output

	val = 0.0d0;
	ntot = eig0%ntot;
	nev1 = min(nev,ntot);

	i1 = 0;
	do p=0,min(m,nsym); ! nsym=nact in mode=1
	 fac = dsqrt((m - p + 1)*1.0d0); ! sqrt(n_photon + 1)
	 i2 = i1 + basis%sec(p)%ntot*basis%sec(nsym-p)%ntot;
	 do j=1,nev1 ! final states... in Nex+1 space...
	  val(j)=val(j)+ fac * sum(eig0%evec(i1+1:i2,1) * 
     .                      eig(i)%evec(i1+1:i2,j) );
	 end do
	 i1 = i2;
	enddo

	open(12,file='absorption.dat',action='write',position='append')
		! output file data order: n, m, A^0, A^1, ...., A^nev1, Einit^LP_0, E^0, E^1,...,E^nev1
		write(12,'(2i10,2x,100000f25.15)') n, m, val(1:nev1),
     .      eig0%eval(1), eig(i)%eval(1:nev1)
	close(12)

	!write(6,*)'fix sqrt(n_photon) issue above.... '

	return
	end subroutine absmatelem
	!.....................................................


	!.....................................................
	! matrix elements for optical emission from a polariton condensate
	!.....................................................
	subroutine plmatelem(i) !(i,ijob,n, nsym,m,m1max, mv, chi, newm)
	implicit none
	integer, intent(in) :: i
	integer :: n1,n2, nev1, i1,i2,p
	double precision, dimension(nev) :: val
	double precision :: fac
	type(Eigensystems) :: eig0

	! calc H at N_ex=m to get its LP_0 eigenstate
	if(newm) then
		call HamParts(nsym,m,mv)
		! make up block for the additional site
		!if(mode>1) call Ham1Parts(nsym,m,mv)
		!newm = .false. because we are going to calc m+1 next anyway.
	endif
	call MakeHhtc(nsym, ijob, mode) ! ijob===> sets wr,delta,lambda,wv values
	call diagonalise(i)

	!save LP_0 at N_ex=m to a new variable eig0
	if(allocated(eig0%evec)) then
		deallocate(eig0%evec,eig0%eval)
	endif
	n1 = eig(i)%n1;
	n2 = eig(i)%n2;
	eig0%ntot = n1;
	!eig0%n1 = n1; 
	!eig0%n2 = n2;
	allocate(eig0%evec(n1,n2), eig0%eval(n2))
	eig0%evec = eig(i)%evec
	eig0%eval = eig(i)%eval

	!write(6,*) '===========>> '
	m1max = min(m-1,nsym); ! or m-1??? ! used in various ham routines...

	! calc H at N_ex=m+1 and final states
	call HamParts(nsym,m-1,mv)
	call MakeHhtc(nsym, ijob, mode)
	call diagonalise(i)


	! calc diploe matrix elements and write output

	val = 0.0d0;
	ntot = eig(i)%ntot;
	nev1 = min(nev,ntot);

	
c	do j=1,nev1 ! final states... in Nex+1 space...
c	 val(j)=val(j)+sum(eig0%evec(1:ntot,1) * 
c     .                      eig(i)%evec(1:ntot,j) );
c	end do

	!write(6,'(a,10000f9.4)')'Nex: LP_0 : ',eig0%evec(:,1)
	!do j=1,nev
	!	write(6,'(a,10000f7.4)')'Nex-1: LP_j : ',eig(i)%evec(:,j)
	!end do

	i1 = 0;
	do p=0,min(m-1,nsym); ! nsym=nact in mode=1
	 fac = dsqrt((m - p)*1.0d0); ! sqrt(n_photon)
	 i2 = i1 + basis%sec(p)%ntot*basis%sec(nsym-p)%ntot;
	 do j=1,nev1 ! final states... in Nex+1 space...
	  val(j)=val(j)+ fac * sum(eig0%evec(i1+1:i2,1) * 
     .                      eig(i)%evec(i1+1:i2,j) );
	 end do
	 i1 = i2;
	enddo



	open(12,file='emission.dat',action='write',position='append')
		! output file data order: n, m, A^0, A^1, ...., A^nev1, Einit^LP_0, E^0, E^1,...,E^nev1
		write(12,'(2i10,2x,100000f25.15)') n, m, val(1:nev1),
     .      eig0%eval(1), eig(i)%eval(1:nev1)
	close(12)


	return
	end subroutine plmatelem
	!.....................................................
	!.....................................................
	! matrix elements for hopmatelem in a polariton condensate, without vibraation case only.
	!.....................................................
	subroutine hopmatelem(i) !(i,ijob,n, nsym,m,m1max, mv, chi, newm)
	implicit none
	integer, intent(in) :: i
	integer :: n1,n2, nev1
	double precision, dimension(nev) :: val


	! parts of hamiltonian
	if(newm) then
		call HamParts(nsym,m,mv)
		! make up block for the additional site
		call Ham1Parts(nsym,m,mv)
		newm = .false.
	endif
	call MakeHhtc(nsym, ijob, mode) ! ijob===> sets wr,delta,lambda,wv values
	call MakeH1(nsym, ijob) ! uses Hhtc 
	! create full hamiltonin
	call glueHblocksD(nsym, ijob, mv) ! uses H1 and Hhtc
	! diagonalise
	call diagonalise(i)

	! hopping matrix elements, both up/down channels.
	call matelem(i,nact,m)

	return
	end subroutine hopmatelem
	!.....................................................

	subroutine getcdms(i, ijob, mode, newm, ij1)
	implicit none
	integer, intent(in) :: i, ijob, mode
	logical, intent(inout) :: newm
	integer, intent(inout) :: ij1

	! calculate Hamiltonian
		if(newm) then
			call HamParts(nsym,m,mv)
			! make up block for the additional site
			if(mode==2) call Ham1Parts(nsym,m,mv)
			newm = .false.
		endif
	! things need to be done for every job
	call MakeHhtc(nsym, ijob, mode) ! ijob===> sets wr,delta,lambda,wv values
	if(mode==2) then
		call MakeH1(nsym, ijob) ! uses Hhtc 
		! create full hamiltonin
		call glueHblocksD(nsym, ijob, mv) ! uses H1 and Hhtc
	endif
	! diagonalise
	call diagonalise(i)

	if(mode==1) then
	 ! calc dms to free mem or wait for more jobs?
	 call checkifgoforcdms(nj,i,goforcdms)
	 if(goforcdms) then
		njl = i-ij1
		 call cdms(ij1,njl,n,m,mv,nev)
		! reset variables for next iteration
		goforcdms = .false.;
		ij1 = i;
	endif
	elseif(mode==2) then
	 !call cdms(i,norig+1,m); ! norig because here n might be very small compared to norig
	 !call cdms1(i, n,m,mv,nev,ij1,njl) ! calc and write dms at each job iteration.
	 call cdms1(i, nact,m,mv,nev) ! calc and write dms at each job iteration.
	 ij1 = i; ! not needed?
	endif

	return
	end 	subroutine getcdms
	!.....................................................

	subroutine setm0variables(i,n,mv,newm)
	implicit none
	integer, intent(in) :: i,n,mv
	logical, intent(inout) :: newm
	
	write(*,*) 'main: m=0 .... ' 
	! .or. param(ijob-1)%m==0 testing absorption at m=0 using chi2
			! allocate space for basis, and calc maps, etc.
	!---------- called once for a given n,m,mv [all m?]---------
	call getmap(n,mv) !,ntot) ! ntot output
	call PermSymBasis(n,mv)

	! set eig etc manually...
	ntot = 1;
	ntotb= 0;
	ntotdn=ntot*(mv+1);
	ntotup=0;
	ntotg=ntotdn;
	if (allocated(eig(i)%eval)) deallocate(eig(i)%eval)
	if (allocated(eig(i)%evec)) deallocate(eig(i)%evec)
	allocate(eig(i)%eval(1))
	allocate(eig(i)%evec(ntotg,1))
	eig(i)%eval(1) = 0.0d0;
	eig(i)%evec(:,1) = 0.0d0;
	eig(i)%evec(1,1) = 1.0d0;

	newm = .true.; ! ? 
	return
	end subroutine setm0variables


	!.....................................................
	subroutine seteig0(i)
	implicit none
	integer, intent(in) :: i
	integer :: n1, n2
	
	!save LP_0 at N_ex=m to a new variable eig0
	if(allocated(eig0%evec)) then
		deallocate(eig0%evec,eig0%eval)
	endif
	n1 = eig(i)%n1;
	n2 = eig(i)%n2;
	eig0%ntot = n1;
	eig0%n1 = n1; 
	eig0%n2 = n2;
	allocate(eig0%evec(n1,n2), eig0%eval(n2))
	eig0%evec = eig(i)%evec
	eig0%eval = eig(i)%eval
	return
	end subroutine seteig0

	!.....................................................


	!.....................................................


	!.....................................................


	!.....................................................





	end program chi2
