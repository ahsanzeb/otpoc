

! Copyright (C) 2018 M. Ahsan Zeb
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

	module modmain
	implicit none

	! to set kind of integers in ksub
	! kind 1 (8bits) can go upto 127, which should be enough 
	! for us if nsites <= 127.
	!	kind 2 (16 bits) max is 32,767
	! selected_int_kind(R) ==> smallest kind, to 10^R exclusive 
	integer, parameter :: isk=selected_int_kind(2);

	! tasks(maxtasks)
	!integer, parameter :: maxtasks = 100; ! some large numeber
	!integer, dimension(maxtasks) :: tasks

	integer :: task, mode ! 1=all N sites sym; 2= N-1 sym & 1 site explicit; 3= N-1 sym & 1 & 1_D.
	integer :: n, m, mv
	double precision :: wr, delta, lambda, wv ! set in main? fro ith job
	integer :: ntot ! size of Hhtc, old hamiltonian for N permut syn sites.
	integer :: ntotb ! size of first min(m-1,n) blocks of N permut sym sites.
	integer :: ntotdn, ntotup, ntotg ! sizes of bigger blocks of new full hamiltonian including one additional site explicitly described apart from the rest of N permut sym ones. 
	!ntotdn: when this new site is down
	!nntotup: when this new site is down
	!ntotg: grand total, size of FULL Hamiltonian
	integer :: ntotgg ! for Hhop = Hf x Hdbl, Hdbl is doubly occupied site...

	integer :: nact, nsym

	! time evolution to get correlation
	logical :: td, fft, melem, hopresp, fixrhoex, amelem
	integer :: chi
	double precision :: dt,w1,w2, kappa2
	integer :: nt,nw, prntstep

	! participation ratios
	logical :: pratios
	! Nex=1: large N? use only low energy subspace [use a small N~5-10 to build permutation sym basis that make up the low energy space for larger N]
	logical :: LargeN
	integer :: ndummy
	logical :: lcdms ! conditional density matrices ?
	
	logical :: debug

	logical :: ddiagOK ! direct diagonalisation OK? calc and saves wavefunctions... 
	integer :: node, num_procs, njobs ! number of processors, number of jobs
	integer :: mynode
	!integer(kind=MPI_OFFSET_KIND), dimension(:,:), allocatable :: iodisp 

	integer :: memlim ! memory limit per node, in bytes

	logical :: writewfs ! write wavefunctions
	logical :: writecdms	! calc and write cond density matrices
	logical :: writebmap ! wirte basis and map files and stop
	logical :: onlyenergy ! write only energy eigenvalues 
	
	type :: njblocksforcmds
		integer :: ij1, nj
	end type njblocksforcmds
	type(njblocksforcmds), dimension(:), allocatable :: njpntr


	
	! to manage independent parallel jobs on num_procs nodes.
	type :: cluster
		integer :: njobs, i1, i2
	end type cluster
	type(cluster), dimension(:), allocatable :: jobs

	type :: parameters
		integer :: m
		double precision :: wr, del, lam, wv, lamd
		integer :: ntot, ntotdn, ntotb, ntotup, ntotg
	end type parameters
	type(parameters), dimension(:), allocatable :: param

	integer :: nm, nwr, ndel, nlam, nwv, nex
	integer, dimension(:), allocatable :: ms
	double precision, dimension(:), allocatable::
     .             wrs, dels, lams, wvs
	double precision :: lamd, lamd2wv    



	type :: darkdata
		integer :: n
		integer, dimension(:), allocatable :: i
		double precision, dimension(:), allocatable :: c
	end type darkdata
!	type :: darkdata0
!		integer :: n
!		integer, dimension(n) :: i
!		double precision, dimension(n) :: c
!	end type darkdata0
	type(darkdata), dimension(:), allocatable :: evecd
	integer :: ndark
	integer, dimension(:), allocatable ::darkind
	integer, dimension(:), allocatable :: dcoefind
	integer :: ncavity, nexciton


	type :: BSectors
		integer :: ntot ! will see later if these ntot are needed?
		integer(kind=isk), allocatable :: sets(:,:)
		double precision, allocatable :: Nv(:), P(:), r(:,:) ! actually integers, but define real to avoid upper limit.; 
		!r for ratios:
		!# avoids large factorials
		!# sec(p)%r[m*,j_{p-1}] = sqrt( P_{p-1}[j]/P_p[k] ); where k=map[m*,j];

		integer, allocatable :: f(:,:) ! to store freq of occupations for all perm symm basis; used in Hb calc.	
	end type BSectors

	type :: BasisSet
		integer :: ntot
		!logical :: xst
		! number of active sites
		!integer :: n
		integer :: maxk ! max k of k-sub it has.
		integer, allocatable :: pntr(:)
		type(BSectors), allocatable :: sec(:)
	end type BasisSet

	type :: Eigensystems
		integer :: ntot, n1,n2 ! size of hilbert space, =nrows, ncols=nevecs
		double precision, allocatable :: eval(:)
		double precision, allocatable :: evec(:,:)
	end type Eigensystems

	type :: HSectors
		integer :: ntot, nnu, nmu, nnz
		integer, allocatable :: row(:)
		integer, allocatable :: col(:)
		!integer, allocatable :: rowpntr(:)
		double precision, allocatable :: vdat(:)
		!double precision :: hg ! fixed number, pth-sec: Sqrt[(m-p)(p+1)(N-p)]
	end type HSectors


	type :: Ham
	!Hg%ntot,Hg%nnz, Hg%coodat, Hg%coo1, Hg%coo2, Hg%dat, Hg%col, Hg%rowpntr
		logical :: xst,dense
		integer :: n,m,m1
		integer :: nev, ncv
		integer :: ntot ! HilbertSpace dimension
		integer :: nnz  ! no of non-zero elements
		integer :: srptr ! size of row pointers, not eq to ntot
		! to hold data before CSR conversion
		integer, allocatable :: coo1(:), coo2(:)
		double precision, allocatable :: coodat(:)
		! final CSR format 
		integer, allocatable :: col(:)
		integer, allocatable :: rowpntr(:)
		double precision, allocatable :: dat(:)
		! direct diagon, dense form
		double precision, allocatable :: h(:,:)
		!integer, allocatable :: row(:)
		 ! pointers to ksub sectors in rowpntr
		integer, allocatable :: spntr(:)
		! to hold data for largest m1 possible for reuse for smaller m1 values
		type(HSectors), allocatable :: sec(:)
	end type Ham


	type(BasisSet):: basis
	! change name of HilbertSpace to something like eigensystem ??
	!	13 hamiltonians and eigensystems
	type(Eigensystems), dimension(:), allocatable :: eig

	type(Eigensystems) :: eig0

!	type(Ham), allocatable:: Hg
	type(Ham) :: Hhtc, Hg, Hb ! off-diagonal terms
	!diagonal terms, simple forms, maybe taken as simple arrays... 
	double precision, dimension(:), allocatable :: Hd, Hv, Hdv
	double precision, dimension(:), allocatable :: Vx

	double precision, dimension(:,:,:,:), allocatable :: dmup, dmdn


	type(Ham) :: Hf, H1, H1b, H1g, Hgup, Hbup
	double precision, dimension(:), allocatable::H1d,H1v,H1dv,Hvup

	type(Ham) :: Hhop
	double precision, dimension(:), allocatable:: Hvd, eigd


	type :: brightdarkstatesNex1
		! brightdarkstatesNex1 defined to be used in case Nex = 1
		! total number of exciton block states coupled to a given photon block state
		integer :: ntot
		! absolute indices of these ntot states: nphoton + *
		integer, dimension(:), allocatable :: indx ! 1:ntot significant
		! gs:couplings of these ntot states to the photon state
		! b: bright state; b(i)=gs(i)  so one array is enough
		double precision, dimension(:), allocatable :: gs !, b ! 1:ntot significant
		double precision :: geff  ! (sum_i gi^2)^1/2	 = g.sqrt(N) always!
		! ds: non-orthogoanl dark states constructed by hand
		! ods: orthogonalised dark states from ds
		double precision, dimension(:,:), allocatable :: ds !, ods ! 1:ntot,1:ntot-1 significant
	end type brightdarkstatesNex1

	type(brightdarkstatesNex1), allocatable:: bds(:) ! array elements - for each photon state


	double precision, dimension(:,:),allocatable :: Ubdd, HbU,UdHbU ! bright-dark basis transformation matrix, dense
	double precision, dimension(:,:),allocatable :: Hdark
	type(Ham) :: Ubd 


	integer :: m1max
	integer, dimension(:), allocatable :: origin
	integer :: nev ! number of eigenvectors to compute

	
	! indexes: map(0:mv, 1,ntotmaxsites), keep in mind when allocate it
	integer, dimension(:,:), allocatable :: map ! stores map of p-1 sites to p sites.

	! above this size, iterative arpack solver will be used
	!	to diagonalise Hamiltonians
	integer:: nmaxddiag ! max size for direct diagonalisation
	integer :: diagmaxitr! key to set it in input: DirectSolverSize
	
	contains
	!===================================================================
	integer (kind=16) function factorial(p) 
	! length of block that would be added
	! to the map if we add a site to p-1 sites
	implicit none
	integer, intent(in):: p
	integer :: i
	factorial = 1;
	do i=2,p
		factorial = factorial*i
	end do
	return
	end function factorial
	!===================================================================
	integer (kind=16) function binomial(n,m) 
	! length of block that would be added
	! to the map if we add a site to p-1 sites
	implicit none
	integer, intent(in):: n,m
	integer :: i, smaller, greater

	if(n<m) then
		binomial = 0;
	elseif(n==m) then
		binomial = 1;		
	else
		!binomial = factorial(n)/(factorial(m)*factorial(n-m))	
		! remove the largest factorial
		if(m .ge. n-m) then
			greater = m
			smaller = n-m
		else
			greater = n-m
			smaller = m
		endif

		binomial = 1
		do i=greater+1,n
				binomial = binomial*i
		end do
		binomial = binomial/factorial(smaller)
			
	endif
	return
	end function binomial
!--------------- timer -----------------
	double precision function clock()
		real*4:: etime, tm(2)
		clock = etime( tm )
	return
	end function clock
!---------------------------------------



	Subroutine timestamp(node)
!--------------------------------------------------------------------
! 	prints welcome message with date and time
!--------------------------------------------------------------------

	implicit none
	integer, intent(in):: node
	integer	:: dt(8), x
	character*10 :: bb(3)
	double precision :: dummy
	logical, save :: start = .true.;
	integer, save	:: dti(8)
!--------------------------------------------------------------------

	call date_and_time(bb(1), bb(2), bb(3), dt)
	! start random sequence with seed = milliseconds in current time
	call srand(dt(8)) !rand(dt(8))
	if(node .ne. 0 ) return


	if (start) then
		dti = dt;
		start = .false.;
		!call srand(0)
		!write(6,*) 'random = ', rand(0)		
		write(6,*)
		!write(6,*)"**************************"//
    ! .		"******************************"
		write(6,'(/,a)')"			TransOC Started! "
		!write(6,'(/,a)')" TransOC: Transport in Organic Cavities"
		!write(6,'(/,a)')
    ! .        '	        . . . Running v 0.0 . . . '
		write(6,'(/,a,i2,a,i2,2a,i2,a,i2.2,a,i4,a)')
     . 	"	   ",
     .                   dt(5)," Baj kar ",dt(6), " min",
     .               " aor tarekh hai ",dt(3),"/",dt(2),"/",dt(1)
		write(6,*)"         -----------------"//
     .		"----------------------------"
	else
		! assume prog will run on the same day
		x= dt(3)*86400+dt(5)*3600 + dt(6)*60 + dt(7) ! sec
		x = x - (dti(3)*86400+dti(5)*3600 + dti(6)*60 + dti(7));
		dt(3) = x/86400; ! days
		x = x - dt(3)*86400;
		dt(5) = x/3600; ! hours
		x = x - dt(5)*3600;
		dt(6) = x/60; ! min
		x = x - dt(6)*60;
		dt(7) = x; ! sec

		write(6,
     . '("Time taken: ",i2,":",i2.2,":",i2.2,":",i2.2,"(d:h:m:s)")')
     .   dt(3),dt(5:7)
	endif

	return
	END Subroutine timestamp
!========================================================	
	subroutine input(node)
	implicit none
	integer, intent(in):: node
	! local
	character(40) :: block
	integer :: iostat
	logical :: file_exists
	logical :: parameters
	double precision :: memory
	integer :: i, nevdm, nexdm, anev
	logical :: fixrhoexdm
	
	debug = .false.
	! set default
	nmaxddiag = 20
	diagmaxitr = 500
	n = 2
	m = 1
	mv = 5
	parameters = .false.
	memory = 1.00 ; ! 1 GB per node
	nev = 1;
	writewfs = .false.
	writecdms = .true.
	writebmap = .false.
	onlyenergy = .false.
	td= .true.; melem= .false.; fixrhoex=.false.;
	amelem= .false.;
	lcdms = .false.; fixrhoexdm = .false.; nexdm=1;
	chi= 1;
	hopresp = .false.; lamd = 1.0d0;
	lamd2wv = 0.0d0; 
	fft=.true.
	nt = 1000;
	nw=500;
	w1=-0.25d0;
	w2=+2.5d0;
	dt = 0.01;
	kappa2=0.05;
	prntstep = 1000;
	pratios = .false.
	LargeN = .false.
	ndummy = -1;
	anev = 1;
	mode = -1;
	task = -1;
	ddiagOK = .true.
	!--------------------------!
	!     read from input.in   !
	!--------------------------!
	INQUIRE(FILE="input.in", EXIST=file_exists)
	if (.not. file_exists)	then
		write(6,*) 'Error: input.in does not exist!'
		stop
	endif

	open(50,file='input.in',action='READ',
     . 					status='OLD',form='FORMATTED',iostat=iostat)
	if (iostat.ne.0) then
		write(*,*)
		write(*,'("Error(readinput): error opening input.in")')
		write(*,*)
		stop
	end if

10		continue

	read(50,*,end=30) block
	! check for a comment
	if ((scan(trim(block),'!').eq.1).or.
     .  (scan(trim(block),'#').eq.1)) goto 10

	select case(trim(block))


	! task= 100,200 and 300 series for 'ResponseFunction', 
	! 'MatrixElements', and 'DensityMatrices' calculations....

	! task: 
	! Response functions: 100 series
	! 101 abs, 102 emission
	! 103 hopping up, 104 hopping down
	! 105 density up, 106 density down
	! Matrix Elements: 200 series
	! 201 abs, 202 emission
	! 203 hopping without vibration
	! Density Matrices: 300 series
	! 301 conditional vibrational dms like objects in symmetric space, mode=1
	! 302 conditional vibrational dms like objects in full space needed for hopping, mode=2
	! 310 photon reduced dm in the condensate state

	case('task') 
	 read(50,*,err=20) task

	case('nstates') 
	 read(50,*,err=20) nev

	case('TimeEvolutionParam')
		read(50,*,err=20) kappa2, dt, nt, nw, w1,w2, fft
		kappa2 = kappa2/2.0d0; ! just set it here....
	case('TimeEvolutionPrintStep')
		read(50,*,err=20) prntstep

	case('FixedRhoex') ! used with MatrixElements/ConditionalDensityMatrices
		read(50,*,err=20) fixrhoex, nex
	 
	case('N')
		read(50,*,err=20) nact !n_original
		if(nact .le. 1) then
			stop 'Error(input): Set N to at least two sites.... '
		endif
	case('M')
		read(50,*,err=20) mv
	case('debug')
		read(50,*,err=20) debug
	case('parameters')
		parameters = .true.
		read(50,*,err=20) nm, nwr, ndel, nlam, nwv
		allocate(ms(nm))
		allocate(wrs(nwr))
		allocate(dels(ndel))
		allocate(lams(nlam))
		allocate(wvs(nwv))
		read(50,*,err=20) ms
		read(50,*,err=20) wrs
		read(50,*,err=20) dels
		read(50,*,err=20) lams
		read(50,*,err=20) wvs
	case('lambdaD')
		read(50,*,err=20) lamd
	case('LargeN')
		read(50,*,err=20) LargeN, n
		if(LargeN .and. n .le. 1) then
			stop 'Error(input): Set Ndummy to at least two sites.... '
		endif
	case('DirectSolverSize')
		read(50,*,err=20) nmaxddiag
	case('IterSolverMaxIter')
		read(50,*,err=20) diagmaxitr
	case('')
		goto 10
	case default
	write(*,*)
	write(*,'("Error(readinput): invalid block: ",A)')trim(block)
	write(*,*)
	stop
	end select
	goto 10
20		continue
	write(*,*)
	write(*,'("Error(readinput): error reading from input.in")')
	write(*,'("Problem occurred in ''",A,"'' block")') trim(block)
	write(*,'("Check input convention in manual")')
	write(*,*)
	stop
30			continue
	close(50)


	! FixedRhoex case overrides the ms(:) array in parameters block
	if(fixrhoex .and. task > 199) then ! disable fixrhoex for response calculations ??
	 nm= 1;
	 deallocate(ms)
	 allocate(ms(nm))
	 ms(1) = nex;
	 write(*,*)'WARNING: FixedRhoex ===> setting Nex = ',nex 
	endif

	! hopping matrix elements, only without vibrations.
	! overrides M: mv; and N_lambda, lambda: nlam, lams(:) (parameters block)
	if(task == 203) then
	 mv = 0;
	 nlam=1;
	 deallocate(lams)
	 allocate(lams(nlam))
	 lams = 0.0d0;
	 if(.not. fixrhoex) then ! override the ms(:) array in parameters block
	  nm= 2*(nact+1);
	   deallocate(ms)
	   allocate(ms(nm))
	   ms = (/ (i,i=1,nm) /);
	   write(*,*)'WARNING:======>>> setting Nex = 1:N ::1' 
		endif
	endif


	! set mode
	! 1: all N symmetric
	! 2: N-1 sym + 1
	! 3: N-1 sym + 1 + 1_D
	! mode 2,3: n/nsym used for perm sym sites; one site added later for which full up/dn basis states are included, so total becomes given N.

	select case(task)
	 case(101,102,201,202,301)
	 	mode = 1;
	 	nsym = nact;
	 case(105,106,203,302)
	 	mode = 2;
	 	nsym = nact - 1; ! make room for the site explicitly described.
	  !n = n - 1; ! safety to ensure n <= nact
	 case(103,104)
	 	mode = 3;
	 	nsym = nact - 1; ! make room for the site explicitly described.
	  !n = n - 1; ! safety to ensure n <= nact
	  hopresp = .true.;
	 case default
	  write(*,*)
	  write(*,'("Error: task ",A10, " not recognised!")')task
	  write(*,*)
	  stop
	end select

	! response calculations need sparse Hf, so never used dense... even if we can in cases where the state is to be calculated for m+1/m-1 (abs/emission)...
	if(task>100 .and. task < 199) then
	 	nmaxddiag = 0; ! to always make sparse Hf used in td.
	 	nev = 1; ! we only need LP state of a given Nex for td.
	endif



	! set n [ndummy] 
	if(.not. LargeN) then
		n = nsym
	elseif(n > nsym) then
			n= nsym
			write(6,*)'input: given n_dummy > n_sym! reset...'
	endif
	! set ndummy global variable so that it can be used in Hamiltonian.f
	! routines there have local variable named n [which is actually given to them = nsym] so we have to use another global variable... 
	ndummy=n;

	! parameters block must be given by user
	if(.not. parameters) then
	write(*,*)
		write(*,'("Error(readinput): parameters not given in input.in")')	
		stop
	endif

	return


99		write(*,*)
	write(*,'("Error(readinput): invalid block: ",A)')trim(block)
	write(*,*)
	stop		

	
	end subroutine input
!==================================================
	end module

