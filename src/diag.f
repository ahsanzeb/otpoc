
	module diag
	implicit none

	public :: diagonalise
	private	:: iterdiag, ddiag
	private	:: matvec

	contains
!------------------------------------------
!	call to diagonalise() will diagonalise
! all newly calculated Hamiltonians
	subroutine diagonalise(ijob)
	use modmain, only: map, Hf,eig,nev
	implicit none
	integer, intent(in) :: ijob
	! local
	integer:: i,ntot,ncv,j

	!Hf%nev = 1
	!Hf%ncv = min(10,nmaxddiag)

	ntot = Hf%ntot;
	! calculate eigenpairs
	if(Hf%dense) then ! direct diagonalisation
		! find all eigenpairs
		eig%ntot = ntot ! dim of final hilbert space	
		eig%n1 = ntot ! dim of hilbert space
		eig%n2 = ntot; ! number of vectors
		!nev = min(nev,ntot);
		write(*,*) ' ddiag.... '
		! upper triangular H stored in evec, & eval allocated
		call ddiag(ntot,Hf%h,ijob,min(nev,ntot))
		!write(*,*) ' eval = ', eig%eval

	else
		! iterdiag will find the storage format of H
		!	and use appropriate matvec routines
		! it will set the global variables 
		!	eig%evec and eig%eval
		!nev = 6; !Hf%nev; Number of eigenpair to be computed
		ncv = nev + 10; !Hf%ncv; a number bigger than nev for iterative diagonalisation blocks, in lobpcg?
		!write(*,*) "=====> ntot, nev,ncv =",ntot,nev,ncv					
		call iterdiag(ntot, nev, ncv,ijob)
		!write(*,*) "diag: iter done.... "
	endif

	return
	end	subroutine diagonalise
!------------------------------------------
	subroutine iterdiag(n,nev,ncv,ijob)
	use modmain, only: diagmaxitr, eig, Hf
	implicit none
	integer, intent(in) :: n ! n x n dimension of Hg
	integer, intent(in) :: nev ! number of eigenpairs to compute
	integer, intent(in) :: ncv ! leading dimensions for all arrays
	integer, intent(in) :: ijob
	integer :: maxn, maxnev, maxncv,ldv
	! arrays
	double precision,allocatable, dimension(:,:) :: v
	double precision,allocatable, dimension(:) ::  workl
	double precision,allocatable, dimension(:) :: workd
	double precision,allocatable, dimension(:,:) :: d
	double precision,allocatable, dimension(:) :: resid

	logical, allocatable, dimension(:) :: selectt
	integer :: iparam(11), ipntr(11)

	! scalars
	character :: bmat*1, which*2
	integer :: ido, lworkl, info, ierr, j
	integer :: nconv, maxitr, mode, ishfts
	double precision :: sigma, zero = 0.0d+0
	double precision :: tol = 1.0d-8

c	dsaupd documentation:
c  NCV     Integer.  (INPUT)
c          Number of columns of the matrix V (less than or equal to N).
c          This will indicate how many Lanczos vectors are generated 
c          at each iteration.  After the startup phase in which NEV 
c          Lanczos vectors are generated, the algorithm generates 
c          NCV-NEV Lanczos vectors at each subsequent update iteration.
c          Most of the cost in generating each Lanczos vector is in the 
c          matrix-vector product OP*x. (See remark 4 there).

	!n = nx*nx
	maxn = n;
	!ncv = 2*nev; ! ncv <= n,  

	maxnev = nev; 
	maxncv= ncv;
	ldv = maxn;

	allocate(v(ldv,maxncv))
	allocate(workl(maxncv*(maxncv+8)))
	allocate(workd(3*maxn))
	allocate(d(maxncv,2))
	allocate(resid(maxn))
	allocate(selectt(maxncv))

	! requirements:
	!	NEV + 1 <= NCV	 <= n 

	which = 'SA'; ! lowest eigenvalues?

	lworkl = ncv*(ncv+8) ! at least ncv*(ncv+8)
	!tol = zero 
	info = 0
	ido = 0

	iparam(1) = 1 ! ishfts 
	iparam(3) = diagmaxitr ! max iterations 
	iparam(7) = 1 ! mode, standard eigenvalue problem

	!write(*,*)" ---- started ---- n,nev,ncv ",n,nev,ncv

		!-------------------------------------
		!Beginning of reverse communication
		!call dsaupd ( ido, bmat, n, which, nev, tol, resid, 
		!&	 		 	 	  ncv, v, ldv, iparam, ipntr, workd, workl,
		!&	 		 	 	  lworkl, info )
10		call dsaupd ( ido, 'I', n, which, nev, tol, resid, 
     .   ncv, v, ldv, iparam, ipntr, workd, workl,
     .   lworkl, info )

		if (ido .eq. -1 .or. ido .eq. 1) then
				call matvec(n, workd(ipntr(1):ipntr(1)+n-1),
     .                    workd(ipntr(2):ipntr(2)+n-1))
			go to 10
		end if 
		!End of Reverse communication
		!-------------------------------------

	!write(*,*)" =======>>>>> info = ",info


	if ( info .lt. 0 ) then ! see documentation of DSAUPD
		write(*,*)'Error(diag): Error with _saupd, info = ', info
		write(*,*)'Error(diag): Check documentation in _saupd '
		stop
	else
		!call dseupd (.true., 'All', selectt, d, v, ldv, sigma, 
		! bmat, n, which, nev, tol, resid, ncv, v, ldv, 
		! iparam, ipntr, workd, workl, lworkl, ierr )
		call dseupd (.true., 'All', selectt, d, v, ldv, sigma, 
     .   'I', n, which, nev, tol, resid, ncv, v, ldv, 
     .   iparam, ipntr, workd, workl, lworkl, ierr )
c        %-------------------  dseupd  -----------------%
c        | Eigenvalues are returned in the first column |
c        | of the two dimensional array D and the       |
c        | corresponding eigenvectors are returned in   |
c        | the first NEV columns of the two dimensional |
c        | array V if requested.  Otherwise, an         |
c        | orthogonal basis for the invariant subspace  |
c        | corresponding to the eigenvalues in D is     |
c        | returned in V.                               |
c        %----------------------------------------------%
	endif

	if ( ierr .lt. 0 ) then ! see documentation of DSEUPD
		write(*,*)'Error(diag): Error with _seupd, info = ', ierr
		write(*,*)'Error(diag): Check documentation in _seupd '
		stop
	endif


	!write(*,*)" diag : -------------- done--------------------"
	! sure Hf%col, dat, rowpntr is allocated!
	! to lower memory usage deallocate Hf before space for eig is allocated
	!deallocate(Hf%col)
	!deallocate(Hf%dat)
	!deallocate(Hf%rowpntr)


	! set eig
	eig(ijob)%ntot = n
	eig(ijob)%n1 = n
	eig(ijob)%n2 = nev
	if(allocated(eig(ijob)%evec))deallocate(eig(ijob)%evec)
	if(allocated(eig(ijob)%eval))deallocate(eig(ijob)%eval)
	allocate(eig(ijob)%evec(n,nev))
	allocate(eig(ijob)%eval(nev))
	eig(ijob)%evec = v(1:n,1:nev) !reshape(v(1:n,1:nev),(/n*nev/))
	eig(ijob)%eval = d(1:nev,1) ! dim of d = (maxncv,2), so first nev values?

	!write(*,*) "diag: n1,n2 [n,nev]= ",n,nev
	!write(*,*) "diag: eval [d] = ",d(1:nev,1)
	!open(17,file='eigval-ijob',action='write',position='append')
	!write(17,*) ijob, d(1:nev,1)
	!close(17)

	!open(17,file='wv-ijob',action='write',position='append')
	!write(17,*) (d(j+1,1)-d(j,1),j=1,nev-1)
	!close(17)

	deallocate(v,workl,workd,d,resid,selectt)

9000		continue

	return
	end subroutine iterdiag

	!-----------------------------------
	subroutine matvec(n,x,y)
	! multiplies Hamiltonian Hf with input vector x
	! gives y as output: y = H.x
	! Hamiltonian is stored in CSR format
	use modmain, only: Hf ! Hf dim = n x n
	implicit none
	integer, intent(in) :: n
	double precision, dimension(n), intent(in):: x
	double precision, dimension(n), intent(out):: y
	integer :: i,j,k

	y = 0.0d0
 	do i=1,Hf%ntot
		do j=Hf%rowpntr(i), Hf%rowpntr(i+1)-1
			k = Hf%col(j); ! icol
			y(i) = y(i) + x(k)*Hf%dat(j)
			y(k) = y(k) + x(i)*Hf%dat(j) ! Hf upper or lower triangular only. 
			! Hb might have lower triangular elements?! 
			! no issues with having mixed upper/lower elements as long as they appear once in either upper or lower.
			! diagonal terms in Hf are halved to balance double counting.
		enddo
 	enddo
 	
	return
	end subroutine matvec
!-----------------------------------
	subroutine ddiag(ntot,H,ijob,nev)
	use modmain, only: eig, clock
	implicit none
	integer, intent(in) :: ntot,ijob,nev
	!double precision, dimension(ntot,ntot), intent(inout):: H	
	!double precision, dimension(ntot), intent(inout):: W
	double precision, dimension(ntot,ntot), intent(in):: H	
	! local
	double precision, dimension(ntot) :: W
	integer :: lwork != 3*ntot ! LWORK >= max(1,3*N-1)
	double precision, dimension(3*ntot):: WORK
	integer :: info,j
	double precision ::starttime, endtime,starttime1, endtime1
	EXTERNAL	 DSYEV ! LAPACK/BLAS

	!write(*,'(/,a)') " diagonal: diagonalising H ... "
	starttime1 = clock()

	lwork = 3*ntot ! LWORK >= max(1,3*N-1)

	if(1==0) then
	write(*,*)'************************************'
	do j=1,ntot
		write(*,'(10000E8.1)') H(j,:)
	end do
	write(*,*)'************************************'
	endif
	
	
	CALL DSYEV('Vectors','Upper',ntot,H,ntot,W,WORK,LWORK,INFO)
	! Check for convergence.
	IF( INFO.GT.0 ) THEN
		write(*,'(/,a)')
     .   'diagonal: The algorithm failed to compute eigenvalues.'
		STOP
	END IF

	! ******** set global variable eig *********
	! set eig, lowest polariton only
	eig(ijob)%ntot = ntot
	eig(ijob)%n2 = nev
	if(allocated(eig(ijob)%evec))deallocate(eig(ijob)%evec)
	if(allocated(eig(ijob)%eval))deallocate(eig(ijob)%eval)
	allocate(eig(ijob)%evec(ntot,nev))
	allocate(eig(ijob)%eval(nev))
	eig(ijob)%evec = H(1:ntot,1:nev) !reshape(H(1:ntot,1:nev),(/ntot*nev/)) ! columns are eigenvectors??? or rows???
	eig(ijob)%eval = W(1:nev)
	
	endtime1 = clock()
	!write(*,'(/,a,E8.2)')" diagonal: total time taken	= ",
  !   . 			endtime1-starttime1

	!open(17,file='eigval-ijob',action='write',position='append')
	!write(17,*) ijob, w(1:min(ntot,5))
	!close(17)

	!open(17,file='wv-ijob',action='write',position='append')
	!write(17,*) (w(j+1)-w(j),j=1,min(ntot,5)-1)
	!close(17)
	
	return
	END subroutine ddiag
!---------------------------------------	
	end 	module diag
