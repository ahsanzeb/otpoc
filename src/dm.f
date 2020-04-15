

	module dmat
	use modmain
	implicit none

	public:: cdms, rwallnodes, cdms1, dmphot
	private:: cdmup, cdmdn 

	contains

	!--------------------------------------------------------------------
	! p-th block cdm for spin up
	!--------------------------------------------------------------------
	subroutine cdmup(ij1,nj,n,p,mv,nev,dm)
	implicit none
	integer, intent(in) :: ij1,nj, n,p, mv,nev ! ij1: now nj can be smaller than actual total nj for a node
	double precision,dimension(mv+1,mv+1,nj,nev),intent(out) :: dm
	integer :: l1,l2, lout, nnp,ntotp
	!double precision, dimension(:,:), allocatable :: v
	integer :: ii,i1,j,j1,i,jj,i11,j11,ij, ijob, m1ij,k,xj
	double precision :: fac, Pii, Pi1,Pj1,fac2
	type(Eigensystems), dimension(:), allocatable :: v
	
	l1 = basis%pntr(p); l2 = basis%pntr(p+1)

	!write(*,*)'p, l1, l2 = ',p, l1,l2
	!lout = basis%sec(p)%ntot;
	nnp = basis%sec(n-p)%ntot;
	!allocate(v(lout))
	! tracemu(l,lout,nnp,v) ! traces over mu type
	!v = tracemu(l2-l1, lout, nnp, eig%evec(l1:l2-1,1)) 
	!dm = trmol(n,p,mv,lout,v) ! function to trace over p-1 other molecules' vibrational states

	ntotp = min(basis%sec(p-1)%ntot,basis%sec(ndummy-1)%ntot); ! states for which p block states are available to map onto... 

	dm =0.0d0

	allocate(v(nj))
	k = 0
	do i=ij1+1,ij1+nj
		k = k + 1
		ijob = jobs(node)%i1 + i -1
		m1ij = min(param(ijob)%m, n)
		if(p .le. m1ij) then
			allocate(v(k)%evec(l2-l1,nev))
			v(k)%evec = eig(i)%evec(l1:l2-1,:)
		endif
	enddo
	!allocate(v(l2-l1,nj)) ! lout*nnp = l2-l1
	!do ij=1,nj
	!	v(:,ij) = eig(ij)%evec(l1:l2-1)
	!end do
	!================================================
	do ii=1, ntotp !basis%sec(p-1)%ntot
		!Pii = basis%sec(p-1)%P(ii)
		do i=0,mv
			i1 = map(i,ii) ! {nu|mu}_(k-1) + i ==> i1-th basis state of {nu|mu}_k
			!Pi1 = basis%sec(p)%P(i1)
			!write(*,*)'i1 = ', i1
			!fac = Pii/dsqrt(Pi1); ! fac due to vib permutations
			fac = basis%sec(p-1)%r(i,ii) ! r(i,ii) = dsqrt(Pii/Pi1)
			i1 = (i1-1)*nnp;
		
			do j=i,mv ! upper triangular only
				j1 = map(j,ii) 
				!Pj1 = basis%sec(p)%P(j1)
				!fac2 = fac/dsqrt(Pj1); ! fac due to vib permutations
				fac2 = fac*basis%sec(p-1)%r(j,ii) ! r(j,ii) = dsqrt(Pii/Pj1)
				!write(*,*)'j1 = ', j1
				j1 = (j1-1)*nnp;
				k = 0
				do ij=ij1+1,ij1+nj,1
					k=k+1
					ijob = jobs(node)%i1 + ij - 1;
					if(p .le. min(param(ijob)%m,n)) then
						do xj=1,nev
						dm(i+1,j+1,k,xj) = dm(i+1,j+1,k,xj)
!     .         + fac2*sum(v(i1+1:i1+nnp,ij)*v(j1+1:j1+nnp,ij))
     . + fac2*sum(v(k)%evec(i1+1:i1+nnp,1)*v(k)%evec(j1+1:j1+nnp,xj))
     				enddo ! xj final state from LP
					endif
				end do
			enddo
		enddo
	end do ! ii
	! from perm sym of exciton-photon states
	dm = dm*p*1.0d0/n ! p out of n are up
	!================================================
	
	deallocate(v)
	return
	end subroutine cdmup
	!--------------------------------------------------------------------
	! p-th block cdm for spin down
	!--------------------------------------------------------------------
	subroutine cdmdn(ij1,nj,n,p,mv,nev,dm)
	implicit none
	integer, intent(in) :: ij1,nj, n, p, mv,nev ! ij1: now nj can be smaller than actual total nj for a node
	double precision,dimension(mv+1,mv+1,nj,nev),intent(out) :: dm
	integer :: l1,l2, lout,nnp,ntotp
	!double precision, dimension(:,:), allocatable :: v
	integer :: ii,i1,j,j1,i,jj,i11,j11,ij, ijob, m1ij,k
	double precision :: fac, Pii, Pi1,Pj1,fac2
	type(Eigensystems), dimension(:), allocatable :: v

	dm =0.0d0
	
	l1 = basis%pntr(p); l2 = basis%pntr(p+1)
	!lout = basis%sec(n-p)%ntot;
	nnp = basis%sec(n-p)%ntot; ! mu_(n-p)
	!allocate(v(lout))
	! tracenu(l,lout,np,v) ! traces over nu type
	!v = tracenu(l2-l1, lout, np, eig%evec(l1:l2-1,1)) 
	!dm = trmol(n,n-p,mv,lout,v) ! function to trace over k-1=n-p-1 other molecules' vibrational states

	!write(*,*)'l1,l2 = ',l1,l2
	ntotp = min(basis%sec(n-p-1)%ntot,basis%sec(ndummy-1)%ntot); ! states for which n-p block states are available to map onto... 



	dm =0.0d0
	!allocate(v(l2-l1,nj))
	!do ij=1,nj
	!	v(:,ij) = eig(ij)%evec(l1:l2-1)
	!end do

	allocate(v(nj))
	k = 0
	do i=ij1+1,ij1+nj
		k = k + 1
		ijob = jobs(node)%i1 + i -1
		m1ij = min(param(ijob)%m, n)
		if(p .le. m1ij) then
			allocate(v(k)%evec(l2-l1,nev))
			v(k)%evec = eig(i)%evec(l1:l2-1,:)
		endif
	enddo



	!allocate(v(l2-l1,nj)
	!do ij=1,nj
	!	v(:,ij) = eig(ij)%evec(l1:l2-1)
	!end do

	!write(*,*)'size(basis%sec(p)%r) =',size(basis%sec(p)%r) 
	!================================================
	do ii=1, ntotp !basis%sec(n-p-1)%ntot
		!Pii = basis%sec(n-p-1)%P(ii)
		do i=0,mv
			i1 = map(i,ii) ! {mu}_(n-p-1) + i ==> i1-th basis state of {mu}_(n-p)
			!Pi1 = basis%sec(n-p)%P(i1)
			!fac = Pii/dsqrt(Pi1); ! fac due to vib permutations
			fac = basis%sec(n-p-1)%r(i,ii) ! r(i,ii) = dsqrt(Pii/Pi1)
			do j=i,mv ! upper triangular only
				j1 = map(j,ii) 
				!Pj1 = basis%sec(n-p)%P(j1)
				!fac2 = fac/dsqrt(Pj1); ! fac due to vib permutations
				fac2 = fac * basis%sec(n-p-1)%r(j,ii); ! r(j,ii) = dsqrt(Pii/Pj1)
				do jj=1,basis%sec(p)%ntot ! trace over nu_{p} 
					i11 = (jj-1)*nnp + i1;
					j11 = (jj-1)*nnp + j1;
					k = 0
					do ij = ij1 + 1,ij1+nj,1
						ijob = jobs(node)%i1 + ij - 1;
						k = k + 1
						if(p .le. min(param(ijob)%m,n)) then
							dm(i+1,j+1,k,:) = dm(i+1,j+1,k,:)
!     .         						+ v(i11,ij)*v(j11,ij)*fac2;
     .         + v(k)%evec(i11,1)*v(k)%evec(j11,:)*fac2;
						endif
					end do
					!dm(i+1,j+1) = dm(i+1,j+1) + v(i11)*v(j11)*fac2;
				end do ! jj
			enddo
		enddo
	end do ! ii
	! from perm sym of exciton-photon states
	dm = dm*(n-p)*1.0d0/n ! n-p out of n down
	!================================================

	deallocate(v)
	return
	end subroutine cdmdn
	!--------------------------------------------------------------------
	! full cdms combine all blocks
	!--------------------------------------------------------------------
	subroutine cdms(ij1, nj,n,m,mv,nev)
	implicit none
	integer, intent(in) :: ij1, nj,n,m,mv,nev
	double precision,dimension(mv+1,mv+1,nj,nev):: dmup,dmdn
	double precision,dimension(mv+1,mv+1,nj,nev) :: dm
	integer :: p,m1,ij,m1max, k
	double precision :: tr

	! total number of jobs on this node
	!nj = jobs(node)%njobs
	!write(*,*)'jobs(node)%i1, ij1 = ', jobs(node)%i1, ij1
	m1max = 0;
	do ij= jobs(node)%i1 + ij1, jobs(node)%i1 + ij1 + nj-1     !jobs(node)%i1, jobs(node)%i2
		m1max = max(m1max, min(param(ij)%m,n) )
	end do
	
	dmup = 0.0d0; dmdn =0.0d0;
	do p=0,m1max
		if(p==0) then ! all spin down
			call cdmdn(ij1,nj,n,p,mv,nev,dm)
			dmdn = dmdn + dm;
		elseif(p==n) then ! all spin up
			call cdmup(ij1,nj,n,p,mv,nev,dm)
			dmup = dmup + dm;
		else ! p up, n-p down
			call cdmdn(ij1,nj,n,p,mv,nev,dm)
			dmdn = dmdn + dm;
			call cdmup(ij1,nj,n,p,mv,nev,dm)
			dmup = dmup + dm;
		endif	
	end do

	! test trace
	if(debug) then
	dm = dmup + dmdn
	do ij=1,nj
		tr = 0.0d0;
		do p=1,mv+1
			tr = tr + dm(p,p,ij,1) ! tr = 1 only for LP to LP scattering, last index=1
			!write(*,*)'p, tr = ', p, tr
		end do
		write(*,*)'node, ij, trace[dmup+dmdn] = ',node, ij, tr
	end do
	endif


	! write output files at each node
	call writeatnode(ij1,nj,dmup,'dmup')
	call writeatnode(ij1,nj,dmdn,'dmdn')
	! deallocate eig%evec for the done jobs to free memory
	!write(*,*)'main: freeing memory ... '
	do k = ij1+1, ij1+nj
	 if(allocated(eig(k)%evec)) deallocate(eig(k)%evec)
	 !? if(allocated(eig(k)%eval)) deallocate(eig(k)%eval)
	enddo



	!call writecdms(mv,nj,nev, dmup, dmdn)


	
	return
	end subroutine cdms
	!--------------------------------------------------------------------

















	subroutine dmphot(ij1,nj,n,p,mv,nev,dm)
	implicit none
	integer, intent(in) :: ij1,nj, n,p, mv,nev ! ij1: now nj can be smaller than actual total nj for a node
	double precision,dimension(mv+1,mv+1,nj,nev),intent(out) :: dm
	integer :: l1,l2, lout, nnp,ntotp
	!double precision, dimension(:,:), allocatable :: v
	integer :: ii,i1,j,j1,i,jj,i11,j11,ij, ijob, m1ij,k,xj
	double precision :: fac, Pii, Pi1,Pj1,fac2
	type(Eigensystems), dimension(:), allocatable :: v
	
	l1 = basis%pntr(p); l2 = basis%pntr(p+1)
	nnp = basis%sec(n-p)%ntot;
	ntotp = min(basis%sec(p-1)%ntot,basis%sec(ndummy-1)%ntot); ! states for which p block states are available to map onto... 

	dm =0.0d0

	allocate(v(nj))
	k = 0
	do i=ij1+1,ij1+nj
		k = k + 1
		ijob = jobs(node)%i1 + i -1
		m1ij = min(param(ijob)%m, n)
		if(p .le. m1ij) then
			allocate(v(k)%evec(l2-l1,nev))
			v(k)%evec = eig(i)%evec(l1:l2-1,:)
		endif
	enddo

	if(p .le. min(param(ijob)%m,n)) then
	do xj=1,nev
		dm(i+1,j+1,k,xj) = dm(i+1,j+1,k,xj)
     . + fac2*sum(v(k)%evec(i1+1:i1+nnp,1)*v(k)%evec(j1+1:j1+nnp,xj))
	enddo ! xj final state from LP
	endif
	
	deallocate(v)
	return
	end subroutine dmphot



































	! mode 2: N-1 sym + 1 site case.
	! conditional reduced vibrational density matrix of the site explicitly described... 
	! so will be diff from the sym site's cdms for the excited states due to different excited states from the non-sym sector.
	subroutine cdms1(i,n,m,mv, nev)
	implicit none
	integer, intent(in) :: i,n,m,mv, nev !, ij1, nj
	double precision, dimension(mv+1,mv+1,nev) :: dmup, dmdn
	integer :: i1,i2,j1,j2, k1,k2,j, k
	double precision, dimension(nev) :: trdn, trup

	!write(*,*)'ntot, ntotb, ntotg, ntotdn = ',
  !   .               ntot, ntotb, ntotg, ntotdn
	
	trup(:) = 0.0d0; trdn(:) = 0.0d0;
	dmdn = 0.0d0;
	dmup = 0.0d0;

	do j=1,nev
		! dn
		do k1=0,mv
			i1 = k1*ntot; 
			i2 = i1 + ntot;
			do k2=0,mv
				j1 = k2*ntot; 
				j2 = j1 + ntot;
	
				dmdn(k1+1,k2+1,j) = sum(eig(i)%evec(i1+1:i2,1) * 
     .                  eig(i)%evec(j1+1:j2,j) )
			end do
		end do

		! up 
		do k1=0,mv
			i1 = ntotdn + k1*ntotb; 
			i2 = i1 + ntotb;
			do k2=0,mv
				j1 = ntotdn + k2*ntotb; 
				j2 = j1 + ntotb;
				dmup(k1+1,k2+1,j) = sum(eig(i)%evec(i1+1:i2,1) * 
     .                  eig(i)%evec(j1+1:j2,j) )
			end do
		end do

	end do ! j=1,nev

	! write output files at each node
	call writeatnode(i,1,dmup,'dmup')
	call writeatnode(i,1,dmdn,'dmdn')
	


	! trace(cdms) to another file useful with FixedRhoex 
	open(12,file='matelemv.dat',action='write',position='append')
		write(12,'(2i10,3x,100000f25.15)') n,m, 
     .      trdn(1)**2 , ( trup(j)**2, j=1,nev ) 
	! LP -> LP (j=1) has diff traces for the two channels... 
	!other transitions have the same amp, just opposite signs...
	close(12)


	! write eigvalues... 
	open(12,file='evals.dat',action='write',position='append')
		write(12,'(10000f25.15)') eig(i)%eval(1:nev)
	close(12)

	! deallocate eig%evec for the done jobs to free memory
	!write(*,*)'main: freeing memory ... '
	 if(allocated(eig(i)%evec)) deallocate(eig(i)%evec)
	 !? if(allocated(eig(i)%eval)) deallocate(eig(i)%eval)

	return
	end 	subroutine cdms1
	!=============================================================	



	subroutine writecdmsfiles(mv,nj,nev, dmup, dmdn)
	implicit none
	integer, intent(in) :: mv, nj,nev
	double precision,dimension(mv+1,mv+1,nj,nev),
     .                         intent(in)::dmup,dmdn
	integer :: i,j, k1, k2


	! write output files
	open(10,file='dmdn.dat',action='write',position='append')
	open(12,file='dmup.dat',action='write',position='append')

	do i=1, nj
	 do j=1, nev
	  ! down
	  do k1=1,mv+1
	   write(10,'(1000f25.15)') (dmdn(k1,k2,i,j), k2=1,mv+1)
	  end do
		! up
	  do k1=1,mv+1
	   write(12,'(1000f25.15)') (dmup(k1,k2,i,j), k2=1,mv+1)
	  end do
	 enddo
	enddo
	
	close(10)
	close(12)

	return
	end subroutine writecdmsfiles




	subroutine writeatnode(ij1,nj,dm,filename)
	implicit none
	integer, intent(in) :: ij1,nj
	double precision, dimension(mv+1,mv+1,nj,nev), intent(in) :: dm
	character(len=*), intent(in) :: filename
	integer, dimension(3) :: dsize
	integer :: thefile
	integer :: i,ij,j,xj
	character :: rank*30, fname*100
	
	write(rank,'(i6.6)') node
	fname = trim(filename)//'-'//trim(rank)
	! write unformatted file
	if(ij1==0) then
		open(1,file=trim(fname), form="formatted", action="write")
	else
		open(1,file=trim(fname), form="formatted", action="write",
     .                                      position="append")
	endif
	do ij=1,nj
		do xj=1,nev
			do i=1,mv+1
				write(1,*) (dm(i,j,ij,xj), j=1,mv+1)
			end do
		enddo ! xj
	end do
	close(1)
	return
	end 	subroutine writeatnode
	!=============================================================
	subroutine rwallnodes(filename,n)
	implicit none
	character(len=*), intent(in) :: filename
	integer, intent(in) :: n
	integer :: i,ij,k,j,shift,xj
	character :: rank*30, fname*100
	!local
	double precision, dimension(mv+1,mv+1,njobs,nev) :: dm
	character :: dir*30
	logical :: ex


	! trim() does not work if you make dir by "write(dir,*) n"
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

	dm =0.0d0

	shift=0;
	do i=0,min(njobs,num_procs)-1
		write(rank,'(i6.6)') i
		fname = trim(filename)//'-'//trim(rank)
		! read unformatted file
		open(1,file=trim(fname), form="formatted", action="read")
		do ij=jobs(i)%i1,jobs(i)%i2!1,jobs(i)%njobs
			do xj=1,nev
				do k=1,mv+1
					read(1,*) (dm(k,j,ij,xj), j=1,mv+1)
				end do
			enddo
			!call checktrace(mv+1,dm(:,:,ij))
		end do
		!close(1)
		close(1, status='delete')
		!shift = jobs(i)%njobs
	end do

	! write output file with data from all nodes
	fname = trim(dir)//'/'//trim(filename)//'.dat'
	! read unformatted file
	open(1,file=trim(fname), form="formatted", action="write",
     .      position="append")
	do ij=1,njobs
		do xj=1,nev
			call symmetrise(mv+1,dm(:,:,ij,xj)) ! complete lower triangular
			do i=1,mv+1
				write(1,*) (dm(i,j,ij,xj), j=1,mv+1)
			end do
		end do
	end do
	close(1)

	return
	end 	subroutine rwallnodes
	!=============================================================

	
	subroutine symmetrise(l,a)
	implicit none
	integer, intent(in) :: l
	double precision, dimension(l,l), intent(inout) :: a
	integer :: i,j
	do i=1,l-1
		do j=i+1,l,1
			a(j,i) = a(i,j)
		end do
	end do
	return
	end 	subroutine symmetrise
	!=============================================================




	end 	module dmat
