
	module hamiltonian1
	use modmain	
	use hamiltonian, only: coocsr, coo2dense
	implicit none

	! H1 is equivalent to Hhtc block but for one less excitation
	! that has been taken up by the newly added site.
	public :: MakeH1, Ham1Parts, glueHblocks, glueHblocksD
	private:: MakeH1g, MakeH1d,MakeH1v, MakeH1b ! in H1's Hilbert space
	private:: MakeHgup, MakeHbup, MakeHvup ! in Full hilbert space
	private:: setsizes
	contains

! ===========================================
	subroutine Ham1Parts(n,m,mv) ! n=norig set in calling prog main.f
	implicit none
	integer, intent(in):: n,m,mv

	! sets sizes of various blocks/subblocks of full hamiltonian
	call setsizes(m,n,mv)
	if(hopresp) then
		write(6,*)'ntot, ntotdn, ntotb, ntotup,  ntotg, ntotgg = ',
     .               ntot, ntotdn, ntotb, ntotup, ntotg, ntotgg
	else
		write(6,*)'ntot, ntotdn, ntotb, ntotup,  ntotg = ',
     .               ntot, ntotdn, ntotb, ntotup, ntotg
	endif

	!write(6,*)'Ham1Parts: ---------- 1'
	! H1's pats in its own space
	call MakeH1g(n,m)

	!write(6,*)'Ham1Parts: ---------- 2'

	call MakeH1d(n,m)
	!write(6,*)'Ham1Parts: ---------- 3'

	call MakeH1v(n,m)
	!write(6,*)'Ham1Parts: ---------- 4'
	if(mv>0) call MakeH1b(n,m,mv)
	!write(6,*)'Ham1Parts: ---------- 5'

	! full hilbert space parts
	call MakeHgup(n, mv)
	if(mv>0) call MakeHbup(mv)
	call MakeHvup(mv)
	

	return
	end subroutine Ham1Parts
! ===========================================
!	subroutine MakeHam(n,m,mv)
!	implicit none
!	integer, intent(in):: n,m,mv
!
!	call MakeHhtc(n, wr,delta,lambda,wv)
!
!	subroutine
!	end 	subroutine MakeHam
! ===========================================


	!----------------------------------------------------------
	subroutine MakeH1g(n,m)
	! constructs sparse Hg for m excitations using Hg%sec
	! saves the Hamiltonians in CSR format
	implicit none
	integer, intent(in) :: n,m
	double precision:: hgxp
	integer :: m1,nnz,ntot,i,p,ind1,ind2,ind
	integer :: nnup, nmunp ! n-nu-p; n-mu-n-minus-p


	m1=min(m-1,n); ! the new site we added takes one excitation, so m-1 here.
	
	! starting ndex for p-th block
	if(allocated(origin))deallocate(origin)
	if(allocated(basis%pntr))deallocate(basis%pntr)
	allocate(origin(0:m1max+1))
	allocate(basis%pntr(0:m1max+1))
	origin(0) = 0;
	do p=1,m1max+1
		nnup = basis%sec(p-1)%ntot;
		nmunp = basis%sec(n-p+1)%ntot;
		origin(p) = origin(p-1) + nnup*nmunp
	end do
	
	basis%pntr = origin + 1 ! info repeat... use pntr instead of origin?!
	!write(*,*) 'origin = ', origin
	!write(*,*) 'basis%pntr = ', basis%pntr
	
	! total nnz
	nnz = sum(Hg%sec(0:m1-1)%nnz)
	H1g%nnz = nnz
	! aux arrays to hold data from all p-blocks in coo format
	if(allocated(H1g%coo1))deallocate(H1g%coo1)
	if(allocated(H1g%coo2))deallocate(H1g%coo2)
	if(allocated(H1g%coodat))deallocate(H1g%coodat)
	allocate(H1g%coo1(nnz))
	allocate(H1g%coo2(nnz))
	allocate(H1g%coodat(nnz))

	!nnz = 0;
	!do p=0,m1-1
	!	nnz += Hg%sec(0:m1-1)%nnz
	!end do ! p
	ntot = sum(Hg%sec(0:m1)%ntot) ! Hg%sec%ntot  set in makeHgSec ??
	H1g%ntot = ntot
	!nrow = Hg%ntot;

	! coo format but with global indexing 
	!shift = 0;
	i = 0;
	do p=0,m1-1
		! one less photon: m-1 photons
		hgxp = dsqrt((m-1-p)*(p+1)*(n-p)*1.0d0); ! exciton-photon matrix element
		do ind=1,Hg%sec(p)%nnz
			!ind1 = Hg%sec(p)%row(ind)
			!ind2 = Hg%sec(p)%col(ind)
			! global position: i,j
			!i = ind + shift;
			i = i + 1;
			H1g%coo1(i) = origin(p) + Hg%sec(p)%row(ind);
			H1g%coo2(i) = origin(p+1) + Hg%sec(p)%col(ind);
			H1g%coodat(i) = hgxp * Hg%sec(p)%vdat(ind); ! full matrix element		
		end do! ind: nnz elements in p-th sector
	end do! p: exciton-photon sectors

	return
	end 	subroutine MakeH1g
	!----------------------------------------------------------

	!----------------------------------------------------------
	subroutine MakeH1d(n,m)
	implicit none
	integer, intent(in) :: n,m
	integer :: i1,m1,ntot,p

	if(allocated(H1d)) deallocate(H1d)
	allocate(H1d(H1g%ntot))
	H1d(:) = Hd(1:H1g%ntot) - 1.0d0; ! one less photon, or one more exciton
	return
	end 	subroutine MakeH1d
	!----------------------------------------------------------
	subroutine MakeH1v(n,m)
	implicit none
	integer, intent(in) :: n,m
	integer :: i1,m1,ntot,p, Nvi, i2,i

	if(allocated(H1v)) deallocate(H1v)
	allocate(H1v(H1g%ntot)) 
	H1v(:) = Hv(1:H1g%ntot) ! same vibrations
	return
	end 	subroutine MakeH1v
	!----------------------------------------------------------



		!----------------------------------------------------------
	subroutine MakeH1b(n,m,mv)
	implicit none
	integer, intent(in) :: n,m,mv
	integer :: m1,nnz

	! Hb already calculated:
	! Hb%nnz: set to only significant size of coo Hb, 
	! Do not reset it to sum of Hb%sec(p)%nnz

	m1=min(m-1,n); ! one less excitation [reduces m1 only if m-1 < n].

	nnz  = Hb%spntr(m1+1) - 1; ! nnz for H1b
	H1b%nnz = nnz;
	if(nnz == 0) then
		!write(6,*)'H1b%nnz = 0; is m=1?' 
		return
	endif 	 
		
	! Hb's array size larger than its significant size nnz.
	!nnz = Hb%nnz; ! signficant size of Hb arrays [already calculated].

	if(allocated(H1b%coo1))deallocate(H1b%coo1)
	if(allocated(H1b%coo2))deallocate(H1b%coo2)
	if(allocated(H1b%coodat))deallocate(H1b%coodat)
	allocate(H1b%coo1(nnz))
	allocate(H1b%coo2(nnz))
	allocate(H1b%coodat(nnz))

	H1b%coo1(:) = Hb%coo1(1:nnz) 
	H1b%coo2(:) = Hb%coo2(1:nnz) 
	H1b%coodat(:) = Hb%coodat(1:nnz) 

	return
	end 	subroutine MakeH1b
	!----------------------------------------------------------

	!----------------------------------------------------------
	! makes H1 block, Hhtc's equivalent with one less excitation when the additional site is up.
	subroutine MakeH1(n,ijob)
	implicit none
	integer, intent(in) :: n,ijob
	!double precision, intent(in) :: wr, delta, lambda, wv

	! local
	double precision :: g, lamwv
	integer :: i,n1,n2,n3,nnz, Hbnnzm

	! parameters for this job
	wr = param(ijob)%wr;
	delta = param(ijob)%del;
	lambda = param(ijob)%lam;
	wv = param(ijob)%wv;

	! hg : upper triangular
	! Hb might have lower triangular elements?! 
	! no issues with having mixed upper/lower elements as long as they appear once in either upper or lower.
	! diagonal: Hdv are halved for matvec()

	!..............................................
	! combine all terms to make full Hamiltnian in coo
	!..............................................

	g = wr/dsqrt(nact*1.0d0);
	lamwv = -lambda*wv;

	nnz = H1g%nnz + H1b%nnz; ! diagonal will be kept separate

	!write(*,*)'Hg%nnz, Hb%nnz, Hg%ntot =',Hg%nnz,Hb%nnz,Hg%ntot
	
	H1%nnz = nnz;
	if(allocated(H1%coo1))deallocate(H1%coo1)
	if(allocated(H1%coo2))deallocate(H1%coo2)
	if(allocated(H1%coodat))deallocate(H1%coodat)
	allocate(H1%coo1(nnz))
	allocate(H1%coo2(nnz))
	allocate(H1%coodat(nnz))

	if(allocated(H1dv)) deallocate(H1dv)
	allocate(H1dv(H1g%ntot))  ! set H1g%ntot first

	H1dv = delta*H1d +  wv*H1v; ! summ arrays, values.

	n1 = 0;
	n3 = H1g%nnz;
	n2 = n3;

	!write(*,*) '-----   Hg%nnz = ',Hg%nnz

	! Hg first
	H1%coo1(n1+1:n2) = H1g%coo1(1:n3)
	H1%coo2(n1+1:n2) = H1g%coo2(1:n3)
	H1%coodat(n1+1:n2) = g * H1g%coodat(1:n3)


	if(H1b%nnz>0) then
		!write(*,*) 'Hg: n1, n2 = ',n1, n2 
		n1 = n2;
		n3 = H1b%nnz; ! Hb%nnz set to spref[= significant no of elements]
								! Do not reset it to sum of Hb%sec(p)%nnz
		n2 = n1 + n3;
		! Hb afterwards
		H1%coo1(n1+1:n2) = H1b%coo1(1:n3)
		H1%coo2(n1+1:n2) = H1b%coo2(1:n3)
		H1%coodat(n1+1:n2) = lamwv * H1b%coodat(1:n3)
	endif
	
	return
	end subroutine MakeH1

	!--------------------------------------------------

	!--------------------------------------------------
	! sets global variables ntot, ntotb, ntotdn, ntotup, ntotg
	subroutine setsizes(m,n,mv)
	implicit none
	integer, intent(in) :: n,m,mv
	integer :: p

	!write(*,*)'basis%sec(p)%ntot=',basis%sec(:)%ntot
	ntot = 0
	do p=0,min(m,n)
		ntot = ntot + basis%sec(p)%ntot * basis%sec(n-p)%ntot
		!write(*,*)'setsizes: p, ntot = ', p, ntot
	end do

	ntotb = 0
	do p=0,min(m-1,n)
		ntotb = ntotb + basis%sec(p)%ntot * basis%sec(n-p)%ntot
	end do

	ntotdn = ntot * (mv+1);
	ntotup = ntotb * (mv+1);
	ntotg = ntotdn + ntotup;

	! in case hopping response is requested
	ntotgg = ntotg * (mv+1)


	return
	end subroutine setsizes

	! coo format: Hhtc%coo1, Hhtc%coo2, Hhtc%coodat
	!--------------------------------------------------
	subroutine MakeHgup(n,mv) ! upper triangular
	implicit none
	integer, intent(in) :: n, mv
	integer :: k, j1,j2,j3,p
	integer, dimension(ntotb) :: iaux

	iaux = (/ (k,k=1,ntotb) /);

	Hgup%nnz = ntotup;
	if(allocated(Hgup%coo1))deallocate(Hgup%coo1)
	if(allocated(Hgup%coo2))deallocate(Hgup%coo2)
	if(allocated(Hgup%coodat))deallocate(Hgup%coodat)

	allocate(Hgup%coo1(ntotup))
	allocate(Hgup%coo2(ntotup))
	allocate(Hgup%coodat(ntotup))
	do k=0,mv
		j1= k*ntotb
		j2= (k+1)*ntotb
		Hgup%coo1(j1+1:j2) = k*ntot + iaux(:); ! k*ntot+1 : k*ntot+ntotb
		Hgup%coo2(j1+1:j2) = ntotdn + k*ntotb + iaux(:);

		! different p sectors have diff photon number... goes in the sqrt.
		do p=0,min(m-1,n)
			j3 = j1 + basis%sec(p)%ntot*basis%sec(n-p)%ntot
			Hgup%coodat(j1+1:j3)= dsqrt((m-p)*1.0d0); ! * g
			j1 = j3;
		end do
		if(j2 .ne. j3) then
			write(*,*)'Hgup: =====> j2 .ne. j3, ERROR?'
			stop
		endif
	end do

	return
	end subroutine 
	!--------------------------------------------------
	subroutine MakeHbup(mv)
	implicit none
	integer, intent(in) :: mv
	integer :: k, j1,j2
	integer, dimension(ntotb) :: iaux

	iaux = (/ (k,k=ntotdn+1,ntotdn+ntotb) /);

	!write(*,*)'===========> iaux='
	!write(*,'(1000i4)') iaux

	
	Hbup%nnz = ntotb * mv; 
	if(allocated(Hbup%coo1))deallocate(Hbup%coo1)
	if(allocated(Hbup%coo2))deallocate(Hbup%coo2)
	if(allocated(Hbup%coodat))deallocate(Hbup%coodat)

	allocate(Hbup%coo1(Hbup%nnz))
	allocate(Hbup%coo2(Hbup%nnz))
	allocate(Hbup%coodat(Hbup%nnz))

	do k=0,mv-1
		j1= k*ntotb
		j2= (k+1)*ntotb
		Hbup%coo1(j1+1:j2) = k*ntotb + iaux(:); ! k*ntot+1 : k*ntot+ntotb
		Hbup%coo2(j1+1:j2) = (k+1)*ntotb + iaux(:);
		Hbup%coodat(j1+1:j2)= dsqrt(1.0d0*(k+1)); ! * lambda*wv
	end do

	!write(*,*)'===========> Hbup%coo1='
	!write(*,'(1000i4)') Hbup%coo1

	return
	end subroutine 
	!--------------------------------------------------
	subroutine MakeHvup(mv) ! a diagonal part for exciton and vibrations on the newly added site. In full hilbert space, just like Hgup and Hbup....
	implicit none
	integer, intent(in) :: mv
	integer :: k, j1,j2

	if(allocated(Hvup)) deallocate(Hvup)
	allocate(Hvup(ntotg)) 

	! dn block
	do k=0,mv
		j1= k*ntot
		j2= (k+1)*ntot
		Hvup(j1+1:j2) = k*1.0d0;
	end do

	! up block
	do k=0,mv
		j1= ntotdn + k*ntotb
		j2= ntotdn + (k+1)*ntotb
		Hvup(j1+1:j2) = k*1.0d0;
	end do

	return
	end subroutine 
	!--------------------------------------------------
	! combine all blocks coo arrays to make Hf%coo1(), Hf%coo2(), Hf%coodat()
	subroutine glueHblocks(n,ijob,mv)
	implicit none
	integer, intent(in) :: n,ijob,mv
	integer :: nnz, nref
	double precision :: g, lamwv
	integer :: i, i1,i2,j1,j2,k, n1,n2

	! parameters for this job
	wr = param(ijob)%wr;
	delta = param(ijob)%del;
	lambda = param(ijob)%lam;
	wv = param(ijob)%wv;
	
	g = wr/dsqrt(nact*1.0d0);
	lamwv = -lambda*wv;

	! size of full hilbert space and Hf
	Hf%ntot = ntotg;

	! calc nnz
	nnz = (Hhtc%nnz + H1%nnz) * (mv+1); ! (mv+1) Hhtc and H1 blocks 
	nnz = nnz + ntotg; ! diagonal terms;
	nnz = nnz + Hgup%nnz; !
	nnz = nnz + Hbup%nnz; ! 
	Hf%nnz = nnz;
	
	if(allocated(Hf%coo1))deallocate(Hf%coo1)
	if(allocated(Hf%coo2))deallocate(Hf%coo2)
	if(allocated(Hf%coodat))deallocate(Hf%coodat)
	allocate(Hf%coo1(nnz))
	allocate(Hf%coo2(nnz))
	allocate(Hf%coodat(nnz))

	! dn block, Hhtc
	nref = 0;
	do k=0,mv
		j1= k*Hhtc%nnz + 1;
		j2= (k+1)*Hhtc%nnz
		Hf%coo1(j1:j2) = k*ntot + Hhtc%coo1(1:Hhtc%nnz)
		Hf%coo2(j1:j2) = k*ntot + Hhtc%coo2(1:Hhtc%nnz)
		Hf%coodat(j1:j2) = Hhtc%coodat(1:Hhtc%nnz);
	end do
	nref = nref + Hhtc%nnz * (mv+1); ! mv+1 blocks with Hhtc%nnz elem

	! up block, H1
	do k=0,mv
		j1= nref + k*H1%nnz + 1;
		j2= nref + (k+1)*H1%nnz
		Hf%coo1(j1:j2) = ntotdn + k*ntotb + H1%coo1(1:H1%nnz)
		Hf%coo2(j1:j2) = ntotdn + k*ntotb + H1%coo2(1:H1%nnz)
		Hf%coodat(j1:j2) = H1%coodat(1:H1%nnz)
	end do
	nref = nref + H1%nnz * (mv+1); ! mv+1 blocks with H1%nnz elem

	! Hgup: upper triangular (or maybe lower, one of the two)
	j1 = nref + 1;
	j2 = nref + Hgup%nnz; ! Hgup%nnz = ntotup
	Hf%coo1(j1:j2) =	 Hgup%coo1(:) ! Hgup already has absolute indices
	Hf%coo2(j1:j2) =	 Hgup%coo2(:) ! Hgup already has absolute indices
	Hf%coodat(j1:j2) = g*Hgup%coodat(:);
	nref = nref + Hgup%nnz;

	if(mv>0) then
		! Hbup: upper triangular
		j1 = nref + 1;
		j2 = nref + Hbup%nnz;
		Hf%coo1(j1:j2) =	 Hbup%coo1(:) ! Hgup already has absolute indices
		Hf%coo2(j1:j2) =	 Hbup%coo2(:) ! Hgup already has absolute indices
		Hf%coodat(j1:j2) = lamwv*Hbup%coodat(:)
		nref = nref + Hbup%nnz;
	endif

	! diagonal parts: half them for matvec routine.
	! Hvup: full space;
	!	Hdv, H1dv: Hhtc and H1 spaces;
	! H1dv includes the energy of one more exciton on the new site. 

	!write(*,*) 'Hdv(1:45)=' 
	!write(*,'(1000f5.2)') Hdv(1:45)
	!write(*,*) 'Hdv(46:66)=' 
	!write(*,'(1000f5.2)')  Hdv(46:66)

	! dn block
	do k=0,mv
		i1 = k*ntot + 1;
		i2 = (k+1)*ntot;
		j1= nref + i1;
		j2= nref + i2;
		Hf%coo1(j1:j2) = (/ (i, i=i1,i2) /);
		Hf%coo2(j1:j2) = Hf%coo1(j1:j2);
		!write(6,*)'i1:i2, size(Hdv) = ',i2-i1+1, size(Hdv)
		Hf%coodat(j1:j2) = (wv*Hvup(i1:i2) + Hdv)*0.5d0; ! half it for matvec
	end do
	nref = nref + ntot * (mv+1) 

	! up block 
	do k=0,mv
		i1 = k*ntotb + 1;
		i2 = (k+1)*ntotb;
		j1= nref + i1;
		j2= nref + i2;
		Hf%coo1(j1:j2) = (/ (i, i=ntotdn+i1,ntotdn+i2) /);
		Hf%coo2(j1:j2) = Hf%coo1(j1:j2);
		Hf%coodat(j1:j2)=(wv*Hvup(ntotdn+i1:ntotdn+i2) +H1dv)*0.5d0;
		! half it for matvec
	end do
	nref = nref + ntotb * (mv+1) 

	!...........................................................
	! convert Hf to CSR... 
	!...........................................................

	!write(*,*)'============ Full hamiltonain =========='
	Hf%dense = .false.
	Hf%ntot = ntotg ! dimension of full hilbert space
	if (Hf%ntot .le. nmaxddiag) then 
		! Hf in dense format and use direct diagonalisation
		Hf%dense = .true.
		n1 =Hf%ntot;
		n2 = Hf%nnz;
		if(allocated(Hf%h))deallocate(Hf%h)
		allocate(Hf%h(n1,n1))
		call coo2dense(n1, n2, Hf%coo1(1:n2),
     .      Hf%coo2(1:n2),Hf%coodat(1:n2), Hf%h)

	if(1==0) then
		do k=1,n1
			write(*,'(10000f4.1)') Hf%h(k,:)
		end do
	endif
     
	else ! iterative diaognalisation
		!..............................................
		! change the format to CSR 
		!..............................................
		if(allocated(Hf%col))deallocate(Hf%col)
		if(allocated(Hf%dat))deallocate(Hf%dat)
		if(allocated(Hf%rowpntr))deallocate(Hf%rowpntr)
		allocate(Hf%col(nnz))
		allocate(Hf%dat(nnz))
		allocate(Hf%rowpntr(Hf%ntot + 1))

		!write(*,*) 'Hf%ntot, Hf%nnz = ',Hf%ntot, Hf%nnz 

		!write(*,'(1000i4)')Hf%coo2
		!write(*,*) any(Hf%coo1 > 261)
		!write(*,*) any(Hf%coo2 > 261)	

			
		call coocsr(Hf%ntot, Hf%nnz, 
     .  Hf%coodat, Hf%coo1, Hf%coo2,  
     .  Hf%dat, Hf%col, Hf%rowpntr)
	endif

	!write(*,*)'============ Full hamiltonain done =========='

	return
	end subroutine glueHblocks
	!--------------------------------------------------

	!--------------------------------------------------
	! combine all blocks coo arrays to make Hf%coo1(), Hf%coo2(), Hf%coodat()
	! 	and Hhop%coo[1/2/dat]
	subroutine glueHblocksD(n,ijob,mv)
	implicit none
	integer, intent(in) :: n,ijob,mv
	integer :: nnz, nref, nrefd, nnz0
	double precision :: g, lamwv,lamwvd
	integer :: i, i1,i2,j1,j2,k, n1,n2, j3

	! parameters for this job
	wr = param(ijob)%wr;
	delta = param(ijob)%del;
	lambda = param(ijob)%lam;
	wv = param(ijob)%wv;
	
	g = wr/dsqrt(nact*1.0d0);
	lamwv = -lambda*wv;
	lamwvd = -param(ijob)%lamd*wv;
	lamd2wv = param(ijob)%lamd**2 * wv;


	! size of full hilbert space and Hf
	Hf%ntot = ntotg;

	! calc nnz
	nnz = (Hhtc%nnz + H1%nnz) * (mv+1); ! (mv+1) Hhtc and H1 blocks 
	nnz = nnz + ntotg; ! diagonal terms;
	nnz = nnz + Hgup%nnz; !
	nnz = nnz + Hbup%nnz; ! 
	Hf%nnz = nnz;
	nnz0 = Hf%nnz - ntotg; ! size of Hf arrays without the diagonal terms... 
	
	if(allocated(Hf%coo1))deallocate(Hf%coo1)
	if(allocated(Hf%coo2))deallocate(Hf%coo2)
	if(allocated(Hf%coodat))deallocate(Hf%coodat)
	allocate(Hf%coo1(nnz))
	allocate(Hf%coo2(nnz))
	allocate(Hf%coodat(nnz))

	! dn block, Hhtc
	nref = 0;
	do k=0,mv
		j1= k*Hhtc%nnz + 1;
		j2= (k+1)*Hhtc%nnz
		Hf%coo1(j1:j2) = k*ntot + Hhtc%coo1(1:Hhtc%nnz)
		Hf%coo2(j1:j2) = k*ntot + Hhtc%coo2(1:Hhtc%nnz)
		Hf%coodat(j1:j2) = Hhtc%coodat(1:Hhtc%nnz);
	end do
	nref = nref + Hhtc%nnz * (mv+1); ! mv+1 blocks with Hhtc%nnz elem

	! up block, H1
	do k=0,mv
		j1= nref + k*H1%nnz + 1;
		j2= nref + (k+1)*H1%nnz
		Hf%coo1(j1:j2) = ntotdn + k*ntotb + H1%coo1(1:H1%nnz)
		Hf%coo2(j1:j2) = ntotdn + k*ntotb + H1%coo2(1:H1%nnz)
		Hf%coodat(j1:j2) = H1%coodat(1:H1%nnz)
	end do
	nref = nref + H1%nnz * (mv+1); ! mv+1 blocks with H1%nnz elem

	! Hgup: upper triangular (or maybe lower, one of the two)
	j1 = nref + 1;
	j2 = nref + Hgup%nnz; ! Hgup%nnz = ntotup
	Hf%coo1(j1:j2) =	 Hgup%coo1(:) ! Hgup already has absolute indices
	Hf%coo2(j1:j2) =	 Hgup%coo2(:) ! Hgup already has absolute indices
	Hf%coodat(j1:j2) = g*Hgup%coodat(:) !* 0.0d0; ! TESTING g=0 here...
	nref = nref + Hgup%nnz;

	if(mv>0) then
		! Hbup: upper triangular
		j1 = nref + 1;
		j2 = nref + Hbup%nnz;
		Hf%coo1(j1:j2) =	 Hbup%coo1(:) ! Hgup already has absolute indices
		Hf%coo2(j1:j2) =	 Hbup%coo2(:) ! Hgup already has absolute indices
		Hf%coodat(j1:j2) = lamwv*Hbup%coodat(:)
		nref = nref + Hbup%nnz;
	endif

	! diagonal parts: half them for matvec routine.
	! Hvup: full space;
	!	Hdv, H1dv: Hhtc and H1 spaces;
	! H1dv includes the energy of one more exciton on the new site. 

	!write(*,*) 'Hdv(1:45)=' 
	!write(*,'(1000f5.2)') Hdv(1:45)
	!write(*,*) 'Hdv(46:66)=' 
	!write(*,'(1000f5.2)')  Hdv(46:66)

	! dn block
	do k=0,mv
		i1 = k*ntot + 1;
		i2 = (k+1)*ntot;
		j1= nref + i1;
		j2= nref + i2;
		Hf%coo1(j1:j2) = (/ (i, i=i1,i2) /);
		Hf%coo2(j1:j2) = Hf%coo1(j1:j2);
		!write(6,*)'i1:i2, size(Hdv) = ',i2-i1+1, size(Hdv)
		Hf%coodat(j1:j2) = (wv*Hvup(i1:i2) + Hdv)*0.5d0; ! half it for matvec
	end do
	nref = nref + ntot * (mv+1) 

	! up block 
	do k=0,mv
		i1 = k*ntotb + 1;
		i2 = (k+1)*ntotb;
		j1= nref + i1;
		j2= nref + i2;
		Hf%coo1(j1:j2) = (/ (i, i=ntotdn+i1,ntotdn+i2) /);
		Hf%coo2(j1:j2) = Hf%coo1(j1:j2);
		Hf%coodat(j1:j2)=(wv*Hvup(ntotdn+i1:ntotdn+i2) +H1dv)*0.5d0;
		! half it for matvec
	end do
	nref = nref + ntotb * (mv+1) 

	! let's use Hf to construct Hhop... 
	! first use non-diagonal terms of Hf only, use its diag terms later
	! along with diag of D site... 
	!call makeHop()


	if(hopresp) then
		Hhop%ntot = ntotgg;
		Hhop%nnz = Hf%nnz * (mv+1) + ntotg * mv ! diag blocks + vib coupling
		
		if(allocated(Hhop%coo1))deallocate(Hhop%coo1)
		if(allocated(Hhop%coo2))deallocate(Hhop%coo2)
		if(allocated(Hhop%coodat))deallocate(Hhop%coodat)
		allocate(Hhop%coo1(Hhop%nnz))
		allocate(Hhop%coo2(Hhop%nnz))
		allocate(Hhop%coodat(Hhop%nnz))


		! D site: vib coupling...
		do k=0,mv-1
			j1= k*ntotg + 1;
			j2= (k+1)*ntotg;
			j3= (k+2)*ntotg;
			Hhop%coo1(j1:j2) = (/ (i, i=j1,j2) /);
			Hhop%coo2(j1:j2) = (/ (i, i=j2+1,j3) /);
			Hhop%coodat(j1:j2) = lamwvd*dsqrt((k+1)*1.0d0);
		end do
		nrefd = ntotg * mv; 

		! Hf blocks (without their diag terms)
		! nnz0: first nnz0 elem of Hf arrays are non-diagonal terms...
		do k=0,mv
			j1= nrefd + k*nnz0 + 1;
			j2= nrefd + (k+1)*nnz0;
			! nnz0 non diag terms in Hf.... already been build into Hf arrays... 
			Hhop%coo1(j1:j2) = k*ntotg + Hf%coo1(1:nnz0) ! k*ntotg = start pos for kth block of Hf in Hhop
			Hhop%coo2(j1:j2) = k*ntotg + Hf%coo2(1:nnz0)
			Hhop%coodat(j1:j2) = Hf%coodat(1:nnz0)
		end do
		nrefd = nrefd + nnz0*(mv+1);
		
		! diagonal term... halved.
		do k=0,mv
			j1= nrefd + k*ntotg + 1;
			j2= nrefd + (k+1)*ntotg;
			! nnz0 non diag terms in Hf.... already been build into Hf arrays... 
			Hhop%coo1(j1:j2) = (/ (i, i=k*ntotg+1, (k+1)*ntotg) /);
			Hhop%coo2(j1:j2) = Hhop%coo1(j1:j2);
			Hhop%coodat(j1:j2) = 0.5d0*wv*k + Hf%coodat(nnz0+1:Hf%nnz) ! Hf diag already halved...
		end do

		! change the format to CSR 
		if(allocated(Hhop%col))deallocate(Hhop%col)
		if(allocated(Hhop%dat))deallocate(Hhop%dat)
		if(allocated(Hhop%rowpntr))deallocate(Hhop%rowpntr)
		allocate(Hhop%col(Hhop%nnz))
		allocate(Hhop%dat(Hhop%nnz))
		allocate(Hhop%rowpntr(Hhop%ntot + 1))
		call coocsr(Hhop%ntot, Hhop%nnz, 
     .  Hhop%coodat, Hhop%coo1, Hhop%coo2,  
     .  Hhop%dat, Hhop%col, Hhop%rowpntr)

		! set eigd
		if(allocated(eigd))deallocate(eigd)
		allocate(eigd(mv+1));
		do k=0,mv
			eigd(k+1) = (-lamd)**k/factorial(k);
		end do
		eigd = eigd * dexp(-0.5d0*lamd**2)

	endif ! hopresp

	!...........................................................
	! convert Hf to CSR... 
	!...........................................................
	!write(*,*)'============ Full hamiltonain =========='
	Hf%dense = .false.
	Hf%ntot = ntotg ! dimension of full hilbert space
	if (Hf%ntot .le. nmaxddiag .and. ddiagOK) then 
		! Hf in dense format and use direct diagonalisation
		Hf%dense = .true.
		n1 =Hf%ntot;
		n2 = Hf%nnz;
		if(allocated(Hf%h))deallocate(Hf%h)
		allocate(Hf%h(n1,n1))
		call coo2dense(n1, n2, Hf%coo1(1:n2),
     .      Hf%coo2(1:n2),Hf%coodat(1:n2), Hf%h)
		!write(*,*) 'Hf%dense = ',Hf%dense

	if(1==0) then
		do k=1,n1
			write(*,'(10000f4.1)') Hf%h(k,:)
		end do
	endif
     
	else ! iterative diaognalisation
		!..............................................
		! change the format to CSR 
		!..............................................
		if(allocated(Hf%col))deallocate(Hf%col)
		if(allocated(Hf%dat))deallocate(Hf%dat)
		if(allocated(Hf%rowpntr))deallocate(Hf%rowpntr)
		allocate(Hf%col(nnz))
		allocate(Hf%dat(nnz))
		allocate(Hf%rowpntr(Hf%ntot + 1))

		!write(*,*) 'Hf%ntot, Hf%nnz = ',Hf%ntot, Hf%nnz 

		!write(*,'(1000i4)')Hf%coo2
		!write(*,*) any(Hf%coo1 > 261)
		!write(*,*) any(Hf%coo2 > 261)	

			
		call coocsr(Hf%ntot, Hf%nnz, 
     .  Hf%coodat, Hf%coo1, Hf%coo2,  
     .  Hf%dat, Hf%col, Hf%rowpntr)
	endif

	!write(*,*)'============ Full hamiltonain done =========='

	return
	end subroutine glueHblocksD
	!--------------------------------------------------

	end !module


	
