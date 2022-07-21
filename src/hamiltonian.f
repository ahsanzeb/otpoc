
	module hamiltonian
	use modmain	
	implicit none
	
	public :: MakeHhtc, HamParts
	private:: MakeHgSec,MakeHg
	private:: MakeHd,MakeHv
	private:: MakeHbSec,MakeHb
	private:: sumdup, rowadddup
	public :: coocsr, coo2dense
	contains

! ===========================================
	subroutine HamParts(n,m,mv) ! arg n is norig here, set in main.f
	implicit none
	integer, intent(in):: n,m,mv

	!write(*,*)'------------------1 '
	call MakeHgSec(n,m,mv)
	!write(*,*)'------------------2 '
	if(.not. writewfs) then
		call MakeHg(n,m)
	else
		call MakeHgAndVx(n,m)
	endif
	!write(*,*)'------------------3 '
	call MakeHd(n,m)
	!write(*,*)'------------------4 '
	call MakeHv(n,m)
	!write(*,*)'------------------5 '
	if(mv > 0) then
		call MakeHbSec(n,m,mv)
		!write(*,*)'------------------6 '
		call MakeHb(n,m,mv)
	else
		Hb%nnz = 0; !use to calc size of sparse Hhtc in MakeHhtc
	endif
	!write(*,*)'------------------7 '
	return
	end subroutine HamParts
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
	subroutine MakeHgSec(n,m,mv)
	implicit none
	integer, intent(in) :: n,m,mv
	integer :: p,nnz,nnp,nnp1,ind,i,j,k,ii,jj,ind1,ind2
	double precision :: pmui,pmuf,pnui,pnuf, ri, rjj
	integer :: ntotp, ntotnp
	
	! basis info: Normalisations of vibrational states
	! define indexes from 0 for sec of Ham and basis so that p sec is p-th elem.
	! similarly, map ind from 0,0
	! vibrational maps that connect perm symm states of p sites with those of p+1 sites... 


	if(allocated(Hg%sec)) deallocate(Hg%sec)
	!allocate(Hg%sec(0:m1max))
	allocate(Hg%sec(0:n))
	!allocate(Hg%nnzs(m1max))	! array for number of nnz in all blocks

	! loop over exciton-photon blocks.
	! separate p=0 case....????
	!write(*,*)'Hg: m1max = ',m1max
	!write(*,*) basis%sec(:)%ntot
	 
	do p=0,m1max-1 ! m1=Min(N,m); m1max= max if calc for multiple m's are desired.
		! nnz = N_{p}*N_{N-p-1}*(M+1) non-zero elements in pth block of Hg: p to p+1 up spins
		!nnz=(basis%sec(p)%ntot)*(basis%sec(n-p-1)%ntot)*(mv+1) 

c		ii = basis%sec(min(m1max-1, ndummy-1))%ntot;
c		ntotp = min(basis%sec(p)%ntot,ii); ! states for which p+1 block states are available for coupling through Hg.... 
c		ntotnp=min(basis%sec(n-p-1)%ntot,ii);! states for which n-p block states are available for coupling through Hg.... 

		ntotp = min(basis%sec(p)%ntot,basis%sec(ndummy-1)%ntot); ! states for which p+1 block states are available for coupling through Hg.... 
		ntotnp=min(basis%sec(n-p-1)%ntot,basis%sec(ndummy-1)%ntot);! states for which n-p block states are available for coupling through Hg.... 


		nnz = ntotp* ntotnp *(mv+1)
		Hg%sec(p)%nnz = nnz
		allocate(Hg%sec(p)%row(nnz))
		allocate(Hg%sec(p)%col(nnz))
		allocate(Hg%sec(p)%vdat(nnz))

		! exciton-photon factor
		!Hg%sec(p)%hgxp = dsqrt((m-p)*(p+1)*(n-p)) ! depends on m, so we can set it later when making full hamiltonian CSR format 
		
		!ntot = basis%sec(p-1)%ntot + 1 ! +1?
		!allocate(Hg%sec(p)%rowpntr(ntot))
		!write(*,*)'Hg: ntotp, ntotnp = ', ntotp, ntotnp
		nnp = basis%sec(n-p)%ntot; ! number of perm symm basis for n-p sites (for n-p initial down spins)
		nnp1 = basis%sec(n-p-1)%ntot;! number of perm symm basis for n-p-1 sites (for n-p-1 final down spins)
		! |i,j> in pth sector to |ii,jj> in p+1 sector; photon emission term.
		! |i,j>: i for p up, j for N-p down;
		ind=0;
		do i=1, ntotp ! basis%sec(n-p)%ntot
		do jj = 1, ntotnp ! basis%sec(n-p-1)%ntot
		do k=0,mv
			ind = ind + 1;
			! maps, order of loops ii,j,k and shape of the maps.... fortran style = row first or col???
			ii = map(k,i);
			j = map(k,jj); 
			ind1 = (i-1)*nnp + j ! initial state index inside p block
			ind2 = (ii-1)*nnp1 + jj ! final state index inside p+1 block
			Hg%sec(p)%row(ind) = ind1 ! Fortran: column first, faster memoery access
			Hg%sec(p)%col(ind) = ind2
			!Pnui=basis%sec(p)%P(i);
			!Pnuf=basis%sec(p+1)%P(ii);
			!Pmui=basis%sec(n-p)%P(j);
			!Pmuf=basis%sec(n-p-1)%P(jj);		
			!Hg%sec(p)%vdat(ind) = dsqrt((Pnui*Pmuf)/(Pmui*Pnuf)) ! smaller in numerator, bigger in denominator
			! we have: 	basis%sec(-)%r(k,i) = dsqrt(Pnui/Pnuf) 	= ri
			! 						basis%sec(-)%r(k,jj) = dsqrt(Pmuf/Pmui) 	= rjj
			ri = basis%sec(p)%r(k,i); ! 
			rjj = basis%sec(n-p-1)%r(k,jj);
			Hg%sec(p)%vdat(ind) = ri*rjj;
			!write(*,*)'Hg: k, i,jj: ri, rjj = ', k, i,jj, ri, rjj
		!write(*,'(a,5i5)')'k,i,j; ii,jj = ',k,i,j, ii,jj 
		!write(*,'(a,5f10.3)')'P,s = ',Pnui,Pmui, Pnuf,Pmuf


		end do ! k 
		end do ! j
		end do ! ii
		!write(*,*)'Hg%sec(p)%vdat(:) = ', Hg%sec(p)%vdat(:)
	end do ! p

	! to use for size of the full basis of p-blocks (i.e., with p up spins)
	do p = 0,n
		Hg%sec(p)%ntot = basis%sec(p)%ntot * basis%sec(n-p)%ntot
	end do

	return
	end 	subroutine MakeHgSec
	!----------------------------------------------------------




	!----------------------------------------------------------
	subroutine MakeHg(n,m)
	! constructs sparse Hg for m excitations using Hg%sec
	! saves the Hamiltonians in CSR format
	implicit none
	integer, intent(in) :: n,m
	double precision:: hgxp
	integer :: m1,nnz,ntot,i,p,ind1,ind2,ind
	integer :: nnup, nmunp ! n-nu-p; n-mu-n-minus-p


	m1=min(m,n);
	
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
	Hg%nnz = nnz
	! aux arrays to hold data from all p-blocks in coo format
	if(allocated(Hg%coo1))deallocate(Hg%coo1)
	if(allocated(Hg%coo2))deallocate(Hg%coo2)
	if(allocated(Hg%coodat))deallocate(Hg%coodat)
	allocate(Hg%coo1(nnz))
	allocate(Hg%coo2(nnz))
	allocate(Hg%coodat(nnz))

	!nnz = 0;
	!do p=0,m1-1
	!	nnz += Hg%sec(0:m1-1)%nnz
	!end do ! p
	ntot = sum(Hg%sec(0:m1)%ntot) ! Hg%sec%ntot  set in makeHgSec ??
	Hg%ntot = ntot
	!nrow = Hg%ntot;

	! coo format but with global indexing 
	!shift = 0;
	i = 0;
	do p=0,m1-1
		hgxp = dsqrt((m-p)*(p+1)*(n-p)*1.0d0); ! exciton-photon matrix element
		!write(*,*)'hd: hgxp = ',hgxp
		do ind=1,Hg%sec(p)%nnz
			!ind1 = Hg%sec(p)%row(ind)
			!ind2 = Hg%sec(p)%col(ind)
			! global position: i,j
			!i = ind + shift;
			i = i + 1;
			Hg%coo1(i) = origin(p) + Hg%sec(p)%row(ind);
			Hg%coo2(i) = origin(p+1) + Hg%sec(p)%col(ind);
			Hg%coodat(i) = hgxp * Hg%sec(p)%vdat(ind); ! full matrix element		
		end do! ind: nnz elements in p-th sector
	end do! p: exciton-photon sectors


	! AHSAN, JULY
	if(allocated(Hg%sec)) deallocate(Hg%sec)


	!-----------------------------------------------
	! allocate mem to hold final Hg in CSR format
	!if(allocated(Hg%rowpntr)) deallocate(Hg%rowpntr)
	!if(allocated(Hg%col)) deallocate(Hg%col)
	!if(allocated(Hg%dat)) deallocate(Hg%dat)
	!allocate(Hg%rowpntr(ntot + 1))
	!allocate(Hg%col(nnz))
	!allocate(Hg%dat(nnz))
	! convert COO to CSR sparse format
	! coocsr( nrow, nnz, a, ir, jc, ao, jao, iao )
	!Inp: Hg%ntot,Hg%nnz, Hg%coodat, Hg%coo1, Hg%coo2
	!Out: Hg%dat, Hg%col, Hg%rowpntr
	!call coocsr(Hg%ntot,Hg%nnz, 
  !   .  Hg%coodat, Hg%coo1, Hg%coo2, Hg%dat, Hg%col, Hg%rowpntr)
	!-----------------------------------------------
	!deallocate(Hg%coodat, Hg%coo1, Hg%coo2) ! we can reuse these ....
	!	by simply scaling sqrt[(m-p)/(m_old-p)] for sec p.
	! and then adding extra blocks for larger m, if needed, and then converting to CSR.

	return
	end 	subroutine MakeHg
	!----------------------------------------------------------








	!----------------------------------------------------------
	subroutine MakeHgAndVx(n,m) ! makes Hg and Vx [exciton decay]
	! constructs sparse Hg for m excitations using Hg%sec
	! saves the Hamiltonians in CSR format
	implicit none
	integer, intent(in) :: n,m
	double precision:: hgxp, vxxd
	integer :: m1,nnz,ntot,i,p,ind1,ind2,ind
	integer :: nnup, nmunp ! n-nu-p; n-mu-n-minus-p


	m1=min(m,n);
	
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
	Hg%nnz = nnz
	! aux arrays to hold data from all p-blocks in coo format
	if(allocated(Hg%coo1))deallocate(Hg%coo1)
	if(allocated(Hg%coo2))deallocate(Hg%coo2)
	if(allocated(Hg%coodat))deallocate(Hg%coodat)
	allocate(Hg%coo1(nnz))
	allocate(Hg%coo2(nnz))
	allocate(Hg%coodat(nnz))

	if(allocated(Vx))deallocate(Vx)
	allocate(Vx(nnz))

	!nnz = 0;
	!do p=0,m1-1
	!	nnz += Hg%sec(0:m1-1)%nnz
	!end do ! p
	ntot = sum(Hg%sec(0:m1)%ntot) ! Hg%sec%ntot  set in makeHgSec ??
	Hg%ntot = ntot
	!nrow = Hg%ntot;

	! coo format but with global indexing 
	!shift = 0;
	i = 0;
	do p=0,m1-1
		! MAZ 12 May 2019: matrix elem for p+1 to p excitons  
		hgxp = dsqrt((m-p)*(p+1)*(n-p)*1.0d0); ! exciton-photon matrix element
		vxxd = dsqrt((p+1)*(n-p)*1.0d0); ! exciton-decay matrix element
		do ind=1,Hg%sec(p)%nnz
			!ind1 = Hg%sec(p)%row(ind)
			!ind2 = Hg%sec(p)%col(ind)
			! global position: i,j
			!i = ind + shift;
			i = i + 1;
			Hg%coo1(i) = origin(p) + Hg%sec(p)%row(ind);
			Hg%coo2(i) = origin(p+1) + Hg%sec(p)%col(ind);
			Hg%coodat(i) = hgxp * Hg%sec(p)%vdat(ind); ! full matrix element		
			Vx(i) = vxxd * Hg%sec(p)%vdat(ind); 
		end do! ind: nnz elements in p-th sector
	end do! p: exciton-photon sectors

	!-----------------------------------------------
	! allocate mem to hold final Hg in CSR format
	!if(allocated(Hg%rowpntr)) deallocate(Hg%rowpntr)
	!if(allocated(Hg%col)) deallocate(Hg%col)
	!if(allocated(Hg%dat)) deallocate(Hg%dat)
	!allocate(Hg%rowpntr(ntot + 1))
	!allocate(Hg%col(nnz))
	!allocate(Hg%dat(nnz))
	! convert COO to CSR sparse format
	! coocsr( nrow, nnz, a, ir, jc, ao, jao, iao )
	!Inp: Hg%ntot,Hg%nnz, Hg%coodat, Hg%coo1, Hg%coo2
	!Out: Hg%dat, Hg%col, Hg%rowpntr
	!call coocsr(Hg%ntot,Hg%nnz, 
  !   .  Hg%coodat, Hg%coo1, Hg%coo2, Hg%dat, Hg%col, Hg%rowpntr)
	!-----------------------------------------------
	!deallocate(Hg%coodat, Hg%coo1, Hg%coo2) ! we can reuse these ....
	!	by simply scaling sqrt[(m-p)/(m_old-p)] for sec p.
	! and then adding extra blocks for larger m, if needed, and then converting to CSR.

	return
	end 	subroutine MakeHgAndVx
	!----------------------------------------------------------




	!----------------------------------------------------------
	subroutine MakeHd(n,m) ! diagonal part of cavity-exciton hamiltonian
	! constructs sparse Hd for m excitations using Hg%sec
	implicit none
	integer, intent(in) :: n,m
	integer :: i1,m1,ntot,p

	if(allocated(Hd)) deallocate(Hd)
	allocate(Hd(Hg%ntot))  ! set Hg%ntot first
	!write(*,*) 'Hg%ntot = ', Hg%ntot 
	m1=min(m,n);
	! coo format but with global indexing 
	i1 = 0;
	do p=0,m1
		ntot = basis%sec(p)%ntot * basis%sec(n-p)%ntot ! N_{nu_p} * N_{mu_p}
		!write(*,*) 'p, pth_ntot = ', p, ntot
		Hd(i1+1:i1+ntot) = (m-p)*1.0d0
		i1 = i1 + ntot;
	end do


	return
	end 	subroutine MakeHd
	!----------------------------------------------------------



	!----------------------------------------------------------
	subroutine MakeHv(n,m) ! diagonal part of vibrational hamiltonian
	! constructs sparse Hv using basis%sec
	implicit none
	integer, intent(in) :: n,m
	integer :: i1,m1,ntot,p, Nvi, i2,i

	if(allocated(Hv)) deallocate(Hv)
	allocate(Hv(Hg%ntot)) ! set Hg%ntot first

	m1 = min(n,m);
	! coo format but with global indexing 
	i1 = 0;
	do p=0,m1
		!ntot = basis%sec(p)%ntot * basis%sec(n-p)%ntot ! N_{nu_p} * N_{mu_p}
		! vib on p excited and n-p unexcited molecules
		do i=1,basis%sec(p)%ntot
			!write(*,*)'i,p, basis%sec(p)%ntot =',i,p, basis%sec(p)%ntot
			Nvi = basis%sec(p)%Nv(i)
			i2 = i1 + basis%sec(n-p)%ntot
			Hv(i1+1:i2) = Nvi + basis%sec(n-p)%Nv(:)
			i1 = i2; 
		end do
	end do

	return
	end 	subroutine MakeHv
	!----------------------------------------------------------

	



	!----------------------------------------------------------
	subroutine MakeHbSec(n,m,mv) ! m for mv==M
	implicit none
	integer, intent(in) ::n,m, mv
	integer :: fi, fi1
	integer :: p, nnz,nnp,ind,k,i,j,ii !,indi,indf 
	double precision :: Pnui,Pnuf
	integer:: ntotp
	
	! basis info: Normalisations of vibrational states
	! define indexes from 0 for sec of Ham and basis so that p sec is p-th elem.
	! similarly, map ind from 0,0
	! vibrational maps that connect perm symm states of p sites with those of p+1 sites... 

	!m1max set in main program

	if(allocated(Hb%sec)) deallocate(Hb%sec)
	allocate(Hb%sec(0:m1max))
	!allocate(Hg%nnzs(m1max))	! array for number of nnz in all blocks
	!write(6,*)'Hbsec: m1max=',m1max
	! loop over exciton-photon blocks.
	! p=0, no hb
	do p=1,m1max ! m1=Min(N,m); m1max= max if calc for multiple m's are desired.
		! nnz = N_{p-1}*M*N_{n-p} non-zero elements in pth block of Hb: some may have the same coordinates, will be combined before converting to csr.
		!nnz=(basis%sec(p-1)%ntot)*mv !*(basis%sec(n-p)%ntot) ! N_{n-p} states with down spins.; Hb is diagonal in these. treated later.
		ntotp = min(basis%sec(p-1)%ntot, basis%sec(ndummy-1)%ntot);
		nnz= ntotp * mv;
		Hb%sec(p)%nnz = nnz
		allocate(Hb%sec(p)%row(nnz))
		allocate(Hb%sec(p)%col(nnz))
		allocate(Hb%sec(p)%vdat(nnz))

		nnp = basis%sec(n-p)%ntot; ! number of perm symm basis for n-p sites [for mu, spin down n-p sies]

		Hb%sec(p)%ntot = basis%sec(p)%ntot ! nu_p basis only, mu_{n-p} treated later

		ind=0;
		do ii = 1, ntotp !basis%sec(p-1)%ntot
		do k=0,mv-1 
			ind = ind + 1; 
			i = map(k,ii);
			j = map(k+1,ii);	

			Hb%sec(p)%row(ind) = i !(l,l=(i-1)*nnp+1,i*nnp)
			Hb%sec(p)%col(ind) = j !(l,l=(j-1)*nnp+1,j*nnp)

			! number of total perm that basis states i, j stand for.
			!Pnui=basis%sec(p)%P(i);
			!Pnuf=basis%sec(p)%P(j);
			! frequency of k in i-th perm symm state, fi = ?
			fi = basis%sec(p)%f(k,i)	 ! from basis .... tabulated
			fi1 = basis%sec(p)%f(k+1,i)	

			!Hb%sec(p)%vdat(ind) = fi*dsqrt((k+1)*Pnui/Pnuf) ! ahsan, 14-01-2020
			!Hb%sec(p)%vdat(ind) = fi * dsqrt((k+1)*1.0d0/(fi+1))
			Hb%sec(p)%vdat(ind) = dsqrt((k+1)*1.0d0*fi*(fi1+1))

			! for details, see notepad starting 13 jan 2020 on chi2 code...
			! Pnui/Pnuf = 1/(fi+1)
		end do ! k 
		end do ! ii
	end do ! p

	return
	end 	subroutine MakeHbSec
	!----------------------------------------------------------







	!----------------------------------------------------------
	subroutine MakeHb(n,m,mv) !,spref)
	! constructs sparse Hb for m excitations using Hb%sec
	! NOTE: only first spref elements of Hb%coo1 , Hb%coo12, Hb%coodat arrays are significant.
	implicit none
	integer, intent(in) :: n,m,mv
	!integer, intent(out) :: spref
	integer :: spref
	integer, dimension(:), allocatable :: iro, jco
	double precision, dimension(:), allocatable :: aoo
	integer :: i1,i2,l1,l2,j1,j2,ii,k,i,j,m1,maxsize
	integer :: nnp,nnp1,nnz,nnzu,p,nrow,ref
	double precision :: val
	
	! find max size of sector to allocate aux arrays
	maxsize = 0; nnz = 0;
	m1=min(m,n);
	do p=1,m1 ! exclude p=0 for Hb 
		!write(*,*) 'p, Hb%sec(p)%nnz = ',p,Hb%sec(p)%nnz
		maxsize=max(maxsize,Hb%sec(p)%nnz) ! maxsize in nu_p space
		nnp = basis%sec(n-p)%ntot ! N_{N-p}; (for unexcited sites)
		nnz = nnz + Hb%sec(p)%nnz * nnp ! nnz elem in full [nu x mu] space
	enddo




	!write(*,*)'maxsize = ',maxsize
	
	! aux arrays:
	allocate(iro(maxsize))
	allocate(jco(maxsize))
	allocate(aoo(maxsize))
	
	! full Hb coo arrays: in FULL [nu x mu] space.
	if(allocated(Hb%coo1))deallocate(Hb%coo1)
	if(allocated(Hb%coo2))deallocate(Hb%coo2)
	if(allocated(Hb%coodat))deallocate(Hb%coodat)
	allocate(Hb%coo1(nnz))
	allocate(Hb%coo2(nnz))
	allocate(Hb%coodat(nnz))
	if(allocated(Hb%spntr))deallocate(Hb%spntr)
	allocate(Hb%spntr(m1+1))
	Hb%spntr(1) = 1;

	! total nnz
	!m1=min(m,n);
	ref = basis%sec(n-0)%ntot !N_{p=0}; all photons, no exciton block: binomial(n+m,m) states.
	spref = 0;
	do p=1,m1
		! sorted coo format & sum over duplicate entries.
		! subroutine sumdup( nrow, nnz, ir, jc, a, nnzu, iro, jo, aoo )
		nrow = basis%sec(p)%ntot ! not full Hg%sec(p)%ntot
		nnz = Hb%sec(p)%nnz;
		
		call sumdup(nrow, nnz, 
     .          Hb%sec(p)%row, Hb%sec(p)%col, Hb%sec(p)%vdat,
     .          nnzu, iro(1:nnz), jco(1:nnz), aoo(1:nnz))
		!...........................................................
		! get the global indices by taking into account diagonal mu states...
		! and combine data from all sectors into a larger coo sparse...
		!...........................................................
		! local indexes in p-th block: diagonal in mu
		!	I_{inu,imu} = (inu-1)*N_{N-p} + imu
		!	J_{jnu,imu} = (jnu-1)*N_{N-p} + imu
		nnp = basis%sec(n-p)%ntot ! N_{N-p}; number of perm sym basis for n-p sites; (for unexcited sites)
		l1 = 0; i1=0; j1=0;
		do ii=1,nnzu ! only nnzu elements of iro,jco,aoo are significant
			! diagonal in mu_{n-p} basis; nnp total states of type mu			
			i = iro(ii);
			j = jco(ii);
			val = aoo(ii);
			!..........................................
			i1 = ref + (i-1)*nnp; ! row indices global
			i2 = i1 + nnp;
			!..........................................	
			j1 = ref + (j-1)*nnp;! col indices global
			j2 = j1 + nnp;
			!..........................................			
			l1 = spref + (ii-1)*nnp; ! sparse coo indices
			l2 = l1 + nnp;
			!..........................................			 
			Hb%coo1(l1+1:l2) = (/(k,k=i1+1,i2)/);
			Hb%coo2(l1+1:l2) = (/(k,k=j1+1,j2)/);				
			Hb%coodat(l1+1:l2) = val
		end do	
		ref = ref + (basis%sec(p)%ntot)*nnp !nnp=(basis%sec(n-p)%ntot) ! global position of p-th block
		spref = spref + nnzu*nnp; ! sparse coo reference for p-th block
		Hb%spntr(p+1) = spref + 1
	end do


	! AHSAN, JULY
	deallocate(iro)
	deallocate(jco)
	deallocate(aoo)
	if(allocated(Hb%sec)) deallocate(Hb%sec)


	! use this to construct full hamiltonian, multiply wv, lambda etc there....
	! NOTE: only first spref elements of Hb%coo1 , Hb%coo12, Hb%coodat arrays are significant.
	Hb%nnz = spref ! only significant size of coo Hb.
	! Hb%nnz set to spref[= significant no of elements]
	! Do not reset it to sum of Hb%sec(p)%nnz
	return
	end 	subroutine MakeHb
	!----------------------------------------------------------

	subroutine MakeHhtc(n,ijob,mode)!, wr,delta,lambda,wv)
	implicit none
	integer, intent(in) :: n,ijob,mode
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

	nnz = Hg%nnz + Hb%nnz ; !+ Hg%ntot ! Hg%ntot for size of diagonal term, Hv+Hd
	if(mode==1) nnz = nnz + Hg%ntot; ! Hg%ntot for size of diagonal term, Hv+Hd 

	!write(*,*)'Hg%nnz, Hb%nnz, Hg%ntot =',Hg%nnz,Hb%nnz,Hg%ntot
	
	Hhtc%nnz = nnz;
	if(allocated(Hhtc%coo1))deallocate(Hhtc%coo1)
	if(allocated(Hhtc%coo2))deallocate(Hhtc%coo2)
	if(allocated(Hhtc%coodat))deallocate(Hhtc%coodat)
	allocate(Hhtc%coo1(nnz))
	allocate(Hhtc%coo2(nnz))
	allocate(Hhtc%coodat(nnz))

	if(allocated(Hdv)) deallocate(Hdv)
	allocate(Hdv(Hg%ntot))  ! set Hg%ntot first

	!write(*,*)'Hg%ntot = ',Hg%ntot

	Hdv = delta*Hd +  wv*Hv; ! summ arrays, values.
	! diagonal will be kept seperate

	n1 = 0;
	n3 = Hg%nnz;
	n2 = n3;

	!write(*,*) '-----   Hg%nnz = ',Hg%nnz

	! Hg first
	Hhtc%coo1(n1+1:n2) = Hg%coo1(1:n3)
	Hhtc%coo2(n1+1:n2) = Hg%coo2(1:n3)
	Hhtc%coodat(n1+1:n2) = g * Hg%coodat(1:n3)

	if (mv > 0) then
		!write(*,*) 'Hg: n1, n2 = ',n1, n2 
		n1 = n2;
		n3 = Hb%nnz; ! Hb%nnz set to spref[= significant no of elements]
								! Do not reset it to sum of Hb%sec(p)%nnz
		n2 = n1 + n3;
		! Hb afterwards
		Hhtc%coo1(n1+1:n2) = Hb%coo1(1:n3)
		Hhtc%coo2(n1+1:n2) = Hb%coo2(1:n3)
		Hhtc%coodat(n1+1:n2) = lamwv * Hb%coodat(1:n3)
	endif

	if(mode==1) then ! make Hf sparse here....

	n1 = n2;
	n3 = Hg%ntot;
	n2 = n1 + n3;
	! now diagonal terms
	Hhtc%coo1(n1+1:n2) = (/ (i,i=1,n3) /)
	Hhtc%coo2(n1+1:n2) = (/ (i,i=1,n3) /) 
	Hhtc%coodat(n1+1:n2) = Hdv(1:n3)*0.5d0 ! half it for matvec() using upper triangular only

	!write(*,'(a,3x,10000f15.10)') 'Hdv = ',Hhtc%coodat(n1+1:n2)

	!write(*,*) 'Hdv: n1, n2 = ',n1, n2 



	! AHSAN, JULY
	if(allocated(Hg%coo1))deallocate(Hg%coo1)
	if(allocated(Hg%coo2))deallocate(Hg%coo2)
	if(allocated(Hg%coodat))deallocate(Hg%coodat)

	if(allocated(Hb%coo1))deallocate(Hb%coo1)
	if(allocated(Hb%coo2))deallocate(Hb%coo2)
	if(allocated(Hb%coodat))deallocate(Hb%coodat)

	if(allocated(Hd)) deallocate(Hd)
	if(allocated(Hv)) deallocate(Hv)







	!write(*,*)'============ htc =========='
	
	Hf%dense = .false.
	Hf%ntot = Hg%ntot ! dimension of full hilbert space
	Hf%nnz = Hhtc%nnz; 
	if (Hf%ntot .le. nmaxddiag .and. ddiagOK) then 
		! Hhtc in dense format and use direct diagonalisation
		Hf%dense = .true.
		n1 =Hf%ntot;
		if(allocated(Hf%h))deallocate(Hf%h)
		allocate(Hf%h(n1,n1))
		call coo2dense(n1, n2, Hhtc%coo1(1:n2),
     .      Hhtc%coo2(1:n2),Hhtc%coodat(1:n2), Hf%h)
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
		call coocsr(Hf%ntot, Hf%nnz, 
     .  Hhtc%coodat, Hhtc%coo1, Hhtc%coo2,  
     .  Hf%dat, Hf%col, Hf%rowpntr)
	endif

	
	endif

	! ASHAN, JULY
	if(allocated(Hhtc%coo1))deallocate(Hhtc%coo1)
	if(allocated(Hhtc%coo2))deallocate(Hhtc%coo2)
	if(allocated(Hhtc%coodat))deallocate(Hhtc%coodat)

	return
	end subroutine MakeHhtc
	!--------------------------------------------------



	subroutine coo2dense(ntot,nnz,row,col,val,a)
	implicit none
	integer, intent(in) :: ntot,nnz
	integer, dimension(nnz), intent(in)  :: row, col
	double precision, dimension(nnz), intent(in)  :: val
	double precision, dimension(ntot,ntot), intent(out) :: a
	integer :: i,j,k
		a =0.0d0
		do k=1,nnz
			a(row(k),col(k)) = val(k)
		end do
		! symmetrise
		a = a + transpose(a); 
		! diagonal terms were halved so
		! no need to minus diagonalmatrix(diagonal(a))
	return
	end 	subroutine coo2dense

	!--------------------------------------------------
	subroutine coocsr( nrow, nnz, a, ir, jc, ao, jao, iao )
	!----------------------------------------------------------
	! MAZ Oct 2018, Modified from sparskit, @ Youcef Saad, 07 January 2004
!! COOCSR converts COO to CSR.
!    This routine converts a matrix that is stored in COO coordinate format
!    a, ir, jc into a CSR row general sparse ao, jao, iao format.
! IN:
!    Input, integer NROW, the row dimension of the matrix.
!    Input, integer NNZ, the number of nonzero elements.
!
! a,
! ir,
! jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz
!         nonzero elements of the matrix with a(k) = actual real value of
!         the elements, ir(k) = its row number and jc(k) = its column
!        number. The order of the elements is arbitrary.
!
! OUT:
!    Output, real AO(*), JAO(*), IAO(NROW+1), the matrix in CSR
!    Compressed Sparse Row format.
!
	implicit none

	! in
	integer, intent(in):: nrow, nnz
	integer, dimension(nnz), intent(in):: ir, jc
	double precision, dimension(nnz), intent(in) :: a
	! out
	integer, dimension(nnz), intent(out):: jao
	integer, dimension(nrow+1), intent(out):: iao
	double precision, dimension(nnz), intent(out) :: ao

	! local
	integer:: i, iad,j,k,k0
	double precision:: x

	iao(1:nrow+1) = 0
!
!  Determine the row lengths.
!
	do k = 1, nnz
		iao(ir(k)) = iao(ir(k)) + 1
	end do
!
!  The starting position of each row.
!
	k = 1
	do j = 1, nrow+1
		k0 = iao(j)
		iao(j) = k
		k = k + k0
	end do
!
!  Go through the structure once more.  Fill in output matrix.
!
	do k = 1, nnz
		i = ir(k)
		j = jc(k)
		x = a(k)
		iad = iao(i)
		ao(iad) = x
		jao(iad) = j
		iao(i) = iad + 1
	end do
!
!  Shift back IAO.
!
	do j = nrow, 1, -1
		iao(j+1) = iao(j)
	end do
	iao(1) = 1
	
	return
	end subroutine coocsr 
	!----------------------------------------------------------








	!--------------------------------------------------
	subroutine sumdup( nrow, nnz, ir, jc, a, nnzu, iro, jo, aoo )
	!----------------------------------------------------------
	! MAZ Oct 2018, sum duplicate entries in coo format.
	! IN:
	!    Input, integer NROW, the row dimension of the matrix.
	!    Input, integer NNZ, the number of nonzero elements.
	!
	! a,
	! ir,
	! jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz
	!         nonzero elements of the matrix with a(k) = actual real value of
	!         the elements, ir(k) = its row number and jc(k) = its column
	!        number. The order of the elements is arbitrary.
	!
	! OUT:
	!    Output, iro, jo, aoo unique entries in  COO
	!
	implicit none

	! in
	integer, intent(in):: nrow, nnz
	integer, dimension(nnz), intent(in):: ir, jc
	double precision, dimension(nnz), intent(in) :: a
	! out
	integer, intent(out) :: nnzu ! nnz-unique
	integer, dimension(nnz), intent(out):: jo, iro ! col and row indexes; only first nnzout significant
	double precision, dimension(nnz), intent(out) :: aoo ! values; only first nnzout significant

	! local
	integer i, iad,j,k,k0, nnzout, i1,i2, ind1, ind2,nnzr
	double precision:: x
	integer, dimension(nrow+1):: iao	
	integer, dimension(nnz):: jao
	double precision, dimension(nnz) :: ao
	double precision, dimension(:), allocatable :: valout
	integer, dimension(:), allocatable :: colout


	!write(*,*)'0    ir = ',ir
	!write(*,*)'0    jc = ',jc
	! ........................................
	!			sort w.r.t row index
	! ........................................
	iao(1:nrow+1) = 0
	!  Determine the row lengths.
	do k = 1, nnz
		iao(ir(k)) = iao(ir(k)) + 1
	end do

	!write(*,*)'0    iao = ',iao


	!  The starting position of each row.
	k = 1
	do j = 1, nrow+1
		k0 = iao(j)
		iao(j) = k
		k = k + k0
	end do
	!  Go through the structure once more.  Fill in output matrix.
	do k = 1, nnz
		i = ir(k)
		j = jc(k)
		x = a(k)
		iad = iao(i)
		ao(iad) = x
		iro(iad) = i
		jao(iad) = j
		iao(i) = iad + 1
	end do
	!write(*,*)'rev iao = ',iao

	! iro, jao, ao: sorted w.r.t row, in coo format
	!  Shift back IAO.
	do j = nrow, 1, -1
		iao(j+1) = iao(j)
	end do
	iao(1) = 1
	!write(*,*)'iao = ',iao
	!write(*,*)'jao = ',jao
	
	! ........................................
	!			sum over duplicate terms.
	! ........................................
	ind1 = 0;
	do i=1, nrow
		! range of entries for this row
		nnzr = iao(i+1)-iao(i); ! number of elem in ith row
		!......................................................
		if (nnzr<1) cycle ! go to next iteration; here, its probably last iteration.
		!......................................................		
		i1 = iao(i); ! starting index
		i2 = iao(i+1) - 1; ! end index
		! rowadddup(nnz,col,val, nnzout, colout,valout)
		allocate(colout(nnzr))
		allocate(valout(nnzr))
		!write(*,*) '----- nnzr = ',nnzr

		call rowadddup(nnzr,jao(i1:i2),ao(i1:i2),nnzout,colout,valout)
		ind2 = ind1 + nnzout;
		iro(ind1+1:ind2) = i; ! row index final
		!iaoo(i+1) = iaoo(i) + nnzout; ! rowpntr final
		jo(ind1+1:ind2) = colout(1:nnzout) ! col index final
		aoo(ind1+1:ind2) = valout(1:nnzout) ! value final
		ind1 = ind2; ! advance index for next iteration
		deallocate(colout,valout)
	enddo
	nnzu = ind2; ! total number of unique nnz values

	return
	end subroutine sumdup
!----------------------------------------------------------------------


!----------------------------------------------------------------------
	subroutine rowadddup(nnz,col,val, nnzout, colout,valout)
	! adds duplicate terms for a given row.
	implicit none
	integer, intent(in):: nnz ! nnz in this row === nnzrow 
	integer, dimension(nnz), intent(in):: col ! col index, with possible repetion
	double precision, dimension(nnz), intent(in):: val ! values
	integer, intent(out):: nnzout ! number of elements with distinct columns indices
	integer, dimension(nnz), intent(out):: colout ! distict col indeces 
	double precision, dimension(nnz), intent(out):: valout ! sum of duplicate values
	! local
	integer:: l,nset,k,j
	logical:: found

	!write(*,*) '----- nnz = ',nnz
	colout=0;
	valout=0.0d0;

	nset=0;
	! set the first element
	nset = nset + 1; 	! set the number of distinct col values set so far.
	colout(nset) = col(nset);
	valout(nset) = val(nset);

	found = .false.
	do l=2,nnz ! 
		j = col(l);
		do k=1,nset
			if(j==colout(k)) then ! found
				valout(k) = valout(k) + val(l) ! add l-th value to this box
				found = .true. 
				exit ! exit the inner do loop, over k.
			endif
		enddo ! k		
		if (.not. found) then
			nset = nset + 1; ! advance index to make space for new distict j value
			colout(nset) = j;
			valout(nset) = val(l);				
		endif
		found = .false. ! reset found for next iteration
	end do ! l

	!set total non-zero elements in out arrays
	nnzout = nset

	if(nnzout .ne. nnz) then
		write(*,*) 'rowadddup is useful: nnz, nnzout = ',nnz, nnzout
	endif
	
	return
	end subroutine rowadddup
!----------------------------------------------------------------------


	end !module


	
