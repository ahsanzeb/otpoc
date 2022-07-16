


	module bases
	use modmain !, only: map, binomial
	
	implicit none

	public :: PermSymBasis, writebasis
	private :: facttable
	
	contains

	!===================================================================
	!..........................................
	! permutation symmetric vibrational basis
	!..........................................
	subroutine PermSymBasis(n,mv)
	!	calculates Nv, Perm, frequiencies of occup etc 
	! can be made more efficients.... by only calc the req sections of basis... 
	!e.g., for m excitations, we only need 0:m and n-m:n sections.... 
	!but, let's leave it for another time, if it really becomes necessary...
	!modmain: nsym
	implicit none

	integer, intent(in) :: n,mv
	! local
	integer, dimension(:,:), allocatable :: seti, setf
	integer :: ntot
	integer :: i,j,k, ntoti,rk,s, p
	double precision :: fac
	
	! with n sites, mv+1 vibrational states per molecule
	ntot = binomial(n+mv,mv)! largest size a sector can get,
	! aux arrays to hold actual sets of occup for all
	! perm sym basis for a given number of sites.
	allocate(seti(n,ntot))	 ! use largest for n-1 sites?
	allocate(setf(n,ntot))


	! allocate sectors
	allocate(basis%sec(0:nsym))


	! ------------ 0 site case: for (all up no down) or (no up, all down) cases
	allocate(basis%sec(0)%Nv(1))
	!allocate(basis%sec(0)%P(1))
	allocate(basis%sec(0)%f(0:mv,1))
	basis%sec(0)%Nv = (/ 0 /)
	!basis%sec(0)%P = (/ 1.0 /)
	basis%sec(0)%f(:,1) = (/(0,i=0,mv)/)
	basis%sec(0)%ntot = 1
	!..........................................
	! manulally set values for first site
	!set = 0
	! one site: n=1; mv+1 states= 0,1,2,...,mv
	! prepend sites one by one, so fist site is at nth position.
	! why prepend?, to use the same map_{p-1}:==> p-sector as in Polariton code.
	seti(n,1:mv+1) = (/(i,i=0,mv)/) ! set nth site occupations, add sites in everse direction.
	!ntoti = mv+1; ! number of states for 1 site [used for p=1 block, and, n-p=1 block]
	ntot = mv+1 !basis%sec(1)%ntot ! set somewhere....
	basis%sec(1)%ntot = ntot
	allocate(basis%sec(1)%Nv(ntot))
	!allocate(basis%sec(1)%P(ntot))
	allocate(basis%sec(1)%f(0:mv,ntot))
	basis%sec(1)%Nv(:) = (/(i,i=0,mv)/)
	!basis%sec(1)%P(:) = 1.0; !(/(1.0,i=0,mv)/)
	basis%sec(1)%f = 0
	do i=1,ntot ! ntot=mv+1
		basis%sec(1)%f(i-1,i) = 1; ! every occup in every site has freq 1
	enddo
	!..........................................
	! use data for 1 site to generate data for 2, and so on....	
	!..........................................
	ntoti = mv+1 !ntot
	do s=1,n-1 ! take a given s sector, make s+1 sector.
		ntot = binomial(s+1+mv,mv)
		!write(*,*)'ntot = binomial(s+1+m,m) =',ntot
		basis%sec(s+1)%ntot = ntot
		allocate(basis%sec(s+1)%Nv(ntot))
		!allocate(basis%sec(s+1)%P(ntot))
		allocate(basis%sec(s+1)%f(0:mv,ntot))
		j = 0;
		do i=1,ntoti
			!write(*,*) 'seti(n-s+1,i), j ,i= ',seti(n-s+1,i), j,i
			do k= seti(n-s+1,i), mv ! prepend only >= occupations 
				j = j + 1 ! advance index of final sets
				!write(*,*) 'k, j = ', k, j
				setf(n-s+1:n,j) = seti(n-s+1:n,i) ! copy all occup in s case 
				setf(n-s,j) = k ! set occup of additional site, prepended. 

				basis%sec(s+1)%Nv(j) = basis%sec(s)%Nv(i) + k ! number of vib quanta
				basis%sec(s+1)%f(:,j) = basis%sec(s)%f(:,i) ! freq of parent set
				basis%sec(s+1)%f(k,j) = basis%sec(s+1)%f(k,j) + 1 ! freq of k increased by 1
				! { P(j) = P(i) * (s+1)*rk!/(rk+1)! }  ! no need to repeat comput.
				!	use table of factorials: fact(i) = i!; i in [0,n] ; 
				!													table computed by a function once for all.
				!rk = basis%sec(s)%f(k,i);
				!fac = (s+1)*1.0d0/(rk+1) !(s+1)*fact(rk)/fact(rk+1);
				!basis%sec(s+1)%P(j) = basis%sec(s)%P(i)*fac
			end do! k
		end do ! i
		ntoti = j ! =ntot! next iteration builds sets based on sets made in this iteration
		seti(n-s:n,1:j) = setf(n-s:n,1:j)
		!setf(n-s:n,1:j) = 0
	enddo ! s
	!..........................................
	deallocate(seti,setf)

	! the rest of the block in case nlarge=T (and so nsym > ndummy [=n here] )
	! set all p > n blocks equal to nth block
	ntot = basis%sec(n)%ntot;
	!write(*,*)'basis: basis%sec(n)%ntot = ', basis%sec(n)%ntot
	do p=n+1,nsym
		basis%sec(p)%ntot = ntot;
		allocate(basis%sec(p)%Nv(ntot))
		allocate(basis%sec(p)%f(0:mv,ntot))
		basis%sec(p)%Nv = basis%sec(n)%Nv;
		basis%sec(p)%f = basis%sec(n)%f;
	end do

	!call writebasisf()
	
	! calc ratios: basis%sec(p)%r(k,i) = sqrt(P_{p}(i)/P_{p+1}(j)); j=map(k,i)
	call PermRatios(n,mv) 
	
	return
	end 	subroutine PermSymBasis
	!===================================================================

	subroutine PermRatios(n,mv)
	!	calculates ratios: basis%sec(:)%r(:,:)
	! basis%sec(p)%r(k,i) = sqrt(P_{p}(i)/P_{p+1}(j)); j=map(k,i)
	! use modmain, only: nsym, m, basis
	implicit none
	integer, intent(in) :: n, mv ! n can be ndummy if given in main to PermSymBasis()
	integer :: p, i, k, ntot
	integer :: gk1, np

	do p=0, nsym
		ntot = basis%sec(p)%ntot;
		
		if(p .le. n) then
			np = 0; ! all p sites used to calc basis%sec%f...
		else
			np = p-n; ! only ndummy used to calc basis%sec%f, add the rest of the sites (p-ndummy)
			! actual (p) minus the dummy (ndummy).
		endif
		
		allocate(basis%sec(p)%r(0:mv,ntot)) ! r = sqrt(P_{p}(i)/P_{p+1}(j)), j=map(k,i)
		do i=1,basis%sec(p)%ntot
			! k=0 case:
			k=0;
			gk1 = basis%sec(p)%f(k,i) + np + 1; ! all np molecules in the vib ground state.... 
			basis%sec(p)%r(k,i) = dsqrt(gk1*1.0d0/(p+1));
			! k > 0 cases:
			do k=1,mv
				gk1 = basis%sec(p)%f(k,i) + 1;
				basis%sec(p)%r(k,i) = dsqrt(gk1*1.0d0/(p+1));
			end do
		end do ! i

	end do ! p

	
	
	return
	end subroutine PermRatios
	!===================================================================
	
	function facttable(n)
	implicit none
	integer, intent(in) :: n
	integer, dimension(0:n) :: facttable ! length = n+1 
	integer :: i
	facttable(0) = 1;
	do i=1,n
		facttable(i) = facttable(i-1)*i
	enddo
	return	
	end function facttable
	!===================================================================
	! only writes basis%sec%P (number of permutations).
	subroutine writebasis()
	implicit none
	integer :: p
		open(1,file='basis.dat', form="unformatted", action="write")
		do p=0,n		
			write(1) basis%sec(p)%ntot
			write(1) basis%sec(p)%P(1:basis%sec(p)%ntot)
		enddo
		close(1)
	return
	end 	subroutine writebasis
	!===================================================================

	!===================================================================
	! writes basis%sec%f
	subroutine writebasisf()
	implicit none
	integer :: i
		open(1,file='basis-f.dat', form="formatted", action="write")
			write(1,*) basis%sec(n)%ntot
			do i=1,basis%sec(n)%ntot
			 write(1,'(1000i10)') basis%sec(n)%f(:,i)
			end do
		close(1)
	return
	end 	subroutine writebasisf
	!===================================================================


	end !module basis
