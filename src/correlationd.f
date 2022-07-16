
	module correlationd
	use modmain, only: Hhop, eig, mv, kappa2,prntstep,fft,
     .								ntotgg, ntotg, ntot, ntotb, ntotdn, eigd,lamd2wv
	use correlation, only: writeatnode
	implicit none

	public :: tcorrd
	private :: rk4step, yprime
	private :: makepsi0up, makepsi0dn
	private :: FourierT, writecorrt, writecorrw
		
	double precision :: dt, dth, dt6, twopi,kappa2r
	double complex :: iotam, phase
	double complex, dimension(:), allocatable:: corrt,corrw
	double precision, dimension(:), allocatable:: ws, ts,wss
	integer:: nt, nw, m
	character :: fname*100

	double complex,dimension(:), allocatable :: psi0, psit
	double precision, dimension(:), allocatable ::psitnorm
	integer :: ij
	
	contains
! Ahsan Zeb, 04 Nov, 2019, modified from poalriton code routine tcorr.f
! This subroutine takes the hamiltonian, initial state and other parameters, evolve the state and compute the two time correlation that is given as output.

	subroutine tcorrd(dt0,w1,w2,nt0,nw0,ij0,task) 
	implicit none
	double precision, intent(in) :: dt0,w1,w2
	integer, intent(in) :: nt0,nw0, ij0
	integer, intent(in):: task
	! local
	integer :: i 
	double precision:: dw, e0

	dt = dt0;	
	nt = nt0;
	if(fft) then
		nw = nt0; !FFTW: nw0;
	else
		nw = nw0;
	endif
	ij = ij0;

	! set variables for RKK4 integration
	dth = dt/2.0d0;
	dt6 = dt/6.0d0;
	iotam = (0.0d0,-1.0d0);

	kappa2r = kappa2 !/m;
	
	allocate(corrt(nt))
	allocate(psitnorm(nt))
	allocate(ws(nw))
	!allocate(wss(nw))
	allocate(ts(nt))
	
	! calc ws: frequencies....
	ws = 0.0d0;
	twopi= 2.0d0*3.14159265359d0;
	if(fft) then
		dw = twopi/(nt*dt); !(w2-w1)*1.0d0/nw;
	else
		dw = (w2-w1)*1.0d0/(nw-1);
	endif

	do i=1,nw
		!ws(i) = w1 + (i-1)*dw 
		ws(i) =  w1 + (i-1)*dw 
	end do
	! calc ts: times
	ts=0.0d0;
	do i=1,nt
		ts(i) = (i-1)*dt 
	end do


	! psi0, psit
	!ntot = ntotgg;
	allocate(psi0(ntotgg))
	allocate(psit(ntotgg))

	! construct init state for time evolution
	if(task == 103) then 	
		call makepsi0up(ij0) ! exciton projector |up><up| for hopping.... 
	elseif(task==104) then
		call makepsi0dn(ij0) ! exciton projector |dn><dn| for hopping.... 
	else
		stop "Error: Only task 103,104 for up/down hopping response!"
	endif

	write(6,*)'norm of psi0 = ', sum(abs(psi0)**2)
	!write(6,*)'kappa2 = ',kappa2r
	
	! set psi0 and corr(t=0)
	psitnorm(1) = real(DOT_PRODUCT(psi0, psi0));
	psi0 = psi0/dsqrt(psitnorm(1));
	psitnorm(1)=1.0d0;
	psit(:) = psi0(:);
	corrt(1) = dcmplx(1.0, 0.0); !DOT_PRODUCT(psi0, psit);  

	write(6,'(a,f15.10)')'ELP= ',eig(ij0)%eval(1);
	do i=2,nt
		call rk4step(psit)
		e0 = eig(ij0)%eval(1)-lamd2wv;
		phase = zexp((-iotam*e0-kappa2r)*(i-1)*dt);
		corrt(i) = DOT_PRODUCT(psi0, psit) * phase

		psitnorm(i) = real(DOT_PRODUCT(psit, psit))

		if (mod(i,prntstep) == 0) then
			write(6,'(a20,2x,i10,a10,i10)') ' time evolution step ',i,
     .  ' out of ',nt
			write(6,'(a,E15.10)')'|psit|/|psi| =', psitnorm(i)/psitnorm(1)
		endif

		if(psitnorm(i) > 10.0*psitnorm(1)) then
			write(*,*)'corr: RKK4 diverging... stopping..., decrease dt'
			write(*,*)'corr: psitnorm = ', psitnorm(i)
			stop
		endif
	end do

	call writecorrt(task)
	if(fft) then
		call FFourierT()
	else
		call FourierT()
	endif
	call writecorrw(task)

	deallocate(psi0,psit,corrt,corrw,psitnorm,ws,ts)

	write(6,'(a)')'tcorr: finished...!'
	return
	end subroutine tcorrd
	!----------------------------------------
	subroutine yprime(vin,vout)
	implicit none
	double complex, intent(in) :: vin(ntotgg)
	double complex, intent(out) :: vout(ntotgg)
	integer :: i,j	, k

	! Hhop upper or lower triangular only; diag elements halved.
	vout = 0.0d0
	do i=1,Hhop%ntot !nrp-1
		do j=Hhop%RowPntr(i),Hhop%RowPntr(i+1)-1
			k = Hhop%Col(j); 
			vout(i) = vout(i) + vin(k)*Hhop%dat(j);
			vout(k) = vout(k) +  vin(i)*Hhop%dat(j);
		end do
	end do
	vout = vout*iotam - kappa2r*vin
	return
	end subroutine yprime
	!----------------------------------------
	subroutine rk4step(y)
	implicit none
	double complex, intent(inout) :: y(ntotgg)
	! local
	double complex, dimension(ntotgg) :: k1,k2,k3
      call yprime(y,k1)
      call yprime(y+k1*dth,k2)
      k1 = k1 + 2*k2
      call yprime(y+k2*dth,k3)
      k1 = k1 + 2*k3
      call yprime(y+k3*dt,k2)
      k1 = k1 + k2
      y = y + k1*dt6
	return
	end subroutine rk4step
	!----------------------------------------

	subroutine makepsi0up(ijob) ! exciton projector... for hopping...
	implicit none
	integer,intent(in) :: ijob
	integer :: i, i1,i2, j, j1, j2

	! swap the ij-elem to ji-th elem:
	! i = active site's vib state, j = D's vib state.
	psi0(:) = 0.0d0;
	do i=0,mv
		i1= ntotdn + i*ntotb
		i2= i1 + ntotb
		do j=0,mv
			! orig pos: j*ntotg + i*ntotb to j*ntotg + (i+1)*ntotb
			! swap pos: i*ntotg + j*ntotb to i*ntotg + (j+1)*ntotb
			j1 = i*ntotg + ntotdn + j*ntotb
			j2 = j1 + ntotb
			psi0(j1+1:j2) = eig(ijob)%evec(i1+1:i2,1) * eigd(j+1) ! RHS orig pos given above
		end do
	end do

	return
	end 	subroutine makepsi0up
	!----------------------------------------
	subroutine makepsi0dn(ijob) ! exciton projector... for hopping...
	implicit none
	integer,intent(in) :: ijob
	integer :: i, i1,i2, j, j1, j2

	! swap the ij-elem to ji-th elem:
	! i = active site's vib state, j = D's vib state.
	psi0(:) = 0.0d0;
	do i=0,mv
		i1= i*ntot
		i2= i1 + ntot
		do j=0,mv
			! orig pos: j*ntotg + i*ntot to j*ntotg + (i+1)*ntot
			! swap pos: i*ntotg + j*ntot to i*ntotg + (j+1)*ntot
			j1 = i*ntotg + j*ntot
			j2 = j1 + ntot
			psi0(j1+1:j2) = eig(ijob)%evec(i1+1:i2,1) * eigd(j+1) ! RHS orig pos given above
		end do
	end do
	! last ntotup*(mv+1) elements of psi0 are zero
	return
	end 	subroutine makepsi0dn
	!----------------------------------------
	subroutine FourierT()
	!use fftpack5
	implicit none
	integer :: i

	allocate(corrw(nw))
	do i=1,nw
		corrw(i) = sum(zexp(-iotam*ts*ws(i))*corrt)
	end do
	corrw = corrw/dsqrt(twopi*nw*1.0d0);
	
	return
	end 	subroutine FourierT
	!----------------------------------------

	subroutine FFourierT()
	!use fftpack5
	implicit none
	integer :: i
	
	call cfftnd(1,nt,+1,corrt)

	! nw = nt ? cant we have more ws?
	allocate(corrw(nt))
	do i=1,nw/2 ! nw power of 2
		corrw(i) = corrt(nw/2+i)
		corrw(nw/2 + i) = corrt(i)
		ws = ws - ws(nw/2)
	end do
	corrw = corrw/dsqrt(twopi*nw*1.0d0);
	
	return
	end 	subroutine FFourierT
	!----------------------------------------
	subroutine writecorrt(task)
	implicit none
	integer, intent(in) :: task
	integer :: i
	open(10,file='hop-vs-t.dat',
     .             action='write',position='append')
     
	if(task==103) then
		write(10,'(a)')'# Exciton |up><up|'
	elseif(task==104) then
		write(10,'(a)')'# Exciton |dn><dn|'
	endif
	do i=1,nt
		write(10,'(f20.10,3x,2E20.10,2x,f15.10)') ts(i), 
     .   real(corrt(i)), aimag(corrt(i)), psitnorm(i)
	end do
	write(10,*)
	write(10,*)
	close(10)
	return
	end 	subroutine writecorrt
	!----------------------------------------
	subroutine writecorrw(task)
	integer, intent(in) :: task
	integer :: i
	open(10,file='hop-vs-w.dat',
     .             action='write',position='append')
	if(task==103) then
		write(10,'(a)')'# Exciton |up><up|'
	elseif(task==104) then
		write(10,'(a)')'# Exciton |dn><dn|'
	endif
	do i=1,nw
		write(10,'(f20.10,3x,2E20.10)') ws(i),
     .             real(corrw(i)), aimag(corrw(i))
	end do
	write(10,*)
	write(10,*)
	close(10)
	
	return
	end 	subroutine writecorrw
	!----------------------------------------
	subroutine mkfname(n)
	implicit none
	integer, intent(in) :: n
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
	!if(.not. ex) stop 'data directory does not exist!'

	fname = trim(dir)//'/wfs.dat'
	write(6,'(a)') fname
	
	return
	end subroutine mkfname
	!-------------------------------------------

	end module correlationd

