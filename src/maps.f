
	!program maps
	module maps
	use modmain, only: map, factorial, binomial

	implicit none
	public :: getmap, writemap

	
	contains

	!===================================================================
	! map for two sites
	subroutine mapn2(m) ! m = vibrational cutoff
	implicit none
	integer, intent(in):: m
	integer :: k,i,j
	k = 0;
	do i=0,m ! vib
		do j=i,m ! basis states for single site
			map(i,j+1) = k ! +1 to col index to make it 1:ntot for basis
			map(j,i+1) = k
			k = k + 1;
		end do
	end do
	
	!write(6,*) map(0:5,:)

	return
	end subroutine mapn2
	!===================================================================
	!integer function arglistlength(p,m) 
	! length of block that would be added
	! to the map if we add a site to p-1 sites
	!implicit none
	!integer, intent(in):: p,m
	!arglistlength = binomial(p,m) - binomial(p-1,m) 
	!return
	!end function arglistlength
	!===================================================================
	! map for more than two sites; make triangles and left side columns, 
	! see color shaded map table in doc of polariton code.

	subroutine getmap(n,m) !,ntot) !n=no of sites; m = Vibratioanl cutoff
	implicit none
	integer, intent(in) :: n,m
	!integer, intent(out) :: ntot
	! local
	integer :: ntot
	integer, dimension(m+1) :: colist
	integer, dimension(:), allocatable :: list1, list2
	integer :: lmax, iarg,ii,j1,l1,l2, karg,nn,nrows,iii,i,x
	integer :: lcol

	!write(*,*)'binomial(n-1+m,m) = ',binomial(n-1+m,m)
	allocate(map(0:m, binomial(max(n,2)-1+m,m)))
	map(:,:) = 0

	!do i=1,30
	!	write(*,*)'n, N_p =',i, binomial(i-1+m,m)
	!end do
	!do i=1,30
	!	write(*,*)'n, n! =',i, factorial(i)
	!end do
	!stop










	! map(0:m,1:m+1) set by mapn2(m)
	call mapn2(m) ! sets map for two sites, the starting square block of m+1 x m+1.
	if (n < 3 ) then
		ntot = (m+1) ! number of states of n-1 sites in the map
		map = map + 1
		return
	endif

	!return
	! max size list1,2 can get:
	lmax = binomial(n+m,m) ! roughly speaking... 
	allocate(list1(lmax))
	allocate(list2(lmax))

	!write(6,*) 'lmax = ',lmax
	
	!karg = mapfull[-1,-1] +1;
	karg = map(m,m+1) + 1; ! argument for next block's triangles.
	list1(1) = m !(/ ( m ) /);
	l1 = 1; ! length of list1
	j1 = m+1; ! last baiss index in the map for two sites; 1 :=> 2 sites

	do nn=3,n ! no of sites
		!write(*,*) 'list1 = ',list1(1:l1)
		do ii=1,l1
			iarg = list1(ii)
			!py: colist = mapfull[-1,0:m-iarg+1]+1; py: x=[0,1,2]; x[0:1] =[0]
			colist(1:m-iarg+1) = map(0:m-iarg,j1) + 1 ! j1 last row so far
			nrows = (iarg+1)*iarg/2;
			lcol = m-iarg+1
			call triangles(nrows,lcol,colist(1:lcol),iarg,karg,m,j1+1) ! karg inout
			j1 = j1 + nrows
		end do 
		!call mkarglist(l1,list1(1:l1),l2,list2(1:l2))
		!..............................
		if(nn<n) then
			l2 = 0;
			do i=1,l1
				x = list1(i);
				list2(l2+1:l2+x) = (/ (ii,ii=x,1,-1) /)
				l2 = l2 + x;
			end do
			!..............................	
			!write(*,*) 'l1, l2= '	,l1,l2,'lst=',list2(1:l2)
			list1(1:l2) = list2(1:l2)!
			l1 = l2;
		endif
	end do

	map = map + 1 ! basis indices 1:*
	ntot = j1 ! col dim of map
	deallocate(list1,list2)
	
	return
	end subroutine getmap
	!===================================================================
	! makes blocks of map for more than two sites using some info of prev blocks
	subroutine triangles(nrows,lcol,colist,iin,k,m, j1)
	implicit none
	integer, intent(in) ::nrows,lcol, iin, m, j1 
	integer, dimension(lcol), intent(in) :: colist
	integer, intent(inout) :: k
	! local
	integer :: i,j,jj,ii,i0,yi,yf,mjj,mj,j2,km,iii,j11

	!nrows = (iin+1)*iin//2;

	! fill all left side columns for this iin
	j2 = j1+nrows
	i = 0;
	do j=1,lcol
		jj = colist(j)
		map(i,j1:j2-1) = (/(ii,ii=jj,jj+nrows-1)/) ! py: range(0:1) =[0]
		i = i + 1;
	end do

	!write(*,*)'--------------------'
	
	! fill the triangles
	yi = 0;
	do i = 0, iin-1 ! iin \in {1,2,...,M}
		i0 = i+m-iin+1;
		yf = (i+1)*iin - i*(i+1)/2;
		do jj = yi, yf-1
			j = jj - yi;
			do mjj = i0+j, m
				mj = mjj - i0;
				map(mjj,j1+jj) = k
				map(i0+j,j1+yi+mj) = k
				if (j == 0 .and. mjj==m) km = k
				k = k + 1;
			end do ! mjj
		end do ! jj
		!write(*,*)'l:183 maps: i0 = ', i0
		!do iii=1,21
		!write(*,'(6i4)') map(:,iii)
		!end do
		!write(*,*) 'i0, j1+yf, j2-1, km+1, k-1 =',i0,j1+yf,j2-1,km+1,k-1
		!py: Map[yf:,i0] = range(km+1,k)
		map(i0,j1+yf:j2-1) = (/(ii,ii=km+1,k-1)/) ! left column below this triangle
		yi = yf;
	end do !i
	! k has advanced.... goes in the output
	return
	end subroutine triangles
	!===================================================================
	subroutine writemap()
	implicit none
		open(1,file='map.dat', form="unformatted", action="write")
			write(1) shape(map)
			write(1) map	
		close(1)
	return
	end 	subroutine writemap
	!===================================================================

	end !program maps
