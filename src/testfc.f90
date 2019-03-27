

! Main program computes the AP likelihood presented in Li et al. 2016 ()
!  in \Omega_m = 0.27 Lambda CDM
program main
use AP_funs
implicit none
  
  integer :: i, iz
  
  real(rt) :: DAs(nzfc), Hs(nzfc), chisqs(nzfc-1), t1,t2, omfc, wfc

  type(omwpar) :: nowpar
  logical :: smutabstds_inited, inited, printinfo
  integer :: n, j
  real(rt), allocatable :: A(:,:), EVals(:), spEvals(:), EVecs(:,:), A2(:,:)
  real(rt) :: tmpx
  
  ! test of spmat ....
  if(.false.) then
	  open(file='bigcov1.txt',unit=1000)
	  n = 0
	  do while (.true.) 
	     read(1000,*,end=100) tmpx
	     n = n+1
	     cycle
	100  exit
	  enddo
	  close(1000)
	  open(file='bigcov1.txt',unit=1000); allocate(A(n,n),EVals(n),EVecs(n,n),A2(n,n),spEvals(n))
	  do i = 1, n
	     read(1000,*) A(i,:)
	  enddo
	  close(1000)
	!  print *, A

         call spmat_linsys_inout(A,n)
         open(file='bigcov_test.txt',unit=2000)
         do i = 1, n
           write(2000,'(<n>(e14.7,1x))') A(i,:)
         enddo
         close(2000)
         stop

	  call dsyev_linsys (A, n, EVals, EVecs )
	  !do i = 1, n
	  !  print *, '#########################'
	  !  print *, EVals(i)
	 !   print *, EVecs(i,:)
	 ! enddo
	  call spvec(Evals,n,spEvals)
	  do i = 1, n
	    print *, '#########################'
	    print *, EVals(i)
	    print *, spEvals(i)
	  enddo
	!  stop
	  !call reconmat_linsys (Evals, EVecs, n, A2)
	  call reconmat_linsys (spEvals, EVecs, n, A2)
	  do i = 1, n
	    do j = 1, n
	      print *, abs((A2(i,j) - A(i,j)) / A(i,j))
	    enddo
	  enddo
  endif
  
!
  !call dsyev_test()
!  stop
  ForeCast = .true.
  AP_inited = .false.
  printinfo = .true.
  
  do i = 1, 1
  nowpar%omegam = 0.16_rt + 0.01_rt * (i-1); nowpar%w = -1.0_rt
  print *, '###################################################'
  print *, '** Compuate AP chisqs for Lambda CDM model: '
  print *, '   omega_matter, w = ', nowpar%omegam, nowpar%w
    
  ! DAs & Hzs
  do iz = 1, nzfc
    DAs(iz) = DAz_wcdm(nowpar,zeffsfc(iz))
    Hs(iz)  = Hz_wcdm(nowpar,zeffsfc(iz))
!    print *, iz, DAs(iz), Hs(iz)
  enddo
      
  ! AP likelihood
  print *, 'omegam = ', nowpar%omegam
  call AP_Likefc(DAs, Hs, 0.26_rt, -1.0_rt, chisqs, printinfo)
  if(i.eq.1)   call cpu_time(t1)
  enddo
  call cpu_time(t2)
  print *, 'Total time used: ', t2-t1


end program
