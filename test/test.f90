

program main
use AP_funs
implicit none
  
  integer :: iz, i
  
  real(rt) :: DAs(nz), Hs(nz), chisqs(nz-1), chisqs_uncored(nz-1)

  type(omwpar) :: nowpar
  logical :: smutabstds_inited, printinfo

  AP_inited = .false.
  printinfo = .false.
  
  open(file='test_lambdacdm.txt', unit=1923847)
  do i = 0, 0
    nowpar%omegam = 0.27_rt+0.005*i; nowpar%w = -1.0_rt
    print *, '###################################################'
    print *, '** Compuate AP chisqs for Lambda CDM model: '
    print *, '   omega_matter, w = ', nowpar%omegam, nowpar%w
    
    ! DAs & Hzs
    do iz = 1, nz
      DAs(iz) = DAz_wcdm(nowpar,zeffs(iz))
      Hs(iz)  = Hz_wcdm(nowpar,zeffs(iz))
      !print *, iz, DAs(iz), Hs(iz)
    enddo
      
  ! AP likelihood
!  call AP_Like(DAs, Hs,  chisqs, printinfo) 
    call AP_Like(DAs, Hs,  chisqs, printinfo, compute_covmat=.false., use_bigcovmat=.true., covmat_suffixstr='', chisqs__uncored=chisqs_uncored)
    write(1923847, *) nowpar%omegam, sum(chisqs), sum(chisqs_uncored)

    print *
    print *
    print *
    print *, '!-------------------'
    print *, ' For Lambda CDM model with omegam = ', real(nowpar%omegam)
    print *, ' Chisq value of AP method, with / without sys correction = ', sum(chisqs), sum(chisqs_uncored)
  enddo
  close(1923847)

end program
