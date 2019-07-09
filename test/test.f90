

program main
use AP_funs
implicit none
  
<<<<<<< HEAD
  integer :: iz
  
  real(rt) :: DAs(nz), Hs(nz), chisqs(nz-1)
=======
  integer :: iz, i
  
  real(rt) :: DAs(nz), Hs(nz), chisqs(nz-1), chisqs_uncored(nz-1)
>>>>>>> 8e1c01b9cc8737f7369102ebb757359593e20ef2

  type(omwpar) :: nowpar
  logical :: smutabstds_inited, printinfo

  AP_inited = .false.
  printinfo = .false.
  
<<<<<<< HEAD
  nowpar%omegam = 0.27_rt; nowpar%w = -1.0_rt
=======
  open(file='test_lambdacdm.txt', unit=1923847)
  do i = 1, 100
  nowpar%omegam = 0.1_rt+0.005*i; nowpar%w = -1.0_rt
>>>>>>> 8e1c01b9cc8737f7369102ebb757359593e20ef2
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
<<<<<<< HEAD
  call AP_Like(DAs, Hs,  chisqs, printinfo) 
=======
!  call AP_Like(DAs, Hs,  chisqs, printinfo) 
  call AP_Like(DAs, Hs,  chisqs, printinfo, compute_covmat=.false., use_bigcovmat=.true., covmat_suffixstr='', chisqs__uncored=chisqs_uncored)
  write(1923847, *) nowpar%omegam, sum(chisqs), sum(chisqs_uncored)
>>>>>>> 8e1c01b9cc8737f7369102ebb757359593e20ef2

  print *
  print *
  print *
  print *, '!-------------------'
  print *, ' For Lambda CDM model with omegam = ', real(nowpar%omegam)
<<<<<<< HEAD
  print *, ' Chisq value of AP method = ', sum(chisqs)
=======
  print *, ' Chisq value of AP method, with / without sys correction = ', sum(chisqs), sum(chisqs_uncored)
  enddo
  close(1923847)
>>>>>>> 8e1c01b9cc8737f7369102ebb757359593e20ef2

end program
