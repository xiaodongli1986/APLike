

program main
use AP_funs
implicit none
  
  integer :: iz
  
  real(rt) :: DAs(nz), Hs(nz), chisqs(nz-1)

  type(omwpar) :: nowpar
  logical :: smutabstds_inited, printinfo

  AP_inited = .false.
  printinfo = .false.
  nowpar%omegam = 0.277302185883641       
  nowpar%w = -1.1635234496361
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
  call AP_Like(DAs, Hs,  chisqs, printinfo) 

  print *
  print *
  print *
  print *, '!-------------------'
  print *, ' For Lambda CDM model with omegam = ', real(nowpar%omegam)
  print *, ' Chisq value of AP method = ', sum(chisqs)

end program
