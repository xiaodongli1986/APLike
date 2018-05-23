

program main
use AP_funs
implicit none
  
  integer :: iz, iom
  
  real(rt) :: DAs(nz), Hs(nz), chisqs(nz-1)

  type(omwpar) :: nowpar
  logical :: smutabstds_inited, printinfo

  AP_inited = .false.
  printinfo = .false.
  print *, '###################################################'
  print *, '** Compuate AP chisqs for Lambda CDM model: '
    
  !open(file='sdssdr12_AP_lcdm_chisq.mucut0.97.s20to25.withsys.table',unit=999)
  open(file='sdssdr12_AP_lcdm_chisq.mucut0.97.s20to25.table',unit=999)
  ! DAs & Hzs
  do iom = 1, 999
    nowpar%omegam = 0.001*iom
    nowpar%w = -1
    do iz = 1, nz
      DAs(iz) = DAz_wcdm(nowpar,zeffs(iz))
      Hs(iz)  = Hz_wcdm(nowpar,zeffs(iz))
    enddo
    call AP_Like(DAs, Hs,  chisqs, printinfo) 
    print *, '   omega_matter, w, chisq = ', real(nowpar%omegam), real(nowpar%w), real(sum(chisqs))
    write(999,*) nowpar%omegam, sum(chisqs)
    !print *, iz, DAs(iz), Hs(iz)
  enddo
  close(999)
      
  ! AP likelihood

  print *
  print *
  print *
  print *, '!-------------------'
  print *, ' For Lambda CDM model with omegam = ', real(nowpar%omegam)
  print *, ' Chisq value of AP method = ', sum(chisqs)

end program
