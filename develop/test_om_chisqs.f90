

program main
use AP_funs
implicit none
  
  integer :: iz, i
  
  real(rt) :: DAs(nz), Hs(nz), chisqs(nz-1)

  type(omwpar) :: nowpar
  logical :: smutabstds_inited, printinfo, compute_covmat, use_bigcovmat
  character(len=1000) :: covmat_suffixstr

  AP_inited = .false.
  printinfo = .true.
  
  !compute_covmat = .false.
  !open(file='result1.txt', unit=100)
  !write(100,*) '# original chisq, 2000 mocks'

  compute_covmat = .false.
  !open(file='result2.txt', unit=101)
  !write(101,*) '# original chisq, 1900 mocks'
  !open(file='result3_LargePCorrection.txt', unit=102)
  open(file='result3_corred.txt', unit=102)
  write(102,*) '# big covmat chisq, 1900 mocks'

  do i = 1, 99
  !nowpar%omegam = 0.27_rt; nowpar%w = -1.0_rt
  nowpar%omegam = 0.2+0.002_rt*i; nowpar%w = -1.0_rt
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
  covmat_suffixstr = ''
  !call AP_Like(DAs, Hs,  chisqs, printinfo, compute_covmat=compute_covmat, use_bigcovmat=.false., covmat_suffixstr=covmat_suffixstr) 
  !!write(100,'(8f14.3,1x)') nowpar%omegam, sum(chisqs), real(chisqs)
  !write(101,'(8f14.3,1x)') nowpar%omegam, sum(chisqs), real(chisqs)

  call AP_Like(DAs, Hs,  chisqs, printinfo, compute_covmat=compute_covmat, use_bigcovmat=.true., covmat_suffixstr=covmat_suffixstr) 
  write(102,'(8f14.3,1x)') nowpar%omegam, sum(chisqs), real(chisqs)
  write(102,'(8f14.3,1x)') nowpar%omegam, sum(chisqs), real(chisqs)

  compute_covmat = .false.

  !print *, '!-------------------'
  !print *, ' For Lambda CDM model with omegam = ', real(nowpar%omegam)
  print *, ' omegam, Chisq value = ', nowpar%omegam, real(chisqs)
  
  enddo
  !close(100); 
  close(101); 
  close(102)
end program
