

program main
use AP_funs
implicit none
  
  integer :: iz, i, mock_id
  
  real(rt) :: DAs(nz), Hs(nz), chisqs(nz-1), chisqs_uncored(nz-1)

  type(omwpar) :: nowpar
  logical :: smutabstds_inited, printinfo
  character(1000) :: tmpstr, filename 

  AP_inited = .false.
  printinfo = .false.

  ! testing ! let PT mock be sysmock 
  if(.false.) then
    mock_IO_test = .true.
    dataself_IO_test = .true.
  else
    mock_IO_test = .true.
    dataself_IO_test = .false.
    gb_syscor_mock = gb_PT_mock
    nsysmocks = 6
  endif
 
  do mock_id = 0, 0
  write(tmpstr,*) mock_id 
  if(mock_id < 10) then
          tmpstr = '000'//trim(adjustl(tmpstr))
  elseif(mock_id < 100) then
          tmpstr = '00'//trim(adjustl(tmpstr))
  elseif(mock_id < 1000) then
          tmpstr = '0'//trim(adjustl(tmpstr))
  endif
  gb_PT_IO_id = trim(adjustl(tmpstr))

  filename = 'test_lambdacdm_mockIO_'//trim(adjustl(tmpstr))//'input_4mocks_as_correction.txt'
  print*, 'Write result to :', trim(adjustl(filename))
  open(file=filename, unit=1923847)
  AP_inited = .false.

  do i = 1, 100
  nowpar%omegam = 0.1_rt+0.005*i; nowpar%w = -1.0_rt
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
  AP_inited = .true.

  print *
  print *
  print *
  print *, '!-------------------'
  print *, ' For Lambda CDM model with omegam = ', real(nowpar%omegam)
  print *, ' Chisq value of AP method, with / without sys correction = ', sum(chisqs), sum(chisqs_uncored)
  enddo
  close(1923847)
  enddo

end program
