

program main
use AP_funs
implicit none
  
  integer :: iz, i, j, i1i2s(N1,N2), i1,i2
  
  real(rt) :: DAs(nz), Hs(nz), chisqs(nz-1), chisqs__uncored(nz-1), sepchisqs_uncored(N1,N2), sepchisqs(N1,N2)

  type(omwpar) :: nowpar
  logical :: smutabstds_inited, printinfo, compute_covmat, use_bigcovmat, debug=.false., only_bigcov
  character(len=100000) :: covmat_suffixstr, tmpstr1, tmpstr2, tmpstr3, tmpstr4, tmpstr5, tmpstr6, suffixstr='', tmpstrID

  debug = .false. 
  AP_inited = .false.
  printinfo = .true.
  only_bigcov = .true.

  !!!#################################################
  !!! suffix str!!! ##############################
  !suffixstr = '_s6to50_polyfitdeg1_dataselfsyscor'
  !suffixstr = '_s6to40_polyfitdeg3'
  !suffixstr = '_s6to40_polyfitdeg3_om0.31_w-1.4'
  !suffixstr = '_s6to40_polyfitdeg3_om0.31_w-0.6'
  !suffixstr = '_s6to40_polyfitdeg3_om0.26_w-0.6'
  !suffixstr = '_s6to40_polyfitdeg3_om0.26_w-1.0_fact4to5_NBComp_Simp'
  !suffixstr = '_s6to40_polyfitdeg3_dataselfsyscor_om0.26_w-1.0_fact4to5_NBComp_Simp'
  !suffixstr = '_s6to40_polyfitdeg1_dataselfsyscor_om0.26_w-1.0_fact4to5_NBComp_Simp'
  gb_mock_IO_test_id=0; write(tmpstrID,*) gb_mock_IO_test_id; suffixstr = '_s6to40_polyfitdeg3_mockIOtest'//trim(adjustl(tmpstrID))//'_om0.26_w-1.0_fact1to1_NBComp_Simp'
  !suffixstr = '_s6to40_polyfitdeg3_om0.26_w-1.0_fact1to1_database150bin_NBComp_Simp'
  if(debug) suffixstr = trim(adjustl(suffixstr))//'_debug'
  compute_covmat = .false.

  write(tmpstr1,*) N1; write(tmpstr2,*) mubins(1); write(tmpstr3,*) mubins(N1)
  write(tmpstr4,*) N2; write(tmpstr5,'(f10.2)') mucuts(1); write(tmpstr6,'(f10.2)') mucuts(N2)
  tmpstr1 = trim(adjustl(tmpstr1))//'mubins'//trim(adjustl(tmpstr2))//'to'//trim(adjustl(tmpstr3))
  tmpstr4 = trim(adjustl(tmpstr4))//'mucuts'//trim(adjustl(tmpstr5))//'to'//trim(adjustl(tmpstr6))
  tmpstr1 = trim(adjustl(tmpstr1))//'_'//trim(adjustl(tmpstr4))//trim(adjustl(suffixstr))
  print *, '* Settings of mubins/mucuts: ', trim(adjustl(tmpstr1))

  tmpstr2 = ''
  do i2=1,N2; do i1=1,N1;
    write(tmpstr3,*) mubins(i1); write(tmpstr4,'(f10.2)') mucuts(i2);
    tmpstr2 = trim(adjustl(tmpstr2))//'  '//trim(adjustl(tmpstr3))//','//trim(adjustl(tmpstr4))
  enddo; enddo
  print *, '** detailed settings: ', trim(adjustl(tmpstr2))
  
  if(.not.only_bigcov) then
    open(file='omw1_'//trim(adjustl(tmpstr1))//'_uncor.txt', unit=101)
    write(101,'(A)') '# original chisq, 1900 mocks, uncored; nbins minor index, mucut major index:  ## omegam; w; chisq; sep_chisqs with mubin/mucuts ='//trim(adjustl(tmpstr2))
    open(file='omw1_'//trim(adjustl(tmpstr1))//'.txt', unit=102)
    write(102,'(A)') '# original chisq, 1900 mocks, cored; nbins minor index, mucut major index:  ## omegam; w; chisq; sep_chisqs with mubin/mucuts ='//trim(adjustl(tmpstr2))
  endif
  open(file='omw2_'//trim(adjustl(tmpstr1))//'_uncor.txt', unit=103)
  write(103,'(A)') '# big covmat chisq, 1900 mocks, uncored; nbins minor index, mucut major index:  ## omegam; w; chisq; sep_chisqs with mubin/mucuts ='//trim(adjustl(tmpstr2))
  open(file='omw2_'//trim(adjustl(tmpstr1))//'.txt', unit=104)
  write(104,'(A)') '# big covmat chisq, 1900 mocks, cored; nbins minor index, mucut major index:  ## omegam; w; chisq; sep_chisqs with mubin/mucuts ='//trim(adjustl(tmpstr2))



  do i = 1, 41
  do j = 1, 51
  !nowpar%omegam = 0.27_rt; nowpar%w = -1.0_rt
  nowpar%omegam = 0.1+0.01_rt*(i-1); nowpar%w = -1.5_rt+(j-1)*0.02
  if(debug) then; nowpar%omegam = 0.26; nowpar%w = -1.0; endif
  write(*,'(A,f8.3,"/",f8.3)') '############ LCDM, om/w = ', nowpar%omegam, nowpar%w
    
  ! DAs & Hzs
  do iz = 1, nz
    DAs(iz) = DAz_wcdm(nowpar,zeffs(iz))
    Hs(iz)  = Hz_wcdm(nowpar,zeffs(iz))
    !print *, iz, DAs(iz), Hs(iz)
  enddo
      
  ! AP likelihood
  covmat_suffixstr = ''
  if(.not.only_bigcov) then
    call AP_Like(DAs, Hs,  chisqs, printinfo, compute_covmat=compute_covmat, use_bigcovmat=.false., covmat_suffixstr=covmat_suffixstr, chisqs__uncored=chisqs__uncored, sepchisqs_uncored_result=sepchisqs_uncored, sepchisqs_result=sepchisqs) 
    write(101,'(<3+N1*N2>f14.3,1x)') nowpar%omegam, nowpar%w, sum(chisqs__uncored),  sepchisqs_uncored
    write(102,'(<3+N1*N2>f14.3,1x)') nowpar%omegam, nowpar%w, sum(chisqs),  sepchisqs
  !print *, ' omegam, Chisq value = ', nowpar%omegam, real(chisqs)
  endif

  if(i.eq.1 .and. j.eq.1)  AP_inited = .false.
  call AP_Like(DAs, Hs,  chisqs, printinfo, compute_covmat=compute_covmat, use_bigcovmat=.true., covmat_suffixstr=covmat_suffixstr, chisqs__uncored=chisqs__uncored, sepchisqs_uncored_result=sepchisqs_uncored, sepchisqs_result=sepchisqs)
  write(103,'(<3+N1*N2>f14.3,1x)') nowpar%omegam, nowpar%w, sum(chisqs__uncored),  sepchisqs_uncored
  write(104,'(<3+N1*N2>f14.3,1x)') nowpar%omegam, nowpar%w, sum(chisqs),  sepchisqs
  write(*,'(A,f4.2,<nz-1>f10.3)') ' omegam, Chisq value = ', nowpar%omegam, real(chisqs)
 
  compute_covmat = .false.

  !print *, '!-------------------'
  !print *, ' For Lambda CDM model with omegam = ', real(nowpar%omegam)
  if(debug) stop
  enddo
  enddo
  if(.not.only_bigcov) close(101); close(102)
  close(103); close(104)
end program
