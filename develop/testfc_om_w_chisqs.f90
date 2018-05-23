

program main
use AP_funs
implicit none
  
  integer :: iz, i, j, i1i2s(N1,N2), i1,i2, nom, nw
  
  real(rt) :: DAs(nz), Hs(nz), DAsfc(nzfc), Hsfc(nzfc), chisqs(nz-1), chisqs__uncored(nz-1), &
    chisqsfc(nzfc-1), chisqs__uncoredfc(nzfc-1), sepchisqs_uncored(N1,N2), sepchisqs(N1,N2), &
    ommin, ommax, wmin, wmax, omfid

  type(omwpar) :: nowpar
  logical :: smutabstds_inited, printinfo, compute_covmat, use_bigcovmat, debug=.false., only_bigcov
  character(len=100000) :: covmat_suffixstr, tmpstr1, tmpstr2, tmpstr3, tmpstr4, tmpstr5, tmpstr6, suffixstr='', tmpstrID, fcstr,nzfcstr,nowfile

!  ommin = 0.20d0; ommax = 0.30d0; wmin = -1.20d0; wmax = -0.80d0;
  ommin = 0.23d0; ommax = 0.29d0; wmin = -1.10d0; wmax = -0.90d0;
  nom = 41; nw = 51;

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
  suffixstr = 's6to40'

  gb_fc_spmat=.true.; 

  ! name fo the run
  if(trim(adjustl(gb_fcname)).eq.'') then
    write(nzfcstr,*) nzfc;  fcstr = 'fc_DESI'//trim(adjustl(nzfcstr))//'bin'
  else
    fcstr = 'fc_'//trim(adjustl(gb_fcname))
  endif
  if(gb_fc_spmat) fcstr = trim(adjustl(fcstr))//'_spmat' 

  omfid = 0.3121_rt
  if(abs(omfid-0.26) .gt. 0.0001) then
    write(tmpstr1,'(f10.4)') omfid; write(tmpstr2,*) nom; write(tmpstr3,*) nw;
    fcstr = trim(adjustl(fcstr))//'_omfid'//trim(adjustl(tmpstr1)) !//'_nom'//trim(adjustl(tmpstr2))//'_nw'//trim(adjustl(tmpstr3))
    ommin = omfid - 0.03_rt*3; ommax = omfid + 0.03_rt*3
  endif


!  gb_mock_IO_test_id=0; write(tmpstrID,*) gb_mock_IO_test_id; suffixstr = '_s6to40_polyfitdeg3_mockIOtest'//trim(adjustl(tmpstrID))//'_om0.26_w-1.0_fact1to1_NBComp_Simp'
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

!###########################
  if(.true.) then
!  open(file='omw2_'//trim(adjustl(tmpstr1))//'_uncor.txt', unit=103)
!  write(103,'(A)') '# big covmat chisq, 1900 mocks, uncored; nbins minor index, mucut major index:  ## omegam; w; chisq; sep_chisqs with mubin/mucuts ='//trim(adjustl(tmpstr2))
  nowfile = 'omw2_'//trim(adjustl(fcstr))//'_'//trim(adjustl(tmpstr1))//'.txt'
  write(*,'(A,A)') 'Begin writing results into file: ', trim(adjustl(nowfile)); !stop
  open(file=nowfile, unit=104)
  write(104,'(A)') '# big covmat chisq, 1900 mocks, cored; nbins minor index, mucut major index:  ## omegam; w; chisq; sep_chisqs with mubin/mucuts ='//trim(adjustl(tmpstr2))


  do i = 1, nom
  do j = 1, nw
  !nowpar%omegam = 0.27_rt; nowpar%w = -1.0_rt
  nowpar%omegam = ommin+(ommax-ommin)/dble(nom-1)*(i-1);
  nowpar%w = wmin+(wmax-wmin)/dble(nw-1)*(j-1);
!  nowpar%omegam = 0.1+0.01_rt*(i-1); nowpar%w = -1.5_rt+(j-1)*0.02
  if(debug) then; nowpar%omegam = 0.26; nowpar%w = -1.0; endif
  write(*,'(A,f8.3,"/",f8.3)') '############ LCDM, om/w = ', nowpar%omegam, nowpar%w
    
  ! DAs & Hzs
  do iz = 1, nzfc
    DAsfc(iz) = DAz_wcdm(nowpar,zeffsfc(iz))
    Hsfc(iz)  = Hz_wcdm(nowpar,zeffsfc(iz))
    !print *, iz, DAs(iz), Hs(iz)
  enddo
      
  ! AP likelihood
  covmat_suffixstr = ''


  if(i.eq.1 .and. j.eq.1)  AP_inited = .false.
  call AP_Likefc(DAsfc, Hsfc, 0.3121_rt, -1.0_rt, chisqsfc, printinfo, compute_covmat=compute_covmat, covmat_suffixstr=covmat_suffixstr, chisqs__uncored=chisqs__uncoredfc, sepchisqs_uncored_result=sepchisqs_uncored, sepchisqs_result=sepchisqs)
! subroutine AP_likefc(DAs, Hs, omfc, wfc, chisqs,  printinfo, compute_covmat,  covmat_suffixstr, chisqs__uncored, sepchisqs_uncored_result, sepchisqs_result)
!  call AP_Like(DAs, Hs,  chisqs, printinfo, compute_covmat=compute_covmat, use_bigcovmat=.true., covmat_suffixstr=covmat_suffixstr, chisqs__uncored=chisqs__uncored, sepchisqs_uncored_result=sepchisqs_uncored, sepchisqs_result=sepchisqs)
!  write(103,'(<3+N1*N2>f14.3,1x)') nowpar%omegam, nowpar%w, sum(chisqs__uncored),  sepchisqs_uncored
  write(104,'(<3+N1*N2>(e14.7,1x))') nowpar%omegam, nowpar%w, sum(chisqs),  sepchisqs
  write(*,'(A,f8.4,<nz-1>f10.1)') ' omegam, Chisq value = ', nowpar%omegam, real(chisqs)

  compute_covmat = .false.

  !print *, '!-------------------'
  !print *, ' For Lambda CDM model with omegam = ', real(nowpar%omegam)
  if(debug) stop
  enddo
  enddo
!  if(.not.only_bigcov) close(101); close(102)
!  close(103); 
  close(104)
  write(*,'(A,A)') 'Finishing writing results into file: ', trim(adjustl(nowfile))
  endif
!###########################

  stop ! not calculating real data ...

! real data result
  dataself_IO_test = .true.
  debug = .false. 
  AP_inited = .false.
  printinfo = .true.
  only_bigcov = .true.

  print *, '* Settings of mubins/mucuts: ', trim(adjustl(tmpstr1))

  tmpstr2 = ''
  do i2=1,N2; do i1=1,N1;
    write(tmpstr3,*) mubins(i1); write(tmpstr4,'(f10.2)') mucuts(i2);
    tmpstr2 = trim(adjustl(tmpstr2))//'  '//trim(adjustl(tmpstr3))//','//trim(adjustl(tmpstr4))
  enddo; enddo
  print *, '** detailed settings: ', trim(adjustl(tmpstr2))
  
!  open(file='omw2_'//trim(adjustl(tmpstr1))//'_uncor.txt', unit=103)
!  write(103,'(A)') '# big covmat chisq, 1900 mocks, uncored; nbins minor index, mucut major index:  ## omegam; w; chisq; sep_chisqs with mubin/mucuts ='//trim(adjustl(tmpstr2))
  open(file='omw2_dr12_'//trim(adjustl(tmpstr1))//'.txt', unit=104)
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


  if(i.eq.1 .and. j.eq.1)  AP_inited = .false.
!  call AP_Like(DAs, Hs, 0.26_rt, -1.0_rt, chisqs, printinfo, compute_covmat=compute_covmat, covmat_suffixstr=covmat_suffixstr, chisqs__uncored=chisqs__uncored, sepchisqs_uncored_result=sepchisqs_uncored, sepchisqs_result=sepchisqs)
! subroutine AP_likefc(DAs, Hs, omfc, wfc, chisqs,  printinfo, compute_covmat,  covmat_suffixstr, chisqs__uncored, sepchisqs_uncored_result, sepchisqs_result)
  call AP_Like(DAs, Hs,  chisqs, printinfo, compute_covmat=compute_covmat, use_bigcovmat=.true., covmat_suffixstr=covmat_suffixstr, chisqs__uncored=chisqs__uncored, sepchisqs_uncored_result=sepchisqs_uncored, sepchisqs_result=sepchisqs)
!  write(103,'(<3+N1*N2>f14.3,1x)') nowpar%omegam, nowpar%w, sum(chisqs__uncored),  sepchisqs_uncored
  write(104,'(<3+N1*N2>e14.7,1x)') nowpar%omegam, nowpar%w, sum(chisqs),  sepchisqs
  write(*,'(A,f4.2,<nz-1>e15.7)') ' omegam, Chisq value = ', nowpar%omegam, real(chisqs)
 
  compute_covmat = .false.

  !print *, '!-------------------'
  !print *, ' For Lambda CDM model with omegam = ', real(nowpar%omegam)
  if(debug) stop
  enddo
  enddo
!  if(.not.only_bigcov) close(101); close(102)
!  close(103); 
  close(104)
end program
