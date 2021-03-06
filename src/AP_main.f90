
!Rhct
! Sys-corred:      16.0076275       26.95636 56       12.7971764       13.1067686       19.3084583       88.176397148907327     
! Sys-uncorred:    16.5017948       23.0967922       12.0079880       11.4178925       14.0837126       77.108180501000902

! LCDM model. om =   0.27   &nbs p;
! Sys-corred:      14.0279408       21.7550316       11.8282537       10.6070719       12.7000446       70.918342580960200     
! param  mean           sddev          lower1         upper1         lower2         upper2         lower3         upper3
!    1  0.2774861E+00  0.2212044E-01  0.2600098E+00  0.2961426E+00  0.2395020E+00  0.3320312E+00  0.2160645E+00  0.3835449E+00 \Omega_m
! Sys-uncorred:    14.5382814       18.5175323       12.0481081       9.10813332       11.9272099       66.1392657 22333093  


!New result obtained near the end of 2017 (some difference found in the sys-corred result )
 !LCDM model. om =   0.27000000402331353     
 !Sys-corred:      14.0265808       21.7533760       11.8316622       10.6083479       12.6783466       70.898314222190251     
 !Sys-uncorred:    14.5382814       18.5175323       12.0481081       9.10813332       11.9272099       66.139265722333093 


!Result for omegam = 0.27 (precisely, double precision, 0.27, @ 20171109):
 !LCDM model. om =   0.27000000402331353     
! Sys-corred:      14.0311451       21.7494164       11.8326349       10.6039591       12.6807947       70.897949381941828     
! Sys-uncorred:    14.5427723       18.5142632       12.0503693       9.10425854       11.9387131       66.150375278426054 

! Time consumed computing 100 sets of parameters:
 ! ifort: 3.852 secs
 ! gfortran: 13.696 secs


! Main program computes the AP likelihood presented in Li et al. 2016 ()
!  in \Omega_m = 0.27 Lambda CDM
program main
use AP_funs
implicit none
  
  integer :: i, iz
  
  real(rt) :: DAs(nz), Hs(nz), chisqs(nz-1), t1,t2, chisqs__uncored(nz-1), APlnL, APlnL2

  type(omwpar) :: nowpar
  logical :: smutabstds_inited, inited, printinfo

  AP_inited = .false.
  printinfo = .true.
  
  open(unit=100001,file='AP_main_bigcov.txt',action='write')
  open(unit=100002,file='AP_main_smallcov.txt',action='write')
  write(100001,'(A)') '# omegam, lnlike, lnlike (uncorre)'
  write(100002,'(A)') '# omegam, lnlike, lnlike (uncorre)'

  do i = 1, 100
  nowpar%omegam = 0.1_rt + 0.004_rt * (i-1); nowpar%w = -1.0_rt
  print *, '###################################################'
  print *, '** Compuate AP chisqs for Lambda CDM model: '
  print *, '   omega_matter, w = ', nowpar%omegam, nowpar%w
    
  ! DAs & Hzs
  do iz = 1, nz
    DAs(iz) = DAz_wcdm(nowpar,zeffs(iz))
    Hs(iz)  = Hz_wcdm(nowpar,zeffs(iz))
    !print *, iz, DAs(iz), Hs(iz)
  enddo
      
!subroutine AP_like(DAs, Hs, chisqs,  printinfo, compute_covmat, use_bigcovmat, covmat_suffixstr, chisqs__uncored, sepchisqs_uncored_result, sepchisqs_result)
!    implicit none
!    real(rt), intent(in) :: DAs(nz), Hs(nz)
!    logical, intent(in) :: printinfo
!    real(rt), intent(out):: chisqs(nz-1)
!    real(rt), intent(out), optional :: chisqs__uncored(nz-1)
!    real(rt), intent(out), optional :: sepchisqs_uncored_result(n1,n2), sepchisqs_result(n1,n2)
!    logical, intent(in), optional :: compute_covmat, use_bigcovmat
!    character(*), intent(in), optional :: covmat_suffixstr


  ! AP likelihood
 ! call AP_Like(DAs, Hs,  chisqs, printinfo) 
  NBComp = .true.; chisqs=0; chisqs__uncored=0
  call AP_Like(DAs,Hs,chisqs,printinfo=printinfo,use_bigcovmat=.true.,chisqs__uncored=chisqs__uncored,compute_covmat=.false.)
  APlnL  = chisqs(1)/2.0;
  APlnL2 = chisqs__uncored(1)/2.0
  write(*,'(A,2f10.3)') '  lnlike (sys-corred, uncorred),    big    covmat: ', real(APlnL), real(APlnL2)
  write(100001,'(3e14.7)') nowpar%omegam, APlnL, APlnL2
  NBComp = .false.; chisqs=0; chisqs__uncored=0
!  if(i.eq.1) AP_inited=.false.
!  call AP_Like(DAs,Hs,chisqs,printinfo=printinfo,use_bigcovmat=.false.,chisqs__uncored=chisqs__uncored,compute_covmat=.false.)
!  APlnL  = sum(chisqs(1:nz-1))/2.0
!  APlnL2 = sum(chisqs__uncored(1:nz-1))/2.0
!  write(*,'(A,2f10.3)') '  lnlike (sys-corred, uncorred), ''small'' covmat: ', real(APlnL), real(APlnL2)
!  write(100002,'(3e14.7)') nowpar%omegam, APlnL, APlnL2



  if(i.eq.1)   call cpu_time(t1)
  enddo
  call cpu_time(t2)
  close(100001); close(100002);
  print *, 'Total time used: ', t2-t1


end program
