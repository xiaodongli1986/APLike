
! ### Finishing writing dataselfsyscor code. check & test.

! 1. Plans

! 1. The file system: can find all files;
! 2. Loadin & ximu: can load in 2pcf files, and compute ximu
! 3. Cosmo convert: convert 2pcf from 1 to another
! 4. Covmat: compute covmat from all covmat mocks
! 5. Chisq: load in files, compute chisqs!


!###############################################################
!## Module: Types, constants, tools
module AP_type_constants
implicit none

  integer, parameter  :: rt =  kind(1.0d0)
  integer, parameter  :: charlen=1000
  real(rt), parameter :: CONST_C = 299792.458 !unit: km/s
  type :: omwpar 
    real(rt) :: omegam, w
    real(rt) :: wa = 0.0_rt
    real(rt) :: omegak = 0.0_rt
  end type

!########## binning in s !for 2-d binning...
  integer, parameter  :: gbnumsbin = 1 ! binning in s

!###########################################################################
!! Big Settings  
  ! Setting: The binning schemes (number of bins in the mu space )
!  integer, parameter :: N1=36
!  integer, parameter :: mubins(N1) = (/ 5,6,7,8,9,10,11,12,13,14,&
!                  15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30, &
!                  31,32,33,34,35,36,37,38,39,40 /)

!!  integer, parameter ::  N2=26
!!  real(rt), parameter :: mucuts(N2) = (/ 0.99_rt, 0.98_rt, 0.97_rt, 0.96_rt, 0.95_rt, 0.94_rt, &
!!     0.93_rt, 0.92_rt, 0.91_rt, 0.90_rt, 0.89_rt, 0.88_rt, 0.87_rt, 0.86_rt, 0.85_rt, 0.84_rt, 0.83_rt, 0.82_rt, &
!!     0.81_rt, 0.80_rt, 0.75_rt, 0.70_rt, 0.65_rt, 0.60_rt, 0.55_rt, 0.50_rt /) !, &
!!     !0.65_rt, 0.60_rt, 0.55_rt, 0.50_rt, 0.45_rt, 0.40_rt, 0.35_rt, 0.30_rt /)

!  integer, parameter ::  N2=35
!  real(rt), parameter :: mucuts(N2) = (/ 0.99_rt, 0.98_rt, 0.97_rt, 0.96_rt, 0.95_rt, 0.94_rt, &
!     0.93_rt, 0.92_rt, 0.91_rt, 0.90_rt, 0.89_rt, 0.88_rt, 0.87_rt, 0.86_rt, 0.85_rt, 0.84_rt, 0.83_rt, 0.82_rt, &
!     0.81_rt, 0.80_rt, 0.79_rt, 0.78_rt, 0.77_rt, 0.76_rt, 0.75_rt, 0.73_rt, 0.70_rt, &
!     0.65_rt, 0.60_rt, 0.55_rt, 0.50_rt, 0.45_rt, 0.40_rt, 0.35_rt, 0.30_rt /)
!###########################################################################

!####################################
! * Standard for old 1st-bin reference method
  integer, parameter :: N1=6, N2=1
!  integer, parameter :: mubins(N1) = (/ 6,7,8,9,10,11,12,13,14,15 /)
!  integer, parameter :: mubins(N1) = (/ 10,11,12,13,14,15,16,17,18,19,20 /)
!  integer, parameter :: mubins(N1) = (/ 15,16,17,18,19,20 /)
!  integer, parameter :: mubins(N1) = (/ 21,22,23,24,25,26 /)
  integer, parameter :: mubins(N1) = (/ 20,21,22,23,24,25 /)
!  integer, parameter :: mubins(N1) = (/ 25,26,27,28,29,30 /)
!  integer, parameter :: mubins(N1) = (/ 30,31,32,33,34,35 /)
  real(rt), parameter :: mucuts(N2) = (/ 0.97_rt/) !,0.96_rt, 0.95_rt, 0.94_rt, 0.93_rt, 0.92_rt, 0.91_rt, 0.90_rt, 0.89_rt, 0.88_rt, 0.87_rt, 0.86_rt, 0.85_rt /) 
!  real(rt), parameter :: mucuts(N2) = (/  0.97_rt,0.96_rt,0.95_rt,0.94_rt,0.93_rt  /)
!  real(rt), parameter :: mucuts(N2) = (/  0.90_rt,0.89_rt,0.88_rt,0.87_rt,0.86_rt  /)
!  real(rt), parameter :: mucuts(N2) = (/  0.97_rt /)
!####################################
!####################################
! * Standard for bigcov
!  integer, parameter :: N1=6, N2=11
!  integer, parameter :: mubins(N1) = (/ 15,16,17,18,19,20 /)
!  real(rt), parameter :: mucuts(N2) = (/  0.95_rt,0.94_rt,0.93_rt,0.92_rt,0.91_rt,0.90_rt,0.89_rt,0.88_rt,0.87_rt,0.86_rt,0.85_rt /) 
!####################################

  !integer, parameter :: N1=11
  !integer, parameter :: mubins(N1) = (/ 10,11,12,13,14,15,16,17,18,19,20 /)

!  integer, parameter :: mubins(N1) = (/ 11 /)
!  integer, parameter :: mubins(N1) = (/ 10,11,12,13,14,15 /)
!  integer, parameter :: mubins(N1) = (/ 10*gbnumsbin,11*gbnumsbin,12*gbnumsbin,13*gbnumsbin,14*gbnumsbin,15*gbnumsbin /)
!  integer, parameter :: mubins(N1) = (/ 35,36,37,38,39,40 /)
  !integer, parameter :: mubins(N1) = 10 !(/ 10,11,12,13,14,15,16,17,18,19,19 /)
!  integer, parameter :: mubins(N1) = (/ 10,11,12,13,14,15 /)
!  integer, parameter :: mubins(N1) = (/ 25,26,27,28,29,30 /)
!  integer, parameter :: mubins(N1) = (/ 10,11,12,13,14,15,16,17,18,19,20 /)

!  integer, parameter :: N1=6
!  integer, parameter :: mubins(N1) = (/ 20,21,22,23,24,25 /)
!  integer, parameter :: mubins(N1) = (/ 26,27,28,29,30 /)



!! Impose a cut mu < 0.97
!!   Fiber collision effect & Strong FoG at mu > 0.97
!!!   Checked that results stable in the range of 0.85-0.99 
  !integer, parameter :: N2=7
  !real(rt), parameter :: mucuts(N2) = (/  0.99_rt, 0.95_rt, 0.90_rt, 0.85_rt, 0.80_rt, 0.75_rt, 0.70_rt  /) 

  ! Compare i-th bin with the (i-1) th bin (rather than choosing one specific bin as the reference...)
  logical :: NBComp = .true., NBComp_Simp = .false.
 
!!! ABANDON !!!
  logical, parameter :: output_sep_schemes = .true.
  logical, parameter :: output_chisq_eachredbin = .false.
! In total we have N1*N2 schemes with different (#-bin, mucut)

!! Path of covmat files
!!  Setting: directory for covmat files, chisq files
  !character(len=charlen), parameter :: covmatdir = '/home/xiaodongli/LSS/2PCF_AP/covmat_files/', &
  !  chisqdir = '/home/xiaodongli/LSS/2PCF_AP/chisqs/'
  character(len=charlen), parameter :: covmatdir = '/home/xiaodongli/software/APLike/covmat_files/', &
    chisqdir = '/home/xiaodongli/software/APLike/chisqs/'
!  character(len=charlen), parameter :: covmatdir = '/home/xue/workspace/APLike/covmat_files/', &
!    chisqdir = '/home/xue/workspace/APLike/chisqs/'
    
! Setting: directory of the 2pCF (2-point correlation funtion) files 
!! Path of 2-ponit correlation function (2PCF) files, for covmat 
  character(len=charlen), parameter :: covmatfiledir = '/home/xiaodongli/LSS/boss2pcf/data/'
!! Path of 2-ponit correlation function (2PCF) files, for systematic correction
  character(len=charlen), parameter :: syscorfiledir = '/home/xiaodongli/software/APLike/'
  character(len=charlen), parameter :: data2pcffiledir = '/home/xiaodongli/software/APLike/'
!  character(len=charlen), parameter :: filedir = '/home/xue/workspace/APLike/'

!! Settings of 2pCF (rmax: maximal s; nbins: number of binning in s; mubins: number of binning in mu)
!!  data: observational data; sys: Horizon Run 4 N-body mock for systematic correction; cov: Multidark-Patchy mocks for covariance estimation
  integer :: &
                        smax_database=150, nbins_database=750, mubins_database=600, & ! 2pCF of observational data, in "baseline" cosmologies
                        !smax_database=150, nbins_database=150, mubins_database=120, & !use this if using mock_IO_test!!!
!  integer, parameter :: smax_database=150, nbins_database=1200, mubins_database=960, & ! 2pCF of observational data, in "baseline" cosmologies
                        smax_data=51, nbins_data=51, mubins_data=120, & ! 2pCF of observational data, in general cosmologies
                        smax_sysmock=150, nbins_sysmock=150, mubins_sysmock=120, & ! 2pCF from Horizon Run 4 mocks, used for estimation of covmat 
                        smax_covmock=51,  nbins_covmock=51, mubins_covmock=120, &  ! 2pCF from Patchy mocks, used for systematic correction
                        ncovmocks = 1900, nsysmocks = 4 ! number of mocks used for covmat-estimation and systematic-correction
                        !ncovmocks = 2000, nsysmocks = 4 ! number of mocks used for covmat-estimation and systematic-correction
!                        ncovmocks = 50, nsysmocks = 4 ! number of mocks used for covmat-estimation and systematic-correction

! Setting: Integration limits of s
  real(rt), parameter :: ints1 = 6.0_rt, ints2 = 40.0_rt
! Settings: polynomical fitting degeree for dintxi_sys
  integer, parameter :: polyfitdeg = 1

  ! Many settings about IO tests...
  logical :: mock_IO_test = .false.
  integer :: gb_mock_IO_test_id = 1 !!! We shall do self IO test...
  integer, parameter :: gb_J08_mock = 1, gb_PT_mock = 2, gb_mock_type = gb_PT_mock
  integer :: gb_syscor_mock = gb_J08_mock
  character(len=4) :: gb_PT_IO_id = '0000'
  real(rt) :: gb_om_fiducial = 0.26_rt, gb_w_fiducial = -1.0_rt
  logical :: mock_IO_test_usemock3=.false.
  logical :: dataself_IO_test = .false. !!! Using real data itself for systematic correction.
                                 ! Obtained result should be consistent with 
                                 !  "the fiducial cosmology in which the data-2pcf is measured".
                                 ! The option which can recover the fiducial cosmology is considered as the 
                                 !  "the optimised option for the set of data used for analysis" 
                                 ! Usually we shall set polyfitdeg as 1
                                 !  or any kind of option would be able to recover the fiducial
                        
! Settings: effectove redshifts of the redshift bins
  integer, parameter :: nz = 6                        
  real(rt), parameter :: zeffs(nz) = (/ 0.2154098242_rt,  &
  0.3156545036_rt, 0.3869159269_rt, &
  0.4785988874_rt, 0.5408209467_rt, 0.6186561000_rt  /)
  integer, parameter :: ngals(nz) = (/ 82746 +37841,82745+37838 ,82745+37844 , 188227+68966,188206+68973,188227+68968 /)
  
  ! forecast settings 
  logical :: ForeCast = .false., gb_fc_spmat=.false.! 

  !!! Settings for DESI

  real(rt), parameter :: zeffsDESI(29) = (/ 0.05d0, 0.15d0, 0.25d0, 0.35d0, 0.45d0, &
        0.65d0, 0.75d0, 0.85d0, 0.95d0, 1.05d0, 1.15d0, 1.25d0, 1.35d0, 1.45d0, 1.55d0, 1.65d0, 1.75d0, 1.85d0, &
        1.96d0, 2.12d0, 2.28d0, 2.43d0, 2.59d0, 2.75d0, 2.91d0, 3.07d0, 3.23d0, 3.39d0, 3.55d0 /)
  real(rt), parameter :: ngalsDESI(29) = (/ 1631000.0d0, 4303600.0d0, 2672600.0d0, 1024800.0d0, 168000.0d0, &
        1663200.0d0, 4634000.0d0, 3704400.0d0, 3406200.0d0, 2189600.0d0, 2024400.0d0, 1983800.0d0, 848400.0d0, 771400.0d0, 582400.0d0, 298200.0d0,&
        121800.0d0, 120400.0d0, 183680.0d0, 154559.0d0, 115010.0d0, 93310.0d0, 82879.0d0, 69440.0d0, 58239.0d0, 47039.0d0, 35840.0d0,  &
        29119.0d0, 20159.0d0 /)


  ! DESI BGS
  character(len=100000) :: gb_fcname = 'DESI_BGS'
  integer, parameter :: nzFC = 5
  real(rt), parameter :: zeffsFC(nzFC) = (/ 0.05d0, 0.15d0, 0.25d0, 0.35d0, 0.45d0/) !&
  real(rt), parameter :: ngalsFC(nzFC) = (/ 1631000.0d0, 4303600.0d0, 2672600.0d0, 1024800.0d0, 168000.0d0 /)

  ! DESI Main
!  character(len=100000) :: gb_fcname = 'DESI_main'
!  integer, parameter :: nzFC = 13
!  real(rt), parameter :: zeffsfc(nzfc) = zeffsDESI(6:18), ngalsfc(nzfc) = ngalsdesi(6:18)


  ! DESI BGS+main
!  character(len=100000) :: gb_fcname = 'DESI_BGS_main'
!  integer, parameter :: nzFC = 18
!  real(rt), parameter :: zeffsfc(nzfc) = zeffsDESI(1:nzfc), ngalsfc(nzfc) = ngalsdesi(1:nzfc)

  ! DESI first nzfc bins
!  character(len=100000) :: gb_fcname = ''
!  integer, parameter :: nzFC = 17
!  real(rt), parameter :: zeffsfc(nzfc) = zeffsDESI(1:nzfc), ngalsfc(nzfc) = ngalsdesi(1:nzfc)


  ! DESI Main, drop low 1
!  character(len=100000) :: gb_fcname = 'DESI_maindroplow1'
!  integer, parameter :: nzFC = 12
!  real(rt), parameter :: zeffsfc(nzfc) = zeffsDESI(7:18), ngalsfc(nzfc) = ngalsdesi(7:18)
  ! DESI Main, drop up 1
!  character(len=100000) :: gb_fcname = 'DESI_maindropup1'
!  integer, parameter :: nzFC = 12
!  real(rt), parameter :: zeffsfc(nzfc) = zeffsDESI(6:17), ngalsfc(nzfc) = ngalsdesi(6:17)
  ! DESI Main, add low 1
!  character(len=100000) :: gb_fcname = 'DESI_mainaddlow1'
!  integer, parameter :: nzFC = 14
!  real(rt), parameter :: zeffsfc(nzfc) = zeffsDESI(5:18), ngalsfc(nzfc) = ngalsdesi(5:18)
  ! DESI Main, add up 1
!  character(len=100000) :: gb_fcname = 'DESI_mainaddup1'
!  integer, parameter :: nzFC = 14
!  real(rt), parameter :: zeffsfc(nzfc) = zeffsDESI(6:19), ngalsfc(nzfc) = ngalsdesi(6:19)


  ! DESI QSO
!  integer, parameter :: nzFC = 11
!  real(rt), parameter :: zeffsfc(nzfc) = zeffsDESI(19:29), ngalsfc(nzfc) = ngalsdesi(19:29)


  ! DESI 6 bin
!  integer, parameter :: nzFC = 6
!  real(rt), parameter :: zeffsFC(nzFC) = (/ 0.05d0, 0.15d0, 0.25d0, 0.35d0, 0.45d0, 0.55d0 /) !&
!  real(rt), parameter :: ngalsFC(nzFC) = (/ 1631000.0d0, 4303600.0d0, 2672600.0d0, 1024800.0d0, 168000.0d0, 1663200.0d0 /)

  ! DESI 7 bin
!  integer, parameter :: nzFC = 7
!  real(rt), parameter :: zeffsFC(nzFC) = (/ 0.05d0, 0.15d0, 0.25d0, 0.35d0, 0.45d0, 0.65d0, 0.75d0 /) 
!  real(rt), parameter :: ngalsFC(nzFC) = (/ 1631000.0d0, 4303600.0d0, 2672600.0d0, 1024800.0d0, 168000.0d0, 1663200.0d0, 4634000.0d0 /)
!, 3704400.0d0, 3406200.0d0, 2189600.0d0, 2024400.0d0, 1983800.0d0, 848400.0d0, 771400.0d0, 582400.0d0, 298200.0d0, 


!  integer, parameter :: nzFC = 29
!  real(rt), parameter :: zeffsFC(nzFC) = (/ 0.05d0, 0.15d0, 0.25d0, 0.35d0, 0.45d0, &
!        0.65d0, 0.75d0, 0.85d0, 0.95d0, 1.05d0, 1.15d0, 1.25d0, 1.35d0, 1.45d0, 1.55d0, 1.65d0, 1.75d0, 1.85d0, &
!        1.96d0, 2.12d0, 2.28d0, 2.43d0, 2.59d0, 2.75d0, 2.91d0, 3.07d0, 3.23d0, 3.39d0, 3.55d0 /)
!  real(rt), parameter :: ngalsFC(nzFC) = (/ 1631000.0d0, 4303600.0d0, 2672600.0d0, 1024800.0d0, 168000.0d0, &
!        1663200.0d0, 4634000.0d0, 3704400.0d0, 3406200.0d0, 2189600.0d0, 2024400.0d0, 1983800.0d0, 848400.0d0, 771400.0d0, 582400.0d0, 298200.0d0, 121800.0d0, 120400.0d0, &
!        183680.0d0, 154559.0d0, 115010.0d0, 93310.0d0, 82879.0d0, 69440.0d0, 58239.0d0, 47039.0d0, 35840.0d0, 29119.0d0, 20159.0d0 /)

  !!! Settings for BOSS
!  integer, parameter :: nzFC = 6                        
!  real(rt), parameter :: zeffsFC(nzFC) = (/ 0.2154098242_rt,  &
!  0.3156545036_rt, 0.3869159269_rt, &
!  0.4785988874_rt, 0.5408209467_rt, 0.6186561000_rt  /)
!  real(rt), parameter :: ngalsFC(nzFC) = (/ 82746 +37841d0,82745+37838d0 ,82745+37844d0 , 188227+68966d0,188206+68973d0,188227+68968d0 /)
!1.6310e+06, 4.3036e+06d0, 2.6726e+06d0, 1.0248e+06d0, 1.6800e+05d0, &
       !1.6632e+06d0, 4.6340e+06d0, 3.7044e+06d0, 3.4062e+06d0, 2.1896e+06d0, 2.0244e+06d0, 1.9838e+06d0, 8.4840e+05d0, 7.7140e+05d0, 5.8240e+05d0, 2.9820e+05d0, 1.2180e+05d0, 1.2040e+05d0, &
 !      1.8368e+05d0, 1.5456e+05d0, 1.1501e+05d0, 9.3310e+04d0, 8.2880e+04d0, 6.9440e+04d0, 5.8240e+04d0, 4.7040e+04d0, 3.5840e+04d0, 2.9120e+04d0, 2.0160e+04d0 /)
  !integer, parameter :: ngalsFC(nzFC) = 



  integer, parameter :: exclude_bin = -1 ! exclude one redshift bin from the analysis!!!
  

contains

      SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA )
      CHARACTER*(*)    DESC
      INTEGER          M, N, LDA
      DOUBLE PRECISION A( LDA, * )
      INTEGER          I, J
      WRITE(*,*)
      WRITE(*,*) DESC
      DO I = 1, M
         WRITE(*,9998) ( A( I, J ), J = 1, N )
      END DO
 9998 FORMAT( 11(:,1X,e14.7) )
      RETURN
      END SUBROUTINE PRINT_MATRIX


  !------------------------------------------
  ! solve a linear system (compute eigen values & vectors)
  !------------------------------------------ 
      subroutine dsyev_linsys (Minput, n, EVals, EVecs )

      integer, intent(in) :: n
      real(rt), intent(in) :: Minput(n,n)
      real(rt), intent(out) :: EVals(n), EVecs(n,n)

      real(rt), allocatable :: Work(:)
!      real(rt), allocatable :: A(:,:),  WORK(:), B(:,:), lams(:,:), C(:,:), AT(:,:)
      real(rt) :: tmpx, normA
      integer, parameter :: LWMAX = 1000
      !W( N ), WORK( LWMAX )
      integer :: i, j, INFO, LWORK

      EVecs = Minput

      LWORK = -1
      allocate(work(LWMAX))
      CALL DSYEV( 'Vectors', 'Upper', n, EVecs, n, EVals, WORK, LWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      CALL DSYEV( 'Vectors', 'Upper', n, EVecs, N, EVals, WORK, LWORK, INFO )


      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'(dsyev_linsys) The algorithm failed to compute eigenvalues.'
         STOP
      END IF

      end subroutine dsyev_linsys

  !------------------------------------------
  ! reconstruct A from its Evals and Evecs
  !------------------------------------------ 
      subroutine reconmat_linsys (EVals, EVecs, n, Moutput )
      integer, intent(in) :: n
      real(rt), intent(in) ::EVals(n), EVecs(n,n)
      real(rt), intent(out) :: Moutput(n,n)

      real(rt) :: lams(n,n), EVecsT(n,n)
      integer :: i,j

      lams = 0
      do i = 1, n
        lams(i,i) = EVals(i)
        do j = 1, n
          EVecsT(j,i) = EVecs(i,j)
        enddo
      enddo

      Moutput = matmul(EVecs, matmul(lams, EVecsT))
      end subroutine reconmat_linsys 

  !------------------------------------------
  ! construct a sp mat by abondoning negative evals
  !------------------------------------------ 
   ! result saved in a new matrix
      subroutine spmat_linsys(Minput, Moutput, n)
      integer, intent(in) :: n
      real(rt), intent(in) :: Minput(n,n)
      real(rt), intent(out) :: Moutput(n,n)
      real(rt) :: Evals(n), EVecs(n,n), spEvals(n)
      call dsyev_linsys (Minput, n, EVals, EVecs )
      call spvec(Evals,n,spEvals)
      call reconmat_linsys (spEvals, EVecs, n, Moutput)
      end subroutine spmat_linsys
   ! result saved to the same matrix
      subroutine spmat_linsys_inout(M, n)
      integer, intent(in) :: n
      real(rt), intent(inout) :: M(n,n)
      real(rt) :: Evals(n), EVecs(n,n), spEvals(n)
      call dsyev_linsys (M, n, EVals, EVecs )
!      call spvec(Evals,n,spEvals)
      spEvals = abs(Evals)
      call reconmat_linsys (spEvals, EVecs, n, M)
!      M(1,1) = M(1,1)
      end subroutine spmat_linsys_inout

  !------------------------------------------
  ! semi-positive vector
  !------------------------------------------ 
      subroutine spvec(V,n,Vout)
        integer, intent(in) :: n
        real(rt), intent(in) :: V(n)
        real(rt), intent(out) :: Vout(n)
        integer :: i
        Vout = V
        do i = 1, n
          if(Vout(i).lt.0) Vout(i) = 0 !abs(V(i))
        enddo
      end subroutine spvec

  !------------------------------------------
  ! test of the subroutine dsyev (linear system)
  !------------------------------------------ 
      subroutine dsyev_test ( )

      real(rt), allocatable :: A(:,:), W(:), WORK(:), B(:,:), lams(:,:), C(:,:), AT(:,:)
      real(rt) :: tmpx, normA
      integer, parameter :: LWMAX = 1000
      !W( N ), WORK( LWMAX )
      integer :: i, j, N
      integer :: INFO, LWORK

      open(unit=1000,file='bigcov1.txt')
      n = 0
      do while(.true.)
      read(1000,*,end=100) tmpx
      n=n+1
      cycle
100   exit
      enddo
      close(1000)

      allocate(A(n,n), B(n,n), lams(n,n), W(n), work(LWMAX), C(n,n), AT(n,n))
      open(unit=1001,file='bigcov1.txt')
      do i = 1, n
      read(1001,*) A(i,:)
      B(i,:) = A(i,:)
      enddo
      close(1001)
      !print *, INFO

      LWORK = -1
      CALL DSYEV( 'Vectors', 'Upper', N, A, N, W, WORK, LWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      CALL DSYEV( 'Vectors', 'Upper', N, A, N, W, WORK, LWORK, INFO )


      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         STOP
      END IF

      CALL PRINT_MATRIX( 'Eigenvalues', 1, N, W, 1 )
      do i = 1, n
      normA = 0
      do j = 1, n
      normA = normA + A(j,i)**2.0
      enddo
      print *, i, normA
      !CALL PRINT_MATRIX( 'Eigenvectors (stored columnwise)', N, 1, A(:,i), N)
      enddo


      lams = 0.0
      do i = 1, n
      lams(i,i) = W(i)
      enddo
      do i = 1, n
      do j = 1, n
      AT(i,j) = A(j,i)
      enddo
      enddo
      C = matmul(AT,matmul(lams,A))

      ! we can get AT . B. A = lams
      C = matmul(AT,matmul(B,A))
      !C = matmul(AT,A)
      do i = 1, n
      print *, i
      print *, C(i,:)
      !print *, B(i,:)
      print *, lams(i,:)
!      stop
      enddo

      ! Sp we can get A . lams. AT = B!
      C = matmul(A,matmul(lams,AT))
      !C = matmul(AT,A)
      do i = 1, n
      print *, i
      print *, C(i,:)
      print *, B(i,:)
      !print *, lams(i,:)
!      stop
      enddo

      ! 
      end subroutine dsyev_test


  !------------------------------------------
  ! Invert a matrix
  !------------------------------------------ 
  subroutine nizhen(aa,b,n)
    ! Arguments
    real(rt), intent(in) :: aa(n,n)
    integer, intent(in) :: n
    real(rt), intent(out) :: b(n,n)
    ! Local
    integer :: i,j,k
    real(rt) :: a(n,n)
    a=aa
    b=0.0_rt
    do i=1,n
      b(i,i)=1
    enddo
    do i=1,n
      b(i,:)=b(i,:)/a(i,i)
      a(i,i:n)=a(i,i:n)/a(i,i) 
      do j=i+1,n     
        do k=1,n   
          b(j,k)=b(j,k)-b(i,k)*a(j,i)
        enddo
        a(j,i:n)=a(j,i:n)-a(i,i:n)*a(j,i)
      enddo
    enddo
    do i=n,1,-1
      do j=i-1,1,-1
        do k=1,n
          b(j,k)=b(j,k)-b(i,k)*a(j,i)
        enddo
      enddo
    enddo
  end subroutine nizhen
  !------------------------------------------
  ! Polynomial Regression:
  !  polynomial fitting to data points
  ! Y(i) = A(1) + A(2)*X(1) + A(3)*X(2)^2 + ... 
  !        + A(n+1)*X(n)^n
  !------------------------------------------          
        subroutine poly_fit(X,Y,A,ndat,n)
                ! Dummy
                real(rt), intent(in) :: X(ndat),Y(ndat)
                integer, intent(in) :: ndat,n
                real(rt), intent(out) :: A(n+1)
                ! Local
                real(rt) :: CapX(ndat,n+1), CapMatA(n+1,n+1), CapMatB(n+1,n+1), CapMatC(n+1,ndat)
                integer :: i,j,k
                ! CapX(i,j) = X(i) ^ (j-1)
                do i = 1, ndat
                do j = 1, n+1
                        CapX(i,j) = X(i)**(j-1)
                enddo
                enddo
                ! CapMatA = ( CapX^T CapX)
                do i = 1,n+1
                do j = 1,n+1
                        CapMatA(i,j) = 0.0_rt
                        do k = 1,ndat
                                CapMatA(i,j) = CapMatA(i,j) + CapX(k,i)*CapX(k,j)
                        enddo
                enddo
                enddo
                ! CapMatB = ( CapX^T CapX )^(-1)
                call nizhen(CapMatA,CapMatB,n+1)
                ! CapMatC = ( CapX^T CapX )^(-1) CapX^T
                do i = 1, n+1
                do j = 1, ndat
                        CapMatC(i,j) = 0.0_rt
                        do k = 1, n+1
                                CapMatC(i,j) = CapMatC(i,j) + CapMatB(i,k)*CapX(j,k)
                        enddo
                enddo
                enddo
                ! A = ( CapX^T CapX )^(-1) CapX^T Y
                do i = 1, n+1
                        A(i) = 0.0_rt
                        do j = 1, ndat
                                A(i) = A(i) + CapMatC(i,j)*Y(j)
                        enddo
                enddo
          end subroutine poly_fit
  !------------------------------------------
  ! Value of a n-th polynomial 
  !  at some value of x
  !------------------------------------------           
          real(rt) function poly(x,A,n)
                  ! Dummy
                  real(rt), intent(in) :: x, A(n+1)
                  integer, intent(in) :: n
                  ! local
                  integer :: i
                  poly = A(1)
                  do i = 1, n
                          poly = poly + A(i+1)*x**(i)
                  enddo
          end function poly
  !------------------------------------------
  ! polynomial regression of curve Y
  !------------------------------------------
   function polyfitY(X,Y,n,polyfitdeg)
     integer, intent(in) :: n, polyfitdeg
     real(rt), intent(in) :: X(n), Y(n)
     real(rt) :: polyfitY(n), coeff(polyfitdeg+1)
     integer :: i
     call poly_fit(X,Y,coeff,n,polyfitdeg)
     !print *, coeff
     do i = 1, n
       polyfitY(i) = poly(X(i),coeff,polyfitdeg)
     enddo
   end function polyfitY
  !------------------------------------------
  ! correction factor for covmat due to 
  !  finite number of bins/mocks
  !------------------------------------------
   subroutine Percival_cs(Ns,Nb,Npar,factD,factm1,factm2)
                !'''
                !Corretion to covariance to account for errors in deriving covmat from simulations
                !Ns: number of simulations
                !Nb: number of 'bands' in power spectrum or 2pCF
                !Npar: nuber of parameters (in likelihood estimation)
                !'''
                integer, intent(in) :: Ns,Nb,Npar
                real(rt), intent(out) :: factD,factm1,factm2
                real(rt) :: factA,factB
                factD = (Nb+1.0_rt)/(Ns-1.0)
                factA = 2.0_rt / (Ns-Nb-1.0) / (Ns-Nb-4.0)
                factB = (Ns-Nb-2.0_rt) / (Ns-Nb-1.0) / (Ns-Nb-4.0)
                factm1 = (1.0_rt+factB*(Nb-Npar)) / (1.0+factA+factB*(Npar+1))
                factm2 = 1.0_rt/(1.0-factD) * factm1
                !return D, m1, m2
  end subroutine Percival_cs
  !------------------------------------------
  ! Find out the four index of the N smallest number
  !------------------------------------------
  subroutine FindNMinIndice(A,nA,imins,nmins)
    integer, intent(in) :: nA, nmins
    real(rt), intent(in) :: A(nA)
    integer, intent(out) :: imins(nmins)
    real :: Bmax, B(nA)
    integer :: imin,i,j
    do i = 1, nA
      B(i) = A(i)
    enddo
    Bmax = maxval(B)
    do i = 1, nmins
      imins(i:i) = minloc(B)
!      print *, minloc(B), B(minloc(B))
      B(imins(i)) = Bmax
    enddo
!    print *, Amax, A(maxloc(A))
!    imins(1:nmins) = 1
  end subroutine FindNMinIndice
  
  
end module AP_type_constants
!###############################################################

!###############################################################
!## Module: Cosmological functions (DA and H)
module AP_cosmo_related

use AP_type_constants

implicit none


contains

!---------------------------------------------------------------
! Simpson integration for some cosmological function 
  real(rt) function CosmoSimpson(nowpar,fun,xleft,xright,N)
    real(rt), EXTERNAL :: fun
    real(rt), intent(in) :: xleft, xright
    integer, intent(in) :: N
    real(rt) :: x1,x2,BC,f1,f2
    integer :: i
    type(omwpar) :: nowpar
    BC=(xright-xleft)/dble(N)
    x1=xleft;x2=x1+BC;
    f1=fun(nowpar,x1);f2=fun(nowpar,x2);CosmoSimpson=(f1+fun(nowpar,(x1+x2)*0.5d0)*4.0d0+f2)*BC/6.0d0;
    do i = 2,N
      x1=x2;f1=f2;x2=x2+BC;
      f2=fun(nowpar,x2);CosmoSimpson=CosmoSimpson+(f1+fun(nowpar,(x1+x2)*0.5d0)*4.0d0+f2)*BC/6.0d0
    enddo
  end function CosmoSimpson
!----------------------------------------------------------------
! inv_ez = 1/E(z); E(z) = H(z)/H0
  real(rt) function inv_ez(nowpar,z)
    type(omwpar) :: nowpar
    real(rt) :: z
    !inv_ez =  1. / dsqrt(dble(nowpar%omegam*(1.+z)**3.  &
    !  + (1.-nowpar%omegam)*(1.0+z)**(3.*(1.+nowpar%w))))
    inv_ez =  1. / dsqrt(dble(nowpar%omegam*(1.+z)**3. + nowpar%omegak*(1.+z)**2. &
      + (1.-nowpar%omegam-nowpar%omegak)*exp(-3.0*nowpar%wa*z/(1.0+z))*(1.0+z)**(3.*(1.+nowpar%w+nowpar%wa))))
  end function inv_ez 
!----------------------------------------------------------------
! rz = \int 1/H(z) dz, in unit of Mpc/h
  real(rt) function rz(nowpar,z)
    type(omwpar) :: nowpar
    real(rt), intent(in) :: z
    ! use at least 128 bins in redshift space; 128 bins for every redshift interval of 0.25
    rz = CosmoSimpson(nowpar,inv_ez,0.0_rt,z,N=max(128,128*ceiling(z/0.25_rt)))*CONST_C/100.d0
  end function rz
!----------------------------------------------------------------
  real(rt) function fk(distance, omegak)
   real(rt) :: distance, sqrtOk0, omegak
   sqrtOk0 = sqrt(abs(omegak))
   if(omegak < 0.0d0) then
     fk = sin(sqrtOk0*distance)/sqrtOk0
   elseif(omegak > 0.0d0) then
     fk = sinh(sqrtOk0*distance)/sqrtOk0
   else
     fk = distance
   endif
  end function fk
! DA(z) in unit of Mpc/h, in w-cdm model
  real(rt) function DAz_wcdm(nowpar,z)
    type(omwpar) :: nowpar
    real(rt), intent(in) :: z
    real(rt) ::  Intinve, H0fact
    H0fact = 100.0 / CONST_C
    Intinve = rz(nowpar,z) * H0fact
    DAz_wcdm = fk(Intinve, nowpar%omegak) / H0fact / (1.0+z)
    !DAz_wcdm = rz(nowpar,z) / (1.0+z)
    !DAz_wcdm = rz(nowpar,z) / (1.0+z)
    !print *, DAz_wcdm, rz(nowpar,z) / (1.0+z), nowpar%omegak
  end function DAz_wcdm
!----------------------------------------------------------------
! H(z), in unit of km/s/(Mpc/h), in w-cdm model
  real(rt) function Hz_wcdm(nowpar,z)
    type(omwpar) :: nowpar
    real(rt), intent(in) :: z
    Hz_wcdm=100./inv_ez(nowpar,z)
  end function Hz_wcdm
end module AP_cosmo_related



module redbin_weights_bossdr12_6bins
use AP_type_constants
implicit none


  real(rt) :: dr12v4_cmass_N_numgal(3,1000),dr12v4_cmass_S_numgal(3,1000),&
            dr12v4_lowz_N_numgal(3,1000),dr12v4_lowz_S_numgal(3,1000)
  character(charlen) :: dr12v4_cmass_N_numgal_file = '../datfiles/nbar-cmass-dr12v4-N-Reid.dat',&
                      dr12v4_cmass_S_numgal_file = '../datfiles/nbar-cmass-dr12v4-S-Reid.dat',&
                      dr12v4_lowz_N_numgal_file = '../datfiles/nbar-lowz-dr12v4-N-Reid.dat',&
                      dr12v4_lowz_S_numgal_file = '../datfiles/nbar-lowz-dr12v4-S-Reid.dat'
  real(rt), parameter :: numgal_zmin=0.0, numgal_zmax=1.0
  integer, parameter  :: numgal_nbin = 200
  real(rt), parameter :: numgal_dz = (numgal_zmax-numgal_zmin) / dble(numgal_nbin)
  real(rt) :: numgals(numgal_nbin), numgal_zcenters(numgal_nbin), redbin_weights(6,numgal_nbin)
  logical :: numgal_inited = .false.
  real(rt) :: redbin_edges(7) = (/ 0.15_rt,  0.2741_rt, 0.3510_rt, 0.43_rt, 0.51_rt,  0.5720_rt, 0.6929_rt  /)

  contains

  subroutine numgal_init()
    integer  :: i,j,iredbin,i1,i2
    character:: tmpchar
    real(rt) :: zcen,zlow,zhigh,nbar,wfkp,shell_vol,numgal
    
    numgals = 0.0_rt
    do i = 1, 4
      if(i.eq.1) open(82785,file=dr12v4_lowz_N_numgal_file,action='read')
      if(i.eq.2) open(82785,file=dr12v4_lowz_S_numgal_file,action='read')
      if(i.eq.3) open(82785,file=dr12v4_cmass_N_numgal_file,action='read')
      if(i.eq.4) open(82785,file=dr12v4_cmass_S_numgal_file,action='read')
      do j = 1, 2
        read(82785,*) tmpchar
      enddo
      do j = 1, numgal_nbin
        read(82785,*)zcen,zlow,zhigh,nbar,wfkp,shell_vol,numgal
        if(zcen .gt. 0.43 .and. (i.eq.1.or.i.eq.2)) cycle
        if(zcen .lt. 0.43 .and. (i.eq.3.or.i.eq.4)) cycle
        numgals(j) = numgals(j) + numgal
      enddo
    enddo

    print *, 'Finishing initialization of #-gal.'
    do j = 1, numgal_nbin
      numgal_zcenters(j) = numgal_zmin+numgal_dz*(j-0.5)
      !print *, real(numgal_zmin+numgal_dz*(j-1)), real(numgal_zmin+numgal_dz*j), numgals(j)
    enddo

    redbin_weights = 0.0_rt
    do iredbin = 1, 6
      zlow=redbin_edges(iredbin); zhigh = redbin_edges(iredbin+1)
      i1=floor(zlow/numgal_dz+0.0001)+1; i2=ceiling(zhigh/numgal_dz-0.00001)
      do i = i1+1,i2-1
        redbin_weights(iredbin,i) = numgals(i)
      enddo
      redbin_weights(iredbin,i1) = (numgal_zmin+(i1)*numgal_dz - zlow) / numgal_dz * numgals(i1)
      redbin_weights(iredbin,i2) = (zhigh - (numgal_zmin+(i2-1)*numgal_dz)) / numgal_dz * numgals(i2)
      print *, 'iredbin = ', iredbin
!      print *, zlow, numgal_zmin+(i1)*numgal_dz,  (numgal_zmin+(i2-1)*numgal_dz), zhigh
      print *, 'weights in redshift bins: '
      do i = 1, numgal_nbin
        if (redbin_weights(iredbin,i) .ge. 0.00001) &
          print *,  numgal_zcenters(i), real(redbin_weights(iredbin,i))
      enddo
    enddo
    numgal_inited = .true.
  end subroutine numgal_init


  

end module redbin_weights_bossdr12_6bins

module AP_2pcf_tools
use AP_cosmo_related
use redbin_weights_bossdr12_6bins
implicit none
  


! Structure storing covmats
  type :: M2
    real(rt), allocatable :: A(:,:)
    integer :: nA=-1
  end type
  type(M2),save :: covmats(N1,N2,nz-1), bigcovmats(N1,N2),  bigcovmatsfc(N1,N2)
  real(rt) :: dintxi_syscor(maxval(mubins(:)), N1,N2,nz-1)
  real(rt), allocatable ::  global_smutabstds(:,:,:,:,:)
     
!----------------------------
! 1. 

contains

  character(25) function omwstr(omegam, w)
    real(rt) :: omegam, w
    character(10) :: str1, str2
    write(str1, '(f8.4)') omegam
    write(str2, '(f8.4)') w
    omwstr = 'om'//trim(adjustl(str1))//'_w'//trim(adjustl(str2))
  end function omwstr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Settings: Name of 2pCF files
!-----------------------------------
! data (baseline)    
  character(len=charlen) function data2pcffile_base(iz, ombase, wbase)
    integer :: iz
    real(rt) :: ombase, wbase
    character(len=charlen) :: filename, filestr, tmpstr1,tmpstr2,tmpstr3,tmpstrID
    
    if(mock_IO_test) then
      if(gb_mock_type .eq. gb_PT_mock) then
              print *, 'WARNING!!!! Using PT as datafile!!!'
              smax_database=100; nbins_database=500; mubins_database=600
              data2pcffile_base = PT2pcffile(iz);    
              gb_om_fiducial =  0.307115_rt; gb_w_fiducial = -1.0_rt
      elseif(gb_mock_type .eq. gb_J08_mock) then
              print *, 'WARNING!!!! Using J08 as datafile!!!'
              gb_om_fiducial = 0.26_rt; gb_w_fiducial = -1.0_rt
              smax_database=150; nbins_database=150; mubins_database=120
              data2pcffile_base = data2pcffile_base_J08(iz);    
      endif
      print *, 'WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      print *, 'WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      print *, 'WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      print *, 'WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      print *, 'WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      print *, 'WARNING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,'(A, f6.2, f6.2)') ' (data2pcffile_base) use om, w fiducial values of ', gb_om_fiducial, gb_w_fiducial
      return
    ! if forecast then just use the first bin!
!    elseif(Forecast) then
!      write(tmpstr1,*) nbins_database; write(tmpstr2,*) mubins_database; write(tmpstr3,*) smax_database
!      filestr = '.rmax'//trim(adjustl(tmpstr3))//'.'//trim(adjustl(tmpstr1))//'rbins.'&
!        //trim(adjustl(tmpstr2))//'mubins.'
!      filename = trim(adjustl(data2pcffiledir))//'/2pcfs/DR12_LOWZ_data.xyzw.1of3.cosmo-converted.'//&
!        trim(adjustl(omwstr(ombase,wbase)))//trim(adjustl(filestr))//'2pcf'
    endif
    
    gb_om_fiducial = 0.26_rt; gb_w_fiducial = -1.0_rt
    write(*,'(A, f6.2, f6.2)') ' (data2pcffile_base) use om, w fiducial values of ', gb_om_fiducial, gb_w_fiducial
    write(tmpstr1,*) nbins_database; write(tmpstr2,*) mubins_database; write(tmpstr3,*) smax_database
    filestr = '.rmax'//trim(adjustl(tmpstr3))//'.'//trim(adjustl(tmpstr1))//'rbins.'&
      //trim(adjustl(tmpstr2))//'mubins.'
    if (iz .eq. 1) then
      filename = trim(adjustl(data2pcffiledir))//'/2pcfs/DR12_LOWZ_data.xyzw.1of3.cosmo-converted.'//&
        trim(adjustl(omwstr(ombase,wbase)))//trim(adjustl(filestr))//'2pcf'
    elseif (iz .eq. 2) then
      filename = trim(adjustl(data2pcffiledir))//'/2pcfs/DR12_LOWZ_data.xyzw.2of3.cosmo-converted.'//&
        trim(adjustl(omwstr(ombase,wbase)))//trim(adjustl(filestr))//'2pcf'
    elseif (iz .eq. 3) then
      filename = trim(adjustl(data2pcffiledir))//'/2pcfs/DR12_LOWZ_data.xyzw.3of3.cosmo-converted.'//&
        trim(adjustl(omwstr(ombase,wbase)))//trim(adjustl(filestr))//'2pcf'
    elseif (iz .eq. 4) then
      filename = trim(adjustl(data2pcffiledir))//'/2pcfs/DR12_CMASS_data.xyzw.1of3.cosmo-converted.'//&
        trim(adjustl(omwstr(ombase,wbase)))//trim(adjustl(filestr))//'2pcf'
    elseif (iz .eq. 5) then
      filename = trim(adjustl(data2pcffiledir))//'/2pcfs/DR12_CMASS_data.xyzw.2of3.cosmo-converted.'//&
        trim(adjustl(omwstr(ombase,wbase)))//trim(adjustl(filestr))//'2pcf'
    elseif (iz .eq. 6) then
      filename = trim(adjustl(data2pcffiledir))//'/2pcfs/DR12_CMASS_data.xyzw.3of3.cosmo-converted.'//&
        trim(adjustl(omwstr(ombase,wbase)))//trim(adjustl(filestr))//'2pcf'
    endif
    data2pcffile_base = filename
  end function data2pcffile_base
  

  character(len=charlen) function PT2pcffile(iz, IO_id)
    integer :: iz
    character(len=4), intent(in), optional :: IO_id 
    character(len=4) :: now_IO_id
    real(rt) :: ombase, wbase
    character(len=charlen) :: filename, tmpstr
    if(.not. present(IO_id)) then
            now_IO_id = gb_PT_IO_id
    else
            now_IO_id = IO_id
    endif
    if (iz .le. 3) then
      write(tmpstr,*) iz
      filename = &
      '/home/xiaodongli/data/DR12/DR12v4-LOWZ-PATCHY_dense_2pcfs/DR12v4-LOWZ/PatchyV6C.RSD.'//trim(adjustl(now_IO_id))//'.xyzw.'//trim(adjustl(tmpstr))//'of3.rmax100.500rbins.600mubins.fastmiddle.2pcf'
    else
      write(tmpstr,*) iz-3
      filename = &
      '/home/xiaodongli/data/DR12/DR12v4-CMASS-PATCHY_dense_2pcfs/DR12v4-CMASS/PatchyV6C.RSD.'//trim(adjustl(now_IO_id))//'.xyzw.'//trim(adjustl(tmpstr))//'of3.rmax100.500rbins.600mubins.fastmiddle.2pcf'
    endif
!    print *, 'PATCHY as data!!!', trim(adjustl(filename))
    PT2pcffile = filename
  end function PT2pcffile

  
  character(len=charlen) function data2pcffile_base_J08(iz)
    integer :: iz
    real(rt) :: ombase, wbase
    character(len=charlen) :: filename, filestr, tmpstr1,tmpstr2,tmpstr3,tmpstrID
    write(tmpstr1,*) nbins_database; write(tmpstr2,*) mubins_database; write(tmpstr3,*) smax_database
    write(tmpstrID,*) gb_mock_IO_test_ID
    filestr = '.rmax'//trim(adjustl(tmpstr3))//'.'//trim(adjustl(tmpstr1))//'rbins.'&
      //trim(adjustl(tmpstr2))//'mubins.'
    if (iz .eq. 1) then
      filename = trim(adjustl(data2pcffiledir))//'/2pcfs/DR12_LOWZ_J08.RSD.00'//trim(adjustl(tmpstrID))//'.xyzw.1of3'//trim(adjustl(filestr))//'2pcf'
    elseif (iz .eq. 2) then
      filename = trim(adjustl(data2pcffiledir))//'/2pcfs/DR12_LOWZ_J08.RSD.00'//trim(adjustl(tmpstrID))//'.xyzw.2of3'//trim(adjustl(filestr))//'2pcf'
    elseif (iz .eq. 3) then
      filename = trim(adjustl(data2pcffiledir))//'/2pcfs/DR12_LOWZ_J08.RSD.00'//trim(adjustl(tmpstrID))//'.xyzw.3of3'//trim(adjustl(filestr))//'2pcf'
    elseif (iz .eq. 4) then
      filename = trim(adjustl(data2pcffiledir))//'/2pcfs/DR12_CMASS_J08.RSD.00'//trim(adjustl(tmpstrID))//'.xyzw.1of3'//trim(adjustl(filestr))//'2pcf'
    elseif (iz .eq. 5) then
      filename = trim(adjustl(data2pcffiledir))//'/2pcfs/DR12_CMASS_J08.RSD.00'//trim(adjustl(tmpstrID))//'.xyzw.2of3'//trim(adjustl(filestr))//'2pcf'
    elseif (iz .eq. 6) then
      filename = trim(adjustl(data2pcffiledir))//'/2pcfs/DR12_CMASS_J08.RSD.00'//trim(adjustl(tmpstrID))//'.xyzw.3of3'//trim(adjustl(filestr))//'2pcf'
    endif
    data2pcffile_base_J08 = filename
    print *, 'J08 as data!!!', trim(adjustl(filename))
  end function data2pcffile_base_J08
!-----------------------------------
! data (in general, i.e. in different cosmologies)
  character(len=charlen) function data2pcffile_gen(iz, nowpar)
    ! args
    integer :: iz
    type(omwpar) :: nowpar
    ! local
    character(20) :: nowomwstr
    character(len=charlen) :: filename
    nowomwstr = omwstr(nowpar%omegam, nowpar%w)
    if (iz .eq. 1) then
      filename = trim(adjustl(data2pcffiledir))//'DR12v4-LOWZ/xyzw.binsplitted/data.xyzw.1of3.cosmo-converted.'&
                   //trim(adjustl(nowomwstr))//'.rmax51.51rbins.120mubins.2pcf'
    elseif (iz .eq. 2) then
      filename = trim(adjustl(data2pcffiledir))//'DR12v4-LOWZ/xyzw.binsplitted/data.xyzw.2of3.cosmo-converted.'&
                   //trim(adjustl(nowomwstr))//'.rmax51.51rbins.120mubins.2pcf'
    elseif (iz .eq. 3) then
      filename = trim(adjustl(data2pcffiledir))//'DR12v4-LOWZ/xyzw.binsplitted/data.xyzw.3of3.cosmo-converted.'&
                   //trim(adjustl(nowomwstr))//'.rmax51.51rbins.120mubins.2pcf'
    elseif (iz .eq. 4) then
      filename = trim(adjustl(data2pcffiledir))//'DR12v4-CMASS/xyzw.binsplitted/data.xyzw.1of3.cosmo-converted.'&
                   //trim(adjustl(nowomwstr))//'.rmax51.51rbins.120mubins.2pcf'
    elseif (iz .eq. 5) then
      filename = trim(adjustl(data2pcffiledir))//'DR12v4-CMASS/xyzw.binsplitted/data.xyzw.2of3.cosmo-converted.'&
                   //trim(adjustl(nowomwstr))//'.rmax51.51rbins.120mubins.2pcf'
    elseif (iz .eq. 6) then
      filename = trim(adjustl(data2pcffiledir))//'DR12v4-CMASS/xyzw.binsplitted/data.xyzw.3of3.cosmo-converted.'&
                   //trim(adjustl(nowomwstr))//'.rmax51.51rbins.120mubins.2pcf'
    endif
    data2pcffile_gen = filename
  end function data2pcffile_gen
!-----------------------------------
! mock (for systematic correction)
  character(len=charlen) function syscor2pcffile(iz, imock)
    integer :: iz, imock
    character(len=charlen) :: filename
    character(len=15) :: str1
    write(str1, '(i3)') imock-1 
    str1 = 'J08.RSD.00'//trim(adjustl(str1))
    if (iz .eq. 1) then
      filename = trim(adjustl(syscorfiledir))//'/2pcfs/DR12_LOWZ_'//trim(adjustl(str1))//&
                               '.xyzw.1of3.rmax150.150rbins.120mubins.2pcf'
    elseif (iz .eq. 2) then
      filename = trim(adjustl(syscorfiledir))//'/2pcfs/DR12_LOWZ_'//trim(adjustl(str1))//&
                               '.xyzw.2of3.rmax150.150rbins.120mubins.2pcf'
    elseif (iz .eq. 3) then
      filename = trim(adjustl(syscorfiledir))//'/2pcfs/DR12_LOWZ_'//trim(adjustl(str1))//&
                               '.xyzw.3of3.rmax150.150rbins.120mubins.2pcf'
    elseif (iz .eq. 4) then
      filename = trim(adjustl(syscorfiledir))//'/2pcfs/DR12_CMASS_'//trim(adjustl(str1))//&
                               '.xyzw.1of3.rmax150.150rbins.120mubins.2pcf'
    elseif (iz .eq. 5) then
      filename = trim(adjustl(syscorfiledir))//'/2pcfs/DR12_CMASS_'//trim(adjustl(str1))//&
                               '.xyzw.2of3.rmax150.150rbins.120mubins.2pcf'
    elseif (iz .eq. 6) then
      filename = trim(adjustl(syscorfiledir))//'/2pcfs/DR12_CMASS_'//trim(adjustl(str1))//&
                               '.xyzw.3of3.rmax150.150rbins.120mubins.2pcf'
    endif
    syscor2pcffile = filename
  end function syscor2pcffile
!-----------------------------------
! mock for covariance matrix estimation
  character(len=charlen) function cov2pcffile(iz, imock)
    integer :: iz, imock
    character(len=charlen) :: filename
    character(30) :: str1
    write(str1, '(i5)') imock-1 
    if (imock .le. 10) then
       str1 = 'PatchyV6C.RSD.000'//trim(adjustl(str1))
    elseif (imock .le. 100) then
       str1 = 'PatchyV6C.RSD.00'//trim(adjustl(str1))
    elseif (imock .le. 1000) then
       str1 = 'PatchyV6C.RSD.0'//trim(adjustl(str1))
    elseif (imock .le. 10000) then
       str1 = 'PatchyV6C.RSD.'//trim(adjustl(str1))
    endif
    if (iz .eq. 1) then
      filename = trim(adjustl(covmatfiledir))//'DR12v4-LOWZ/xyzw.binsplitted/'//trim(adjustl(str1))//&
                               '.xyzw.1of3.rmax51.51rbins.120mubins.2pcf'
    elseif (iz .eq. 2) then
      filename = trim(adjustl(covmatfiledir))//'DR12v4-LOWZ/xyzw.binsplitted/'//trim(adjustl(str1))//&
                               '.xyzw.2of3.rmax51.51rbins.120mubins.2pcf'
    elseif (iz .eq. 3) then
      filename = trim(adjustl(covmatfiledir))//'DR12v4-LOWZ/xyzw.binsplitted/'//trim(adjustl(str1))//&
                               '.xyzw.3of3.rmax51.51rbins.120mubins.2pcf'
    elseif (iz .eq. 4) then
      filename = trim(adjustl(covmatfiledir))//'DR12v4-CMASS/xyzw.binsplitted/'//trim(adjustl(str1))//&
                               '.xyzw.1of3.rmax51.51rbins.120mubins.2pcf'
    elseif (iz .eq. 5) then
      filename = trim(adjustl(covmatfiledir))//'DR12v4-CMASS/xyzw.binsplitted/'//trim(adjustl(str1))//&
                               '.xyzw.2of3.rmax51.51rbins.120mubins.2pcf'
    elseif (iz .eq. 6) then
      filename = trim(adjustl(covmatfiledir))//'DR12v4-CMASS/xyzw.binsplitted/'//trim(adjustl(str1))//&
                               '.xyzw.3of3.rmax51.51rbins.120mubins.2pcf'
    endif
    cov2pcffile = filename
  end function cov2pcffile
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Load in in the xi(s,mu) table
!   All files are in format of : smin, smax, minimal value of 1-mu, maximal value of 1-mu, DD, DR, RR, three estimators of xi
!   First row in comment
  subroutine ximu_loadsmufile(filename, smutab, smax, numnbins, nummubins, printflag)
    ! argument
    character(len=charlen), intent(in) :: filename
    integer, intent(in) :: smax, numnbins, nummubins
    real(rt), intent(out) :: smutab(numnbins,nummubins,3) ! Four values at each s/mu: DD, DR, RR, xi
    logical, intent(in), optional :: printflag
    ! variables
    character(len=charlen) :: nowstr
    integer :: i, j
    real :: tmpx

    smutab = 0.0_rt
    
    open(unit=44817,file=filename)
    if (.not.present(printflag)) then
            write(*,'(A)') ' (ximu_loadsmufile) Load in 2pCF file: ', trim(adjustl(filename))
    elseif(printflag) then
            write(*,'(A)') ' (ximu_loadsmufile) Load in 2pCF file: ', trim(adjustl(filename))
    endif
    read(44817,*) nowstr
!    print *, 'numnbins, nummubins = ', numnbins, nummubins
!    print *, nowstr
    do i = 1, numnbins
    do j = 1, nummubins
!        print *, i,j
        read(44817,*) tmpx, tmpx, tmpx, tmpx, smutab(i,j,1:3) !, tmpx, tmpx, smutab(i,j,4)
!        print *, smutab(i,j,1:3)
    enddo
    enddo
!    print *, 'End of load in all lines...'
    close(44817)
  end subroutine ximu_loadsmufile

! Compute covmats
  subroutine calc_covmats() !DEBUG c
    integer :: i,j,i1,i2,n, iz,imock, maxmubin
    real(rt), allocatable :: intxis(:,:,:,:,:), tmpX(:), dintxi(:,:)
    real(rt) :: smutab_covmock(nbins_covmock,mubins_covmock,3), deltas, tmpx1, tmpx2, tmpx3
    character(charlen) :: nowfile, tmpstr
    logical :: logvar, printflag, debug_calc_covmats=.false.
    ! 1. Initialize the structure "covmats"
    do i = 1, N1 !DEBUG c
    do j = 1, N2 !DEBUG c
    do iz= 2, nz !DEBUG c
      if (allocated(covmats(i,j,iz-1)%A)) deallocate(covmats(i,j,iz-1)%A) !DEBUG c
      n=mubins(i); covmats(i,j,iz-1)%nA = n-1 !DEBUG c
      allocate(covmats(i,j,iz-1)%A(n-1,n-1)) !DEBUG c
      covmats(i,j,iz-1)%A = 0.0_rt !DEBUG c
    enddo
    enddo
    enddo
    ! 2. Compute intxi, for nz redshift bins, all 2000 mocks
    ! 2.1. Initialize the structure intxis
    maxmubin = maxval(mubins)
    print *, 'maxmubin = ', maxmubin
    allocate(intxis(maxmubin,ncovmocks,nz,N1,N2))
    allocate(dintxi(maxmubin,ncovmocks),tmpX(maxmubin))
    intxis = 0.0_rt; dintxi = 0.0_rt
    ! 2.2 Compute intxis
    do iz = 1, nz
      do imock = 1, ncovmocks
        ! 2.2.1 Check existence of 2pCF file
        nowfile=cov2pcffile(iz,imock)
        inquire(file=nowfile,exist=logvar)
        if (.not.logvar) then 
            print *, ' (calc_covmats) file not found: ', trim(adjustl(nowfile)); stop
        endif
        if(mod(imock,50).eq.0) then
                printflag=.true.
        else
                printflag=.false.
        endif
        ! 2.2.2 Load in 2pCF file
        call ximu_loadsmufile(nowfile, smutab_covmock, smax_covmock, nbins_covmock, mubins_covmock,printflag)
        ! 2.2.3 Compute intxi in the N1*N2 binning schemes
        deltas = smax_covmock / dble(nbins_covmock)
        do i = 1, N1
        do j = 1, N2
          call XiFun(smutab_covmock, deltas, nbins_covmock, mubins_covmock, &
            anglemin=1.0_rt-mucuts(j), anglemax=1.0_rt, &
            smin=ints1, smax=ints2, &
            nummuedge=mubins(i)+1, &
            intxi=tmpX(1:mubins(i)))
          call normfun(tmpX(1:mubins(i)),mubins(i),intxis(1:mubins(i)-1,imock,iz,i,j)) ! normalise the amplitude
          if (debug_calc_covmats.and.mubins(i).eq.25 .and. j.eq.1 .and. iz.eq.1 .and. imock.eq.8) then
             print *, 'nowfile = ', nowfile
             print *, 'intxi  = ', real(tmpX(1:mubins(i)))
             print *, ' python rlt =  [15.089811191950611, 13.224163227234758, ... 16.289372112713142, 16.269104458682122]'
             print *, 'intxi, normed = ', real(intxis(1:mubins(i)-1,imock,iz,i,j))
             print *, ' python rlt =  [1.0191467989709562, 0.89314329057311137, ... 1.1001636292688675, 1.0987947775009022]'
          endif
        enddo
        enddo
      enddo
    enddo
    ! 2.3 Compute covariance matrices
    ! 2.3.1 Loop of schemes & redshifts
    do i = 1, N1
    do j = 1, N2
    do iz= 2, nz
      print *, 'mubin, mucut, iz = ', mubins(i), real(mucuts(j)), iz
      ! 2.3.2 dintxi = intxi@high_z - intxi@low_z; compute it for all covmocks
      n = mubins(i)
      do imock= 1, ncovmocks
       if(.not.NBComp) then
         dintxi(1:n-1,imock) = intxis(1:n-1,imock,iz,i,j)-intxis(1:n-1,imock,1,i,j)
       else
         dintxi(1:n-1,imock) = intxis(1:n-1,imock,iz,i,j)-intxis(1:n-1,imock,iz-1,i,j)
       endif
      enddo
      ! 2.3.3 compute covmat using dintxi
      covmats(i,j,iz-1)%A = 0.0_rt
      do i1= 1,n-1
      do i2 = i1,n-1
        ! mean(X*Y), mean(X), mean(Y)
        tmpx1 = 0.0_rt; tmpx2 = 0.0_rt; tmpx3 = 0.0_rt
        do imock = 1, ncovmocks
          tmpx1 = tmpx1 + dintxi(i1,imock)*dintxi(i2,imock)
          tmpx2 = tmpx2 + dintxi(i1,imock)
          tmpx3 = tmpx3 + dintxi(i2,imock)
        enddo
        tmpx1 = tmpx1 / dble(ncovmocks); tmpx2 = tmpx2 / dble(ncovmocks); tmpx3 = tmpx3 / dble(ncovmocks)
        ! covariance = mean(X*Y) - mean(X)*mean(Y)
        covmats(i,j,iz-1)%A(i1,i2) = (tmpx1 - tmpx2*tmpx3) * dble(ncovmocks) / dble(ncovmocks-1)
      enddo
      enddo
      ! symmetric matrix
      do i1=1,n-1
      do i2=1,i1-1
        covmats(i,j,iz-1)%A(i1,i2) = covmats(i,j,iz-1)%A(i2,i1)
      enddo 
      enddo
      if(i.eq.1 .and. j.eq.1) print *, covmats(i,j,iz-1)%A(1:n-1,1:n-1)
    enddo
    enddo
    enddo
  end subroutine calc_covmats
! Compute covmats: using a total covmat for all redshift bins
  subroutine calc_bigcovmats() !DEBUG c
    integer :: i,j,i1,i2,n, iz,imock, maxmubin
    integer :: row1, row2, minrow2, maxrow2, izbin
    real(rt), allocatable :: intxis(:,:,:,:,:), tmpX(:), dintxi(:,:)
    real(rt) :: smutab_covmock(nbins_covmock,mubins_covmock,3), deltas, tmpx1, tmpx2, tmpx3
    character(charlen) :: nowfile, tmpstr
    logical :: logvar, printflag, debug_calc_covmats=.false.
    ! 1. Initialize the structure "covmats"
    do i = 1, N1 !DEBUG c
    do j = 1, N2 !DEBUG c
      if (allocated(bigcovmats(i,j)%A)) deallocate(bigcovmats(i,j)%A) !DEBUG c
      n=mubins(i); bigcovmats(i,j)%nA = (n-1)*(nz-1) !DEBUG c
      allocate(bigcovmats(i,j)%A((n-1)*(nz-1),(n-1)*(nz-1))) !DEBUG c
      bigcovmats(i,j)%A = 0.0_rt !DEBUG c
    enddo
    enddo
    ! 2. Compute intxi, for nz redshift bins, all 2000 mocks
    ! 2.1. Initialize the structure intxis
    maxmubin = maxval(mubins)
    print *, 'maxmubin = ', maxmubin
    allocate(intxis(maxmubin,ncovmocks,nz,N1,N2))
    !allocate(dintxi(maxmubin,ncovmocks),tmpX(maxmubin))
    allocate(dintxi(maxmubin*(nz-1),ncovmocks),tmpX(maxmubin)) ! bc
    intxis = 0.0_rt; dintxi = 0.0_rt
    ! 2.2 Compute intxis
    do iz = 1, nz
      do imock = 1, ncovmocks
        ! 2.2.1 Check existence of 2pCF file
        nowfile=cov2pcffile(iz,imock)
        inquire(file=nowfile,exist=logvar)
        if (.not.logvar) then 
            print *, ' (calc_covmats) file not found: ', trim(adjustl(nowfile)); stop
        endif
        if(mod(imock,50).eq.0) then
                printflag=.true.
        else
                printflag=.false.
        endif
        ! 2.2.2 Load in 2pCF file
        call ximu_loadsmufile(nowfile, smutab_covmock, smax_covmock, nbins_covmock, mubins_covmock,printflag)
        ! 2.2.3 Compute intxi in the N1*N2 binning schemes
        deltas = smax_covmock / dble(nbins_covmock)
        do i = 1, N1
        do j = 1, N2
          call XiFun(smutab_covmock, deltas, nbins_covmock, mubins_covmock, &
            anglemin=1.0_rt-mucuts(j), anglemax=1.0_rt, &
            smin=ints1, smax=ints2, &
            nummuedge=mubins(i)+1, &
            intxi=tmpX(1:mubins(i)))
          call normfun(tmpX(1:mubins(i)),mubins(i),intxis(1:mubins(i)-1,imock,iz,i,j)) ! normalise the amplitude
          if (debug_calc_covmats.and.mubins(i).eq.25 .and. j.eq.1 .and. iz.eq.1 .and. imock.eq.8) then
             print *, 'nowfile = ', nowfile
             print *, 'intxi  = ', real(tmpX(1:mubins(i)))
             print *, ' python rlt =  [15.089811191950611, 13.224163227234758, ... 16.289372112713142, 16.269104458682122]'
             print *, 'intxi, normed = ', real(intxis(1:mubins(i)-1,imock,iz,i,j))
             print *, ' python rlt =  [1.0191467989709562, 0.89314329057311137, ... 1.1001636292688675, 1.0987947775009022]'
          endif
        enddo
        enddo
      enddo
    enddo
    ! 2.3 Compute covariance matrices
    ! 2.3.1 Loop of schemes & redshifts
    do i = 1, N1
    do j = 1, N2
      print *, '(calc_bigcovmat) compute bigcov for: mubin, mucut, iz = ', mubins(i), real(mucuts(j)), iz
      ! 2.3.2 dintxi = intxi@high_z - intxi@low_z; compute it for all covmocks
      n = mubins(i)
      do imock= 1, ncovmocks
       do iz= 2, nz
         if(.not.NBComp) then
           dintxi((iz-2)*(n-1)+1:(iz-1)*(n-1),imock) = intxis(1:n-1,imock,iz,i,j)-intxis(1:n-1,imock,1,i,j) !bc
         else
           dintxi((iz-2)*(n-1)+1:(iz-1)*(n-1),imock) = intxis(1:n-1,imock,iz,i,j)-intxis(1:n-1,imock,iz-1,i,j) !bc
         endif
       enddo
      enddo
      ! 2.3.3 compute covmat using dintxi
      !covmats(i,j,iz-1)%A = 0.0_rt
      bigcovmats(i,j)%A = 0.0_rt
      do i1= 1, (n-1)*(nz-1) !bc
      do i2 = i1,(n-1)*(nz-1) !bc
        ! mean(X*Y), mean(X), mean(Y)
        tmpx1 = 0.0_rt; tmpx2 = 0.0_rt; tmpx3 = 0.0_rt
        do imock = 1, ncovmocks
          tmpx1 = tmpx1 + dintxi(i1,imock)*dintxi(i2,imock)
          tmpx2 = tmpx2 + dintxi(i1,imock)
          tmpx3 = tmpx3 + dintxi(i2,imock)
        enddo
        tmpx1 = tmpx1 / dble(ncovmocks); tmpx2 = tmpx2 / dble(ncovmocks); tmpx3 = tmpx3 / dble(ncovmocks)
        ! covariance = mean(X*Y) - mean(X)*mean(Y)
        bigcovmats(i,j)%A(i1,i2) = (tmpx1 - tmpx2*tmpx3) * dble(ncovmocks) / dble(ncovmocks-1) !bigcovmat
      enddo
      enddo
      ! symmetric matrix
      do i1=1,(n-1)*(nz-1)
      do i2=1,i1-1
        bigcovmats(i,j)%A(i1,i2) = bigcovmats(i,j)%A(i2,i1)
      enddo 
      enddo
      if(i.eq.1 .and. j.eq.1) print *, bigcovmats(i,j)%A(1:(n-1)*(nz-1),1:(n-1)*(nz-1))

      ! enforce the non-box-diagnol elements being zero!
!for row1 in range(len(data1)):
!    izbin = int(row1/nmu)
!    for row2 in range(len(data1)):
!        maxrow2 = (izbin+2)*nmu -1
!        minrow2 = maxrow2 - 3*nmu + 1
!        #print minrow2, maxrow2
!        if row2 < minrow2 or row2 > maxrow2:
!            data1[row1][row2] = 0
!
      if(NBComp.and.NBComp_Simp) then
       do row1=0,(n-1)*(nz-1)-1
        izbin = int(row1/(n-1))
        do row2=0,(n-1)*(nz-1)-1
          maxrow2 = (izbin+2)*(n-1) -1
          minrow2 = maxrow2 - 3*(n-1) + 1
          if(row2<minrow2.or.row2>maxrow2) bigcovmats(i,j)%A(row1+1,row2+1) = 0
        enddo
       enddo
      endif

    enddo
    enddo
  end subroutine calc_bigcovmats
! Name of files storing covariance matrix
  character(charlen) function covmatfilename(mubin,mucut,iz,suffixstr)
    integer, intent(in) :: mubin, iz
    real(rt), intent(in) :: mucut
    character(*), intent(in), optional :: suffixstr  
    character(charlen) :: tmpstr, tmpstr1, tmpstr2, tmpstr3, tmpstr4, tmpstr5, tmpstr6, tmpstr7
    if(present(suffixstr)) then
      tmpstr = suffixstr
    else
      tmpstr = ''
    endif
    write(tmpstr1,*) mubin
    write(tmpstr2,'(f5.2)') mucut
    write(tmpstr3,*) iz
    write(tmpstr4,*) ncovmocks
    write(tmpstr5,'(f5.1)') ints1
    write(tmpstr6,'(f5.1)') ints2
    if(NBComp) then 
      tmpstr7 = '.NBComp'
    else
      tmpstr7 = ''
    endif
    covmatfilename = trim(adjustl(covmatdir))//''//trim(adjustl(tmpstr1)) &
        //'mubins.mumax' // trim(adjustl(tmpstr2))//'.iz'//trim(adjustl(tmpstr3)) &
        //'.CovMock_'//trim(adjustl(tmpstr4)) &
        //'.s'//trim(adjustl(tmpstr5))//'to'//trim(adjustl(tmpstr6))//trim(adjustl(tmpstr))//trim(adjustl(tmpstr7))//'.covmat'
  end function covmatfilename
  character(charlen) function bigcovmatfilename(mubin,mucut,suffixstr)
    integer, intent(in) :: mubin
    real(rt), intent(in) :: mucut
    character(*), intent(in), optional :: suffixstr  
    character(charlen) :: tmpstr, tmpstr1, tmpstr2, tmpstr3, tmpstr4, tmpstr5, tmpstr6, tmpstr7
    if(present(suffixstr)) then
      tmpstr = suffixstr
    else
      tmpstr = ''
    endif
    write(tmpstr1,*) mubin
    write(tmpstr2,'(f5.2)') mucut
    !write(tmpstr3,*) iz
    write(tmpstr4,*) ncovmocks
    write(tmpstr5,'(f5.1)') ints1
    write(tmpstr6,'(f5.1)') ints2
    if(NBComp) then 
      if(NBComp_Simp) then
        tmpstr7 = '.NBComp_Simp'
      else
        tmpstr7 = '.NBComp'
      endif
    else
      tmpstr7 = ''
    endif
    bigcovmatfilename = trim(adjustl(covmatdir))//''//trim(adjustl(tmpstr1)) &
        //'mubins.mumax' // trim(adjustl(tmpstr2)) &
        //'.CovMock_'//trim(adjustl(tmpstr4)) &
        //'.s'//trim(adjustl(tmpstr5))//'to'//trim(adjustl(tmpstr6))//trim(adjustl(tmpstr))//trim(adjustl(tmpstr7))//'.bigcovmat'
  end function bigcovmatfilename
! Output covmats to files
  subroutine output_covmats(suffixstr)
    character(*), intent(in), optional :: suffixstr
    integer :: i,j,iz,k,l,nA
    character(len=charlen) :: nowfile, tmpstr, tmpstr1 
    character(15) :: tmpstr2
    if(present(suffixstr)) then
      tmpstr = suffixstr
    else
      tmpstr = ''
    endif
    do i = 1, n1
    do j = 1, n2
    do iz= 2, nz
      nowfile = covmatfilename(mubins(i),mucuts(j),iz,tmpstr)
      write(*,'(A,A)'), '  Write covmat to : ', trim(adjustl(nowfile))
      open(unit=7733,file=nowfile)
      nA = covmats(i,j,iz-1)%nA
      do k = 1,nA
       write(tmpstr1, '(e15.7)') covmats(i,j,iz-1)%A(k,1)
       do l = 2,nA
         write(tmpstr2,'(e15.7)') covmats(i,j,iz-1)%A(k,l)
         tmpstr1 = trim(adjustl(tmpstr1))//' '//tmpstr2
       enddo
       write(7733,'(A)') trim(adjustl(tmpstr1))
      enddo
      close(7733)
    enddo
    enddo
    enddo
  end subroutine output_covmats
  subroutine output_bigcovmats(suffixstr)
    character(*), intent(in), optional :: suffixstr
    integer :: i,j,iz,k,l,nA
    character(len=charlen) :: nowfile, tmpstr 
    character(len=charlen*100) :: tmpstr1 
    character(15) :: tmpstr2
    if(present(suffixstr)) then
      tmpstr = suffixstr
    else
      tmpstr = ''
    endif
    do i = 1, n1
    do j = 1, n2
      nowfile = bigcovmatfilename(mubins(i),mucuts(j),tmpstr)
      print *, 'Write covmat to : ', trim(adjustl(nowfile))
      open(unit=7733,file=nowfile)
      nA = bigcovmats(i,j)%nA
      do k = 1,nA
       write(tmpstr1, '(e15.7)') bigcovmats(i,j)%A(k,1)
       do l = 2,nA
         write(tmpstr2,'(e15.7)') bigcovmats(i,j)%A(k,l)
         tmpstr1 = trim(adjustl(tmpstr1))//' '//tmpstr2
       enddo
       write(7733,'(A)') trim(adjustl(tmpstr1))
      enddo
      close(7733)
    enddo
    enddo
  end subroutine output_bigcovmats
! Load covmats from files
  subroutine load_covmats(suffixstr)
    character(*), intent(in), optional :: suffixstr
    integer :: i,j,iz,k,l,nA
    character(len=charlen) :: nowfile, tmpstr, tmpstr1 
    character(15) :: tmpstr2
    if(present(suffixstr)) then
      tmpstr = suffixstr
    else
      tmpstr = ''
    endif
    do i = 1, n1
    do j = 1, n2
    do iz= 2, nz
      nowfile = covmatfilename(mubins(i),mucuts(j),iz,tmpstr)
      write(*,'(A,/,A)'), ' (load_covmats) Load covmat from : ', trim(adjustl(nowfile))
      open(unit=7733,file=nowfile,action='read')
      nA = mubins(i)-1
      if (allocated(covmats(i,j,iz-1)%A) .and. covmats(i,j,iz-1)%nA .ne. nA) deallocate(covmats(i,j,iz-1)%A)
      if (.not. allocated(covmats(i,j,iz-1)%A)) allocate(covmats(i,j,iz-1)%A(nA,nA))
      covmats(i,j,iz-1)%nA = nA
      do k = 1, nA
       read(7733,*) covmats(i,j,iz-1)%A(:,k)
      enddo
      close(7733)
    enddo
    enddo
    enddo
  end subroutine load_covmats
  subroutine load_bigcovmats(suffixstr)
    character(*), intent(in), optional :: suffixstr
    integer :: i,j,k,l,nA
    character(len=charlen) :: nowfile, tmpstr, tmpstr1 
    character(15) :: tmpstr2
    if(present(suffixstr)) then
      tmpstr = suffixstr
    else
      tmpstr = ''
    endif
    do i = 1, n1
    do j = 1, n2
      nowfile = bigcovmatfilename(mubins(i),mucuts(j),tmpstr)
      write(*,'(A,/,A)'), ' (load_bigcovmats) Load bigcovmat from : ', trim(adjustl(nowfile))
      open(unit=7733,file=nowfile,action='read')
      nA = (mubins(i)-1) * (nz-1)
      if (allocated(bigcovmats(i,j)%A) .and. bigcovmats(i,j)%nA .ne. nA) deallocate(bigcovmats(i,j)%A)
      if (.not. allocated(bigcovmats(i,j)%A)) allocate(bigcovmats(i,j)%A(nA,nA))
      bigcovmats(i,j)%nA = nA
      do k = 1, nA
       read(7733,*) bigcovmats(i,j)%A(:,k)
      enddo
      close(7733)
    enddo
    enddo
  end subroutine load_bigcovmats
! Invert all covmats
  subroutine invert_covmats()
    integer :: i,j,iz, n
    real(rt), allocatable :: B(:,:)
    do i = 1, n1
      allocate(B(mubins(i)-1,mubins(i)-1))
      do j = 1, n2
        do iz= 2, nz
          call nizhen(covmats(i,j,iz-1)%A, B, mubins(i)-1)
          covmats(i,j,iz-1)%A = B
        enddo
      enddo
      deallocate(B)
    enddo
  end subroutine invert_covmats
! Invert all covmats
  subroutine invert_bigcovmats()
    integer :: i,j, n,iz
    real(rt), allocatable :: B(:,:), tmpB(:,:)
    do i = 1, n1
      n = mubins(i)-1
      allocate(B((mubins(i)-1)*(nz-1),(mubins(i)-1)*(nz-1)))
      allocate(tmpB(n,n))
      do j = 1, n2
          call nizhen(bigcovmats(i,j)%A, B, (mubins(i)-1)*(nz-1))
          bigcovmats(i,j)%A = B
!        !This will ignore the correlation between different redshift bins!
        !do iz = 2, nz
        !  call nizhen(bigcovmats(i,j)%A((iz-2)*n+1:(iz-1)*n,(iz-2)*n+1:(iz-1)*n), tmpB, n)
        !  bigcovmats(i,j)%A((iz-2)*n+1:(iz-1)*n,(iz-2)*n+1:(iz-1)*n) = tmpB
        !enddo
      enddo
      deallocate(B); deallocate(tmpB)
    enddo
  end subroutine invert_bigcovmats
  subroutine spmat_bigcovmatsfc()
    integer :: i,j
    do i = 1, n1
      do j = 1, n2
        call spmat_linsys_inout(bigcovmatsfc(i,j)%A, (mubins(i)-1)*(nzfc-1))
      enddo
    enddo
  end subroutine spmat_bigcovmatsfc
  subroutine invert_bigcovmatsfc()
    integer :: i,j, n,iz
    real(rt), allocatable :: B(:,:), tmpB(:,:)
    do i = 1, n1
      n = mubins(i)-1
      allocate(B((mubins(i)-1)*(nzfc-1),(mubins(i)-1)*(nzfc-1)))
      allocate(tmpB(n,n))
      do j = 1, n2
          call nizhen(bigcovmatsfc(i,j)%A, B, (mubins(i)-1)*(nzfc-1))
          bigcovmatsfc(i,j)%A = B
!        !This will ignore the correlation between different redshift bins!
        !do iz = 2, nz
        !  call nizhen(bigcovmats(i,j)%A((iz-2)*n+1:(iz-1)*n,(iz-2)*n+1:(iz-1)*n), tmpB, n)
        !  bigcovmats(i,j)%A((iz-2)*n+1:(iz-1)*n,(iz-2)*n+1:(iz-1)*n) = tmpB
        !enddo
      enddo
      deallocate(B); deallocate(tmpB)
    enddo
  end subroutine invert_bigcovmatsfc

  subroutine calc_syscorfc()
    dintxi_syscor = 0.0_rt
  end subroutine calc_syscorfc

  subroutine calc_syscor()
    integer :: i,j,i1,i2,n, iz,imock, maxmubin
    real(rt) :: smutab_sysmock(nbins_sysmock,mubins_sysmock,3), &
      intxis(maxval(mubins),nsysmocks,nz,N1,N2),&
      intxi0(maxval(mubins)), intxi1(maxval(mubins)),&
      tmpX(maxval(mubins)),X(maxval(mubins)),mumids(maxval(mubins),N1,N2),&
      deltas, tmpx1, tmpx2, tmpx3
    real(rt) :: smutab_PTmock(500,600,3), smtab_PTmock_simple(100,120,3)
    character(charlen) :: nowfile, tmpstr
    character(len=4) :: IO_id
    logical :: logvar, debug_calc_syscor = .false.
    ! variables used when using data self as systematic correction
    type(omwpar) :: tmpomw
    real(rt) :: smin_mapping=min(1.0_rt,ints1*0.8_rt), smax_mapping=ints2*1.2_rt, &
      DAstd,Hstd, deltas1,deltas2
    intxis = 0.0_rt
    ! 1. Compute intxi, for nz z-bins, all syscor mocks, N1*N2 schemes

    if(dataself_IO_test) then
     write(*,'(A)')  'WARNING (calc_syscor)! Using observational data itself as systematic correction!!!! '
     !if(mock_IO_test) then
     !  write(*,'(A)')  'ERROR (calc_syscor)! mock_IO_test and dataself_IO_test both be true: ', mock_IO_test, dataself_IO_test
     !  stop
     !endif
     if(.not.allocated(global_smutabstds)) then
       write(*,'(A)')  'ERROR (calc_syscor)! global_smutabstds not allocated; not ready for self-sys computation.'
       stop
     endif
     if(size(global_smutabstds,dim=5).gt. 1) then
       write(*,'(A)')  'ERROR (calc_syscor)! find more than one standard cosmology when dataself_IO_test=.true. : ', size(global_smutabstds,dim=5)
       stop
     else
       write(*,'(A,i3)') '(calc_syscor) dataself_IO_test: check size of global_smutabstds dim-5 (should be one): ', size(global_smutabstds,dim=5)
     endif
     deltas1 = smax_database / float(nbins_database)
     deltas2 = smax_sysmock / float(nbins_sysmock)
     do iz = 1, nz
        tmpomw%omegam = 0.26_rt; tmpomw%w = -1.0_rt; ! as long as we use same DA/H to convert, we do not care whether 
                                                     !  standard omegam/w is fiducial or not
        DAstd = DAz_wcdm(tmpomw,zeffs(iz)); Hstd = Hz_wcdm(tmpomw,zeffs(iz))
        ! Using same DA/H to convert; so there is no cosmological mapping; only simple binning. 
        call DSMapping(global_smutabstds(:,:,:,iz,1), nbins_database, mubins_database, smutab_sysmock, &
            nbins_sysmock, mubins_sysmock, DAstd, DAstd, Hstd, Hstd, deltas1,  deltas2,  smin_mapping, smax_mapping)
            !nbins_data, mubins_data, DAstd, DAnew, Hstd, Hnew, deltas1,  deltas2,  smin_mapping, smax_mapping)
        do i = 1, N1
        do j = 1, N2
          tmpX = 0.0_rt
          call XiFun(smutab_sysmock, deltas2, nbins_sysmock, mubins_sysmock, &
            anglemin=1.0_rt-mucuts(j), anglemax=1.0_rt, &
            smin=ints1, smax=ints2, &
            nummuedge=mubins(i)+1, &
            intxi=tmpX(1:mubins(i)),&
            mumids=mumids(1:mubins(i),i,j))
          do imock = 1, nsysmocks
            call normfun(tmpX(1:mubins(i)),mubins(i),intxis(1:mubins(i)-1,imock,iz,i,j)) ! normalise the amplitude
          enddo
        enddo
        enddo
     enddo 

   elseif(gb_syscor_mock .eq. gb_J08_mock) then


     do iz = 1, nz
     do imock = 1, nsysmocks
        ! 1.1 Check existence of 2pCF file
        nowfile=syscor2pcffile(iz,imock)
        if(mock_IO_test_usemock3) then
          print *, 'WARNING!!!!!! Using a particular mock to do systematic-cor!!!'
          nowfile=syscor2pcffile(iz,3)
        endif
        inquire(file=nowfile,exist=logvar)
        if (.not.logvar) then 
            print *, ' (calc_syscor) file not found: ', trim(adjustl(nowfile)); stop
        endif
        ! 1.2 Load in 2pCF file
        !write(*,'(A,/,A)'), ' (calc_syscor) load in syscor file: ', trim(adjustl(nowfile));
        call ximu_loadsmufile(nowfile, smutab_sysmock, smax_sysmock, nbins_sysmock, mubins_sysmock)
        ! 1.3 Compute intxi in the N1*N2 binning schemes
        deltas = smax_sysmock / dble(nbins_sysmock)
        do i = 1, N1
        do j = 1, N2
          tmpX = 0.0_rt
          call XiFun(smutab_sysmock, deltas, nbins_sysmock, mubins_sysmock, &
            anglemin=1.0_rt-mucuts(j), anglemax=1.0_rt, &
            smin=ints1, smax=ints2, &
            nummuedge=mubins(i)+1, &
            intxi=tmpX(1:mubins(i)),&
            mumids=mumids(1:mubins(i),i,j))
          call normfun(tmpX(1:mubins(i)),mubins(i),intxis(1:mubins(i)-1,imock,iz,i,j)) ! normalise the amplitude
          !intxis(maxval(mubins),nsysmocks,nz,N1,N2)
!          print *, real(mumids(1:mubins(i)))
          if (debug_calc_syscor .and. mubins(i).eq.25 .and. j .eq. 1 .and.iz.eq.1.and.imock.eq.1) then
                print *, '##############################'
                print *, 'Check syscor: '
                  print *, ' * iz, imock = ', iz, imock
                  print *, '   intxi = ', real(tmpX(1:mubins(i)))
                  print *, '   python result = [14.579422652657449, 12.460830746410542,  ... 18.923018560735642]'
                  print *, '   intxi, normed = ', real(intxis(1:mubins(i)-1,imock,iz,i,j))
                print *, '   python result = [0.91530223976179914, 0.78229615556164356, ... 1.1965964165392333, 1.18...]'
                  print *
          endif
        enddo
        enddo
!        stop
     enddo
     enddo

   elseif(gb_syscor_mock .eq. gb_PT_mock) then

     write(*,'(A)')  'WARNING (calc_syscor)! Using PATCHY mock as systematic correction!!!! '

     ! again, we have to convert table from one resolution to anther!
     deltas1 = 100. / 500.
     deltas2 = smax_sysmock / float(nbins_sysmock)

     do iz = 1, nz
     do imock = 1, nsysmocks
        ! 1.1 Check existence of 2pCF file
	write(tmpstr,*) imock-1 !imock
	if(imock < 10 ) then
		IO_id = '000'//trim(adjustl(tmpstr))
	else
		print *, 'ERROR! (calc_syscor) not supporting >=10 PT mock for sys correction! imock, nsysmocks = ', imock, nsysmocks
		stop
	endif
	nowfile = PT2pcffile(iz,IO_id)

        inquire(file=nowfile,exist=logvar)
        if (.not.logvar) then 
            print *, ' (calc_syscor) file not found: ', trim(adjustl(nowfile)); stop
        endif
        ! 1.2 Load in 2pCF file
        call ximu_loadsmufile(nowfile, smutab_PTmock, 100, 500,600)

        tmpomw%omegam = 0.26_rt; tmpomw%w = -1.0_rt; ! as long as we use same DA/H to convert, we do not care whether 
                                                     !  standard omegam/w is fiducial or not
        DAstd = DAz_wcdm(tmpomw,zeffs(iz)); Hstd = Hz_wcdm(tmpomw,zeffs(iz))
	call DSMapping(smutab_PTmock, 500, 600, smtab_PTmock_simple, &
            100, 120, DAstd, DAstd, Hstd, Hstd, deltas1,  deltas2,  smin_mapping, smax_mapping)
!        write(*,'(A,/,A)'), ' (calc_syscor) Finish dsmapping.'
        ! 1.3 Compute intxi in the N1*N2 binning schemes
        deltas = 1.0 !smax_sysmock / dble(nbins_sysmock)
        do i = 1, N1
        do j = 1, N2
          tmpX = 0.0_rt
          !call XiFun(smutab_sysmock, deltas, nbins_sysmock, mubins_sysmock, &
          call XiFun(smtab_PTmock_simple, deltas, 100, 120, &
            anglemin=1.0_rt-mucuts(j), anglemax=1.0_rt, &
            smin=ints1, smax=ints2, &
            nummuedge=mubins(i)+1, &
            intxi=tmpX(1:mubins(i)),&
            mumids=mumids(1:mubins(i),i,j))
          call normfun(tmpX(1:mubins(i)),mubins(i),intxis(1:mubins(i)-1,imock,iz,i,j)) ! normalise the amplitude
          !intxis(maxval(mubins),nsysmocks,nz,N1,N2)
!          print *, real(mumids(1:mubins(i)))
          if (debug_calc_syscor .and. mubins(i).eq.25 .and. j .eq. 1 .and.iz.eq.1.and.imock.eq.1) then
!          if (i.eq.1 .and. j .eq. 1 .and.iz.eq.1) then
!          if (.true.) then
                print *, '##############################'
                print *, 'Check syscor: '
                  print *, ' * iz, imock = ', iz, imock
                  print *, '   intxi = ', real(tmpX(1:mubins(i)))
                  print *, '   python result = [14.579422652657449, 12.460830746410542,  ... 18.923018560735642]'
                  print *, '   intxi, normed = ', real(intxis(1:mubins(i)-1,imock,iz,i,j))
                  print *, '   python result = [0.91530223976179914, 0.78229615556164356, ... 1.1965964165392333, 1.18...]'
                  print *
          endif
        enddo
        enddo
!        stop
     enddo
     enddo
    endif

    ! 2. Redshift evolution of intxi, for nz-1 z-bins, all syscor mocks, N1*N2 schemes
    do i = 1, N1
    do j = 1, N2
      ! 2.1 reference intxi at first bin
      intxi0 = 0.0_rt
      do imock = 1, nsysmocks
        intxi0(1:mubins(i)-1) = intxi0(1:mubins(i)-1) + intxis(1:mubins(i)-1,imock,1,i,j)
      enddo
      intxi0 = intxi0 / dble(nsysmocks)
      ! 2.2 intxi at a higher redshift bin
      do iz = 2, nz
        intxi1 = 0.0_rt
        do imock = 1, nsysmocks
          intxi1(1:mubins(i)-1) = intxi1(1:mubins(i)-1) + intxis(1:mubins(i)-1,imock,iz,i,j)
        enddo
        intxi1 = intxi1 / dble(nsysmocks)
        ! 2.3 Redshift evolution
         dintxi_syscor(1:mubins(i)-1,i,j,iz-1) = &
          intxi1(1:mubins(i)-1) - intxi0(1:mubins(i)-1)
!        if(i.eq.1.and.j.eq.1) write(*,*) 'calc_syscor debug dintxi: ',iz,dintxi_syscor(1:mubins(i)-1,i,j,iz-1) ! selfsyscor_debug
        if(NBComp) intxi0 = intxi1
      enddo
    enddo
    enddo
    
    ! polynomial fit 
    if (polyfitdeg .ge. 1) then
      do i = 1, maxval(mubins)
        X(i) = i-1
      enddo
      do i = 1, N1
        do j = 1, N2
          do iz = 2, nz
            dintxi_syscor(1:mubins(i)-1,i,j,iz-1) = &
              polyfitY(X(1:mubins(i)-1),dintxi_syscor(1:mubins(i)-1,i,j,iz-1),mubins(i)-1,polyfitdeg)
!              polyfitY(mumids(1:mubins(i)-1,i,j),dintxi_syscor(1:mubins(i)-1,i,j,iz-1),mubins(i)-1,polyfitdeg) ! using exact values of mumids for polyfit; almost no effect on the resulted chisq values
          enddo
        enddo
      enddo
   endif
  end subroutine calc_syscor
  
! Coordiante transformation of (s,mu) : from one cosmology to another
  subroutine smu__CosmoConvert(s,mu,DA1,DA2,H1,H2,s_prime,mu_prime)
    real(rt), intent(in) :: s,mu,DA1,DA2,H1,H2
    real(rt), intent(out) :: s_prime,mu_prime
    real(rt) :: alpha1, alpha2, s1, s2  ! s1: angular direction; s2: LOS direction 
    s2 = s*mu;
    s1 = sqrt(s*s - s2*s2)
    alpha1 = DA2 / DA1
    alpha2 = H1 / H2
    s_prime  =  sqrt((alpha1*s1)**2 + (alpha2*s2)**2)
    mu_prime =  alpha2*s2 / s_prime
  end subroutine smu__CosmoConvert


! Mapping xi(s,mu) from baseline cosmology to another 
! Dense grid in baseline cosmology, sparse grid in another
  subroutine DSMapping(smutabstd, nums1, nummu1, smutab2, &
          nums2, nummu2, DAstd, DAnew, Hstd, Hnew, deltas1,  deltas2,  &
          smin_mapping, smax_mapping ) ! range of s considered in the coordinate transformation:  
                                       ! basically, smin_mapping < s1 < s2 < smax_mapping.
    ! argument
    real(rt), intent(in) :: smutabstd(nums1,nummu1,3), DAstd, DAnew, Hstd, Hnew, &
      deltas1, deltas2, smin_mapping, smax_mapping
    integer, intent(in) :: nums1, nummu1, nums2, nummu2
    real(rt), intent(out):: smutab2(nums2,nummu2,3)
    ! variables
    real(rt) :: deltamu1,deltamu2, mubound1,mubound2,sbound1,sbound2, &
      scenter,anglecenter,mucenter,scenter2,anglecenter2,mucenter2, &
      scenter2_bound,anglecenter2_bound,anglecenter3,scenter3, countDD,countDR,countRR,&
      s,angle,ds,dangle,d1,rat,rat_bound,rats,ratangle,rat1,rat2,rat3,rat4
    real(rt) :: smutabstd_centers(nums1,nummu1,2)
    integer :: maxs1, mins1, is1,iangle1,is2,iangle2, is_nearby,imu_nearby, &
      is2_bound, is2_b, iangle2_b, iangle2_bound, &
      is1_bound, iangle1_bound, iangle1_nearby, is1_nearby, &
      smutabstd_ismus(nums1,nummu1,2)
    logical :: sboundflag, muboundflag
    
    deltamu1 = 1.0d0/real(nummu1); deltamu2 = 1.0d0/real(nummu2)
    smutab2=0.0d0; smutabstd_ismus=10000
    ! range of s for smutab1, smutab2
    call smu__CosmoConvert(smax_mapping,0.0_rt,DAnew,DAstd,Hnew,Hstd,sbound1,mubound1)
    call smu__CosmoConvert(smax_mapping,1.0_rt,DAnew,DAstd,Hnew,Hstd,sbound2,mubound2)
    maxs1 = floor( max(sbound1,sbound2) / deltas1 + 0.5 )
    if(maxs1>nums1) print *, ' WARNING (mapping_smudata_to_another_cosmology_DenseToSparse)!',&
      ' Outflow of s: nums1, maxs1 = ', nums1, maxs1
    call smu__CosmoConvert(smin_mapping,0.0_rt,DAnew,DAstd,Hnew,Hstd,sbound1,mubound1)
    call smu__CosmoConvert(smin_mapping,1.0_rt,DAnew,DAstd,Hnew,Hstd,sbound2,mubound2)
    mins1 = floor( min(sbound1,sbound2) / deltas1 + 0.5 )
    mins1=max(mins1,0); maxs1=min(maxs1,nums1-1)
    
!    print *, 'deltas1, deltamu1 ', deltas1, deltamu1
!    print *, 'deltas2, deltamu2 ', deltas2, deltamu2
!    print *, 'sbound1, sbound2  ', sbound1, sbound2
    !print *, mins1, maxs1
    
!    maxs2 = min(floor(smax_mapping / deltas2 + 0.5), nums2)
  !      deltas1=0.2, deltamu1=1.0/600.0, deltas2=1.0, deltamu2=1.0/120.0, smin_mapping=1,smax_mapping=51,

!    ''' method can be simple_bin or divided_pixel'''
!    nummu1 = int(1.0/deltamu1 + 0.5)
!    sbound1, mubound1 = smu__CosmoConvert(smax_mapping,0,DAnew,DAstd,Hnew,Hstd)
!    sbound2, mubound2 = smu__CosmoConvert(smax_mapping,1,DAnew,DAstd,Hnew,Hstd)
!    sbound = max(sbound1, sbound2); nums1 = int(sbound / deltas1 + 0.5)
!    
!    sbound1, mubound1 = smu__CosmoConvert(smin_mapping,0,DAnew,DAstd,Hnew,Hstd)
!    sbound2, mubound2 = smu__CosmoConvert(smin_mapping,1,DAnew,DAstd,Hnew,Hstd)
!    sbound = min(sbound1, sbound2); mins1 = int(sbound / deltas1)
!    
!    nums2, nummu2 = int(smax_mapping / deltas2 + 0.5), int(1.0/deltamu2 + 0.5)
!    numrow3=len(smutabstd[0][0])
!    smutab2 = [[[0 for row3 in range(numrow3+1)] for row2 in range(nummu2)] for row1 in range(nums2)]

!    open(unit=9231489,file='smutab_database.txt')! haha

    do is1 = mins1, maxs1-1
      do iangle1 = 0, nummu1-1
        scenter=(is1+0.5_rt)*deltas1; anglecenter=(iangle1+0.5)*deltamu1
        mucenter=1.0_rt-anglecenter
        call smu__CosmoConvert(scenter,mucenter,DAstd,DAnew,Hstd,Hnew,scenter2,mucenter2)
        anglecenter2 = 1.0_rt - mucenter2
        is2 = floor(scenter2 / deltas2); iangle2 = floor(anglecenter2/deltamu2)
        smutabstd_centers(is1+1,iangle1+1,1) = scenter2
        smutabstd_centers(is1+1,iangle1+1,2) = anglecenter2
        smutabstd_ismus(is1+1,iangle1+1,1) = is2
        smutabstd_ismus(is1+1,iangle1+1,2) = iangle2
!        write(9231489,'(i4,i4,3(e14.7))') is1, iangle1, smutabstd(is1+1,iangle1+1,1:3)! haha
      enddo
    enddo
!    close(9231489) ! haha
    
    
    !print *, 'Initialization of tabs done!'

!    
!    #print mins1, nums1, nummu1, nums2, nummu2
!    if method == 'divided_pixel':
!        smutabstd_centers = [[0 for row2 in range(nummu1)] for row1 in range(nums1)]
!    for is1 in range(mins1, nums1):
!        for iangle1 in range(nummu1):
!            scenter, anglecenter = (is1+0.5)*deltas1, (iangle1+0.5)*deltamu1
!            mucenter = 1.0 - anglecenter
!            scenter2, mucenter2 = smu__CosmoConvert(scenter,mucenter,DAstd,DAnew,Hstd,Hnew,)
!            anglecenter2 = 1.0 - mucenter2
!            is2 = int(scenter2  / deltas2  )
!            iangle2 = int(anglecenter2 / deltamu2)
!            if method == 'simple_bin':
!                if is2 < nums2 and iangle2 < nummu2:
!                    for row3 in compute_rows:
!                        smutab2[is2][iangle2][row3] += smutabstd[is1][iangle1][row3]
!                    if save_counts_row != None:
!                        smutab2[is2][iangle2][save_counts_row] += 1
!            elif method == 'divided_pixel':
!                smutabstd_centers[is1][iangle1] = [scenter2, anglecenter2, is2, iangle2]
!                
    do is1 = mins1, maxs1-1
      do iangle1 = 0, nummu1-1
        scenter2=smutabstd_centers(is1+1,iangle1+1,1)
        anglecenter2=smutabstd_centers(is1+1,iangle1+1,2)
        is2=smutabstd_ismus(is1+1,iangle1+1,1)
        iangle2=smutabstd_ismus(is1+1,iangle1+1,2)
        if (is2.ge.nums2 .or. iangle2.ge.nummu2) cycle
        
 !    if method == 'divided_pixel':
!        for is1 in range(mins1, nums1):
!            for iangle1 in range(nummu1):
!                scenter2, anglecenter2, is2, iangle2 = smutabstd_centers[is1][iangle1]
!                if not (is2 < nums2 and iangle2 < nummu2):
!                    continue       
        
        sboundflag=.false.; muboundflag=.false.
        do is1_nearby = max(is1-1,mins1), min(is1+1,maxs1-1)
          is2_b = smutabstd_ismus(is1_nearby+1,iangle1+1,1)
          iangle2_b = smutabstd_ismus(is1_nearby+1,iangle1+1,2)
          if (is2_b .ne. is2) then 
            sboundflag=.true.; is1_bound=is1_nearby; is2_bound=is2_b;
            scenter2_bound=smutabstd_centers(is1_nearby+1,iangle1+1,1)
          endif
        enddo
        do iangle1_nearby = max(iangle1-1,0), min(iangle1+1,nummu1-1)        
          is2_b = smutabstd_ismus(is1+1,iangle1_nearby+1,1)
          iangle2_b = smutabstd_ismus(is1+1,iangle1_nearby+1,2)
          if (iangle2_b .ne. iangle2) then 
            muboundflag=.true.; iangle1_bound=iangle1_nearby; iangle2_bound=iangle2_b;
            anglecenter2_bound=smutabstd_centers(is1+1,iangle1_nearby+1,2)
          endif
        enddo

            !seriously possible bug found in the python code! Decide to stop and checking the python code for a while.
          
            
!                ### firstly, check boundary:
!                sboundflag, muboundflag = False, False
!                for is1_nearby in [max(is1-1,mins1), min(is1+1,nums1-1)]:
!                    is2_b, iangle2_b = smutabstd_centers[is1_nearby][iangle1][2], smutabstd_centers[is1_nearby][iangle1][3]
!                    if is2_b != is2:
!                        sboundflag=True; is1_bound = is1_nearby; is2_bound = is2_b; 
!                        scenter2_bound = smutabstd_centers[is1_nearby][iangle1][0]
!                for iangle1_nearby in [max(iangle1-1,0), min(iangle1+1,nummu1-1)]:
!                    is2_b, iangle2_b = smutabstd_centers[is1][iangle1_nearby][2], smutabstd_centers[is1][iangle1_nearby][3]
!                    if iangle2_b != iangle2:
!                        muboundflag=True; iangle1_bound = iangle1_nearby; iangle2_bound = iangle2_b
!                        anglecenter2_bound = smutabstd_centers[is1][iangle1_nearby][1]
         
!                    
          countDD = smutabstd(is1+1,iangle1+1,1); 
          countDR = smutabstd(is1+1,iangle1+1,2); 
          countRR = smutabstd(is1+1,iangle1+1,3)
          if ( .not.sboundflag .and. .not. muboundflag) then  
            smutab2(is2+1,iangle2+1,1) = smutab2(is2+1,iangle2+1,1)+countDD
            smutab2(is2+1,iangle2+1,2) = smutab2(is2+1,iangle2+1,2)+countDR
            smutab2(is2+1,iangle2+1,3) = smutab2(is2+1,iangle2+1,3)+countRR
          endif

!                
!                ### Then, treat them case by case...
!                ## s, mu are all not near the boundary of tab 2
!                if ((not sboundflag)and(not muboundflag)):
!                    for row3 in compute_rows:
!                        smutab2[is2][iangle2][row3] += smutabstd[is1][iangle1][row3]
!                    if save_counts_row != None: smutab2[is2][iangle2][save_counts_row] += 1


          if ( sboundflag .and. .not. muboundflag) then 
            s = (is2 + is2_bound + 1) * 0.5 * deltas2
            scenter3 = (scenter2+scenter2_bound) / 2.0
            ds = scenter3 -scenter2
            d1 = s-scenter2+ds
            rat = min(d1 / (2.0_rt*ds),1.0_rt)
            rat_bound = 1.0-rat

            smutab2(is2+1,iangle2+1,1) = smutab2(is2+1,iangle2+1,1)+countDD*rat
            smutab2(is2+1,iangle2+1,2) = smutab2(is2+1,iangle2+1,2)+countDR*rat
            smutab2(is2+1,iangle2+1,3) = smutab2(is2+1,iangle2+1,3)+countRR*rat
            
            if (is2_bound < nums2) then
             smutab2(is2_bound+1,iangle2+1,1) = smutab2(is2_bound+1,iangle2+1,1)+countDD*rat_bound
             smutab2(is2_bound+1,iangle2+1,2) = smutab2(is2_bound+1,iangle2+1,2)+countDR*rat_bound
             smutab2(is2_bound+1,iangle2+1,3) = smutab2(is2_bound+1,iangle2+1,3)+countRR*rat_bound
            endif
          endif

!                            
!                ## s is near the boundary of tab2
!                if sboundflag and (not muboundflag):
!                    s = (is2 + is2_bound+1) * 0.5 * deltas2
!                    if False:
!                        rat = (s-scenter2) / (scenter2_bound-scenter2)
!                    else:
!                        scenter3=(scenter2+scenter2_bound) * 0.5
!                        ds = scenter3-scenter2
!                        d1 = s-scenter2+ds
!                        rat = d1 / (2*ds)
!                        rat = min(rat, 1)
!                        
!                    rat_bound = 1-rat
!                    for row3 in compute_rows:
!                        smutab2[is2][iangle2][row3] += smutabstd[is1][iangle1][row3]*rat
!                    if save_counts_row != None: smutab2[is2][iangle2][save_counts_row] += rat
!                    if is2_bound < nums2:
!                      for row3 in compute_rows:
!                        smutab2[is2_bound][iangle2][row3] += smutabstd[is1][iangle1][row3]*rat_bound
!                      if save_counts_row != None: smutab2[is2_bound][iangle2][save_counts_row] += rat_bound
!              
          if ( muboundflag .and. .not. sboundflag) then 
            angle = (iangle2 + iangle2_bound + 1) * 0.5 * deltamu2
            anglecenter3 = (anglecenter2+anglecenter2_bound) / 2.0
            dangle = anglecenter3 -anglecenter2
            d1 = angle-anglecenter2+dangle
            rat = min(d1 / (2.0_rt*dangle),1.0_rt)
            rat_bound = 1.0-rat

            smutab2(is2+1,iangle2+1,1) = smutab2(is2+1,iangle2+1,1)+countDD*rat
            smutab2(is2+1,iangle2+1,2) = smutab2(is2+1,iangle2+1,2)+countDR*rat
            smutab2(is2+1,iangle2+1,3) = smutab2(is2+1,iangle2+1,3)+countRR*rat
            
            if (iangle2_bound < nummu2) then
             smutab2(is2+1,iangle2_bound+1,1) = smutab2(is2+1,iangle2_bound+1,1)+countDD*rat_bound
             smutab2(is2+1,iangle2_bound+1,2) = smutab2(is2+1,iangle2_bound+1,2)+countDR*rat_bound
             smutab2(is2+1,iangle2_bound+1,3) = smutab2(is2+1,iangle2_bound+1,3)+countRR*rat_bound
            endif
          endif      


!                ## mu is near the boundary of tab2
!                if muboundflag and (not sboundflag):
!                    angle = (iangle2 + iangle2_bound+1) * 0.5 * deltamu2
!!                    if False:
 !                       rat = (angle-anglecenter2) / (anglecenter2_bound-anglecenter2)
 !!                   else:
!!                        anglecenter3=(anglecenter2+anglecenter2_bound) * 0.5
!                        dangle = anglecenter3-anglecenter2
!                        d1 = angle-anglecenter2+dangle
!                        rat = d1 / (2*dangle)
!!                        rat = min(rat, 1)
                        
!                    rat_bound = 1-rat
!                    for row3 in compute_rows:
!                        smutab2[is2][iangle2][row3] += smutabstd[is1][iangle1][row3]*rat
!                    if save_counts_row != None: smutab2[is2][iangle2][save_counts_row] += rat
!                    if iangle2_bound < nummu2:
!                      for row3 in compute_rows:
!                        smutab2[is2][iangle2_bound][row3] += smutabstd[is1][iangle1][row3]*rat_bound
 !                     if save_counts_row != None: smutab2[is2][iangle2_bound][save_counts_row] += rat_bound
 !               

          if ( muboundflag .and. sboundflag) then 
            s = (is2 + is2_bound + 1) * 0.5 * deltas2
            angle = (iangle2 + iangle2_bound + 1) * 0.5 * deltamu2
            scenter3 = (scenter2+scenter2_bound) / 2.0
            ds = scenter3 -scenter2
            d1 = s-scenter2+ds
            rats = min(d1 / (2.0_rt*ds),1.0_rt)

            anglecenter3 = (anglecenter2+anglecenter2_bound) / 2.0
            dangle = anglecenter3 -anglecenter2
            d1 = angle-anglecenter2+dangle
            ratangle = min(d1 / (2.0_rt*dangle),1.0_rt)
  !              ## both s, mu are near the boundy...
!                if muboundflag and sboundflag:
!                    s = (is2 + is2_bound+1) * 0.5 * deltas2
!                    angle = (iangle2 + iangle2_bound+1) * 0.5 * deltamu2
!                    if False:
!                        rats = (s-scenter2) / (scenter2_bound-scenter2)
!                        ratangle = (angle-anglecenter2) / (anglecenter2_bound-anglecenter2)
!                    else:
!                        scenter3=(scenter2+scenter2_bound) * 0.5
!                        ds = scenter3-scenter2
!                        d1 = s-scenter2+ds
!                        rats = d1 / (2*ds)
!                        rats = min(rats, 1)
!!                        anglecenter3=(anglecenter2+anglecenter2_bound) * 0.5
 !                       dangle = anglecenter3-anglecenter2
 !                       d1 = angle-anglecenter2+dangle
 !                       ratangle = d1 / (2*dangle)
 !                       ratangle = min(ratangle, 1)            
            rat1 = rats*ratangle
            smutab2(is2+1,iangle2+1,1) = smutab2(is2+1,iangle2+1,1)+countDD*rat1
            smutab2(is2+1,iangle2+1,2) = smutab2(is2+1,iangle2+1,2)+countDR*rat1
            smutab2(is2+1,iangle2+1,3) = smutab2(is2+1,iangle2+1,3)+countRR*rat1

 !                   # original pixel
 !                   rat1 = rats*ratangle
 !                   for row3 in compute_rows:
 !                       smutab2[is2][iangle2][row3] += smutabstd[is1][iangle1][row3]*rat1
 !                   if save_counts_row != None: smutab2[is2][iangle2][save_counts_row] += rat1

            rat2 = (1-rats) * ratangle
            if (is2_bound < nums2) then
             smutab2(is2_bound+1,iangle2+1,1) = smutab2(is2_bound+1,iangle2+1,1)+countDD*rat2
             smutab2(is2_bound+1,iangle2+1,2) = smutab2(is2_bound+1,iangle2+1,2)+countDR*rat2
             smutab2(is2_bound+1,iangle2+1,3) = smutab2(is2_bound+1,iangle2+1,3)+countRR*rat2
            endif
            
            rat3=rats*(1-ratangle)
            if (iangle2_bound < nummu2) then
             smutab2(is2+1,iangle2_bound+1,1) = smutab2(is2+1,iangle2_bound+1,1)+countDD*rat3
             smutab2(is2+1,iangle2_bound+1,2) = smutab2(is2+1,iangle2_bound+1,2)+countDR*rat3
             smutab2(is2+1,iangle2_bound+1,3) = smutab2(is2+1,iangle2_bound+1,3)+countRR*rat3
            endif
            
            rat4=(1-rats)*(1-ratangle)
            if (iangle2_bound < nummu2 .and. is2_bound<nums2) then
             smutab2(is2_bound+1,iangle2_bound+1,1) = smutab2(is2_bound+1,iangle2_bound+1,1)+countDD*rat4
             smutab2(is2_bound+1,iangle2_bound+1,2) = smutab2(is2_bound+1,iangle2_bound+1,2)+countDR*rat4
             smutab2(is2_bound+1,iangle2_bound+1,3) = smutab2(is2_bound+1,iangle2_bound+1,3)+countRR*rat4
            endif
          endif              
          cycle ! haha        
!                    # diff s
 !                   rat2 = (1-rats)*ratangle
!                    if is2_bound < nums2:
!                      for row3 in compute_rows:
!                        smutab2[is2_bound][iangle2][row3] += smutabstd[is1][iangle1][row3]*rat2
!                      if save_counts_row != None: smutab2[is2_bound][iangle2][save_counts_row] += rat2
!                    # diff angle
!                    rat3 = rats*(1-ratangle)
!                    if iangle2_bound < nummu2:
!                      for row3 in compute_rows:
!                        smutab2[is2][iangle2_bound][row3] += smutabstd[is1][iangle1][row3]*rat3
!                      if save_counts_row != None: smutab2[is2][iangle2_bound][save_counts_row] += rat3
!                    # diff s and diff angle
!                    rat4 = (1-rats)*(1-ratangle)
!                    if iangle2_bound < nummu2 and is2_bound < nums2:
!                      for row3 in compute_rows:
!                        smutab2[is2_bound][iangle2_bound][row3] += smutabstd[is1][iangle1][row3]*rat4
!                      if save_counts_row != None: smutab2[is2_bound][iangle2_bound][save_counts_row] += rat4

      enddo
    enddo

!                            
!    if div_counts_rows != [] and save_counts_row != None:
!        for is2 in range(nums2):
!            for iangle2 in range(nummu2):
!                for row3 in div_counts_rows:
!                    smutab2[is2][iangle2][row3] /= smutab2[is2][iangle2][save_counts_row]
!    return smutab2
  end subroutine DSMapping

  ! standard function of \int xi ds as a function of mu
  subroutine XiFun_std(smutab, deltas, numnbins, nummubins, anglemin, anglemax, smin, smax, nummuedge, intxi, mumids)
    ! argument
    integer, intent(in)  :: numnbins, nummubins, nummuedge
    real(rt), intent(in) :: smutab(numnbins,nummubins,3), deltas,  &
      anglemin,anglemax, smin,smax ! range of s,mu to be integrated
    real(rt), intent(out) :: intxi(nummuedge-1)
    real(rt), intent(out), optional :: mumids(nummuedge-1) ! middle value of mu in each mu-bin 
    ! variable
    real(rt) :: deltamu, angleindices(nummuedge), dangle, dangleindex, &
      DD,DR,RR,xi, dk,k1,k2
    integer :: i,j,k,sindex1, sindex2
    logical :: testprint=.false.

    intxi = 0.0_rt
    deltamu = 1.0 / dble(nummubins)
    
    ! first of all, decide the fractional bins in angular space
    dangle = (anglemax-anglemin) / dble(nummuedge-1)
    dangleindex = dangle/deltamu
    angleindices(1) = anglemin / deltamu
    do j = 2, nummuedge
      angleindices(j) = angleindices(j-1) + dangleindex
    enddo
    if(testprint) print *, ' (XiFun_std) real(angleindices) = ', real(angleindices)
    do i = 1, nummuedge
      if(testprint) print *, ' (XiFun_std) i, angleindics(i)*deltamu = ', i, angleindices(i) * deltamu
    enddo
    
    ! This only holds for deltas eq 1!
    sindex1 = int(smin / deltas+1.5)
    sindex2 = int(smax / deltas+0.5)
    if(testprint) print *,  ' (XiFun_std) smin, smax, deltas, sindex1, sindex2 = ', &
       smin, smax, deltas, sindex1, sindex2
    if(abs(deltas-1.0).gt.1.0e-5) then
      print *, ' (XiFun) Warning! deltas not equal to 1. May have significant error in the range of integral!'
      print *, ' deltas = ' , deltas
    endif
    if(testprint) print *, 'sindex1, sindex2 = ', sindex1, sindex2
        
    if(sindex2 > numnbins) then
      print *, ' (XiFun) WARNING!! index of s outflow: sindex, numnbins = ', sindex2, numnbins, &
        '; forcing sindex2 equal to numnbins'
      stop
      sindex2 = numnbins
    endif

    if(present(mumids)) then
      do j = 1,nummuedge-1
        k1 = floor(angleindices(j)+1+0.00001) 
        k2 = floor(angleindices(j+1)+0.00001)
        mumids(j) = dble((k1-1) + k2)/2.0_rt*deltamu
      enddo
    endif

    do i = sindex1, sindex2    
      if(testprint) print *, 'Doing integration for s = ', i
    if (i.eq.sindex1.and.testprint) then
      do j =1,nummuedge-1
        k1 = angleindices(j)
        k2 = angleindices(j+1)
        print *,  floor(k1+1+0.00001), floor(k2+0.00001)
      enddo
    endif
    do j = 1, nummuedge-1
      DD=0.0_rt; DR=0.0_rt; RR=0.0_rt;
      k1 = angleindices(j)
      k2 = angleindices(j+1)
      if(i.eq.sindex1.and.testprint) print *, j
      ! summation of counts: the whole bins
      if (.false.) then
       if(i.eq.sindex1 ) print *, 'k1,k2 = ', k1,k2
       if(i.eq.sindex1 ) print *, 'ceiling(k1),floor(k2) = ', ceiling(k1),floor(k2)
       do k = ceiling(k1),floor(k2)
         if(i.eq.sindex1.and.testprint) print *, 'Integrating: ', k
         DD = DD + smutab(i,k,1)
         DR = DR + smutab(i,k,2)
         RR = RR + smutab(i,k,3)
       enddo
       ! summatin of counts: the fractional bins (tails)
       dk = (ceiling(k1)-k1)
       if(i.eq.sindex1.and.testprint ) print *, 'ceiling(k1)-k1 = ', dk
       DD = DD + smutab(i,ceiling(k1)-1,1)*dk
       DR = DR + smutab(i,ceiling(k1)-1,2)*dk
       RR = RR + smutab(i,ceiling(k1)-1,3)*dk
       dk = (k2-floor(k2))
       if(i.eq.sindex1.and. testprint) print *, 'k2-floor(k2) = ', dk
       DD = DD + smutab(i,floor(k2)+1,1)*dk
       DR = DR + smutab(i,floor(k2)+1,2)*dk
       RR = RR + smutab(i,floor(k2)+1,3)*dk
      else
      ! suppose k1 = 3.5; k2=4.5; then floor(k1)+1 = 4; floor(k-1)+1 = 4
      ! suppose k = 3; then floor(k)+1=4; floor(k-1)+1=3
       do k =  floor(k1+1+0.00001),floor(k2+0.00001)
         DD = DD + smutab(i,k,1)
         DR = DR + smutab(i,k,2)
         RR = RR + smutab(i,k,3)
       enddo
      endif
              
      xi = (DD-2.0_rt*DR) / RR + 1.0_rt
      intxi(j) = intxi(j) + xi 
    enddo
    enddo
  end subroutine XiFun_std
  
  ! \int xi ds as a function of mu, generalized form (for 2D contour analysis)
  ! 2bin
  subroutine XiFun(smutab, deltas, numnbins, nummubins, anglemin, anglemax, smin, smax, nummuedge, intxi, mumids)
  !!!!! To be check!!!! Very , very , very == wrong!!! Check the many integers & binnumbers !!!
    ! argument
    integer, intent(in)  :: numnbins, nummubins, nummuedge!(nummubins / gbnumsbin)
    real(rt), intent(in) :: smutab(numnbins,nummubins,3), deltas,  &
      anglemin,anglemax, smin,smax ! range of s,mu to be integrated
    real(rt), intent(out) :: intxi(nummuedge-1)
    real(rt), intent(out), optional :: mumids(nummuedge-1) ! middle value of mu in each mu-bin 
    ! variable
    real(rt) :: deltamu, angleindices(nummuedge), dangle, dangleindex, mids, k1,k2, sedges(1000)
!      DD,DR,RR,xi, dk,k1,k2, mids, tmpintxis(nummuedge-1,2)
    integer :: actual_nummuedge, i,j,i1,i2
    logical :: testprint=.false.


!!    print *, 'Number of bins: numnbins, mumubins, nummuedge = ', numnbins, nummubins, nummuedge
!    print *, 'angle range: ', anglemin, anglemax

    intxi = 0.0_rt
    deltamu = 1.0 / dble(nummubins)
    
    ! first of all, decide the fractional bins in angular space
    actual_nummuedge = (nummuedge-1) / gbnumsbin +1
!    print *, 'actual_nummuedge = ', actual_nummuedge
    dangle = (anglemax-anglemin) / dble(actual_nummuedge-1)
    dangleindex = dangle/deltamu
    angleindices(1) = anglemin / deltamu
    do j = 2, actual_nummuedge
      angleindices(j) = angleindices(j-1) + dangleindex
    enddo
    if(testprint) print *, 'real(angleindices) = ', real(angleindices)
    do i = 1, nummuedge
      if(testprint) print *, i, angleindices(i) * deltamu
    enddo    
    
    if(present(mumids)) then
      do j = 1,actual_nummuedge-1
        k1 = floor(angleindices(j)+1+0.00001) 
        k2 = floor(angleindices(j+1)+0.00001)
        mumids(j) = dble((k1-1) + k2)/2.0_rt*deltamu
      enddo
!      print *, 'mumids = ', mumids
    endif

!!    call XiFun_std(smutab, deltas, numnbins, nummubins, anglemin, anglemax, smin, smax, nummuedge, intxi)
!!    print *, 'intxi from standard XiFun: ', intxi
!    return

    
    if(mod((nummuedge-1), gbnumsbin) .ne.0) then
      print *, 'Error (XiFun)! numsbin incompatable with total bins: ', (nummuedge-1), gbnumsbin
      stop
    endif
    
    ! 1. Check: is that 
    
    sedges(1)  = smin 
    sedges(2)  = smax
    sedges(3)  = smax 
    
    intxi = 0.0_rt
    do i = 1, gbnumsbin
      i1 = (actual_nummuedge-1) * (i-1) + 1
      i2 = i1 + actual_nummuedge -2
!      print *, i1, i2, sedges(i), sedges(i+1)
      call XiFun_std(smutab, deltas, numnbins, nummubins, anglemin, anglemax, &
              sedges(i), sedges(i+1), actual_nummuedge, intxi(i1:i2))
!!      print *, 'loop = ', i
!!      print *, 'i2-i1+1, actual_nummuedge = ', i2-i1+1, actual_nummuedge
!!      print *, 'intxi(i1:i2)  = ', intxi(i1:i2)
    enddo
!    print *, '**Finaly, intxi = ', intxi
!    stop
!11.102650094483595        10.485374549922113        11.259850842068577        12.053443999700688        12.886048021969620        13.666540539278243        14.144060785526161        14.471828005713531        14.742092741015004        15.039291591991976        1.4832499704521103        1.6759450472752717        2.0541168689675517        2.3899646754337582        2.9534279613576135        3.3712476554309392        3.7423458309859017        3.9567703564102974        3.9376723411126049        3.9704822516105036      
!11.102650094483595        10.485374549922113        11.259850842068577        12.053443999700688        12.886048021969620        13.666540539278243        14.144060785526161        14.471828005713531        14.742092741015004        15.039291591991976        1.4832499704521103        1.6759450472752717        2.0541168689675517        2.3899646754337582        2.9534279613576135        3.3712476554309392        3.7423458309859017        3.9567703564102974        3.9376723411126049        3.9704822516105036 
!12.585900064935709        12.161319597197386        13.313967711036131        14.443408675134448        15.839475983327233        17.037788194709179        17.886406616512065        18.428598362123829        18.679765082127609        19.009773843602481
  end subroutine XiFun

! normalize an array: amplitude shifted to 1; skip the LAST element
  subroutine normfun(A,nA,Anormed)
    integer, intent(in) :: nA
    real(rt), intent(in) :: A(nA)
    real(rt), intent(out) :: Anormed(nA-1)
    real(rt) :: avg
    integer :: imid, i
    avg = sum(A) / dble(nA)
!    i = nA / 2
!    avg = sum(A(1:i)) / dble(i)
    imid = nA / 2 + 1
!    do i = 1, imid-1
    do i = 1, nA-1
      Anormed(i) = A(i) / avg
    enddo
!    do i = imid, nA
!      Anormed(i-1) = A(i) / avg
!    enddo
  end subroutine normfun
  
  subroutine smu_ximu_CalcOmWChisqs(& 
    omlist, numom, wlist, numw, & ! List of omegam, w
    outputdir, baseoutputfile, & ! Basic name of the outputfile
    omstds, wstds, numomwstds, & ! "standard" values of omegam, w: baseline cosmology for (s,mu) mapping. Stored in data2pcffile_base
    weightedstds, avg_counts & ! assign different weights to the "standard" cosmologies
    )
    !arguments
    integer, intent(in) :: numom, numw, numomwstds
    real(rt), intent(in) :: omlist(numom), wlist(numw),  &
       omstds(numomwstds), wstds(numomwstds)
    character(*), intent(in) :: outputdir, baseoutputfile
    logical, intent(in) :: weightedstds, avg_counts
    
    ! variables
    real(rt) :: omwlist(2,numw,numom), &
      smutabstds(nbins_database,mubins_database,3,nz,numomwstds), & !MF
      smutab_data(nbins_data,mubins_data,3),  smutab_data1(nbins_data,mubins_data,3) !MF
    real(rt) :: smin_mapping=1.0_rt, smax_mapping=50_rt, DAstd,DAnew,Hstd,Hnew, deltas1,deltas2,t0,t1,t2,dt,&
      intxis(maxval(mubins),N1,N2,nz,numomwstds), intxi(maxval(mubins)), dintxi(maxval(mubins)), & !MF: intxis
      sumDAHweights,DAHweights(numomwstds,nz), &
      factD,factM1,factM2, chisqs_nosyscor(n1,n2,nz-1), chisqs_syscor(n1,n2,nz-1), &
      chisqs_nosyscor_all(nz-1), chisqs_syscor_all(nz-1)
    type(omwpar) :: parstd, parnew
    character(len=1000) :: tmpstr,tmpstr1,tmpstr2,tmpstr3,tmpstr4,tmpstr5,tmpstr6,tmpstr7,tmpstr8,tmpstr9,&
      tmpstr10,tmpstr11,tmpstr12,tmpstr13,tmpstr14,&
      nowchisqstr,filenames(N1,N2), nowfile, filename_allsch
    integer :: nowfileunit, fileunits(N1,N2), fileunit_allsch, i,j,k,iz,iomwstds,iomwstds2,i1,i2, n,iom,iw,iomw
    integer, parameter :: basefileunit = 58483
    logical :: debug_calcchisq = .false.
    
 
    if(NBComp) then
      print *, 'ERROR! (smu_ximu_CalcOmW3hisqs) NBComp not supported! Found NBComp = ', NBComp 
      stop
    endif
    if(Forecast) then
      print *, 'ERROR! (smu_ximu_CalcOmWChisqs) Forecast not supported! Found Forecast = ', Forecast 
      stop
    endif
    ! 1. Initialization
    ! 1.1 Check covmat files ready or not
    do i=1,N1
    do j=1,N2
    do iz=2,nz
      if (.not.allocated(covmats(i,j,iz-1)%A).or.covmats(i,j,iz-1)%nA.ne.mubins(i)-1) then
        print *, ' (smu_ximu_CalcOmWChisqs) ERROR! Covmats not ready: iz,i,j,mubin,mucut = ', iz,i,j,mubins(i),mucuts(j)
        stop
      endif
    enddo
    enddo
    enddo
    ! 1.2 Load in 2pCF computed in standard cosmology
    do iomwstds = 1, numomwstds
    do iz = 1, nz
      call ximu_loadsmufile(data2pcffile_base(iz,omstds(iomwstds),wstds(iomwstds)), smutabstds(:,:,:,iz,iomwstds), &
              smax_database, nbins_database, mubins_database)
    enddo
    enddo

    ! 2. Names of files storing chisq values
    write(tmpstr1,*) ncovmocks
    nowchisqstr = trim(adjustl(baseoutputfile)) // '.CovMock_'//trim(adjustl(tmpstr1)) 
    nowfileunit = basefileunit
    do i=1,N1
    do j=1,N2
      tmpstr1=''; write(tmpstr1, *) mubins(i) ! how many bins in mu space
      write(tmpstr2, '(f5.2)') mucuts(j) ! maximal value of mu
      write(tmpstr3, '(f5.1)') ints1
      write(tmpstr4, '(f5.1)') ints2
      write(tmpstr5, *) ncovmocks
      tmpstr6=''; if(polyfitdeg.ge.1) write(tmpstr6,'(A,i1)') '.polyfitdeg',polyfitdeg
      write(tmpstr7,*) nbins_database 
      write(tmpstr8,*) mubins_database; 
      write(tmpstr11,*) numomwstds; 
      tmpstr12 = '';
      do iomwstds = 1, 1
        tmpstr12 = trim(adjustl(tmpstr12))//'_'//trim(adjustl(omwstr(omstds(iomwstds),wstds(iomwstds))))
      enddo
      tmpstr13=''; if(weightedstds) tmpstr13 = '.weightedstds'
      tmpstr9='.baseline_'//trim(adjustl(tmpstr11))//'omws_'//trim(adjustl(tmpstr12))// &
              '__nbins'//trim(adjustl(tmpstr7))//'_mubins'//trim(adjustl(tmpstr8))
      tmpstr14=''; if(avg_counts) write(tmpstr14,*) '.avg_counts'
      nowfile = trim(adjustl(outputdir))//'/'//trim(adjustl(nowchisqstr))//'__'//trim(adjustl(tmpstr1)) &
        //'mubins.mumax' // trim(adjustl(tmpstr2))  &
        //'.CovMock_'//trim(adjustl(tmpstr5)) &
        //'.s'//trim(adjustl(tmpstr3))//'to'//trim(adjustl(tmpstr4))//trim(adjustl(tmpstr6)) &
        //trim(adjustl(tmpstr9))//trim(adjustl(tmpstr13))//trim(adjustl(tmpstr14)) &
        //'.txt'

      tmpstr12 = '';
      do iomwstds = 1, numomwstds
        tmpstr12 = trim(adjustl(tmpstr12))//' '//trim(adjustl(omwstr(omstds(iomwstds),wstds(iomwstds))))
      enddo
        
      filenames(i,j) = nowfile; fileunits(i,j) = nowfileunit
      if(output_sep_schemes) then
        print *, ' (smu_ximu_CalcOmWChisqs) Opening file for output: ', &
         trim(adjustl(nowfile))
        open(unit=fileunits(i,j),file=filenames(i,j))
        write(fileunits(i,j),'(A)') '### mumin  omw   chisq_nosyscor  chisq_syscor   chisqs_nosyscor   chisqs_syscor'&
          //';  base line cosmologies: '//trim(adjustl(tmpstr12))
      endif
      nowfileunit = nowfileunit+1
    enddo
    enddo
    tmpstr3=''; write(tmpstr3,*) N1*N2
    fileunit_allsch=fileunits(N1,N2)+1; 
    filename_allsch = trim(adjustl(filenames(1,1)))//'.'//trim(adjustl(tmpstr3))//'schemes'
    open(unit=fileunit_allsch,file=filename_allsch)
    print *, ' (smu_ximu_CalcOmWChisqs) Opening file for output: ', trim(adjustl(filename_allsch))
    write(fileunit_allsch,'(A)') '### mumin  omw   chisq_nosyscor  chisq_syscor   chisqs_nosyscor   chisqs_syscor'    
    
    ! 3. Compute chisq values
    call cpu_time(t0); t1=t0
    dt = 10.0; iomw=0;    
    do iom = 1, numom
    do iw  = 1, numw
      
      ! 3.1 Compute intxi
      iomw=iomw+1; call cpu_time(t2)
      if (t2-t1.gt.dt) then
        write(*,'(f10.1,A,i5,A,f4.1,A)') t2-t0, ' seconds passed.   #-(om,w) = ', &
           iomw, ' (',100*float(iomw)/float(numom*numw),'%)'
        t1=t2
      endif


      if(weightedstds) then
        do iomwstds = 1, numomwstds
          do iz  = 1, nz
            parstd%omegam = omstds(iomwstds); parstd%w = wstds(iomwstds)
            parnew%omegam = omlist(iom); parnew%w=wlist(iw)
            DAstd = DAz_wcdm(parstd,zeffs(iz)); Hstd = Hz_wcdm(parstd,zeffs(iz))
            DAnew = DAz_wcdm(parnew,zeffs(iz)); Hnew = Hz_wcdm(parnew,zeffs(iz))
            DAHweights(iomwstds,iz) = 1.0 / ((DAnew/DAstd - 1.0)**2.0 + (Hnew/Hstd - 1.0)**2.0 + 0.0001)
          enddo
        enddo
        do iz = 1, nz
          sumDAHweights = sum(DAHweights(1:numomwstds,iz))
          do iomwstds = 1, numomwstds
            DAHweights(iomwstds,iz) = DAHweights(iomwstds,iz) * numomwstds / sumDAHweights
          enddo
        enddo
      else
        DAHweights(:,:) = 1.0_rt
      endif

      do iomwstds = 1, numomwstds
      if(avg_counts.and.iomwstds.gt.1) cycle
      do iz  = 1, nz
        ! 3.1 smutabdata obtained from DSMapping
        deltas1 = smax_database / float(nbins_database) 
        deltas2 = smax_data / float(nbins_data)
!       print *, ' (smu_ximu_CalcOmWChisqs) Do DSMapping... iz = ', iz
        if(.not.avg_counts) then
          parstd%omegam = omstds(iomwstds); parstd%w = wstds(iomwstds)
          parnew%omegam = omlist(iom); parnew%w=wlist(iw)
          DAstd = DAz_wcdm(parstd,zeffs(iz)); Hstd = Hz_wcdm(parstd,zeffs(iz))
          DAnew = DAz_wcdm(parnew,zeffs(iz)); Hnew = Hz_wcdm(parnew,zeffs(iz))
          call DSMapping(smutabstds(:,:,:,iz,iomwstds), nbins_database, mubins_database, smutab_data, &
            nbins_data, mubins_data, DAstd, DAnew, Hstd, Hnew, deltas1,  deltas2,  smin_mapping, smax_mapping)
        else
          ! in case of avg_counts, sum up all smutabstds to build one
          parnew%omegam = omlist(iom); parnew%w=wlist(iw)
          DAnew = DAz_wcdm(parnew,zeffs(iz)); Hnew = Hz_wcdm(parnew,zeffs(iz))
          smutab_data = 0.0_rt
          do iomwstds2 = 1, numomwstds
            parstd%omegam = omstds(iomwstds2); parstd%w = wstds(iomwstds2)
            DAstd = DAz_wcdm(parstd,zeffs(iz)); Hstd = Hz_wcdm(parstd,zeffs(iz))
            call DSMapping(smutabstds(:,:,:,iz,iomwstds2), nbins_database, mubins_database, smutab_data1, &
              nbins_data, mubins_data, DAstd, DAnew, Hstd, Hnew, deltas1,  deltas2,  smin_mapping, smax_mapping)
            smutab_data(:,:,:) = smutab_data(:,:,:) + &
              smutab_data1(:,:,:) * DAHweights(iomwstds2,iz) / dble(numomwstds)
          enddo
        endif
          
!       print *, ' (smu_ximu_CalcOmWChisqs) DSMapping done: iz = ', iz
        if(debug_calcchisq ) then !.and. mubins(i).eq.25 .and.j.eq.1
                print *, '####################################'
                print *, 'iz = ',iz
                print *, "Checking: DAstd, DAnew, Hstd, Hnew = ", DAstd, DAnew, Hstd, Hnew
                do i = 30, 31
                  j = 1
                  print *, i, smutab_data(i,j,1:3), (smutab_data(i,j,1)-2*smutab_data(i,j,2))/smutab_data(i,j,3)+1.0_rt
                enddo        
        endif
        
        do i1 = 1, N1
        do i2 = 1, N2
!         print *, ' (smu_ximu_CalcOmWChisqs) Compute \int xi...'
          call XiFun(smutab_data, deltas2, nbins_data, mubins_data, &
            anglemin=1.0_rt-mucuts(i2), anglemax=1.0_rt, &
            smin=ints1, smax=ints2, &
            nummuedge=mubins(i1)+1, &
            intxi=intxi(1:mubins(i1)))
!         print *, ' (smu_ximu_CalcOmWChisqs) Normalize \int xi...'
          call normfun(intxi(1:mubins(i1)),mubins(i1),intxis(1:mubins(i)-1,i1,i2,iz,iomwstds)) 
          if(debug_calcchisq .and. mubins(i1).eq.25 .and.i2.eq.1) then
                  print *, 'intxi = ', real(intxi(1:25))
                print *
                print *, 'intxi (normed) = ', real(intxis(1:24,i1,i2,iz,iomwstds))
          endif
        enddo!        do i1 = 1, N1
        enddo!        do i2 = 1, N2
!       print *, ' (smu_ximu_CalcOmWChisqs) Compute intxi done: iz = ', iz
      enddo!      do iz  = 1, nz
      enddo!do iomwstds = 1, numomwstds
      
      chisqs_nosyscor_all=0.0_rt; chisqs_syscor_all = 0.0_rt;
      chisqs_nosyscor = 0.0_rt; chisqs_syscor = 0.0_rt

      ! 3.2 Compute chisqs
      do i1 = 1,N1
      do i2 = 1,N2
        call Percival_cs(ncovmocks,mubins(i1),2,factD,factm1,factm2)
        do iz = 2, nz
          n = mubins(i1)
          ! Loop for omstd/wstd
          do iomwstds = 1, numomwstds
            if(avg_counts.and.iomwstds.gt.1) cycle
            ! chisq before sys cor
            dintxi(1:n-1) = intxis(1:n-1,i1,i2,iz,iomwstds) - intxis(1:n-1,i1,i2,1,iomwstds)
            if(debug_calcchisq .and. mubins(i1).eq.25 .and.i2.eq.1) then
                    print *, 'iz = ', iz
                  print *
                  print *, 'dintxi = ', dintxi(1:n-1)
            endif
            if(.not.avg_counts) then
              chisqs_nosyscor(i1,i2,iz-1) = chisqs_nosyscor(i1,i2,iz-1) + &
                  chisq_cov_xbar(dintxi(1:n-1), covmats(i1,i2,iz-1)%A, n-1)/dble(numomwstds)*DAHweights(iomwstds,iz)/factm2
            else
              chisqs_nosyscor(i1,i2,iz-1) = &
                  chisq_cov_xbar(dintxi(1:n-1), covmats(i1,i2,iz-1)%A, n-1)/factm2
            endif
            chisqs_nosyscor_all(iz-1) = chisqs_nosyscor_all(iz-1) + chisqs_nosyscor(i1,i2,iz-1)/ dble(N1*N2)
            ! chisq after sys cor
            dintxi(1:n-1) = dintxi(1:n-1) - dintxi_syscor(1:n-1,i1,i2,iz-1)
            if(.not.avg_counts) then
              chisqs_syscor(i1,i2,iz-1) = chisqs_syscor(i1,i2,iz-1) + &
                    chisq_cov_xbar(dintxi(1:n-1), covmats(i1,i2,iz-1)%A, n-1)/dble(numomwstds)*DAHweights(iomwstds,iz)/factm2
            else
              chisqs_syscor(i1,i2,iz-1) =  &
                    chisq_cov_xbar(dintxi(1:n-1), covmats(i1,i2,iz-1)%A, n-1)/factm2
            endif
            chisqs_syscor_all(iz-1)   = chisqs_syscor_all(iz-1) + chisqs_syscor(i1,i2,iz-1) / dble(N1*N2)
            if(debug_calcchisq .and. mubins(i1).eq.25 .and.i2.eq.1) then
                  print *, 'dintxi_sys = ', dintxi_syscor(1:n-1,i1,i2,iz-1)
                  print *, 'dintxi (after syscor) = ', dintxi(1:n-1)
                  !print *, 'invcov = ', real(covmats(i1,i2,iz-1)%A)
                  print *, 'chisq / chisq_syscor = ', chisqs_nosyscor(i1,i2,iz-1), chisqs_syscor(i1,i2,iz-1)
            endif
          enddo!do iomwstds = 1, numomwstds
        enddo!do iz = 2, nz

        ! write chisq values to files...
        if(output_sep_schemes) then
          tmpstr=''; tmpstr1=''
          write(tmpstr,'(f5.2,1x,A,e15.7,e15.7,5x)') 1.0_rt-mucuts(i2), trim(adjustl(omwstr(omlist(iom),wlist(iw)))), &
            sum(chisqs_nosyscor(i1,i2,1:nz-1)), sum(chisqs_syscor(i1,i2,1:nz-1))
          if(output_chisq_eachredbin) then
            do k = 1,nz-1
              write(tmpstr1,'(e15.7)') chisqs_nosyscor(i1,i2,k)
              tmpstr = trim(adjustl(tmpstr))//' '//trim(adjustl(tmpstr1))
            enddo
            do k = 1,nz-1
              write(tmpstr1,'(e15.7)') chisqs_syscor(i1,i2,k)
              tmpstr = trim(adjustl(tmpstr))//' '//trim(adjustl(tmpstr1))
            enddo
          endif
          write(fileunits(i1,i2),'(A)') trim(adjustl(tmpstr))
        endif
      enddo
      enddo
      ! output allscheme chisq value
      tmpstr=''; tmpstr1=''
      write(tmpstr,'(f5.2,1x,A,e15.7,e15.7,5x)') 1.0_rt-mucuts(i2), trim(adjustl(omwstr(omlist(iom),wlist(iw)))), &
        sum(chisqs_nosyscor_all(1:nz-1)), sum(chisqs_syscor_all(1:nz-1))
!     print *, tmpstr
      if(output_chisq_eachredbin) then
        do k = 1,nz-1
          write(tmpstr1,'(e15.7)') chisqs_nosyscor_all(k)
          tmpstr = trim(adjustl(tmpstr))//' '//trim(adjustl(tmpstr1))
        enddo
        do k = 1,nz-1
          write(tmpstr1,'(e15.7)') chisqs_syscor_all(k)
          tmpstr = trim(adjustl(tmpstr))//' '//trim(adjustl(tmpstr1))
        enddo
      endif
      write(fileunit_allsch,'(A)') trim(adjustl(tmpstr))
    enddo
    enddo
    
    do i1 = 1,n1
    do i2 = 1,n2
      if(output_sep_schemes) close(fileunits(i1,i2))
    enddo
    enddo
    close(fileunit_allsch)
!    chisqs_nosyscor(N1,N2)
!    chisqs_syscor(A1,N2)
  end subroutine smu_ximu_CalcOmWChisqs

  subroutine smu_ximu_CalcDAHChisqs(& 
    DAs, Hs, & ! List of omegam, w
    !outputdir, baseoutputfile, & ! Basic name of the outputfile
    omstds, wstds, numomwstds, & ! "standard" values of omegam, w: baseline cosmology for (s,mu) mapping. Stored in data2pcffile_base
    smutabstds, smutabstds_inited, & ! xi(s,mu) table of baseline cosmologies
    chisqs_nosyscor, chisqs_syscor, & ! values of chisqs, separate schemes
    chisqs_nosyscor_all, chisqs_syscor_all, &! values of chisqs, averaged over all schemes, correction factor for covmat (D, m1, m2) considered
    weightedstds, avg_counts, VolumeWeightedDAH & ! assign different weights to the "standard" cosmologies
    )
    !arguments
    integer, intent(in) ::  numomwstds
    real(rt), intent(in) :: DAs(nz), Hs(nz), & ! values of DA, H at the effective redshifts
       omstds(numomwstds), wstds(numomwstds)
!    character(*), intent(in) :: outputdir, baseoutputfile
    real(rt), intent(inout) :: smutabstds(nbins_database,mubins_database,3,nz,numomwstds)
    logical, intent(inout) :: smutabstds_inited
    real(rt), intent(out) :: chisqs_nosyscor(n1,n2,nz-1), chisqs_syscor(n1,n2,nz-1), &
       chisqs_nosyscor_all(nz-1), chisqs_syscor_all(nz-1)
    logical, intent(in) :: weightedstds, avg_counts
    logical, intent(in), optional :: VolumeWeightedDAH
    ! variables
    real(rt) :: smutab_data(nbins_data,mubins_data,3), smutab_data1(nbins_data,mubins_data,3)
    real(rt) :: smin_mapping=min(1.0_rt,ints1*0.8_rt), smax_mapping=ints2*1.2_rt, &
      DAstd,DAnew,Hstd,Hnew, deltas1,deltas2, sumDAHweights,DAHweights(numomwstds,nz), &
      DAvalues(numgal_nbin),Hvalues(numgal_nbin),DAarray(numgal_nbin),Harray(numgal_nbin), &
      DAstds(nz,numomwstds),Hstds(nz,numomwstds),&
      intxis(maxval(mubins),N1,N2,nz,numomwstds), intxi(maxval(mubins)), dintxi(maxval(mubins)), & 
      factD,factM1,factM2
    type(omwpar) :: parstd
!    character(len=1000) :: tmpstr,tmpstr1,tmpstr2,tmpstr3,tmpstr4,tmpstr5,tmpstr6,tmpstr7,tmpstr8,tmpstr9,&
!      tmpstr10,tmpstr11,tmpstr12,&
!      nowchisqstr,filenames(N1,N2), nowfile, filename_allsch
    integer ::  i,j,k,iz,iomwstds,iomwstds2,i1,i2,i_numgal, n
!    integer, parameter :: basefileunit = 58483
    logical :: numgal_weighted=.false., debug_calcchisq = .false.
    
    if(present(VolumeWeightedDAH)) then
      numgal_weighted=VolumeWeightedDAH
    else
      numgal_weighted=.false.
    endif
    
    if(ForeCast) then
      print *, '(smu_ximu_CalcDAHChisqs) ERROR! ForeCast not suppoted! ForeCast = ', Forecast
      stop
      if(numomwstds .ne. 1) then
        print *, '(smu_ximu_CalcDAHChisqs) ERROR! For forecast standard-omw should be 1: numomwstds = ', numomwstds
        stop
      endif
      if(weightedstds) then
        print *, '(smu_ximu_CalcDAHChisqs) ERROR! For forecast weightedstds not supported: weightedstds = ', weightedstds
        stop
      endif
      if(avg_counts) then
        print *, '(smu_ximu_CalcDAHChisqs) ERROR! For forecast avg_counts not supported: avg_counts = ', avg_counts
        stop
      endif
    endif
        

     do iomwstds = 1, numomwstds
       parstd%omegam = omstds(iomwstds); parstd%w = wstds(iomwstds)     
       if(numgal_weighted) then
         if(.not. numgal_inited) call numgal_init()
         do i = 1, numgal_nbin
           DAvalues(i) = DAz_wcdm(parstd,numgal_zcenters(i))
           Hvalues(i)  = Hz_wcdm (parstd,numgal_zcenters(i))
         enddo
        endif
        do iz = 1, nz
            if(.not.numgal_weighted) then
                    DAstds(iz,iomwstds) = DAz_wcdm(parstd,zeffs(iz)); Hstds(iz,iomwstds) = Hz_wcdm(parstd,zeffs(iz))
            else
                      DAarray = DAvalues(1:numgal_nbin)*redbin_weights(iz,1:numgal_nbin)
                      Harray  = Hvalues(1:numgal_nbin) *redbin_weights(iz,1:numgal_nbin)
                      DAstds(iz,iomwstds) = sum(DAarray(1:numgal_nbin)) / sum(redbin_weights(iz,1:numgal_nbin))
                      Hstds(iz,iomwstds)  = sum(Harray(1:numgal_nbin))  / sum(redbin_weights(iz,1:numgal_nbin))
                      !print *, parstd%omegam, parstd%w
                      !print *, 'DAs:', DAstds(iz,iomwstds), DAz_wcdm(parstd,zeffs(iz))
                      !print *, 'Hs:',  Hstds(iz,iomwstds),  Hz_wcdm(parstd,zeffs(iz))
            endif
       enddo
     enddo

    ! 1. Initialization
    ! 1.1 Check covmat files ready or not
    do i=1,N1
    do j=1,N2
    do iz=2,nz
      if (.not.allocated(covmats(i,j,iz-1)%A).or.covmats(i,j,iz-1)%nA.ne.mubins(i)-1) then
        print *, ' (smu_ximu_CalcDAHChisqs) ERROR! Covmats not ready: iz,i,j,mubin,mucut = ', iz,i,j,mubins(i),mucuts(j)
        stop
      endif
    enddo
    enddo
    enddo
    
    ! 1.2 Load in 2pCF computed in standard cosmology
    if (.not.smutabstds_inited) then
      write(*,'(A)') ' (smu_ximu_CalcDAHChisqs) Load in the xi(s,mu) of baseline cosmologies...'
      do iomwstds = 1, numomwstds
      do iz = 1, nz
        call ximu_loadsmufile(data2pcffile_base(iz,omstds(iomwstds),wstds(iomwstds)), smutabstds(:,:,:,iz,iomwstds), &
                smax_database, nbins_database, mubins_database)
      enddo
      enddo
      smutabstds_inited =.true.
    endif

    if(weightedstds) then
        do iomwstds = 1, numomwstds
          do iz  = 1, nz
            parstd%omegam = omstds(iomwstds); parstd%w = wstds(iomwstds)
            DAstd=DAstds(iz,iomwstds); Hstd=Hstds(iz,iomwstds)
            DAnew = DAs(iz); Hnew = Hs(iz)
            DAHweights(iomwstds,iz) = 1.0 / ((DAnew/DAstd - 1.0)**2.0 + (Hnew/Hstd - 1.0)**2.0 + 0.0001)
          enddo
        enddo
        do iz = 1, nz
          sumDAHweights = sum(DAHweights(1:numomwstds,iz))
          do iomwstds = 1, numomwstds
            DAHweights(iomwstds,iz) = DAHweights(iomwstds,iz) * numomwstds / sumDAHweights
          enddo
        enddo
    else
        DAHweights(:,:) = 1.0_rt
    endif


    do iomwstds = 1, numomwstds
      if(avg_counts.and.iomwstds.gt.1) cycle
      do iz  = 1, nz
        ! 3.1 smutabdata obtained from DSMapping
        deltas1 = smax_database / float(nbins_database) 
        deltas2 = smax_data / float(nbins_data)
!       print *, ' (smu_ximu_CalcOmWChisqs) Do DSMapping... iz = ', iz
        if(.not.avg_counts) then
          parstd%omegam = omstds(iomwstds); parstd%w = wstds(iomwstds)
          !DAstd = DAz_wcdm(parstd,zeffs(iz)); Hstd = Hz_wcdm(parstd,zeffs(iz))
          DAstd=DAstds(iz,iomwstds); Hstd=Hstds(iz,iomwstds)
          DAnew = DAs(iz); Hnew = Hs(iz)
          call DSMapping(smutabstds(:,:,:,iz,iomwstds), nbins_database, mubins_database, smutab_data, &
            nbins_data, mubins_data, DAstd, DAnew, Hstd, Hnew, deltas1,  deltas2,  smin_mapping, smax_mapping)
        else
          ! in case of avg_counts, sum up all smutabstds to build one
          DAnew = DAs(iz); Hnew = Hs(iz)
          smutab_data = 0.0_rt
          do iomwstds2 = 1, numomwstds
            parstd%omegam = omstds(iomwstds2); parstd%w = wstds(iomwstds2)
            !DAstd = DAz_wcdm(parstd,zeffs(iz)); Hstd = Hz_wcdm(parstd,zeffs(iz))
            DAstd=DAstds(iz,iomwstds); Hstd=Hstds(iz,iomwstds)
            call DSMapping(smutabstds(:,:,:,iz,iomwstds2), nbins_database, mubins_database, smutab_data1, &
              nbins_data, mubins_data, DAstd, DAnew, Hstd, Hnew, deltas1,  deltas2,  smin_mapping, smax_mapping)
            smutab_data(:,:,:) = smutab_data(:,:,:) + &
              smutab_data1(:,:,:) * DAHweights(iomwstds2,iz) / dble(numomwstds)
          enddo
        endif
!       print *, ' (smu_ximu_CalcDAHChisqs) DSMapping done: iz = ', iz
        if(debug_calcchisq ) then !.and. mubins(i).eq.25 .and.j.eq.1
                print *, '####################################'
                print *, 'iz = ',iz
                print *, "Checking: DAstd, DAnew, Hstd, Hnew = ", DAstd, DAnew, Hstd, Hnew
                do i = 1, 15
                  j = 1
                  print *, i, smutab_data(i,j,1:3), (smutab_data(i,j,1)-2*smutab_data(i,j,2))/smutab_data(i,j,3)+1.0_rt
                enddo        
        endif
        
        do i1 = 1, N1
        do i2 = 1, N2
!         print *, ' (smu_ximu_CalcDAHChisqs) Compute \int xi...'
          call XiFun(smutab_data, deltas2, nbins_data, mubins_data, &
            anglemin=1.0_rt-mucuts(i2), anglemax=1.0_rt, &
            smin=ints1, smax=ints2, &
            nummuedge=mubins(i1)+1, &
            intxi=intxi(1:mubins(i1)))
!         print *, ' (smu_ximu_CalcDAHChisqs) Normalize \int xi...'
          call normfun(intxi(1:mubins(i1)),mubins(i1),intxis(1:mubins(i1)-1,i1,i2,iz,iomwstds)) 
          if(debug_calcchisq .and. i1.eq.1 .and.i2.eq.1) then
                  print *, 'intxi = ', real(intxi(1:mubins(i1)))
                print *
                print *, 'intxi (normed) = ', real(intxis(1:mubins(i1)-1,i1,i2,iz,iomwstds))
          endif
        enddo!        do i1 = 1, N1
        enddo!        do i2 = 1, N2
!       print *, ' (smu_ximu_CalcDAHChisqs) Compute intxi done: iz = ', iz
      enddo!      do iz  = 1, nz
    enddo!do iomwstds = 1, numomwstds

    chisqs_nosyscor_all=0.0_rt; chisqs_syscor_all = 0.0_rt;
    chisqs_nosyscor = 0.0_rt; chisqs_syscor = 0.0_rt

    ! 3.2 Compute chisqs
    do i1 = 1,N1
      do i2 = 1,N2
        call Percival_cs(ncovmocks,mubins(i1),2,factD,factm1,factm2)
        do iz = 2, nz
          n = mubins(i1)
          ! Loop for omstd/wstd
          do iomwstds = 1, numomwstds
            if(avg_counts.and.iomwstds.gt.1) cycle
            ! chisq before sys cor
            if(.not.NBComp) then
              dintxi(1:n-1) = intxis(1:n-1,i1,i2,iz,iomwstds) - intxis(1:n-1,i1,i2,1,iomwstds)
            else
              dintxi(1:n-1) = intxis(1:n-1,i1,i2,iz,iomwstds) - intxis(1:n-1,i1,i2,iz-1,iomwstds)
            endif
            if(debug_calcchisq .and. mubins(i1).eq.25 .and.i2.eq.1) then
                    print *, 'iz = ', iz
                  print *
                  print *, 'dintxi = ', dintxi(1:n-1)
            endif
            if(.not.avg_counts) then
              chisqs_nosyscor(i1,i2,iz-1) = chisqs_nosyscor(i1,i2,iz-1) + &
                  chisq_cov_xbar(dintxi(1:n-1), covmats(i1,i2,iz-1)%A, n-1)/dble(numomwstds)*DAHweights(iomwstds,iz)/factm2
            else
              chisqs_nosyscor(i1,i2,iz-1) = &
                  chisq_cov_xbar(dintxi(1:n-1), covmats(i1,i2,iz-1)%A, n-1)/factm2
            endif
            chisqs_nosyscor_all(iz-1) = chisqs_nosyscor_all(iz-1) + chisqs_nosyscor(i1,i2,iz-1)/ dble(N1*N2)
            ! chisq after sys cor
            dintxi(1:n-1) = dintxi(1:n-1) - dintxi_syscor(1:n-1,i1,i2,iz-1)
            if(.not.avg_counts) then
              chisqs_syscor(i1,i2,iz-1) = chisqs_syscor(i1,i2,iz-1) + &
                    chisq_cov_xbar(dintxi(1:n-1), covmats(i1,i2,iz-1)%A, n-1)/dble(numomwstds)*DAHweights(iomwstds,iz)/factm2
            else
              chisqs_syscor(i1,i2,iz-1) =  &
                    chisq_cov_xbar(dintxi(1:n-1), covmats(i1,i2,iz-1)%A, n-1)/factm2
            endif
            chisqs_syscor_all(iz-1)   = chisqs_syscor_all(iz-1) + chisqs_syscor(i1,i2,iz-1) / dble(N1*N2)
            if(debug_calcchisq .and. mubins(i1).eq.25 .and.i2.eq.1) then
                  print *, 'dintxi_sys = ', dintxi_syscor(1:n-1,i1,i2,iz-1)
                  print *, 'dintxi (after syscor) = ', dintxi(1:n-1)
                  !print *, 'invcov = ', real(covmats(i1,i2,iz-1)%A)
                  print *, 'chisq / chisq_syscor = ', chisqs_nosyscor(i1,i2,iz-1), chisqs_syscor(i1,i2,iz-1)
            endif
          enddo!do iomwstds = 1, numomwstds
        enddo!do iz = 2, nz

 !       ! write chisq values to files...
 !       if(output_sep_schemes) then
 !         tmpstr=''; tmpstr1=''
 !         write(tmpstr,'(f5.2,1x,A,e15.7,e15.7,5x)') 1.0_rt-mucuts(i2), trim(adjustl(omwstr(omlist(iom),wlist(iw)))), &
 !           sum(chisqs_nosyscor(i1,i2,1:nz-1)), sum(chisqs_syscor(i1,i2,1:nz-1))
 !         if(output_chisq_eachredbin) then
 !           do k = 1,nz-1
 !             write(tmpstr1,'(e15.7)') chisqs_nosyscor(i1,i2,k)
 !             tmpstr = trim(adjustl(tmpstr))//' '//trim(adjustl(tmpstr1))
 !           enddo
 !           do k = 1,nz-1
 !             write(tmpstr1,'(e15.7)') chisqs_syscor(i1,i2,k)
 !             tmpstr = trim(adjustl(tmpstr))//' '//trim(adjustl(tmpstr1))
 !           enddo
 !         endif
 !         write(fileunits(i1,i2),'(A)') trim(adjustl(tmpstr))
 !       endif
      enddo!do i1 = 1,N1
    enddo!do i2 = 1,N2
      ! output allscheme chisq value
!      tmpstr=''; tmpstr1=''
!      write(tmpstr,'(f5.2,1x,A,e15.7,e15.7,5x)') 1.0_rt-mucuts(i2), trim(adjustl(omwstr(omlist(iom),wlist(iw)))), &
!        sum(chisqs_nosyscor_all(1:nz-1)), sum(chisqs_syscor_all(1:nz-1))
!     print *, tmpstr
!      if(output_chisq_eachredbin) then
!        do k = 1,nz-1
!          write(tmpstr1,'(e15.7)') chisqs_nosyscor_all(k)
!          tmpstr = trim(adjustl(tmpstr))//' '//trim(adjustl(tmpstr1))
!        enddo
!        do k = 1,nz-1
!          write(tmpstr1,'(e15.7)') chisqs_syscor_all(k)
!          tmpstr = trim(adjustl(tmpstr))//' '//trim(adjustl(tmpstr1))
!        enddo
!      endif
!      write(fileunit_allsch,'(A)') trim(adjustl(tmpstr))
!    enddo
!    enddo
    
!    do i1 = 1,n1
!    do i2 = 1,n2
!      if(output_sep_schemes) close(fileunits(i1,i2))
!    enddo
!    enddo
!    close(fileunit_allsch)
!    chisqs_nosyscor(N1,N2)
!    chisqs_syscor(A1,N2)
  end subroutine smu_ximu_CalcDAHChisqs


  subroutine smu_ximu_CalcDAHChisqs_bigcov(& 
    DAs, Hs, & ! List of omegam, w
    !outputdir, baseoutputfile, & ! Basic name of the outputfile
    omstds, wstds, numomwstds, & ! "standard" values of omegam, w: baseline cosmology for (s,mu) mapping. Stored in data2pcffile_base
    smutabstds, smutabstds_inited, & ! xi(s,mu) table of baseline cosmologies
    chisqs_nosyscor, chisqs_syscor, & ! values of chisqs, separate schemes
    chisqs_nosyscor_all, chisqs_syscor_all, & ! values of chisqs, averaged over all schemes, correction factor for covmat (D, m1, m2) considered
    weightedstds, avg_counts, VolumeWeightedDAH & ! assign different weights to the "standard" cosmologies
    )
    !arguments
    integer, intent(in) ::  numomwstds
    real(rt), intent(in) :: DAs(nz), Hs(nz), & ! values of DA, H at the effective redshifts
       omstds(numomwstds), wstds(numomwstds)
!    character(*), intent(in) :: outputdir, baseoutputfile
    real(rt), intent(inout) :: smutabstds(nbins_database,mubins_database,3,nz,numomwstds)
    logical, intent(inout) :: smutabstds_inited
    real(rt), intent(out) :: chisqs_nosyscor(n1,n2), chisqs_syscor(n1,n2), &
       chisqs_nosyscor_all, chisqs_syscor_all
    logical, intent(in) :: weightedstds, avg_counts
    logical, intent(in), optional :: VolumeWeightedDAH
    
    ! variables
    real(rt) :: smutab_data(nbins_data,mubins_data,3), smutab_data1(nbins_data,mubins_data,3)
    real(rt) :: smin_mapping=min(1.0_rt,ints1*0.8_rt), smax_mapping=ints2*1.2_rt, &
      DAstd,DAnew,Hstd,Hnew, deltas1,deltas2, sumDAHweights,DAHweights(numomwstds,nz), &
      DAvalues(numgal_nbin),Hvalues(numgal_nbin),DAarray(numgal_nbin),Harray(numgal_nbin), &
      DAstds(nz,numomwstds),Hstds(nz,numomwstds),&
      intxis(maxval(mubins),N1,N2,nz,numomwstds), intxi(maxval(mubins)), dintxi(maxval(mubins)), & 
      dintxi_big(maxval(mubins)*(nz-1)),factD,factM1,factM2
    type(omwpar) :: parstd
!    character(len=1000) :: tmpstr,tmpstr1,tmpstr2,tmpstr3,tmpstr4,tmpstr5,tmpstr6,tmpstr7,tmpstr8,tmpstr9,&
!      tmpstr10,tmpstr11,tmpstr12,&
!      nowchisqstr,filenames(N1,N2), nowfile, filename_allsch
    integer ::  i,j,k,iz,iomwstds,iomwstds2,i1,i2,i_numgal, n
!    integer, parameter :: basefileunit = 58483
    logical :: numgal_weighted=.false., debug_calcchisq = .false.
    ! forecast related
    real(rt) :: DAfcs(nzFC), Hfcs(nzFC)! forecast related 

!    if(ForeCast) then
!      if(numomwstds .ne. 1) then
!        print *, '(smu_ximu_CalcDAHChisqs) ERROR! For forecast standard-omw should be 1: numomwstds = ', numomwstds
!        stop
!      endif
!      if(weightedstds) then
!        print *, '(smu_ximu_CalcDAHChisqs) ERROR! For forecast weightedstds not supported: weightedstds = ', weightedstds
!        stop
!      endif
!      if(avg_counts) then
!        print *, '(smu_ximu_CalcDAHChisqs) ERROR! For forecast avg_counts not supported: avg_counts = ', avg_counts
!        stop
!      endif
!      ! init for forecast
!      do iz = 1, nzFC
!        DAfcs(iz) = DAz_wcdm(fcpar,zeffsFC(iz))
!        Hfcs(iz)  = Hz_wcdm(fcpar,zeffsFC(iz))
!      enddo
!    endif
        


    !! test variables 
    
    if(present(VolumeWeightedDAH)) then
      numgal_weighted=VolumeWeightedDAH
    else
      numgal_weighted=.false.
    endif
    

    do iomwstds = 1, numomwstds
       parstd%omegam = omstds(iomwstds); parstd%w = wstds(iomwstds)     
       if(numgal_weighted) then
         if(.not. numgal_inited) call numgal_init()
         do i = 1, numgal_nbin
           DAvalues(i) = DAz_wcdm(parstd,numgal_zcenters(i))
           Hvalues(i)  = Hz_wcdm (parstd,numgal_zcenters(i))
         enddo
        endif
        do iz = 1, nz
            if(.not.numgal_weighted) then
                    DAstds(iz,iomwstds) = DAz_wcdm(parstd,zeffs(iz)); Hstds(iz,iomwstds) = Hz_wcdm(parstd,zeffs(iz))
            else
                      DAarray = DAvalues(1:numgal_nbin)*redbin_weights(iz,1:numgal_nbin)
                      Harray  = Hvalues(1:numgal_nbin) *redbin_weights(iz,1:numgal_nbin)
                      DAstds(iz,iomwstds) = sum(DAarray(1:numgal_nbin)) / sum(redbin_weights(iz,1:numgal_nbin))
                      Hstds(iz,iomwstds)  = sum(Harray(1:numgal_nbin))  / sum(redbin_weights(iz,1:numgal_nbin))
                      !print *, parstd%omegam, parstd%w
                      !print *, 'DAs:', DAstds(iz,iomwstds), DAz_wcdm(parstd,zeffs(iz))
                      !print *, 'Hs:',  Hstds(iz,iomwstds),  Hz_wcdm(parstd,zeffs(iz))
            endif
       enddo
    enddo

    ! 1. Initialization
    ! 1.1 Check covmat files ready or not
    do i=1,N1
    do j=1,N2
      if (.not.allocated(bigcovmats(i,j)%A).or.bigcovmats(i,j)%nA.ne.(mubins(i)-1)*(nz-1)) then
        print *, ' (smu_ximu_CalcDAHChisqs_bigcov) ERROR! Covmats not ready: i,j,mubin,mucut = ', i,j,mubins(i),mucuts(j)
        stop
      endif
    enddo
    enddo
    
    ! 1.2 Load in 2pCF computed in standard cosmology
    if (.not.smutabstds_inited) then
      print *, ' (smu_ximu_CalcDAHChisqs_bigcov) Load in the xi(s,mu) of baseline cosmologies...'
      do iomwstds = 1, numomwstds
      do iz = 1, nz
        call ximu_loadsmufile(data2pcffile_base(iz,omstds(iomwstds),wstds(iomwstds)), smutabstds(:,:,:,iz,iomwstds), &
                smax_database, nbins_database, mubins_database)
      enddo
      enddo
      smutabstds_inited =.true.
    endif

    if(weightedstds) then
        do iomwstds = 1, numomwstds
          do iz  = 1, nz
            parstd%omegam = omstds(iomwstds); parstd%w = wstds(iomwstds)
            DAstd=DAstds(iz,iomwstds); Hstd=Hstds(iz,iomwstds)
            DAnew = DAs(iz); Hnew = Hs(iz)
            DAHweights(iomwstds,iz) = 1.0 / ((DAnew/DAstd - 1.0)**2.0 + (Hnew/Hstd - 1.0)**2.0 + 0.0001)
          enddo
        enddo
        do iz = 1, nz
          sumDAHweights = sum(DAHweights(1:numomwstds,iz))
          do iomwstds = 1, numomwstds
            DAHweights(iomwstds,iz) = DAHweights(iomwstds,iz) * numomwstds / sumDAHweights
          enddo
        enddo
        print *, 'ERROR (smu_ximu_CalcDAHChisqs_bigcov):  DAHweights not supported in bigcovmat!'; stop
    else
        DAHweights = 1.0_rt
    endif

    if(.true.) then
      do iomwstds = 1, numomwstds
        if(avg_counts.and.iomwstds.gt.1) cycle
        do iz  = 1, nz
          ! 3.1 smutabdata obtained from DSMapping
          deltas1 = smax_database / float(nbins_database) 
          deltas2 = smax_data / float(nbins_data)
!         print *, ' (smu_ximu_CalcOmWChisqs) Do DSMapping... iz = ', iz
          if(.not.avg_counts) then
            parstd%omegam = omstds(iomwstds); parstd%w = wstds(iomwstds)
            !DAstd = DAz_wcdm(parstd,zeffs(iz)); Hstd = Hz_wcdm(parstd,zeffs(iz))
            DAstd=DAstds(iz,iomwstds); Hstd=Hstds(iz,iomwstds)
            DAnew = DAs(iz); Hnew = Hs(iz)
            call DSMapping(smutabstds(:,:,:,iz,iomwstds), nbins_database, mubins_database, smutab_data, &
              nbins_data, mubins_data, DAstd, DAnew, Hstd, Hnew, deltas1,  deltas2,  smin_mapping, smax_mapping)
          else
            ! in case of avg_counts, sum up all smutabstds to build one
            DAnew = DAs(iz); Hnew = Hs(iz)
            smutab_data = 0.0_rt
            do iomwstds2 = 1, numomwstds
              parstd%omegam = omstds(iomwstds2); parstd%w = wstds(iomwstds2)
              !DAstd = DAz_wcdm(parstd,zeffs(iz)); Hstd = Hz_wcdm(parstd,zeffs(iz))
              DAstd=DAstds(iz,iomwstds); Hstd=Hstds(iz,iomwstds)
              call DSMapping(smutabstds(:,:,:,iz,iomwstds2), nbins_database, mubins_database, smutab_data1, &
                nbins_data, mubins_data, DAstd, DAnew, Hstd, Hnew, deltas1,  deltas2,  smin_mapping, smax_mapping)
              smutab_data(:,:,:) = smutab_data(:,:,:) + &
                smutab_data1(:,:,:) * DAHweights(iomwstds2,iz) / dble(numomwstds)
            enddo
          endif
!         print *, ' (smu_ximu_CalcDAHChisqs) DSMapping done: iz = ', iz
          if(debug_calcchisq ) then !.and. mubins(i).eq.25 .and.j.eq.1
                print *, '####################################'
                print *, 'iz = ',iz
                print *, "Checking: DAstd, DAnew, Hstd, Hnew = ", DAstd, DAnew, Hstd, Hnew
                do i = 1, 15
                  j = 1
                  print *, i, smutab_data(i,j,1:3), (smutab_data(i,j,1)-2*smutab_data(i,j,2))/smutab_data(i,j,3)+1.0_rt
                enddo        
          endif
        
          do i1 = 1, N1
          do i2 = 1, N2
!           print *, ' (smu_ximu_CalcDAHChisqs) Compute \int xi...'
            call XiFun(smutab_data, deltas2, nbins_data, mubins_data, &
              anglemin=1.0_rt-mucuts(i2), anglemax=1.0_rt, &
              smin=ints1, smax=ints2, &
              nummuedge=mubins(i1)+1, &
              intxi=intxi(1:mubins(i1)))
!           print *, ' (smu_ximu_CalcDAHChisqs) Normalize \int xi...'
            call normfun(intxi(1:mubins(i1)),mubins(i1),intxis(1:mubins(i1)-1,i1,i2,iz,iomwstds)) 
            if(debug_calcchisq .and. i1.eq.1 .and.i2.eq.1) then
                   print *, 'intxi = ', real(intxi(1:mubins(i1)))
                print *
                print *, 'intxi (normed) = ', real(intxis(1:mubins(i1)-1,i1,i2,iz,iomwstds))
            endif
          enddo!        do i1 = 1, N1
          enddo!        do i2 = 1, N2
!         print *, ' (smu_ximu_CalcDAHChisqs) Compute intxi done: iz = ', iz
        enddo!      do iz  = 1, nz
      enddo!do iomwstds = 1, numomwstds

      chisqs_nosyscor = 0.0_rt; chisqs_syscor = 0.0_rt
      chisqs_nosyscor_all=0.0_rt; chisqs_syscor_all = 0.0_rt;

      if(exclude_bin > 0) print *, 'Warning! Exclude bin: ', exclude_bin
    !   3.2 Compute chisqs
    !  TBD To Be modified.... chisq computation
      do i1 = 1,N1
        do i2 = 1,N2
            dintxi_big = 0.0
            if(.not.NBComp) then
              call Percival_cs(ncovmocks,mubins(i1)*(nz-1),2,factD,factm1,factm2)
            else
              call Percival_cs(ncovmocks,mubins(i1),2,factD,factm1,factm2)
            endif
            n = mubins(i1)
            ! Loop for omstd/wstd
            do iomwstds = 1, numomwstds
              if(avg_counts.and.iomwstds.gt.1) cycle
              ! chisq before sys cor
              do iz = 2, nz
                if(.not.NBComp) then
                  dintxi(1:n-1) = intxis(1:n-1,i1,i2,iz,iomwstds) - intxis(1:n-1,i1,i2,1,iomwstds)
                else 
                  dintxi(1:n-1) = intxis(1:n-1,i1,i2,iz,iomwstds) - intxis(1:n-1,i1,i2,iz-1,iomwstds)
                endif
                if(iz.eq.exclude_bin) dintxi = 0.0 ! xiaodong: exclude last bin!!!
                dintxi_big((iz-2)*(n-1)+1:(iz-1)*(n-1)) = dintxi(1:n-1)
              enddo
              if(.not.avg_counts) then
                chisqs_nosyscor(i1,i2) = chisqs_nosyscor(i1,i2) + &
                    chisq_cov_xbar(dintxi_big(1:(nz-1)*(n-1)), bigcovmats(i1,i2)%A, (nz-1)*(n-1))/dble(numomwstds)/factm2
              else
                chisqs_nosyscor(i1,i2) = &
                    chisq_cov_xbar(dintxi_big(1:(nz-1)*(n-1)), bigcovmats(i1,i2)%A, (nz-1)*(n-1))/factm2
              endif
              chisqs_nosyscor_all = chisqs_nosyscor_all + chisqs_nosyscor(i1,i2)/ dble(N1*N2)
              !print *
              !print *, 'dintxi = ', dintxi_big(1:(nz-1)*(n-1))
              !open(unit=100,file='dintxi.txt'); write(100,*) dintxi_big(1:(nz-1)*(n-1)); close(100); 
              !print *, '1,nz,n=',1,nz,n
              !print *, '(nz-1)*(n-1) = ', (nz-1)*(n-1)
              !print *, 'bigcovmats(i1,i2)%A = ', bigcovmats(i1,i2)%A
              !print *, 'chisqs_nosyscor(i1,i2) = ', chisq_cov_xbar(dintxi_big(1:(nz-1)*(n-1)), bigcovmats(i1,i2)%A, (nz-1)*(n-1))
              !print *, 'n1, n2, factm2 = ', n1,n2,factm2
              !print *, 'iomwstds,iz = ', iomwstds,iz
              !print *, 'DAHweights(iomwstds,iz) = ', DAHweights(iomwstds,iz)
              !stop
              ! chisq after sys cor
              do iz = 2, nz
                if(.not.NBComp) then
                  dintxi(1:n-1) = intxis(1:n-1,i1,i2,iz,iomwstds) - intxis(1:n-1,i1,i2,1,iomwstds)
                else
                  dintxi(1:n-1) = intxis(1:n-1,i1,i2,iz,iomwstds) - intxis(1:n-1,i1,i2,iz-1,iomwstds)
                endif
                !if(i1.eq.1.and.i2.eq.1) write(*,*) '(CalcDAHChisqs, dintxi):',iz, dintxi(1:n-1) !selfsyscor_debug
                dintxi(1:n-1) = dintxi(1:n-1) - dintxi_syscor(1:n-1,i1,i2,iz-1)
                if(iz.eq.exclude_bin) dintxi = 0.0 ! xiaodong: exclude last bin!!!
                dintxi_big((iz-2)*(n-1)+1:(iz-1)*(n-1)) = dintxi(1:n-1)
              enddo
              if(.not.avg_counts) then
                chisqs_syscor(i1,i2) = chisqs_syscor(i1,i2) + &
                      chisq_cov_xbar(dintxi_big(1:(nz-1)*(n-1)), bigcovmats(i1,i2)%A, (nz-1)*(n-1))/dble(numomwstds)/factm2!*DAHweights(iomwstds,iz)
              else
                chisqs_syscor(i1,i2) =  &
                      chisq_cov_xbar(dintxi_big(1:(nz-1)*(n-1)), bigcovmats(i1,i2)%A, (nz-1)*(n-1))/factm2
              endif
              chisqs_syscor_all   = chisqs_syscor_all + chisqs_syscor(i1,i2) / dble(N1*N2)

              if(debug_calcchisq .and. mubins(i1).eq.25 .and.i2.eq.1) then
                do iz = 2, nz
                      print *, 'iz = ', iz
                    print *, 'dintxi_sys = ', dintxi_syscor(1:n-1,i1,i2,iz-1)
                    print *, 'dintxi (after syscor) = ', dintxi(1:n-1)
                  !print *, 'invcov = ', real(covmats(i1,i2,iz-1)%A)
                enddo
                  print *, 'chisq / chisq_syscor = ', chisqs_nosyscor(i1,i2), chisqs_syscor(i1,i2)        
              endif
            enddo!do iomwstds = 1, numomwstds
 !         enddo!do iz = 2, nz

 !         ! write chisq values to files...
 !         if(output_sep_schemes) then
 !           tmpstr=''; tmpstr1=''
 !           write(tmpstr,'(f5.2,1x,A,e15.7,e15.7,5x)') 1.0_rt-mucuts(i2), trim(adjustl(omwstr(omlist(iom),wlist(iw)))), &
 !             sum(chisqs_nosyscor(i1,i2,1:nz-1)), sum(chisqs_syscor(i1,i2,1:nz-1))
 !           if(output_chisq_eachredbin) then
 !             do k = 1,nz-1
 !               write(tmpstr1,'(e15.7)') chisqs_nosyscor(i1,i2,k)
 !               tmpstr = trim(adjustl(tmpstr))//' '//trim(adjustl(tmpstr1))
 !             enddo
 !             do k = 1,nz-1
 !               write(tmpstr1,'(e15.7)') chisqs_syscor(i1,i2,k)
 !               tmpstr = trim(adjustl(tmpstr))//' '//trim(adjustl(tmpstr1))
 !             enddo
 !           endif
 !           write(fileunits(i1,i2),'(A)') trim(adjustl(tmpstr))
 !         endif
        enddo!do i1 = 1,N1
      enddo!do i2 = 1,N2
      ! output allscheme chisq value
!      tmpstr=''; tmpstr1=''
!      write(tmpstr,'(f5.2,1x,A,e15.7,e15.7,5x)') 1.0_rt-mucuts(i2), trim(adjustl(omwstr(omlist(iom),wlist(iw)))), &
!        sum(chisqs_nosyscor_all(1:nz-1)), sum(chisqs_syscor_all(1:nz-1))
!     print *, tmpstr
!      if(output_chisq_eachredbin) then
!        do k = 1,nz-1
!          write(tmpstr1,'(e15.7)') chisqs_nosyscor_all(k)
!          tmpstr = trim(adjustl(tmpstr))//' '//trim(adjustl(tmpstr1))
!        enddo
!        do k = 1,nz-1
!          write(tmpstr1,'(e15.7)') chisqs_syscor_all(k)
!          tmpstr = trim(adjustl(tmpstr))//' '//trim(adjustl(tmpstr1))
!        enddo
!      endif
!      write(fileunit_allsch,'(A)') trim(adjustl(tmpstr))
!    enddo
!    enddo
    
!    do i1 = 1,n1
!    do i2 = 1,n2
!      if(output_sep_schemes) close(fileunits(i1,i2))
!    enddo
!    enddo
!    close(fileunit_allsch)
!    chisqs_nosyscor(N1,N2)
!    chisqs_syscor(A1,N2)
    endif
  end subroutine smu_ximu_CalcDAHChisqs_bigcov




  subroutine smu_ximu_CalcDAHChisqs_bigcovfc(& 
    DAs, Hs, & ! List of omegam, w
    omstd, wstd,  & ! "standard" values of omegam, w. Stored in data2pcffile_base
    omfc, wfc, & ! "standard" "forecast values" of omegam, w.
    smutabstds, smutabstds_inited, & ! xi(s,mu) table of baseline cosmologies
    chisqs_nosyscor, chisqs_syscor, & ! values of chisqs, separate schemes
    chisqs_nosyscor_all, chisqs_syscor_all & ! values of chisqs, averaged over all schemes, correction factor for covmat (D, m1, m2) considered
    )
    !arguments
    real(rt), intent(in) :: DAs(nzfc), Hs(nzfc), & ! values of DA, H at the effective redshifts
       omstd, wstd, omfc, wfc
!    character(*), intent(in) :: outputdir, baseoutputfile
    real(rt), intent(inout) :: smutabstds(nbins_database,mubins_database,3,nz,1)
    logical, intent(inout) :: smutabstds_inited
    real(rt), intent(out) :: chisqs_nosyscor(n1,n2), chisqs_syscor(n1,n2), &
       chisqs_nosyscor_all, chisqs_syscor_all
    
    ! variables
    real(rt) :: smutab_data(nbins_data,mubins_data,3), smutab_data1(nbins_data,mubins_data,3)
    real(rt) :: smin_mapping=min(1.0_rt,ints1*0.8_rt), smax_mapping=ints2*1.2_rt, &
      DAstd,DAnew,Hstd,Hnew, deltas1,deltas2, sumDAHweights,DAHweights(nzfc), &
      DAvalues(numgal_nbin),Hvalues(numgal_nbin),DAarray(numgal_nbin),Harray(numgal_nbin), &
      intxis(maxval(mubins),N1,N2,nzfc), intxi(maxval(mubins)), dintxi(maxval(mubins)), & 
      dintxi_big(maxval(mubins)*(nzfc-1)),factD,factM1,factM2

    type(omwpar) :: parstd, parfc
!    character(len=1000) :: tmpstr,tmpstr1,tmpstr2,tmpstr3,tmpstr4,tmpstr5,tmpstr6,tmpstr7,tmpstr8,tmpstr9,&
!      tmpstr10,tmpstr11,tmpstr12,&
!      nowchisqstr,filenames(N1,N2), nowfile, filename_allsch
    integer ::  i,j,k,iz,i1,i2,i_numgal, n,nbox,nA, k1,k2
!    integer, parameter :: basefileunit = 58483
    logical :: numgal_weighted=.false., debug_calcchisq = .false.
    ! forecast related
    real(rt) :: DAfcs(nzFC), Hfcs(nzFC), ngal1,ngal2, ngalfc1,ngalfc2, sigsqrat1,sigsqrat2! forecast related 

    if(.not.NBComp) then
      print *, '(smu_ximu_CalcDAHChisqs_bigcovfc) ERROR! Must use NBComp for forecast! NBComp = ', NBComp
      stop
    endif

    parfc%omegam = omfc; parfc%w = wfc; 
    do iz = 1, nzFC
        DAfcs(iz) = DAz_wcdm(parfc,zeffsFC(iz))
        Hfcs(iz)  = Hz_wcdm(parfc,zeffsFC(iz))
    enddo
        

    
    ! 1.2 Load in 2pCF computed in standard cosmology
    if (.not.smutabstds_inited) then
      print *, ' (smu_ximu_CalcDAHChisqs_bigcov) Load in the xi(s,mu) of baseline cosmologies...'
      do iz = 1, nz
        call ximu_loadsmufile(data2pcffile_base(iz,omstd,wstd), smutabstds(:,:,:,iz,1), &
                smax_database, nbins_database, mubins_database)
      enddo
      smutabstds_inited =.true.
    endif

!    print *, '(smu_ximu_CalcDAHChisqs_bigcovfc) End of load in smutabstds.'
    DAHweights = 1.0_rt

    do iz  = 1, nzfc
!          print *, 'iz = ', iz
          ! 3.1 smutabdata obtained from DSMapping
          deltas1 = smax_database / float(nbins_database) 
          deltas2 = smax_data / float(nbins_data)
!         print *, ' (smu_ximu_CalcOmWChisqs) Do DSMapping... iz = ', iz
          DAstd = DAfcs(iz); Hstd = Hfcs(iz)
          DAnew = DAs(iz);   Hnew = Hs(iz)
!          print *, 'DAstd, Hstd, DAnew, Hnew = ', DAstd, Hstd, DAnew, Hnew
!          call DSMapping(smutabstds(:,:,:,iz,1), nbins_database, mubins_database, smutab_data, &
!              nbins_data, mubins_data, DAstd, DAnew, Hstd, Hnew, deltas1,  deltas2,  smin_mapping, smax_mapping)
          ! Very important! always using first bin as reference...
          call DSMapping(smutabstds(:,:,:,1,1), nbins_database, mubins_database, smutab_data, &
              nbins_data, mubins_data, DAstd, DAnew, Hstd, Hnew, deltas1,  deltas2,  smin_mapping, smax_mapping)
!         print *, ' (smu_ximu_CalcDAHChisqs) DSMapping done: iz = ', iz
          if(debug_calcchisq ) then !.and. mubins(i).eq.25 .and.j.eq.1
                print *, '####################################'
                print *, 'iz = ',iz
                print *, "Checking: DAstd, DAnew, Hstd, Hnew = ", DAstd, DAnew, Hstd, Hnew
                do i = 1, 15
                  j = 1
                  print *, i, smutab_data(i,j,1:3), (smutab_data(i,j,1)-2*smutab_data(i,j,2))/smutab_data(i,j,3)+1.0_rt
                enddo        
          endif
        
          do i1 = 1, N1
          do i2 = 1, N2
!           print *, ' (smu_ximu_CalcDAHChisqs) Compute \int xi...'
            call XiFun(smutab_data, deltas2, nbins_data, mubins_data, &
              anglemin=1.0_rt-mucuts(i2), anglemax=1.0_rt, &
              smin=ints1, smax=ints2, &
              nummuedge=mubins(i1)+1, &
              intxi=intxi(1:mubins(i1)))
!           print *, ' (smu_ximu_CalcDAHChisqs) Normalize \int xi...'
            call normfun(intxi(1:mubins(i1)),mubins(i1),intxis(1:mubins(i1)-1,i1,i2,iz)) 
            if(debug_calcchisq .and. i1.eq.1 .and.i2.eq.1) then
                   print *, 'intxi = ', real(intxi(1:mubins(i1)))
                print *
                print *, 'intxi (normed) = ', real(intxis(1:mubins(i1)-1,i1,i2,iz))
            endif
          enddo!        do i1 = 1, N1
          enddo!        do i2 = 1, N2
!         print *, ' (smu_ximu_CalcDAHChisqs) Compute intxi done: iz = ', iz
    enddo!      do iz  = 1, nzfc

    chisqs_nosyscor = 0.0_rt; chisqs_syscor = 0.0_rt
    chisqs_nosyscor_all=0.0_rt; chisqs_syscor_all = 0.0_rt;

    if(exclude_bin > 0) print *, 'Warning! Exclude bin: ', exclude_bin
    !   3.2 Compute chisqs
    !  TBD To Be modified.... chisq computation
    do i1 = 1,N1
        do i2 = 1,N2
            dintxi_big = 0.0
!            if(.not.NBComp) then
!              call Percival_cs(ncovmocks,mubins(i1)*(nz-1),2,factD,factm1,factm2)
!            else
            call Percival_cs(ncovmocks,mubins(i1),2,factD,factm1,factm2)
!            endif
            n = mubins(i1)
            ! Loop for omstd/wstd
              ! chisq before sys cor
            do iz = 2, nzfc
!                if(.not.NBComp) then
!                  dintxi(1:n-1) = intxis(1:n-1,i1,i2,iz) - intxis(1:n-1,i1,i2,1)
!                else 
                dintxi(1:n-1) = intxis(1:n-1,i1,i2,iz) - intxis(1:n-1,i1,i2,iz-1)
!                endif
                if(iz.eq.exclude_bin) dintxi = 0.0 ! xiaodong: exclude last bin!!!
                dintxi_big((iz-2)*(n-1)+1:(iz-1)*(n-1)) = dintxi(1:n-1)
!                if(i1.eq.1.and.i2.eq.1) print *, 'iz, dintxi = ', iz, real(dintxi)
            enddo
            chisqs_nosyscor(i1,i2) = chisqs_nosyscor(i1,i2) + &
                    chisq_cov_xbar(dintxi_big(1:(nzfc-1)*(n-1)), bigcovmatsfc(i1,i2)%A, (nzfc-1)*(n-1))/factm2
            print *, 'i1,i2,iz, chisq = ', i1,i2,iz, chisq_cov_xbar(dintxi_big(1:(nzfc-1)*(n-1)), bigcovmatsfc(i1,i2)%A, (nzfc-1)*(n-1))/factm2
            chisqs_nosyscor_all = chisqs_nosyscor_all + chisqs_nosyscor(i1,i2)/ dble(N1*N2)
            do iz = 2, nzfc
!                if(.not.NBComp) then
!                  dintxi(1:n-1) = intxis(1:n-1,i1,i2,iz) - intxis(1:n-1,i1,i2,1)
!                else
                dintxi(1:n-1) = intxis(1:n-1,i1,i2,iz) - intxis(1:n-1,i1,i2,iz-1)
!                endif
                !if(i1.eq.1.and.i2.eq.1) write(*,*) '(CalcDAHChisqs, dintxi):',iz, dintxi(1:n-1) !selfsyscor_debug
                dintxi(1:n-1) = dintxi(1:n-1) - dintxi_syscor(1:n-1,i1,i2,iz-1)
                if(iz.eq.exclude_bin) dintxi = 0.0 ! xiaodong: exclude last bin!!!
                dintxi_big((iz-2)*(n-1)+1:(iz-1)*(n-1)) = dintxi(1:n-1)
            enddo
            chisqs_syscor(i1,i2) = chisqs_syscor(i1,i2) + &
                      chisq_cov_xbar(dintxi_big(1:(nzfc-1)*(n-1)), bigcovmatsfc(i1,i2)%A, (nzfc-1)*(n-1))/factm2
            chisqs_syscor_all   = chisqs_syscor_all + chisqs_syscor(i1,i2) / dble(N1*N2)

        enddo!do i1 = 1,N1
    enddo!do i2 = 1,N2
      ! output allscheme chisq value
!      tmpstr=''; tmpstr1=''
!      write(tmpstr,'(f5.2,1x,A,e15.7,e15.7,5x)') 1.0_rt-mucuts(i2), trim(adjustl(omwstr(omlist(iom),wlist(iw)))), &
!        sum(chisqs_nosyscor_all(1:nz-1)), sum(chisqs_syscor_all(1:nz-1))
!     print *, tmpstr
!      if(output_chisq_eachredbin) then
!        do k = 1,nz-1
!          write(tmpstr1,'(e15.7)') chisqs_nosyscor_all(k)
!          tmpstr = trim(adjustl(tmpstr))//' '//trim(adjustl(tmpstr1))
!        enddo
!        do k = 1,nz-1
!          write(tmpstr1,'(e15.7)') chisqs_syscor_all(k)
!          tmpstr = trim(adjustl(tmpstr))//' '//trim(adjustl(tmpstr1))
!        enddo
!      endif
!      write(fileunit_allsch,'(A)') trim(adjustl(tmpstr))
!    enddo
!    enddo
    
!    do i1 = 1,n1
!    do i2 = 1,n2
!      if(output_sep_schemes) close(fileunits(i1,i2))
!    enddo
!    enddo
!    close(fileunit_allsch)
!    chisqs_nosyscor(N1,N2)
!    chisqs_syscor(A1,N2)
  end subroutine smu_ximu_CalcDAHChisqs_bigcovfc
    
  real(rt) function chisq_cov_xbar(xbar, invcov, n)
    integer, intent(in) :: n
    real(rt), intent(in) :: xbar(n), invcov(n,n)
    integer :: i,j
    chisq_cov_xbar = 0.0_rt
    do i = 1, n
    do j = 1, n
      chisq_cov_xbar = chisq_cov_xbar+xbar(i)*invcov(i,j)*xbar(j)
    enddo
    enddo
  end function chisq_cov_xbar
  
  character(charlen) function AP_MCMCstr(numomwstds, omstds, wstds)
    integer, intent(in) :: numomwstds
    real(rt), intent(in) :: omstds(numomwstds), wstds(numomwstds)
    character(charlen) :: tmpstr1,tmpstr2,tmpstr3,tmpstr4,tmpstr5,tmpstr6,tmpstr7,tmpstr8,&
      tmpstr9,tmpstr10,tmpstr11,tmpstr12,tmpstr13,tmpstr14,tmpstr15
    integer :: iomwstds
    
    write(tmpstr1, *) N1 ! how many bins in mu space
    write(tmpstr2, *) mubins(1) ! how many bins in mu space
    write(tmpstr3, *) mubins(N1) ! how many bins in mu space
    write(tmpstr4, *) N2
    write(tmpstr5, '(f4.2)') mucuts(1)
    write(tmpstr6, '(f4.2)') mucuts(N2)
    write(tmpstr7, '(f5.1)') ints1
    write(tmpstr8, '(f5.1)') ints2
    write(tmpstr9, *) ncovmocks
    tmpstr10=''; if(polyfitdeg.ge.1) write(tmpstr10,'(A,i1)') '.polyfitdeg',polyfitdeg
    write(tmpstr11,*) nbins_database 
    write(tmpstr12,*) mubins_database; 
    write(tmpstr13,*) numomwstds; 
    tmpstr14 = '';
    do iomwstds = 1, numomwstds
        tmpstr14 = trim(adjustl(tmpstr14))//'_'//trim(adjustl(omwstr(omstds(iomwstds),wstds(iomwstds))))
    enddo
    tmpstr15='.baseline_'//trim(adjustl(tmpstr13))//'omws_'//trim(adjustl(tmpstr14))// &
              '__sbin'//trim(adjustl(tmpstr11))//'_mubin'//trim(adjustl(tmpstr12))
    AP_MCMCstr = trim(adjustl(tmpstr1)) &
        //'mu'//trim(adjustl(tmpstr2))//'to'//trim(adjustl(tmpstr3))//'_'//trim(adjustl(tmpstr4))&
        //'mumax'// trim(adjustl(tmpstr5))//'to'//trim(adjustl(tmpstr6))  &
        //'_Ints'//trim(adjustl(tmpstr7))//'to'//trim(adjustl(tmpstr8)) &
        //'_'//trim(adjustl(tmpstr9))//'covmocks'//trim(adjustl(tmpstr10)) &
        //'_base'//trim(adjustl(tmpstr13))//'omws'//trim(adjustl(tmpstr14)) &
        //'_'//trim(adjustl(tmpstr11))//'sbin'//'_'//trim(adjustl(tmpstr12))//'mubin'  
  end function AP_MCMCstr
    
end module AP_2pcf_tools


module AP_funs
use AP_2pcf_tools
implicit none

logical :: AP_inited = .false.
contains

 subroutine AP_Like(DAs, Hs, chisqs,  printinfo, compute_covmat, use_bigcovmat, covmat_suffixstr, chisqs__uncored, sepchisqs_uncored_result, sepchisqs_result)
    implicit none
    real(rt), intent(in) :: DAs(nz), Hs(nz)
    logical, intent(in) :: printinfo
    real(rt), intent(out):: chisqs(nz-1)
    real(rt), intent(out), optional :: chisqs__uncored(nz-1)
    real(rt), intent(out), optional :: sepchisqs_uncored_result(n1,n2), sepchisqs_result(n1,n2)
    logical, intent(in), optional :: compute_covmat, use_bigcovmat
    character(*), intent(in), optional :: covmat_suffixstr
    real(rt) :: sepchisqs_uncored(n1,n2,nz-1), sepchisqs(n1,n2,nz-1), chisqs_uncored(nz-1), &
                sepchisqs_uncored_bigcov(n1,n2), sepchisqs_bigcov(n1,n2), fact
    integer :: num_fidcos
    real(rt) :: om_fiducial(100), w_fiducial(100)
    logical :: usebigcovmat, computecovmat
    character(len=1000) :: covmatsuffixstr
    integer :: i1, i2, iz, i


    ! Set mode: whether use bigcovmat
    usebigcovmat = .false.
    if(present(use_bigcovmat)) then
        if(use_bigcovmat) usebigcovmat=.true.
    endif

    computecovmat = .false.
    if(present(compute_covmat)) then
        if(compute_covmat) computecovmat=.true.
    endif
   
    if(present(covmat_suffixstr)) then
      covmatsuffixstr = covmat_suffixstr
    else
      covmatsuffixstr = ''
    endif


    !--------------------------------
    ! Preparation for the compute of AP likelihood
  
    ! Using 2-point CF measured in one fiducial model to 
    !   infer the 2-point CF in non-fiducial models
    ! parameters of the fiducial model

    num_fidcos = 1 
    om_fiducial(1)  = gb_om_fiducial;  w_fiducial(1)  = gb_w_fiducial
    !om_fiducial(1)  = 0.26_rt;  w_fiducial(1)  = -1.50_rt 
    if(.not.AP_inited) then
      print *, '(Begin) Load in necessary files.'
      if(usebigcovmat) then
        if(computecovmat) then
          print *, 'Compute/output big covmats...'; call calc_bigcovmats(); call output_bigcovmats(covmatsuffixstr)
        endif
        print *, '* Load in big covmats:'; call load_bigcovmats(covmatsuffixstr)

        print *, '* Invert big covmats:'; call invert_bigcovmats()
        !open(unit=1087,file='invbigcov.txt')
        !do i1 = 1, bigcovmats(1,1)%nA
        !do i2 = 1, bigcovmats(1,1)%nA
        ! write(1087,*) bigcovmats(1,1)%A(i1,i2)
        !enddo
        !enddo
        !close(1087)
        !stop
      else
        if(computecovmat) then
          print *, 'Compute/output covmats...'; call calc_covmats(); call output_covmats(covmatsuffixstr)
        endif
        print *, '* Load in covmats:'; call load_covmats(covmatsuffixstr)
        print *, '* Invert covmats:'; call invert_covmats()
      endif
!     call system('sleep 0'); print *, 'Compute/output covmats...';call calc_covmats();call output_covmats()
      print *, '* Compute systematic correction:'
      if(.not.allocated(global_smutabstds)) allocate(global_smutabstds(nbins_database,mubins_database,3,nz,num_fidcos))
      if(dataself_IO_test.or.mock_IO_test) then
       do i = 1, 10
         write(*,'(A)'), 'WARNING (AP_like)!!! Found dataself_IO_test OR mock_IO_test set as .true.; pre-load global_smutabstds for self-syscor!!!'
       enddo
       do iz = 1, nz
        call ximu_loadsmufile(data2pcffile_base(iz,om_fiducial(1),w_fiducial(1)), global_smutabstds(:,:,:,iz,1), &
         smax_database, nbins_database, mubins_database)
       enddo
      endif
      call calc_syscor()
    endif
    ! End 
    !--------------------------------

    chisqs_uncored = 0.0; chisqs = 0.0;
    sepchisqs_uncored_bigcov = 0.0; sepchisqs_bigcov = 0.0;
    sepchisqs_uncored = 0.0; sepchisqs = 0.0;
    if(usebigcovmat) then
      call smu_ximu_CalcDAHChisqs_bigcov(& 
        DAs, Hs, & ! Values of DA, H in six cosmologies
        om_fiducial(1:num_fidcos), w_fiducial(1:num_fidcos), num_fidcos, & ! "standard" values of omegam, w: baseline cosmology for (s,mu) mapping. 
        global_smutabstds, AP_inited, & ! xi(s,mu) table of baseline cosmologies
        sepchisqs_uncored_bigcov, sepchisqs_bigcov, & ! values of chisqs, separate schemes
        chisqs_uncored(1), chisqs(1), & ! values of chisqs, averaged over all schemes, correction factor for covmat (D, m1, m2) considered
        weightedstds = .false., avg_counts = .false. &
        ) 
        fact = (1.0/1.0) ! Correction factor due to finite size of Horizon Run 4 simulation (used as correction of systematics)
       !print *, sepchisqs_uncored_bigcov
       !print *, 'sepchisqs_uncored_bigcov'
       !stop
    else
      call smu_ximu_CalcDAHChisqs(& 
        DAs, Hs, & ! Values of DA, H in six cosmologies
        om_fiducial(1:num_fidcos), w_fiducial(1:num_fidcos), num_fidcos, & ! "standard" values of omegam, w: baseline cosmology for (s,mu) mapping. 
        global_smutabstds, AP_inited, & ! xi(s,mu) table of baseline cosmologies
        sepchisqs_uncored, sepchisqs, & ! values of chisqs, separate schemes
        chisqs_uncored, chisqs, & ! values of chisqs, averaged over all schemes, correction factor for covmat (D, m1, m2) considered
        weightedstds = .false., avg_counts = .false. &
        ) 
      fact = (1.0/1.0) ! Correction factor, due to finite size of HR4 & correlation between delta xis
    endif

    sepchisqs_uncored_bigcov = sepchisqs_uncored_bigcov * fact
    sepchisqs_uncored = sepchisqs_uncored * fact

    sepchisqs_bigcov = sepchisqs_bigcov * fact
    sepchisqs = sepchisqs * fact

    chisqs_uncored = chisqs_uncored * fact
    chisqs = chisqs * fact

    if(present(chisqs__uncored)) chisqs__uncored = chisqs_uncored
!    chisqs = chisqs_uncored
!
!    print *, 'Sys-corred:   ', real(sepchisqs_uncored_bigcov)
!    print *, 'Sys-uncorred: ', real(sepchisqs_bigcov)
!    stop
    if(printinfo) then
      write(*,'(A,A)') '  * use ', trim(adjustl(AP_MCMCstr(num_fidcos,om_fiducial(1:num_fidcos), w_fiducial(1:num_fidcos))))
      if(usebigcovmat) then
        write(*,'(A,f10.3,A,f10.3)') '(AP_like) bigcov chisq = ', real(chisqs(1)), ';  chisq_uncored = ', real(chisqs_uncored(1))
      else
        write(*,'(A,<nz>(f10.3),A,<nz>(f10.3))') '(AP_like) 1stref chisq =', real(sum(chisqs(1:nz-1))), real(chisqs(1:nz-1)), ';  chisq_uncored =', real(sum(chisqs_uncored(1:nz-1))),  real(chisqs_uncored(1:nz-1))
      endif
    endif
    if(present(sepchisqs_uncored_result)) then
      if(usebigcovmat) then 
        sepchisqs_uncored_result = sepchisqs_uncored_bigcov
      else
        do i1 = 1, N1; do i2 = 1, N2
         sepchisqs_uncored_result(i1,i2) = sum(sepchisqs_uncored(i1,i2,1:nz-1))
        enddo; enddo
      endif
    endif

    if(present(sepchisqs_result)) then
      if(usebigcovmat) then 
        sepchisqs_result = sepchisqs_bigcov
      else
        do i1 = 1, N1; do i2 = 1, N2
         sepchisqs_result(i1,i2) = sum(sepchisqs(i1,i2,1:nz-1))
        enddo; enddo
      endif
    endif

  end subroutine AP_like


 subroutine init_bigcovmatsfc()    ! 1. Initialization
    integer :: i,j,nbox, nA,iz,k1,k2
    real(rt) :: sigsqrat2, sigsqrat1
    ! 1.1 Check covmat files ready or not
    do i=1,N1
    do j=1,N2
      if (.not.allocated(bigcovmats(i,j)%A).or.bigcovmats(i,j)%nA.ne.(mubins(i)-1)*(nz-1)) then
        print *, ' (smu_ximu_CalcDAHChisqs_bigcov) ERROR! Covmats not ready: i,j,mubin,mucut = ', i,j,mubins(i),mucuts(j)
        stop
      endif

      nbox = mubins(i)-1; nA = nbox*(nzfc-1)
      bigcovmatsfc(i,j)%nA = nA
      if(allocated(bigcovmatsfc(i,j)%A)) deallocate(bigcovmatsfc(i,j)%A)
      allocate(bigcovmatsfc(i,j)%A(nA,nA)); bigcovmatsfc(i,j)%A=0.0_rt
!      print *, 'i,j = ', i,j
      ! for bigcovmats(i,j)%A: 
      !   1, nbox * 1,nbox: sigsq ~ 1/ngal1 + 1/ngal2
      !   nbox+1, 2*nbox * 1,nbox: sigsq ~/ngal2

      ! use the 2*nbox * 2*nbox unit as reference 
      ! infer covmat for forecast
      !compute bigcovmats fc from bigcovmats... !!! very important !! write code, execute, and check result!!!!!!
      do iz = 1, nzfc-1
!        print *, 'iz = ', iz
      !   (1, nbox) * (1,nbox): sigsq ~ 1/ngal1 + 1/ngal2
        sigsqrat1 = 1./dble(ngals(1)) + 1./dble(ngals(2))
        sigsqrat2 = 1./dble(ngalsfc(iz)) + 1./dble(ngalsfc(iz+1))
        do k1 = 1, nbox
        do k2 = 1, nbox
          bigcovmatsfc(i,j)%A(nbox*(iz-1)+k1,nbox*(iz-1)+k2) = &
             sigsqrat2 / sigsqrat1 * bigcovmats(i,j)%A(k1,k2)
        enddo
        enddo
      !   (nbox+1, 2*nbox) * (1,nbox): sigsq ~/ngal2
        if(iz .lt. nzfc-1.and..true.) then
          sigsqrat1 = 1./dble(ngals(2))
          sigsqrat2 = 1./dble(ngalsfc(iz+1)) 
          do k1 = nbox+1, 2*nbox
          do k2 = 1, nbox
            bigcovmatsfc(i,j)%A(nbox*(iz-1)+k1,nbox*(iz-1)+k2) = &
               sigsqrat2 / sigsqrat1 * bigcovmats(i,j)%A(k1,k2)
          enddo
          enddo
          do k1 = 1, nbox
          do k2 = nbox+1, 2*nbox
            bigcovmatsfc(i,j)%A(nbox*(iz-1)+k1,nbox*(iz-1)+k2) = &
               sigsqrat2 / sigsqrat1 * bigcovmats(i,j)%A(k1,k2)
          enddo
          enddo
        endif
!        print *, 'iz / nzfc = ', iz, nzfc
      enddo
      ! test region
      if(.false.) then
!        print *, 'i,j=', i,j
        open(unit=7273399,file='bigcov1.txt')
        open(unit=7273400,file='bigcov2.txt')
        nA = bigcovmats(i,j)%nA
        do k1 = 1, nA
          write(7273399,'(<nA>(e15.7,1x))') bigcovmats(i,j)%A(k1,1:nA)
        enddo
        nA = bigcovmatsfc(i,j)%nA
        do k1 = 1, nA
          write(7273400,'(<nA>(e15.7,1x))') bigcovmatsfc(i,j)%A(k1,1:nA)
        enddo
        close(7273399)
        close(7273400)
      endif

    enddo
    enddo
 end subroutine init_bigcovmatsfc
 subroutine AP_likefc(DAs, Hs, omfc, wfc, chisqs,  printinfo, compute_covmat,  covmat_suffixstr, chisqs__uncored, sepchisqs_uncored_result, sepchisqs_result)
!            AP_like(DAs, Hs, chisqs,               printinfo, compute_covmat, use_bigcovmat, covmat_suffixstr, chisqs__uncored, sepchisqs_uncored_result, sepchisqs_result)
    implicit none
    real(rt), intent(in) :: DAs(nzfc), Hs(nzfc), omfc, wfc
    logical, intent(in) :: printinfo
    real(rt), intent(out):: chisqs(nzfc-1)
    real(rt), intent(out), optional :: chisqs__uncored(nzfc-1)
    real(rt), intent(out), optional :: sepchisqs_uncored_result(n1,n2), sepchisqs_result(n1,n2)
    logical, intent(in), optional :: compute_covmat
    character(*), intent(in), optional :: covmat_suffixstr
    real(rt) :: sepchisqs_uncored(n1,n2,nzfc-1), sepchisqs(n1,n2,nzfc-1), chisqs_uncored(nzfc-1), &
                sepchisqs_uncored_bigcov(n1,n2), sepchisqs_bigcov(n1,n2), fact
    integer :: num_fidcos
    real(rt) :: om_fiducial(100), w_fiducial(100)
    logical :: computecovmat
    character(len=1000) :: covmatsuffixstr
    integer :: i1, i2, iz, i

    computecovmat = .false.
    if(present(compute_covmat)) then
        if(compute_covmat) computecovmat=.true.
    endif
   
    if(present(covmat_suffixstr)) then
      covmatsuffixstr = covmat_suffixstr
    else
      covmatsuffixstr = ''
    endif

    !--------------------------------
    ! Preparation for the compute of AP likelihood
  
    ! Using 2-point CF measured in one fiducial model to 
    !   infer the 2-point CF in non-fiducial models
    ! parameters of the fiducial model

    num_fidcos = 1 
    om_fiducial(1)  = 0.26_rt;  w_fiducial(1)  = -1.00_rt 
    if(.not.AP_inited) then
      print *, '(Begin) Load in necessary files.'
      if(computecovmat) then
        print *, 'Compute/output big covmats...'; call calc_bigcovmats(); call output_bigcovmats(covmatsuffixstr)
      endif
      print *, '* Load in big covmats:'; call load_bigcovmats(covmatsuffixstr)
      print *, '* Init bigcovmatsfc:'; call init_bigcovmatsfc();
      print *, '* Invert big covmats:'; call invert_bigcovmats();
      if(gb_fc_spmat) then
        print *, '* Semi-p big covmats (fc):'; call spmat_bigcovmatsfc();
      endif
      print *, '* Invert big covmats (fc):'; call invert_bigcovmatsfc();
      !open(unit=1087,file='invbigcov.txt')
      !do i1 = 1, bigcovmats(1,1)%nA
      !do i2 = 1, bigcovmats(1,1)%nA
      ! write(1087,*) bigcovmats(1,1)%A(i1,i2)
      !enddo
      !enddo
      !close(1087)
      !stop
!     call system('sleep 0'); print *, 'Compute/output covmats...';call calc_covmats();call output_covmats()
      print *, '* Compute systematic correction:'
      if(.not.allocated(global_smutabstds)) allocate(global_smutabstds(nbins_database,mubins_database,3,nz,num_fidcos))
      if(dataself_IO_test) then
         write(*,*), '(AP_likefc) ERROR: dataself_IO_test not supported for forecast: ', dataself_IO_test
         stop
      endif
      call calc_syscorfc()
    endif
    ! End 
    !--------------------------------

    chisqs_uncored = 0.0; chisqs = 0.0;
    sepchisqs_uncored_bigcov = 0.0; sepchisqs_bigcov = 0.0;
    sepchisqs_uncored = 0.0; sepchisqs = 0.0;
!    print *, 'Beg of calling smu_ximu_CalcDAHChisqs_bigcovfc.'
    call smu_ximu_CalcDAHChisqs_bigcovfc(& 
      DAs, Hs, & ! List of omegam, w
      om_fiducial(1), w_fiducial(1),  & ! "standard" values of omegam, w. Stored in data2pcffile_base
      omfc, wfc, & ! "standard" "forecast values" of omegam, w.
      global_smutabstds, AP_inited, & ! xi(s,mu) table of baseline cosmologies
      sepchisqs_uncored_bigcov, sepchisqs_bigcov, & ! values of chisqs, separate schemes
      chisqs_uncored(1), chisqs(1) & ! values of chisqs, averaged over all schemes, correction factor for covmat (D, m1, m2) considered
    )
    fact = (1.0/1.0) 
!    print *, 'End of calling smu_ximu_CalcDAHChisqs_bigcovfc.'
!  subroutine smu_ximu_CalcDAHChisqs_bigcovfc(& 
!    DAs, Hs, & ! List of omegam, w
!    omstd, wstd,  & ! "standard" values of omegam, w. Stored in data2pcffile_base
!    omfc, wfc, & ! "standard" "forecast values" of omegam, w.
!    smutabstds, smutabstds_inited, & ! xi(s,mu) table of baseline cosmologies
!    chisqs_nosyscor, chisqs_syscor, & ! values of chisqs, separate schemes
!    chisqs_nosyscor_all, chisqs_syscor_all & ! values of chisqs, averaged over all schemes, correction factor for covmat (D, m1, m2) considered
!    )

    sepchisqs_uncored_bigcov = sepchisqs_uncored_bigcov * fact
    sepchisqs_uncored = sepchisqs_uncored * fact

    sepchisqs_bigcov = sepchisqs_bigcov * fact
    sepchisqs = sepchisqs * fact

    chisqs_uncored = chisqs_uncored * fact
    chisqs = chisqs * fact

    if(present(chisqs__uncored)) chisqs__uncored = chisqs_uncored
!    chisqs = chisqs_uncored
!
!    print *, 'Sys-corred:   ', real(sepchisqs_uncored_bigcov)
!    print *, 'Sys-uncorred: ', real(sepchisqs_bigcov)
!    stop
    if(printinfo) then
      write(*,'(A,A)') '  * use ', trim(adjustl(AP_MCMCstr(num_fidcos,om_fiducial(1:num_fidcos), w_fiducial(1:num_fidcos))))
!      print *, '(AP_likefc) bigcov chisq = ', real(chisqs(1)), ';  chisq_uncored = ', real(chisqs_uncored(1))
      write(*,'(A,f10.3,A,f10.3)') '(AP_likefc) bigcov chisq = ', real(chisqs(1)), ';  chisq_uncored = ', real(chisqs_uncored(1))
    endif
    if(present(sepchisqs_uncored_result)) then
        sepchisqs_uncored_result = sepchisqs_uncored_bigcov
    endif

    if(present(sepchisqs_result)) then
        sepchisqs_result = sepchisqs_bigcov
    endif

  end subroutine AP_likefc

subroutine check_load_files()
    character(len=charlen) :: tmpstr1, tmpstr2, inputfile, outputfile, printstr, nowfile
    integer :: i, j, iz, imock
    logical :: logvar
    type(omwpar) :: parstd, parnew, nowpar
    real(rt) :: smutab_database(nbins_database, mubins_database, 3), &
                smutab_sysmock(nbins_sysmock, mubins_sysmock, 3), &
                smutab_covmock(nbins_covmock, mubins_covmock, 3)

    real(rt) :: smutab_data(nbins_data, mubins_data, 3)
    real(rt) :: smin_mapping=1.0_rt, smax_mapping=50_rt, DAstd, DAnew, Hstd, Hnew, deltas1, deltas2, intxi(1000), avg, &
      tmpx1,tmpx2,tmpx3,tmpx4,tmpx5,tmpx6,xi
    integer, parameter :: numom = 1, numw = 1
    real(rt) :: omlist(numom), wlist(numw), omstds(1), wstds(1)



    printstr = "Now it is empty!"
    intxi = 0.0_rt
    parstd%omegam = 0.26_rt; parstd%w = -1.0_rt
    parnew%omegam = 0.26_rt; parnew%w = -1.5_rt

    do iz = 1, nz
        print *, 'Redshift = ', zeffs(iz), '...'
        print *
    enddo

    print *, '  (LSS_ximu) checking existence of data base files...'
    do iz = 1, nz
        nowfile=data2pcffile_base(iz, 0.26_rt, -1.0_rt)
        inquire(file=nowfile,exist=logvar)
        if (.not.logvar) then 
            print *, iz, imock, logvar
            print *, 'file not found: ', trim(adjustl(nowfile))
            stop
        endif

        call ximu_loadsmufile(nowfile, smutab_database, smax_database, nbins_database, mubins_database)
        print *, trim(adjustl(nowfile))
        !call DSMapping(smutabstd, nums1, nummu1, smutab2, &
        !  nums2, nummu2, DAstd, DAnew, Hstd, Hnew, deltas1,  deltas2,  smin_mapping, smax_mapping)
        DAstd=DAz_wcdm(parstd,zeffs(iz))
        DAnew=DAz_wcdm(parnew,zeffs(iz))
        Hstd=Hz_wcdm(parstd,zeffs(iz))
        Hnew=Hz_wcdm(parnew,zeffs(iz))
        deltas1 = smax_database / float(nbins_database) 
        deltas2 = smax_data / float(nbins_data)
        !!############################################
        !!############################################
        !!############################################
        !! Begin checking 
        !!############################################
        !!############################################
        !!############################################

        ! 1. smutab_database
        print *, '####################################'
        print *, 'Checking of smutab_database:'
        do i = 30, 31
          j = 100
          print *, i, smutab_database(i,j,1:3), (smutab_database(i,j,1)-2*smutab_database(i,j,2))/smutab_database(i,j,3)+1.0_rt
        enddo
        print *, ' python result is '
        print *, ' 30 [  1.39243300e-09   4.03922200e-10   4.01837100e-10] 2.45479001317'
        print *, ' 31 [  1.18845100e-09   4.70914000e-10   4.29971000e-10] 1.5735805438'

        ! 2. DAs, Hs
        print *, '####################################'
        print *, "Checking: DAstd, DAnew, Hstd, Hnew = ", DAstd, DAnew, Hstd, Hnew
        print *, ' python result is '
        print *, ' ((507.9209000267865, 532.6121518748575, 109.85497398929201, 100.9493258918237)'

        ! DSMapping
        print *, '####################################'
        print *, 'Checking of smutab_database:'
        call DSMapping(smutab_database, nbins_database, mubins_database, smutab_data, &
          nbins_data, mubins_data, DAstd, DAnew, Hstd, Hnew, deltas1,  deltas2,  smin_mapping, smax_mapping)
        do i = 30, 31
          j = 1
          print *, i, smutab_data(i,j,1:3), (smutab_data(i,j,1)-2*smutab_data(i,j,2))/smutab_data(i,j,3)+1.0_rt
        enddo
        print *, ' python result is'
        print *, '30 [2.1285943755830989e-07, 1.9963138302423669e-07, 1.9887429462126441e-07] 0.0627077830991'
        print *, '31 [2.2624904340654802e-07, 2.1190150929401287e-07, 2.1199168691170416e-07] 0.0681050844047'


        
        !XiFun(smutab, deltas, nbins, mubins, anglemin, anglemax, smin, smax, nummuedge, intxi)
        print *, '####################################'
        print *, 'Checking of smutab_database:'
        call XiFun(smutab_data, deltas2, nbins_data, mubins_data, 0.05_rt, 1.0_rt, 6.0_rt, 40.0_rt, 26, intxi(1:25))
        print *, 'intxi = ', real(intxi(1:25))
        print *, 'python result is '
        print *, '[13.53979955829141, 12.888214311950264, 13.051479126899061, ...'
        avg = sum(intxi(1:25)) / 25.0
        intxi(1:25) = intxi(1:25) / avg
        print *
        print *, 'intxi (normed) = ', real(intxi(1:25))
        print *, 'python result is '
        print *, '[0.83193377117529554, 0.791898040299293, 0.80192961518456618,...'
        
        !systematic correction
        call calc_syscor()
        print *, '####################################'
        print *, 'Checking of systematic correction:'
        do i = 1, N1
        do j = 1, N2
          if(mubins(i).eq.25 .and. j.eq.1) then
            print *, ' #-mubin, mucut = ', mubins(i), mucuts(j)
            print *, ' systematic correction @ 1,2-rd bins:'
            print *, real(dintxi_syscor(1:mubins(i)-1,i,j,1))!TMPTEST
            print *, ' python rlt = ', '2.267e-03 1.026e-02 -1.985e-02 ... -1.401e-03 -3.181e-03 '
            print *, ' systematic correction @ 1,6-rd bins:'
            print *, real(dintxi_syscor(1:mubins(i)-1,i,j,5))!TMPTEST
            print *, ' python rlt = ', '1.210e-01 7.030e-02 3.038e-02 ... -5.863e-02 -3.604e-02 '
          endif
        enddo
        enddo
        
        ! Check of covmats...
        if(.false.) then
                call calc_covmats()
                print *, '######################################'
                print *, 'Checking of covmats: '
                print *, 'Python result (Please run code of /home/xiaodongli/LSS/2PCF_AP/2.2.CheckFortran.ipynb): '
!                print *, 'Python result (Please run code of /home/xue/workspace/APLike/2.2.CheckFortran.ipynb): '
                print *, 'Maximal difference between the two matrice:  4.5039261036e-10'
                print *, 'Maximal difference between the two matrice:  4.99956751346e-11'
                print *, 'Maximal difference between the two matrice:  3.09641352055e-10'
                print *, 'Maximal difference between the two matrice:  4.99967656793e-11'
                print *, 'Maximal difference between the two matrice:  1.53330606241e-10'
                !call output_covmats()
        endif
        
        omlist(1) = 0.26_rt
        wlist(1)  = -1.5_rt
        omstds(1) = 0.26_rt
        wstds(1) = -1.0_rt
        ! Check of calcchisqs...
        call load_covmats()
        call invert_covmats()
        call smu_ximu_CalcOmWChisqs(&
            omlist=omlist, numom=numom, wlist=wlist, numw=numw, & ! List of omegam, w
            !omlist=omlist(1:1), numom=1, wlist=wlist(1:1), numw=1, & ! List of omegam, w
            outputdir='/home/xiaodongli/LSS/2PCF_AP/chisqs', &
!            outputdir='/home/xue/workspace/APLike/chisqs', &
            baseoutputfile='Debug_170418', & ! Basic name of the outputfile
            omstds=omstds, wstds=wstds, numomwstds=1, & ! "standard" values of omegam, w: baseline cosmology for (s,mu) mapping. 
        !    omstd=0.11_rt, wstd=-2.0_rt & ! "standard" values of omegam, w: baseline cosmology for (s,mu) mapping. 
            weightedstds=.false., avg_counts = .false. &
            )
        stop
    enddo
    
    
    

    nowpar%omegam=0.26; nowpar%w=0.0
    do iz = 1, nz
        nowfile=data2pcffile_gen(iz,nowpar)
        inquire(file=nowfile, exist=logvar)
        if (.not.logvar) then 
            print *, iz, imock, logvar
            print *, 'file not found: ', trim(adjustl(nowfile))
            stop
        endif
    enddo

    print *, '  (LSS_ximu) checking existence of mock files (for systematic correction)...'
    do iz = 1, nz
     do imock = 1, nsysmocks
        nowfile=syscor2pcffile(iz,imock)
        inquire(file=nowfile,exist=logvar)
        if (.not.logvar) then 
            print *, iz, imock, logvar
            print *, 'file not found: ', trim(adjustl(nowfile))
            stop
        endif
        call ximu_loadsmufile(nowfile, smutab_sysmock, smax_sysmock, nbins_sysmock, mubins_sysmock)
     enddo
    enddo

    print *, '  (LSS_ximu) checking existence of mock files (for covmat estimation)...'
    do iz = 1, nz
     do imock = 1, ncovmocks
        nowfile=cov2pcffile(iz,imock)
        inquire(file=nowfile,exist=logvar)
        if (.not.logvar) then 
            print *, iz, imock, logvar
            print *, 'file not found: ', trim(adjustl(nowfile))
            stop
        endif
        call ximu_loadsmufile(nowfile, smutab_covmock, smax_covmock, nbins_covmock, mubins_covmock)
     enddo
    enddo


    if(iargc().le.1) then
        print *, printstr
        stop
    endif

    outputfile = ""
    do i = 1, iargc()
        if(mod(i,2).eq.0) cycle
        call getarg(i,tmpstr1)
        call getarg(i+1,tmpstr2)
        if(trim(adjustl(tmpstr1)).eq."-inputfile") then
            read(tmpstr2,"(A)") inputfile
        elseif(trim(adjustl(tmpstr1)).eq."-outputfile") then
            read(tmpstr2,"(A)") outputfile
        else
            print *, "Unkown argument: ", trim(adjustl(tmpstr1))
            write(*,"(A)") trim(adjustl(printstr))
            stop
        endif
    enddo

    if(trim(adjustl(outputfile)).eq."") then
        
    endif

    print *, 'This is an empty program LSS_ximu!!'
end subroutine check_load_files
end module AP_funs 

