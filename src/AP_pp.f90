


! post-processing cosmomc files
!  in \Omega_m = 0.27 Lambda CDM
program main
use AP_funs
implicit none
  
  !integer, parameter :: nbasename = 4, numchain = 3
  !integer, parameter :: nbasename = 12, numchain = 4
  integer, parameter :: nbasename = 2, numchain = 4
  character(len=10000) :: basenames(nbasename), suffix, file1,file2,file3,fullfile2,fullfile3,filep0,filep,fullfilep,filerange0,filerange1,filerange2
  character(len=1000) :: tmpstr, tmpstr1, tmpstr2
  integer :: i,ibase,ichain, iz, icol,omcol,wcol,wacol, H0col, maxcol,omegakcol, f1,f2,f3,fullf2,fullf3,f0,fp,fullfp,ifp, iline,nline, npar
  real(rt) :: DAs(nz), Hs(nz), chisqs(nz-1),chisqs__uncored(nz-1), t1,t2, tmpA(10000), wei,lnL,w,wa,H0,omegak,om,minlnL,APlnL,APlnL2,APlnL0
  type(omwpar) :: nowpar
  logical :: smutabstds_inited, inited, printinfo, use_bigcovmat, compute_covmat=.false., output_fullfile=.true.

  !suffix = '_AP001_bin1ref'
  !use_bigcovmat = .false.
  !suffix = '_AP001_bigcov'
  !use_bigcovmat = .true.
  APlnL0 = 40 ! "typical" APlnL value
  !suffix = '_AP001_bigcov'
  !use_bigcovmat = .true.

  basenames(1) = '/home/xiaodongli/software/cosmomc/PLA/base_w/plikHM_TTTEEE_lowTEB/base_w_plikHM_TTTEEE_lowTEB'
  basenames(2) = '/home/xiaodongli/software/cosmomc/PLA/base/plikHM_TTTEEE_lowEB/base_plikHM_TTTEEE_lowEB'
  
  !basenames(1) = '/home/xiaodongli/software/cosmomc/chains/base_wwa_xue18/PLC_w'
  !basenames(1) = '/home/xiaodongli/software/cosmomc/chains/base_wwa_xue18/PLC_BAO_w'
  !basenames(3) = '/home/xiaodongli/software/cosmomc/chains/base_wwa_xue18/PLC_wa'
  !basenames(4) = '/home/xiaodongli/software/cosmomc/chains/base_wwa_xue18/PLC_BAO_wa'

  !basenames(1) = '/home/xiaodongli/software/cosmomc/chains/plcbaos/PLC_BAO_base'
  !basenames(2) = '/home/xiaodongli/software/cosmomc/chains/plcbaos/PLC_BAO_base_w'
  !basenames(3) = '/home/xiaodongli/software/cosmomc/chains/plcbaos/PLC_BAO_base_omegak'
  !basenames(4) = '/home/xiaodongli/software/cosmomc/chains/plcbaos/PLC_BAO_base_mnu'
  !basenames(5) = '/home/xiaodongli/software/cosmomc/chains/plcbaos/PLC_BAO_base_nnu'
  !basenames(6) = '/home/xiaodongli/software/cosmomc/chains/plcbaos/PLC_BAO_base_nrun'
  !basenames(7) = '/home/xiaodongli/software/cosmomc/chains/plcbaos/PLC_BAO_base_r'
  !basenames(8) = '/home/xiaodongli/software/cosmomc/chains/plcbaos/PLC_BAO_base_w_omegak'
  !basenames(9) = '/home/xiaodongli/software/cosmomc/chains/plcbaos/PLC_BAO_base_w_mnu'
  !basenames(10) = '/home/xiaodongli/software/cosmomc/chains/plcbaos/PLC_BAO_base_w_nnu'
  !basenames(11) = '/home/xiaodongli/software/cosmomc/chains/plcbaos/PLC_BAO_base_w_nrun'
  !basenames(12) = '/home/xiaodongli/software/cosmomc/chains/plcbaos/PLC_BAO_base_w_r'

  !use_bigcovmat = .false.; NBComp = .false.
  !suffix = '_AP001_bin1ref_mubin20to25_mucut0.97_fact1to1'
  !suffix = '_AP001_bin1ref_mubin25to30_mucut0.97_fact1to1'
  !suffix = '_AP001_bin1ref_mubin30to35_mucut0.97_fact1to1'
  !suffix = '_AP001_bin1ref_mubin20to25_mucut0.97_fact2to3'
  !suffix = '_AP001_bin1ref_mubin20to25_mucut0.97_fact1to1'

  !-------------------------------------
  ! standard options
  use_bigcovmat = .true.; NBComp = .true.
  !suffix = '_AP001_bigcov_mubin15to20_mucut0.97_fact1to1_NBComp'
  NBComp_Simp=.false.; suffix = '_AP001_bigcov_mubin15to20_mucut0.97_fact1to1_NBCompFull'

  !suffix = '_AP001_bigcov_mubin20to25_mucut0.97_fact1to1_NBComp'
  !NBComp_Simp=.false.; suffix = '_AP001_bigcov_mubin20to25_mucut0.97_fact1to1_NBCompFull'

  !suffix = '_AP001_bigcov_mubin20to25_mucut0.97_fact1to1_NBComp'
  !suffix = '_AP001_bigcov_mubin20to25_mucut0.97_fact1to1_NBCompFull'
  !suffix = '_AP001_bigcov_mubin25to30_mucut0.97_fact1to1_NBComp'
  !suffix = '_AP001_bigcov_mubin30to35_mucut0.97_fact1to1_NBComp'
  !suffix = '_AP001_bigcov_mubin15to20_mucut0.97_fact1to1_NBComp'
  !suffix = '_AP001_bigcov_mubin20to25_mucut0.95_fact1to1_NBComp'
  !suffix = '_AP001_bigcov_mubin20to25_mucut0.90_fact1to1_NBComp'
  !suffix = '_AP001_bigcov_mubin20to25_mucut0.97_s5to40_fact1to1_NBComp'; compute_covmat=.true.
  !suffix = '_AP001_bigcov_mubin20to25_mucut0.97_s4to40_fact1to1_NBComp'; compute_covmat=.true.
  !suffix = '_AP001_bigcov_mubin20to25_mucut0.97_s7to40_fact1to1_NBComp'; compute_covmat=.true.
  !suffix = '_AP001_bigcov_mubin20to25_mucut0.97_s8to40_fact1to1_NBComp'; compute_covmat=.true.
  !suffix = '_AP001_bigcov_mubin20to25_mucut0.95to0.85_fact1to1_NBCompomegak,'
  !suffix = '_AP001_bigcov_mubin15to20_mucut0.95to0.85_fact1to1_NBComp'
  !suffix = '_AP001_bigcov_mubin20to25_mucut0.97_fact4to5_NBComp'

  AP_inited = .false.
  printinfo = .false.
  do ibase = 1, nbasename
    ! columns of om/H0/w/wa
    print *, '#############################'  
    filep0 =  trim(adjustl(basenames(ibase)))//'.paramnames'
    filerange0 =  trim(adjustl(basenames(ibase)))//'.ranges'
    icol = 1;   wcol=-1;wacol=-1;omcol=-1;H0col=-1;omegakcol=-1
    f0 = 9802;open(unit=f0,file=trim(adjustl(filep0)),action="read")
    do while(.true.)
      read(f0,*,end=100) tmpstr1, tmpstr2
      if(trim(adjustl(tmpstr1)).eq.'w') wcol = icol
      if(trim(adjustl(tmpstr1)).eq.'wa') wacol = icol
      if(trim(adjustl(tmpstr1)).eq.'omegam*') omcol = icol
      if(trim(adjustl(tmpstr1)).eq.'H0*') H0col = icol
      if(trim(adjustl(tmpstr1)).eq.'omegak') omegakcol = icol
      icol = icol+1
      cycle
100   exit
    enddo
    maxcol = max(omcol,wcol,wacol,H0col);
    npar = icol - 1
    write(*,'(A,5i4,A,2i4,3A)') 'columns of om,w,wa,H0,omegak = ', omcol, wcol, wacol, H0col, omegakcol, &
       ' /  maxcol, #-par =', maxcol,npar, &
      ' (from ', trim(adjustl(filep0)), ')'
    close(f0)
 
    ! new paramname file for om/H0/w/wa/omegak
    do ifp = 1, 2
      if(ifp.eq.1) then
         filep =  trim(adjustl(basenames(ibase)))//trim(adjustl(suffix))//'.paramnames'
         fullfilep = trim(adjustl(basenames(ibase)))//trim(adjustl(suffix))//'_full.paramnames'
         filerange1 = trim(adjustl(basenames(ibase)))//trim(adjustl(suffix))//'.ranges'
         filerange2 = trim(adjustl(basenames(ibase)))//trim(adjustl(suffix))//'_full.ranges'
      else
         filep =  trim(adjustl(basenames(ibase)))//trim(adjustl(suffix))//'_nosyscor.paramnames'
         fullfilep = trim(adjustl(basenames(ibase)))//trim(adjustl(suffix))//'_nosyscor_full.paramnames'
         filerange1 = trim(adjustl(basenames(ibase)))//trim(adjustl(suffix))//'_nosyscor.ranges'
         filerange2 = trim(adjustl(basenames(ibase)))//trim(adjustl(suffix))//'_nosyscor_full.ranges'
      endif
      fp = 2113259;open(unit=fp,file=trim(adjustl(filep)),action="write")
      write(fp,'(A)') 'omegam    \Omega_m'
      write(fp,'(A)') 'H0        H_0'
      write(fp,'(A)') 'w         w'
      write(fp,'(A)') 'wa        w_a'
      write(fp,'(A)') 'omegak    \Omega_k'
      close(fp)
      if(output_fullfile) then
        call system('cp '//trim(adjustl(filep0))//' '//trim(adjustl(fullfilep)))
        call system('cp '//trim(adjustl(filerange0))//' '//trim(adjustl(filerange1)))
        call system('cp '//trim(adjustl(filerange0))//' '//trim(adjustl(filerange2)))
      endif
    enddo
    !cycle
    do ichain = 1, numchain
      write(tmpstr, *) ichain
      file1 = trim(adjustl(basenames(ibase)))//'_'//trim(adjustl(tmpstr))//'.txt'
      file2 = trim(adjustl(basenames(ibase)))//trim(adjustl(suffix))//'_'//trim(adjustl(tmpstr))//'.txt'
      file3 = trim(adjustl(basenames(ibase)))//trim(adjustl(suffix))//'_nosyscor_'//trim(adjustl(tmpstr))//'.txt'
      fullfile2 = trim(adjustl(basenames(ibase)))//trim(adjustl(suffix))//'_full_'//trim(adjustl(tmpstr))//'.txt'
      fullfile3 = trim(adjustl(basenames(ibase)))//trim(adjustl(suffix))//'_nosyscor_full_'//trim(adjustl(tmpstr))//'.txt'
      f1 = 1328; f2 = 1026; f3=920374; fullf2 = 913408; fullf3 = 7123859;

      ! get minlnL
      open(unit=f1,file=file1,action='read'); nline=1;
      minlnL = 1.0e30
      do while(.true.)
       read(f1,*,end=1000) tmpA(1:2)
       minlnL = min(tmpA(2),minlnL); nline=nline+1
       cycle
1000   exit
      enddo
      close(f1)

      ! post process
      iline=1
      open(unit=f1,file=file1,action='read'); 
      open(unit=f2,file=file2,action='write'); 
      open(unit=f3,file=file3,action='write'); 
      write(*,'(A,A)') ' * read from: ', trim(adjustl(file1))
      write(*,'(A,A)') ' * write to : ', trim(adjustl(file2))
      write(*,'(A,A)') ' * write to : ', trim(adjustl(file3))
      if(output_fullfile) then
        open(unit=fullf2,file=fullfile2,action='write'); 
        open(unit=fullf3,file=fullfile3,action='write'); 
        write(*,'(A,A)') ' * write to : ', trim(adjustl(fullfile2))
        write(*,'(A,A)') ' * write to : ', trim(adjustl(fullfile3))
      endif
      do while(.true.)
        read(f1,*,end=101) tmpA(1:npar+2);
        wei=tmpA(1); lnL=tmpA(2); om=tmpA(omcol+2); H0=tmpA(H0col+2); 
        w=-1; if(wcol.ne.-1) w=tmpA(wcol+2); 
        wa=0; if(wacol.ne.-1) wa=tmpA(wacol+2); 
        omegak=0; if(omegakcol.ne.-1) omegak=tmpA(omegakcol+2); 

        ! AP like
        nowpar%omegam=om; nowpar%w=w; nowpar%wa=wa; nowpar%omegak=omegak;
        do iz = 1, nz
          DAs(iz) = DAz_wcdm(nowpar,zeffs(iz))
          Hs(iz)  = Hz_wcdm(nowpar,zeffs(iz))
        enddo
        if(iline.eq.1) then; printinfo=.true.; else; printinfo=.false.; endif;
        call AP_Like(DAs,Hs,chisqs,printinfo=printinfo,use_bigcovmat=use_bigcovmat,chisqs__uncored=chisqs__uncored,compute_covmat=compute_covmat) 
        compute_covmat=.false.
        if(use_bigcovmat) then
          APlnL  = chisqs(1)/2.0; 
          APlnL2 = chisqs__uncored(1)/2.0
        else
          APlnL  = sum(chisqs(1:nz-1))/2.0
          APlnL2 = sum(chisqs__uncored(1:nz-1))/2.0
        endif
        write(f2,'(7e15.7)') wei*exp(APlnL0-APlnL),  lnL+APlnL,  om,H0,w,wa, omegak
        write(f3,'(7e15.7)') wei*exp(APlnL0-APlnL2), lnL+APlnL2, om,H0,w,wa, omegak
        if(output_fullfile) then
          write(fullf2,'(<npar+2>e15.7)') wei*exp(APlnL0-APlnL),  lnL+APlnL,  tmpA(3:npar+2)
          write(fullf3,'(<npar+2>e15.7)') wei*exp(APlnL0-APlnL2), lnL+APlnL2, tmpA(3:npar+2)
        endif
        if(mod(iline,200).eq.1) then
          write(*,'(i7,A,i7,A,5f7.3,A,e10.3,e12.5, A,e10.3,A,2f8.2)') iline,'/',nline, ' * om/H0/w/wa/omegak =', om,H0,w,wa,omegak, &
                '; wei/lnL=',wei,lnL, '  (minlnL=',minlnL,') ** APlnL = ', APlnL,APlnL2
        endif
        iline=iline+1
        cycle
101     exit
      enddo
      close(f1); close(f2); close(f3);
      if(output_fullfile) then
        close(fullf2); close(fullf3);
      endif
    enddo
      
  enddo


  stop

  !! test region
  AP_inited = .false.
  printinfo = .true.
  do i = 1, 10
  nowpar%omegam = 0.2_rt + 0.02_rt * (i-1); nowpar%w = -1.0_rt
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
  print *, nowpar%omegam
  call AP_Like(DAs, Hs,  chisqs, printinfo, use_bigcovmat=use_bigcovmat) 
  if(i.eq.1)   call cpu_time(t1)
  enddo
  call cpu_time(t2)
  print *, 'Total time used: ', t2-t1
  !! end of test region

end program
