!==============================================================
!     diffrad_vm.f90 -- combined vector meson generator
!     phi (ivec=3, exp t-form) + J/psi (ivec=4, dipole t-form)
!     with TUNED defaults (June 2026); parameters via input file
!
!     MC Generator for diffractive vector meson electroproduction
!     with QED radiative corrections (collinear approximation)
!     Based on DIFFRAD by I.Akushevich (1998)
!
!     LUND output format (px py pz E mass):
!       1: beam e-      (initial, shifted if ISR)
!       2: target p     (initial, at rest)
!       3: scattered e- (final, shifted if FSR)
!       4: pi+          (from rho decay)
!       5: pi-          (from rho decay)
!       6: recoil p     (final)
!       7: gamma        (final, hard events only)
!==============================================================

!==============================================================
!  input_mod: keyword-based input reader
!  Format:  key  value  [value2]      (! = comment)
!  Range parameters:  key  min  max
!==============================================================
      module input_mod
      implicit none

      integer, parameter :: npars = 29

      type param_t
        character(len=32) :: name
        integer  :: itype    ! 1=real*8, 2=integer
        integer  :: nvals    ! 1=single value, 2=range (min max)
        real*8   :: rval     ! value or range min
        real*8   :: rval2    ! range max (nvals==2 only)
        integer  :: ival
        real*8   :: minval, maxval
        logical  :: is_set
      end type param_t

      contains

!--------------------------------------------------------------
      subroutine read_input(fname, pars, np)
      implicit none
      character(len=*),  intent(in)    :: fname
      integer,           intent(in)    :: np
      type(param_t),     intent(inout) :: pars(np)
      character(len=256) :: line, clean, key, rest
      integer :: ios, ios2, i, ic, p1
      logical :: found
      open(9, file=trim(fname), status='old', iostat=ios)
      if(ios /= 0) then
        write(*,*) 'ERROR: cannot open: '//trim(fname);  stop
      endif
      do
        read(9,'(A)',iostat=ios) line
        if(ios /= 0) exit
        ic = index(line,'!')
        if(ic > 0) line = line(1:ic-1)
        clean = adjustl(line)
        if(len_trim(clean) == 0) cycle
        p1 = scan(clean,' '//char(9))
        if(p1 <= 1) cycle
        key  = clean(1:p1-1)
        rest = adjustl(clean(p1:))
        if(len_trim(rest) == 0) cycle
        found = .false.
        do i = 1, np
          if(trim(key) == trim(pars(i)%name)) then
            if(pars(i)%itype == 1) then
              if(pars(i)%nvals == 2) then
                read(rest,*,iostat=ios2) pars(i)%rval, pars(i)%rval2
              else
                read(rest,*,iostat=ios2) pars(i)%rval
              endif
            else
              read(rest,*,iostat=ios2) pars(i)%ival
            endif
            if(ios2 /= 0) then
              write(*,*) 'ERROR: bad value for: '//trim(key);  stop
            endif
            pars(i)%is_set = .true.
            found = .true.
          endif
        end do
        if(.not.found) write(*,*) 'WARNING: unknown parameter: '//trim(key)
      end do
      close(9)
      end subroutine read_input

!--------------------------------------------------------------
      subroutine validate(pars, np)
      implicit none
      integer,       intent(in) :: np
      type(param_t), intent(in) :: pars(np)
      integer :: i
      do i = 1, np
        if(.not.pars(i)%is_set) then
          write(*,*) 'ERROR: missing parameter: '//trim(pars(i)%name); stop
        endif
        if(pars(i)%itype == 1) then
          if(pars(i)%rval < pars(i)%minval .or. pars(i)%rval > pars(i)%maxval) then
            write(*,*) 'ERROR: out of range: ', trim(pars(i)%name), &
                       '  val=', pars(i)%rval;  stop
          endif
          if(pars(i)%nvals == 2) then
            if(pars(i)%rval2 < pars(i)%minval .or. &
               pars(i)%rval2 > pars(i)%maxval) then
              write(*,*) 'ERROR: out of range: ', trim(pars(i)%name), &
                         '  max=', pars(i)%rval2;  stop
            endif
            if(pars(i)%rval > pars(i)%rval2) then
              write(*,*) 'ERROR: min > max for: '//trim(pars(i)%name); stop
            endif
          endif
        else
          if(dble(pars(i)%ival) < pars(i)%minval .or. &
             dble(pars(i)%ival) > pars(i)%maxval) then
            write(*,*) 'ERROR: out of range: ', trim(pars(i)%name), &
                       '  val=', pars(i)%ival;  stop
          endif
        endif
      end do
      end subroutine validate

!--------------------------------------------------------------
      subroutine print_params(pars, np)
      implicit none
      integer,       intent(in) :: np
      type(param_t), intent(in) :: pars(np)
      integer :: i
      write(*,'(a)') '=== DIFFRAD input parameters ==='
      do i = 1, np
        if(pars(i)%itype == 2) then
          write(*,'(2X,A12,I10)') trim(pars(i)%name), pars(i)%ival
        else if(pars(i)%nvals == 2) then
          write(*,'(2X,A12,F10.4,A4,F10.4)') &
            trim(pars(i)%name), pars(i)%rval,' .. ',pars(i)%rval2
        else
          write(*,'(2X,A12,F12.5)') trim(pars(i)%name), pars(i)%rval
        endif
      end do
      write(*,'(a)') '================================'
      end subroutine print_params

      end module input_mod

!==============================================================
      program diffrad_gen
      use input_mod
      implicit real*8(a-h,o-z)
      real*4 urand

      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/bwpar/amv0,gamv
      common/sxy/s,x,sx,sxp,q2,w2,aly,anu,sqly,an,tamin,tamax,xs,ys
      common/phi/phirad,tdif,phidif,tq,vmax,ivec
      common/tpar/tslope
      common/vv1vv2/aa1,aa2,bb,sib
      common/ivv/vcurr,cutv
      common/cuts/wmin2,ymin,ymax
      common/amf2/taa,atm(8,6),sfm0(8)
      common/sigsig/sigmat0,sigmal0
      common/pri/ipri
      common/phimodel/phi_alf1,phi_alf2,phi_alf3,phi_nuT,phi_bt,phi_cR
      common/jpsimodel/pj_alf1,pj_alf2,pj_alf3,pj_nuT,pj_mg2,pj_cR
      real*8 pj_alf1,pj_alf2,pj_alf3,pj_nuT,pj_mg2,pj_cR

      real*8 k1(4),k2(4),ptar(4),ph(4),pp(4),kgam(4),pip(4),pim(4)
      integer*4 iy
      integer*8 ntry, maxtry
      character(len=256) :: input_file, lund_file, stat_file, vdist_file
      character(len=256) :: cl_arg
      integer :: narg_tot, iarg_cur, llen
      type(param_t) :: pars(npars)

!     Parse command-line arguments: -input <file>  -lund <file>
      input_file = 'gen_input.dat'
      lund_file  = ''
      narg_tot   = command_argument_count()
      iarg_cur   = 1
      do while(iarg_cur .le. narg_tot)
        call get_command_argument(iarg_cur, cl_arg)
        if(trim(cl_arg) .eq. '-input') then
          iarg_cur = iarg_cur + 1
          call get_command_argument(iarg_cur, input_file)
        elseif(trim(cl_arg) .eq. '-lund') then
          iarg_cur = iarg_cur + 1
          call get_command_argument(iarg_cur, lund_file)
        endif
        iarg_cur = iarg_cur + 1
      enddo

!     Initialize parameter table
!                   name    itype nvals rval  rval2 ival   minval   maxval  is_set
      pars(1)  = param_t('bmom',  1,1, 0d0,0d0, 0,    1d0, 1000d0, .false.)
      pars(2)  = param_t('tmom',  1,1, 0d0,0d0, 0,    0d0,10000d0, .false.)
      pars(3)  = param_t('lepton',2,1, 0d0,0d0, 0,    1d0,    2d0, .false.)
      pars(4)  = param_t('ivec',  2,1, 0d0,0d0, 0,    1d0,    4d0, .false.)
      pars(5)  = param_t('cutv',  1,1, 0d0,0d0, 0,  -10d0,   10d0, .false.)
      pars(6)  = param_t('nev',   2,1, 0d0,0d0, 0,    1d0,   1d9,  .false.)
      pars(7)  = param_t('iy',    2,1, 0d0,0d0, 0,    1d0, 2.15d9, .false.)
      pars(8)  = param_t('Q2',    1,2, 0d0,0d0, 0,    0d0,  100d0, .false.)
      pars(9)  = param_t('y',     1,2, 0d0,0d0, 0,    0d0,    1d0, .false.)
      pars(10) = param_t('t',     1,2, 0d0,0d0, 0,    0d0,   10d0, .false.)
      pars(11) = param_t('tslope',1,1, 0d0,0d0, 0,    0d0,  100d0, .false.)
      pars(12) = param_t('iborn', 2,1, 0d0,0d0, 0,    0d0,    1d0, .false.)
      pars(13) = param_t('W',     1,2, 0d0,0d0, 0,    0d0, 1000d0, .false.)
      pars(14) = param_t('xB',          1,2, 0d0,  0d0, 0,    0d0,    1d0, .false.)
!     Optional acceptance cuts (not required; defaults = no cut)
      pars(15) = param_t('momentum_electron',1,2, 0d0,1d4,0, 0d0, 1d4, .false.)
      pars(16) = param_t('theta_electron',   1,2, 0d0,180d0,0, 0d0, 180d0, .false.)
      pars(17) = param_t('momentum_proton',  1,2, 0d0,1d4,0, 0d0, 1d4, .false.)
      pars(18) = param_t('theta_proton',     1,2, 0d0,180d0,0, 0d0, 180d0, .false.)
      pars(19) = param_t('momentum_hplus',   1,2, 0d0,1d4,0, 0d0, 1d4, .false.)
      pars(20) = param_t('theta_hplus',      1,2, 0d0,180d0,0, 0d0, 180d0, .false.)
      pars(21) = param_t('momentum_hminus',  1,2, 0d0,1d4,0, 0d0, 1d4, .false.)
      pars(22) = param_t('theta_hminus',     1,2, 0d0,180d0,0, 0d0, 180d0, .false.)
!     --- Cross-section model parameters (optional, have defaults) ---
!     Defaults are the TUNED values (June 2026):
!       phi  (ivec=3, exponential t):  alf2=-1.245 alf3=0.762 nuT=2.344
!                                      bt=1.284  cR=1.0
!       jpsi (ivec=4, dipole t):       alf2=4.122  alf3=0.32  nuT=3.0
!                                      mg2=3.112 cR=0.4
!     Table defaults below are the PHI values; for ivec=4 any parameter
!     not set in the input file is replaced by the J/psi default in code.
      pars(23) = param_t('alf1',  1,1, 400d0,0d0, 0,  0d0, 1d6,  .false.)
      pars(24) = param_t('alf2',  1,1,-1.245d0,0d0,0, -50d0, 50d0, .false.)
      pars(25) = param_t('alf3',  1,1, 0.762d0,0d0,0, -50d0, 50d0, .false.)
      pars(26) = param_t('nuT',   1,1, 2.344d0,0d0,0,  0d0, 50d0, .false.)
      pars(27) = param_t('bt',    1,1, 1.284d0,0d0,0,  0d0,100d0, .false.)
      pars(28) = param_t('cR',    1,1, 1.000d0,0d0,0,  0d0, 50d0, .false.)
!     --- J/psi only: dipole mass^2 (tuned 2026-06) ---
      pars(29) = param_t('mg2',   1,1, 3.112d0,0d0,0, 0.1d0,100d0, .false.)

!     Read, validate, print
      call read_input(trim(input_file), pars, npars)
      call validate(pars, 14)          ! only required parameters
      call print_params(pars, npars)

!     Extract values
      bmom   = pars(1)%rval
      tmom   = pars(2)%rval
      lepton = pars(3)%ival
      ivec   = pars(4)%ival
      cutv   = pars(5)%rval
      nev    = pars(6)%ival
      iy     = pars(7)%ival
      q2min  = pars(8)%rval;  q2max  = pars(8)%rval2
      ymin   = pars(9)%rval;  ymax   = pars(9)%rval2
      tmin   = pars(10)%rval; tmax   = pars(10)%rval2   ! |t| min max (positive)
      tslope = pars(11)%rval
      iborn  = pars(12)%ival
      wmin   = pars(13)%rval
      xbmin  = pars(14)%rval; xbmax  = pars(14)%rval2
      wmin2  = wmin*wmin
!     Optional cuts — if not set in input file, defaults mean no cut
      if(pars(15)%is_set) then
        pe_min = pars(15)%rval;  pe_max = pars(15)%rval2
      else
        pe_min = 0d0;            pe_max = 1d4
      endif
      if(pars(16)%is_set) then
        the_min = pars(16)%rval; the_max = pars(16)%rval2
      else
        the_min = 0d0;           the_max = 180d0
      endif
      if(pars(17)%is_set) then
        pp_min = pars(17)%rval;  pp_max = pars(17)%rval2
      else
        pp_min = 0d0;            pp_max = 1d4
      endif
      if(pars(18)%is_set) then
        thp_min = pars(18)%rval; thp_max = pars(18)%rval2
      else
        thp_min = 0d0;           thp_max = 180d0
      endif
      if(pars(19)%is_set) then
        php_min = pars(19)%rval;  php_max = pars(19)%rval2
      else
        php_min = 0d0;            php_max = 1d4
      endif
      if(pars(20)%is_set) then
        thhp_min = pars(20)%rval; thhp_max = pars(20)%rval2
      else
        thhp_min = 0d0;           thhp_max = 180d0
      endif
      if(pars(21)%is_set) then
        phm_min = pars(21)%rval;  phm_max = pars(21)%rval2
      else
        phm_min = 0d0;            phm_max = 1d4
      endif
      if(pars(22)%is_set) then
        thhm_min = pars(22)%rval; thhm_max = pars(22)%rval2
      else
        thhm_min = 0d0;           thhm_max = 180d0
      endif
      if(pars(15)%is_set .or. pars(16)%is_set) then
        write(*,'(a)') ' Electron acceptance cuts active:'
        write(*,'(a,f8.3,a,f8.3,a)') '   |p|_e = [',pe_min,' ,',pe_max,' ] GeV/c'
        write(*,'(a,f8.3,a,f8.3,a)') '   the_e = [',the_min,' ,',the_max,' ] deg'
      endif
      if(pars(17)%is_set .or. pars(18)%is_set) then
        write(*,'(a)') ' Proton acceptance cuts active:'
        write(*,'(a,f8.3,a,f8.3,a)') '   |p|_p = [',pp_min,' ,',pp_max,' ] GeV/c'
        write(*,'(a,f8.3,a,f8.3,a)') '   the_p = [',thp_min,' ,',thp_max,' ] deg'
      endif
      if(pars(19)%is_set .or. pars(20)%is_set) then
        write(*,'(a)') ' h+ acceptance cuts active:'
        write(*,'(a,f8.3,a,f8.3,a)') '   |p|h+ = [',php_min,' ,',php_max,' ] GeV/c'
        write(*,'(a,f8.3,a,f8.3,a)') '   th_h+ = [',thhp_min,' ,',thhp_max,' ] deg'
      endif
      if(pars(21)%is_set .or. pars(22)%is_set) then
        write(*,'(a)') ' h- acceptance cuts active:'
        write(*,'(a,f8.3,a,f8.3,a)') '   |p|h- = [',phm_min,' ,',phm_max,' ] GeV/c'
        write(*,'(a,f8.3,a,f8.3,a)') '   th_h- = [',thhm_min,' ,',thhm_max,' ] deg'
      endif

!     --- phi model parameters (table defaults = tuned phi values) ---
      phi_alf1 = pars(23)%rval
      phi_alf2 = pars(24)%rval
      phi_alf3 = pars(25)%rval
      phi_nuT  = pars(26)%rval
      phi_bt   = pars(27)%rval
      phi_cR   = pars(28)%rval
!     --- J/psi model parameters: input value if set, else tuned default ---
      pj_alf1 = 400.d0
      if(pars(23)%is_set) pj_alf1 = pars(23)%rval
      pj_alf2 = 4.122d0
      if(pars(24)%is_set) pj_alf2 = pars(24)%rval
      pj_alf3 = 0.320d0
      if(pars(25)%is_set) pj_alf3 = pars(25)%rval
      pj_nuT  = 3.000d0
      if(pars(26)%is_set) pj_nuT  = pars(26)%rval
      pj_cR   = 0.400d0
      if(pars(28)%is_set) pj_cR   = pars(28)%rval
      pj_mg2  = pars(29)%rval
      if(ivec.eq.4) then
        write(*,'(a)')       ' J/psi model parameters (dipole t):'
        write(*,'(a,f12.4)') '   alf1 = ', pj_alf1
        write(*,'(a,f12.4)') '   alf2 = ', pj_alf2
        write(*,'(a,f12.4)') '   alf3 = ', pj_alf3
        write(*,'(a,f12.4)') '   nuT  = ', pj_nuT
        write(*,'(a,f12.4)') '   mg2  = ', pj_mg2
        write(*,'(a,f12.4)') '   cR   = ', pj_cR
      else
        write(*,'(a)')       ' Phi model parameters (exp t):'
        write(*,'(a,f12.4)') '   alf1 = ', phi_alf1
        write(*,'(a,f12.4)') '   alf2 = ', phi_alf2
        write(*,'(a,f12.4)') '   alf3 = ', phi_alf3
        write(*,'(a,f12.4)') '   nuT  = ', phi_nuT
        write(*,'(a,f12.4)') '   bt   = ', phi_bt
        write(*,'(a,f12.4)') '   cR   = ', phi_cR
      endif

      call setcon(ivec,lepton)
      s = 2d0*(sqrt(tmom**2+amp**2)*sqrt(bmom**2+aml2)+bmom*tmom)
      ebeam = bmom

      write(*,'(a)') '=============================='
      write(*,'(a,f8.3,a)') ' Ebeam  = ',ebeam,' GeV'
      write(*,'(a,f8.3,a)') ' sqrt(S)= ',sqrt(s),' GeV'
      write(*,'(a,i2)')     ' ivec   = ',ivec
      write(*,'(a,f6.3,a)') ' cutv   = ',cutv,' GeV^2'
      write(*,'(a,f6.3,a)') ' wmin   = ',wmin,' GeV'
      write(*,'(a,2f7.4)')  ' xB     = ',xbmin,xbmax
      write(*,'(a,i8)')     ' nev    = ',nev
      if(iborn.eq.1) write(*,'(a)') ' Mode: BORN ONLY (no RC)'
      write(*,'(a)') '=============================='

!     Default lund filename if not given on command line
      if(len_trim(lund_file) .eq. 0) then
        if(iborn.eq.1) then
          lund_file = 'born_events.lund'
        else
          lund_file = 'rc_events.lund'
        endif
      endif
!     Derive stat and vdist filenames: strip .lund, append _stat.dat / _vdist.dat
      llen = len_trim(lund_file)
      if(llen .gt. 5 .and. lund_file(llen-4:llen) .eq. '.lund') then
        stat_file  = lund_file(1:llen-5)//'_stat.dat'
        vdist_file = lund_file(1:llen-5)//'_vdist.dat'
      else
        stat_file  = trim(lund_file)//'_stat.dat'
        vdist_file = trim(lund_file)//'_vdist.dat'
      endif
      write(*,'(a)') ' Output: '//trim(lund_file)
      write(*,'(a)') ' Stats:  '//trim(stat_file)
      open(unit=10, file=trim(lund_file))
      open(unit=11, file=trim(stat_file))
      open(unit=12, file=trim(vdist_file))

      ngen=0; nphys=0; nsoft=0; nhard=0; nisr=0; nfsr=0; ntry=0
      nf_born=0; nf_thresh=0; nf_tkin=0; nf_vmax=0
      nf_sib=0; nf_sshxxh=0; nf_sxqv=0; nf_sigtot=0
      wsum=0d0; wsum2=0d0   ! for cross section integration
      vcurr=0d0; ipri=0      ! for qqt/podinl
      wmax=0d0              ! maximum weight for accept/reject

!     ── Warm-up pass to find wmax ──────────────────────────────
      write(*,'(a)') ' Finding wmax (warm-up)...'
      nwarm = 10000
      do iwarm = 1, nwarm
        call sample_born(q2min,q2max,xbmin,xbmax,tmin,tmax,iy,sg_born,iacc)
        if(iacc.eq.0) cycle
        call conkin(s)
        call sample_bw(amv0,gamv,ivec,iy,sqrt(max(0d0,w2))-amp,amv)
        if(sqrt(w2).lt.amp+amv) cycle
        tt1=w2-q2-amp2; tt2=w2-amp2+amv**2
        disc1=tt1**2+4d0*q2*w2; disc2=tt2**2-4d0*amv**2*w2
        if(disc1.lt.0d0.or.disc2.lt.0d0) cycle
        tdmink=-q2+amv**2-.5d0/w2*(tt1*tt2+sqrt(disc1)*sqrt(disc2))
        tdmaxk=-q2+amv**2-.5d0/w2*(tt1*tt2-sqrt(disc1)*sqrt(disc2))
        if(tdif.lt.tdmink.or.tdif.gt.tdmaxk) cycle
        call bornin(sib)
        if(sib.le.0d0) cycle
!       Warm-up: use sg_born*sib for both Born and RC runs
!       For RC run, multiply by factor 2 as safety margin
!       (exact qqt can exceed Born but typically by <2x)
          wtrial = sg_born * sib
        if(wtrial.gt.wmax) wmax = wtrial
      enddo
      if(iborn.eq.1)then
        wmax = wmax * 1.5d0   ! 50% margin for Born
      else
        wmax = wmax * 3.0d0   ! 3x margin for RC (covers RC enhancement)
      endif
      write(*,'(a,g12.4)') ' wmax = ',wmax
      write(*,'(a)') ' Starting main generation...'

!     ── Main event loop ──────────────────────────────────────────
      do while (ngen .lt. nev)

        ntry = ntry + 1
!       Attempt limit: large for J/psi (low efficiency near threshold),
!       moderate for rho/phi/omega
        maxtry = 1000_8*nev
        if(ivec.eq.4) maxtry = 100000_8*nev
        if(ntry .gt. maxtry) then
          write(*,*) 'ERROR: too many attempts'
          goto 999
        endif

!       Step 1: Sample Born kinematics into common blocks
        call sample_born(q2min,q2max,xbmin,xbmax,tmin,tmax, &
                         iy,sg_born,iacc)
        if(iacc.eq.0)then; nf_born=nf_born+1; goto 100; endif

!       Step 2: Compute derived kinematics
        call conkin(s)

!       Sample vector meson mass from Breit-Wigner
        call sample_bw(amv0,gamv,ivec,iy,sqrt(max(0d0,w2))-amp,amv)

        if(sqrt(w2).lt.amp+amv)then
          nf_thresh=nf_thresh+1; goto 100
        endif

        tt1 = w2-q2-amp2
        tt2 = w2-amp2+amv**2
        disc1 = tt1**2+4d0*q2*w2
        disc2 = tt2**2-4d0*amv**2*w2
        if(disc1.lt.0d0.or.disc2.lt.0d0)then
          nf_born=nf_born+1; goto 100
        endif
        tdmink = -q2+amv**2-.5d0/w2*(tt1*tt2+sqrt(disc1)*sqrt(disc2))
        tdmaxk = -q2+amv**2-.5d0/w2*(tt1*tt2-sqrt(disc1)*sqrt(disc2))
        if(tdif.lt.tdmink.or.tdif.gt.tdmaxk)then
          nf_tkin=nf_tkin+1; goto 100
        endif

!       Step 3: RC quantities
        sxt = sx+tdif
        tq  = q2+tdif-amv**2
        aa1 = (q2*sxp*sxt-(s*sx+2d0*amp2*q2)*tq)/2d0/aly
        aa2 = (q2*sxp*sxt-(x*sx-2d0*amp2*q2)*tq)/2d0/aly
        sqbb1 = sqrt(max(0d0,q2*sxt**2-sxt*sx*tq-amp2*tq**2-amv**2*aly))
        sqbb2 = sqrt(max(0d0,q2*(s*x-amp2*q2)-aml2*aly))
        bb = sqbb1*sqbb2/aly

        vmax_kin = tt2+.5d0/q2*(-tt1*tq &
          -sqrt(max(0d0,tt1**2+4d0*q2*w2)) &
          *sqrt(max(0d0,tq**2+4d0*amv**2*q2))) &
          - 1d-8
        if(cutv.gt.1d-12)then
          vmax = min(vmax_kin,cutv)
        else
          vmax = vmax_kin
        endif
        if(vmax.le.0d0)then; nf_vmax=nf_vmax+1; goto 100; endif

        call bornin(sib)
        if(sib.le.0d0)then; nf_sib=nf_sib+1; goto 100; endif

!       Jacobian dQ²/dy = Q²/y (since Q² = s·xB·y), inverted to give
!       rconv = y/Q²: converts qqt output (Akushevich dx_B dy dt) to
!       Diehl convention (dx_B dQ² dt).  No bornin_ak call needed.
        rconv = ys/q2

        vv1 = aa1/2d0
        vv2 = aa2/2d0
        ssh = x+q2-vv2
        xxh = s-q2-vv1
        if(ssh.le.0d0.or.xxh.le.0d0)then
          nf_sshxxh=nf_sshxxh+1; goto 100
        endif

        dlm = log(q2/aml2)
        delinf_val = (dlm-1d0)*log(vmax**2/ssh/xxh)
        deltavr = (1.5d0*dlm-2d0-.5d0*log(xxh/ssh)**2 &
                   +fspen(1d0-amp2*q2/ssh/xxh)-pi**2/6d0)
        delta_vac = vacpol(q2)
        extai1 = exp(alpha/pi*delinf_val)
        sig_soft = sib*extai1*(1d0+alpha/pi*(deltavr+delta_vac))

!       Exact hard cross section via qqt (integrates podinl)
!       Only needed for RC run; skip for Born-only to avoid slow integration
!
!       sig_hard_fix scheme:
!        - sig_total (the RATE) keeps the validated Bardin-Shumeiko form
!          sig_soft + sig_F with the SIGNED finite remainder sig_F
!          (no clamp -> removes the O(0.1%) clamp bias).
!        - The hard-event SAMPLING cross section is the physical
!          sigma_hard(v>vcut) = sig_F + sib*(alpha/pi)(dlm-1)ln(vmax^2/vcut^2)
!          (the analytic soft-kernel integral over [vcut,vmax] added back).
!          Positive by construction up to O(alpha^2); clamped as numerical
!          safety. vcut_ir = 1d-2 GeV^2 (omega ~ 5 MeV) chosen at the
!          detector-insensitivity scale; must match vcut_use below.
        if(iborn.eq.0) then
          phidif = phirad
          call difflt(q2,w2,tdif,sigmal0,sigmat0)
          call qqt(sig_hard_exact)
!         Scale to Diehl convention: rconv * Akushevich
          sig_F = rconv * sig_hard_exact
          vcut_ir = 1d-2
          ana_soft = sib*alpha/pi*(dlm-1d0)*log(vmax**2/vcut_ir**2)
          sig_hard = max(0d0, sig_F + ana_soft)
          sig_total = sig_soft + sig_F
        else
          sig_hard = 0d0
          sig_total = sig_soft
        endif
        if(sig_total.le.0d0)then; nf_sigtot=nf_sigtot+1; goto 100; endif
        sig_hard = min(sig_hard, sig_total)

        sg_born_full = sg_born * sib
        if(sg_born_full.le.0d0)then
          nf_sigtot=nf_sigtot+1
          goto 100
        endif

!       ── Accept/reject ─────────────────────────────────────
!       Born run: weight = sg_born * sib
!       RC run:   weight = sg_born * sig_total (includes RC)
        r_ar = dble(urand(iy))
        if(iborn.eq.1)then
          ar_weight = sg_born_full
        else
          ar_weight = sg_born * sig_total
        endif
        if(r_ar .gt. ar_weight/wmax) goto 100

!       ── Scattered electron acceptance cuts ────────────────────────
!       |p|_e' and theta from exact kinematics
        Eprime_cut = ebeam*(1d0-ys)
        ap1_cut = sqrt(max(0d0, ebeam**2  - aml2))
        ap2_cut = sqrt(max(0d0, Eprime_cut**2 - aml2))
        if(ap2_cut.lt.pe_min .or. ap2_cut.gt.pe_max) goto 100
        if(ap1_cut.gt.0d0 .and. ap2_cut.gt.0d0) then
          costhe_cut = (2d0*ebeam*Eprime_cut - 2d0*aml2 - q2) &
                       / (2d0*ap1_cut*ap2_cut)
          costhe_cut = min(1d0, max(-1d0, costhe_cut))
          the_cut_deg = acos(costhe_cut) * 180d0 / pi
          if(the_cut_deg.lt.the_min .or. the_cut_deg.gt.the_max) goto 100
        endif
!       ──────────────────────────────────────────────────────────────

        nphys = nphys + 1

!       Step 4: Born only or full RC?
        if(iborn.eq.1)then
!         BORN ONLY - skip all RC
          nsoft = nsoft + 1
          call build_4vectors(ebeam,xs,ys,tdif,phirad,k1,ptar,k2,ph,pp)
          call rotz_event(k2,ph,pp,kgam,pi,dble(urand(iy)))
          call decay_rho(ph,pip,pim,iy,ivec)
          call check_acceptance(pp,pip,pim, &
               pp_min,pp_max,thp_min,thp_max, &
               php_min,php_max,thhp_min,thhp_max, &
               phm_min,phm_max,thhm_min,thhm_max, &
               pars(17)%is_set,pars(18)%is_set, &
               pars(19)%is_set,pars(20)%is_set, &
               pars(21)%is_set,pars(22)%is_set, ipass)
          if(ipass.eq.0) goto 100
          ngen = ngen + 1
          wsum  = wsum  + ar_weight
          wsum2 = wsum2 + ar_weight**2
          call write_lund(10,ngen,k1,ptar,k2,ph,pp,kgam, &
                          .false.,ivec,iy,pip,pim,ebeam,ar_weight)
          goto 100
        endif

!       Full RC: Soft or hard?
        prob_hard = sig_hard/sig_total
        r = dble(urand(iy))

        if(r .gt. prob_hard) then
!         NON-RADIATED EVENT
          nsoft = nsoft + 1
          call build_4vectors(ebeam,xs,ys,tdif,phirad,k1,ptar,k2,ph,pp)
          call rotz_event(k2,ph,pp,kgam,pi,dble(urand(iy)))
          call decay_rho(ph,pip,pim,iy,ivec)
          call check_acceptance(pp,pip,pim, &
               pp_min,pp_max,thp_min,thp_max, &
               php_min,php_max,thhp_min,thhp_max, &
               phm_min,phm_max,thhm_min,thhm_max, &
               pars(17)%is_set,pars(18)%is_set, &
               pars(19)%is_set,pars(20)%is_set, &
               pars(21)%is_set,pars(22)%is_set, ipass)
          if(ipass.eq.0) goto 100
          ngen = ngen + 1
          wsum  = wsum  + ar_weight
          wsum2 = wsum2 + ar_weight**2
          call write_lund(10,ngen,k1,ptar,k2,ph,pp,kgam, &
                          .false.,ivec,iy,pip,pim,ebeam,ar_weight)

        else
!         HARD RADIATED EVENT
          nhard = nhard + 1

!         Sample v from the collinear soft kernel (1/v) weighted by
!         the Born suppression at shifted W'^2 = W^2 - v (sig_hard_fix:
!         replaces podinl-based sampling; podinl is the IR-subtracted
!         remainder, not a sampling density). vcut_use = vcut_ir.
          vcut_use = 1d-2
          rvlogmin = log(vcut_use)
          rvlogmax = log(vmax)
          if(rvlogmax.le.rvlogmin)then
            nf_born=nf_born+1
            goto 100
          endif
          call sample_v_soft(vcut_use,vmax,q2,w2,tdif,iy,vrad, &
                           iok_vrad)
          if(iok_vrad.eq.0)then
            nf_born=nf_born+1
            goto 100
          endif
          write(12,'(f12.6)') vrad

!         ISR vs FSR
          Ee1_val = ebeam
          Ee2_val = ebeam*(1d0-ys)
          p_isr   = (Ee2_val**2)/(Ee1_val**2+Ee2_val**2)

          r2 = dble(urand(iy))
          if(r2 .lt. p_isr) then
!           ISR
            nisr = nisr + 1
            omega  = vrad/(2d0*amp)      ! Bug #1 fix: v=2*Mp*omega
            Ee1_rc = Ee1_val - omega
            if(Ee1_rc.le.0d0)then
              nf_thresh=nf_thresh+1
              goto 100
            endif
            call build_4vectors(Ee1_rc,xs,ys,tdif,phirad, &
                                    k1,ptar,k2,ph,pp)
            kgam(1) = omega
            kgam(2) = 0d0
            kgam(3) = 0d0
            kgam(4) = omega
!           Restore original beam in LUND: pp was computed with Ee1_rc so
!           pp = k1_orig + ptar - kgam - k2 - rho (4-momentum conserved)
            k1(1) = Ee1_val
            k1(2) = 0d0
            k1(3) = 0d0
            k1(4) = sqrt(max(0d0,Ee1_val**2-aml2))
          else
!           FSR
            nfsr = nfsr + 1
            omega  = vrad/(2d0*amp)      ! Bug #1 fix: v=2*Mp*omega
            Ee2_rc = Ee2_val - omega
            if(Ee2_rc.le.0d0)then
              nf_thresh=nf_thresh+1
              goto 100
            endif
            call build_4vectors(ebeam,xs,ys,tdif,phirad, &
                                    k1,ptar,k2,ph,pp)
!           pp from build_4vectors is Born-level: pp = k1+ptar-k2_born-rho
!           Do NOT recompute pp with k2_fsr: k2_fsr+kgam=k2_born exactly,
!           so 4-momentum is conserved with pp_Born.
            ascale = Ee2_rc/Ee2_val
            k2(1) = k2(1)*ascale
            k2(2) = k2(2)*ascale
            k2(3) = k2(3)*ascale
            k2(4) = k2(4)*ascale
            kgam(1) = omega
            kgam(2) = (k2(2)/Ee2_rc)*omega
            kgam(3) = (k2(3)/Ee2_rc)*omega
            kgam(4) = (k2(4)/Ee2_rc)*omega
          endif

          call rotz_event(k2,ph,pp,kgam,pi,dble(urand(iy)))
          call decay_rho(ph,pip,pim,iy,ivec)
          call check_acceptance(pp,pip,pim, &
               pp_min,pp_max,thp_min,thp_max, &
               php_min,php_max,thhp_min,thhp_max, &
               phm_min,phm_max,thhm_min,thhm_max, &
               pars(17)%is_set,pars(18)%is_set, &
               pars(19)%is_set,pars(20)%is_set, &
               pars(21)%is_set,pars(22)%is_set, ipass)
          if(ipass.eq.0) goto 100
          ngen = ngen + 1
          wsum  = wsum  + ar_weight
          wsum2 = wsum2 + ar_weight**2
          call write_lund(10,ngen,k1,ptar,k2,ph,pp,kgam, &
                          .true.,ivec,iy,pip,pim,ebeam,ar_weight)
        endif

  100   continue
      enddo

  999 continue

      write(*,'(a)') ''
      write(*,'(a)') '====== Debug failure counts ======'
      write(*,'(a,i8)') ' sample_born fail : ',nf_born
      write(*,'(a,i8)') ' W threshold fail : ',nf_thresh
      write(*,'(a,i8)') ' t kinemat  fail  : ',nf_tkin
      write(*,'(a,i8)') ' vmax<=0    fail  : ',nf_vmax
      write(*,'(a,i8)') ' sib<=0     fail  : ',nf_sib
      write(*,'(a,i8)') ' ssh/xxh    fail  : ',nf_sshxxh
      write(*,'(a,i8)') ' sx-qv      fail  : ',nf_sxqv
      write(*,'(a,i8)') ' sig_total  fail  : ',nf_sigtot
      write(*,'(a)') '=================================='
      write(*,'(a)') '====== Generator Statistics ======'
      anorm = 1d0/max(1_8,ntry)
!     With accept/reject: sigma = wmax * efficiency (using nphys, not ngen)
      axsec = wmax * dble(nphys)/max(1_8,ntry)
      axsec_err = axsec/sqrt(dble(max(1,nphys)))
      write(*,'(a)') '====== Cross Section ======'
      write(*,'(a,g12.4,a)') ' sigma_Born   = ',axsec,' nb'
      write(*,'(a,g12.4,a)') ' stat error   = ',axsec_err,' nb'
      write(*,'(a,g12.4)')   ' efficiency   = ',dble(nphys)/max(1_8,ntry)
      write(*,'(a)') '=================================='
      write(*,'(a,i8)') ' Events written   : ',ngen
      write(*,'(a,i8)') ' Events physics   : ',nphys
      write(*,'(a,i12)') ' Total attempts   : ',ntry
      write(*,'(a,i8)') ' Non-radiated     : ',nsoft
      write(*,'(a,i8)') ' Hard radiated    : ',nhard
      write(*,'(a,i8)') '   ISR events     : ',nisr
      write(*,'(a,i8)') '   FSR events     : ',nfsr
      write(*,'(a,f8.4)') ' Hard fraction    : ', &
                           dble(nhard)/max(1,ngen)
      write(*,'(a)') '=================================='

      write(11,'(a,i8)') 'ngen      = ',ngen
      write(11,'(a,i8)') 'nphys     = ',nphys
      write(11,'(a,i8)') 'nsoft     = ',nsoft
      write(11,'(a,i8)') 'nhard     = ',nhard
      write(11,'(a,f8.4)') 'hard_frac = ',dble(nhard)/max(1,nphys)
      write(11,'(a,g12.4)') 'sigma_nb  = ',axsec
      write(11,'(a,g12.4)') 'sigma_err = ',axsec_err

      close(10); close(11); close(12)
      end


!==============================================================
!     sample_born
!==============================================================
      subroutine sample_born(q2min,q2max,xbmin,xbmax,tmin,tmax, &
                             iy,sg_born,iacc)
      implicit real*8(a-h,o-z)
      real*4 urand
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/sxy/s,x,sx,sxp,q2,w2,aly,anu,sqly,an,tamin,tamax,xs,ys
      common/phi/phirad,tdif,phidif,tq,vmax,ivec
      common/tpar/tslope
      common/cuts/wmin2,ymin,ymax
      common/jpsimodel/pj_alf1,pj_alf2,pj_alf3,pj_nuT,pj_mg2,pj_cR
      integer*4 iy

      iacc = 0

!     Always draw all 4 random numbers before any early returns.
!     This prevents the LCG from falling into a periodic orbit where
!     early-return events (ys>1, W cuts) consume fewer RNs than full
!     events, causing the generator to lock onto a single kinematic point.

      rq2min = log(q2min)
      rq2max = log(q2max)
      r1 = dble(urand(iy))
      r2 = dble(urand(iy))
      r3 = dble(urand(iy))
      r4 = dble(urand(iy))

      q2 = exp(rq2min + r1*(rq2max-rq2min))
      wjacq2 = q2*(rq2max-rq2min)

!      r1 = dble(urand(iy))
!      r2 = dble(urand(iy))
!      r3 = dble(urand(iy))
!      r4 = dble(urand(iy))

!      q2 = q2min + r1*(q2max-q2min)
!      wjacq2 = (q2max-q2min)
      
!     sample xB with log-sampling: tames the 1/xB flux factor,
!     critical for J/psi where xB spans 3+ decades
      rxbmin = log(xbmin)
      rxbmax = log(xbmax)
      xs = exp(rxbmin + r2*(rxbmax-rxbmin))
      wjacxb = xs*(rxbmax-rxbmin)
      if(xs.ge.1d0.or.xs.le.0d0) return

!     derive y from xB and Q2
      ys = q2/(s*xs)
      if(ys.ge.1d0.or.ys.le.0d0) return
!     user analysis cut: ymin <= y <= ymax
      if(ys.lt.ymin .or. ys.gt.ymax) return

!     W cuts: W^2 = M_p^2 + s*y - Q^2
      w2loc = amp2 + s*ys - q2
!     (1) hard physics threshold: W > M_p + M_V
      w2thr = (amp + amv)**2
      if(w2loc.lt.w2thr) return
!     (2) user analysis cut: W > wmin
      if(w2loc.lt.wmin2) return

!     ── t sampling: shape depends on the meson ──
!     ivec=1,2,3 rho/omega/phi : exponential proposal exp(tslope*t)
!     ivec=4 jpsi              : dipole proposal 1/(amg2-t)^4
!     tslope MUST match bt in sigma_T_phi for phi.
      tdmin_u = -tmax
      tdmax_u = -tmin
      if(ivec.le.3) then
!       Rho/omega/phi: exponential exp(tslope*t)
        abslope = tslope
        abt1 = exp(abslope*tdmin_u)
        abt2 = exp(abslope*tdmax_u)
        abst = abt2 - abt1
        if(abst.le.0d0) return
        ranexp = abt1 + r3*abst
        if(ranexp.le.0d0) return
        tdif  = log(ranexp)/abslope
        wjact = abst/abslope/exp(abslope*tdif)
      else
!       Jpsi: dipole 1/(amg2-t)^4 — matches sigma_T_jpsi
        amg2 = pj_mg2      ! mg2 from input file (common /jpsimodel/)
        absA = (amg2 - tdmin_u)**(-3)    ! at t=-tmax (smaller)
        absB = (amg2 - tdmax_u)**(-3)    ! at t=-tmin (larger)
        abst = absB - absA
        if(abst.le.0d0) return
        ranexp = absA + r3*abst
        if(ranexp.le.0d0) return
        tdif  = amg2 - ranexp**(-1d0/3d0)
        wjact = abst * (amg2 - tdif)**4 / 3d0
      endif

      phirad = 2d0*pi*r4

      sg_born = wjacq2 * wjacxb * wjact * 2d0*pi
      iacc = 1
      return
      end


!==============================================================
!     build_4vectors
!     NOTE: all mass variables use 'am' prefix to avoid implicit integer
!     (m,M,k,j,l,n,i are integer by default in Fortran)
!==============================================================
      subroutine build_4vectors(ebeam,xs,ys,tdif,phirad, &
                                k1,ptar,k2,ph,pp)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/sxy/s,x,sx,sxp,q2,w2,aly,anu,sqly,an,tamin,tamax,xs_c,ys_c
      real*8 k1(4),ptar(4),k2(4),ph(4),pp(4)
      real*8 qvec(3),qhat(3),e1v(3),e2v(3),ph3(3),tmp(3)

!     Beam electron (on-shell with mass)
      ap1_e = sqrt(max(0d0,ebeam**2-aml2))
      k1(1) = ebeam
      k1(2) = 0d0
      k1(3) = 0d0
      k1(4) = ap1_e

!     Target proton at rest
      ptar(1) = amp
      ptar(2) = 0d0
      ptar(3) = 0d0
      ptar(4) = 0d0

!     Scattered electron (on-shell with mass)
      aEe2 = ebeam*(1d0-ys)
      if(aEe2.le.0d0) return
      ap2_e = sqrt(max(0d0,aEe2**2-aml2))
      if(ap2_e.le.0d0) return
!     Q2 = 2*E1*E2 - 2*me^2 - 2*|p1|*|p2|*cos(theta)
      aq2_e = s*xs*ys
      acosthe = (2d0*ebeam*aEe2 - 2d0*aml2 - aq2_e) &
                / (2d0*ap1_e*ap2_e)
      if(abs(acosthe).gt.1d0) return
      asinthe = sqrt(max(0d0,1d0-acosthe**2))

      k2(1) = aEe2
      k2(2) = ap2_e*asinthe
      k2(3) = 0d0
      k2(4) = ap2_e*acosthe

!     Virtual photon
      aEq     = k1(1)-k2(1)
      qvec(1) = k1(2)-k2(2)
      qvec(2) = k1(3)-k2(3)
      qvec(3) = k1(4)-k2(4)
      aqmag   = sqrt(qvec(1)**2+qvec(2)**2+qvec(3)**2)

!     Rho meson 4-vector
!     Use 'amrho2' not 'Mrho2' -- M->m is IMPLICIT INTEGER!
      amrho2 = amv**2
      anu_loc = ebeam - aEe2
      aq2_loc = s*xs*ys

!     q.ph = (-Q2 + Mrho2 - t) / 2
      aqdotph = (-aq2_loc + amrho2 - tdif)/2d0

!     Erho from 4-momentum conservation
      aErho = (2d0*amp*anu_loc - aq2_loc + amrho2 - 2d0*aqdotph) &
              /(2d0*amp)
      if(aErho.lt.amv) return

      aphrho = sqrt(max(0d0,aErho**2-amrho2))

!     Angle of rho wrt virtual photon
      if(aqmag*aphrho.le.0d0) return
      acosalpha = (aEq*aErho - aqdotph)/(aqmag*aphrho)
      if(abs(acosalpha).gt.1d0) acosalpha=sign(1d0,acosalpha)
      asinalpha = sqrt(max(0d0,1d0-acosalpha**2))

!     Orthonormal basis around q
      qhat(1) = qvec(1)/aqmag
      qhat(2) = qvec(2)/aqmag
      qhat(3) = qvec(3)/aqmag

      tmp(1) = 0d0; tmp(2) = 1d0; tmp(3) = 0d0
      call cross3(qhat,tmp,e1v)
      aemag = sqrt(e1v(1)**2+e1v(2)**2+e1v(3)**2)
      if(aemag.lt.1d-10) then
        tmp(1)=0d0; tmp(2)=0d0; tmp(3)=1d0
        call cross3(qhat,tmp,e1v)
        aemag=sqrt(e1v(1)**2+e1v(2)**2+e1v(3)**2)
      endif
      e1v(1)=e1v(1)/aemag; e1v(2)=e1v(2)/aemag; e1v(3)=e1v(3)/aemag
      call cross3(qhat,e1v,e2v)

!     Rho 3-momentum direction
      ph3(1) = aphrho*(acosalpha*qhat(1) &
              + asinalpha*cos(phirad)*e1v(1) &
              + asinalpha*sin(phirad)*e2v(1))
      ph3(2) = aphrho*(acosalpha*qhat(2) &
              + asinalpha*cos(phirad)*e1v(2) &
              + asinalpha*sin(phirad)*e2v(2))
      ph3(3) = aphrho*(acosalpha*qhat(3) &
              + asinalpha*cos(phirad)*e1v(3) &
              + asinalpha*sin(phirad)*e2v(3))

!     Put rho on mass shell: E = sqrt(|p|^2 + Mrho^2)
      ph(2) = ph3(1)
      ph(3) = ph3(2)
      ph(4) = ph3(3)
      ph(1) = sqrt(ph(2)**2+ph(3)**2+ph(4)**2+amrho2)

!     Recoil proton from 4-momentum conservation
      pp(1) = k1(1)+ptar(1)-k2(1)-ph(1)
      pp(2) = k1(2)+ptar(2)-k2(2)-ph(2)
      pp(3) = k1(3)+ptar(3)-k2(3)-ph(3)
      pp(4) = k1(4)+ptar(4)-k2(4)-ph(4)

      return
      end


!==============================================================
!     cross3
!==============================================================
      subroutine cross3(a,b,c)
      implicit real*8(a-h,o-z)
      dimension a(3),b(3),c(3)
      c(1) = a(2)*b(3)-a(3)*b(2)
      c(2) = a(3)*b(1)-a(1)*b(3)
      c(3) = a(1)*b(2)-a(2)*b(1)
      return
      end


!==============================================================
!     decay_rho: vector meson -> h+ h- isotropic in rest frame
!       ivec=1 (rho): pi+ pi-  (M_pi = 0.13957 GeV)
!       ivec=3 (phi): K+  K-   (M_K  = 0.49368 GeV)
!       others: pi+ pi- (default)
!==============================================================
      subroutine decay_rho(ph,pip,pim,iy,ivec)
      implicit real*8(a-h,o-z)
      real*4 urand
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      real*8 ph(4),pip(4),pim(4)
      real*8 abeta(3),pip_rf(4),pim_rf(4)
      integer*4 iy, ivec

      if(ivec.eq.3) then
        ampi = 0.493677d0   ! K+/K- mass
      elseif(ivec.eq.4) then
        ampi = 0.000511d0   ! e+/e- mass
      else
        ampi = 0.13957d0    ! pi+/pi- mass
      endif
      amrho = amv

!     Pion momentum in rho rest frame
      appion = sqrt(max(0d0,(amrho/2d0)**2-ampi**2))

!     Random decay direction
      acosth = 2d0*dble(urand(iy))-1d0
      asinth = sqrt(max(0d0,1d0-acosth**2))
      aphid  = 2d0*pi*dble(urand(iy))

      pip_rf(1) = amrho/2d0
      pip_rf(2) = appion*asinth*cos(aphid)
      pip_rf(3) = appion*asinth*sin(aphid)
      pip_rf(4) = appion*acosth

      pim_rf(1) = amrho/2d0
      pim_rf(2) = -pip_rf(2)
      pim_rf(3) = -pip_rf(3)
      pim_rf(4) = -pip_rf(4)

!     Boost to lab (rho velocity)
      aErho   = ph(1)
      abeta(1) = ph(2)/aErho
      abeta(2) = ph(3)/aErho
      abeta(3) = ph(4)/aErho
      abeta2   = abeta(1)**2+abeta(2)**2+abeta(3)**2

!     Use actual rho mass (ph is on mass shell from build_4vectors)
      agamma = aErho/amrho

      call lorentz_boost(pip_rf,abeta,agamma,abeta2,pip)
      call lorentz_boost(pim_rf,abeta,agamma,abeta2,pim)

      return
      end


!==============================================================
!     lorentz_boost: p_out = boost(p_in) by velocity abeta
!     Standard formula: E'=g(E+b.p), p'=p+(g-1)/b2*(b.p)*b+g*E*b
!==============================================================
      subroutine lorentz_boost(p_in,abeta,agamma,abeta2,p_out)
      implicit real*8(a-h,o-z)
      real*8 p_in(4),abeta(3),p_out(4)

      if(abeta2.lt.1d-20) then
        p_out(1)=p_in(1); p_out(2)=p_in(2)
        p_out(3)=p_in(3); p_out(4)=p_in(4)
        return
      endif

      abdotp = abeta(1)*p_in(2)+abeta(2)*p_in(3)+abeta(3)*p_in(4)
      acoef  = (agamma-1d0)/abeta2*abdotp + agamma*p_in(1)

      p_out(1) = agamma*(p_in(1)+abdotp)
      p_out(2) = p_in(2) + acoef*abeta(1)
      p_out(3) = p_in(3) + acoef*abeta(2)
      p_out(4) = p_in(4) + acoef*abeta(3)

      return
      end


!==============================================================
!     check_acceptance -- proton and decay daughter |p|/theta cuts
!     4-vectors: p(1)=E, p(2)=px, p(3)=py, p(4)=pz
!     pip = h+ daughter, pim = h- daughter
!     ipass=1 if all active cuts pass, ipass=0 otherwise
!==============================================================
      subroutine check_acceptance(pp,pip,pim, &
           pp_min,pp_max,thp_min,thp_max, &
           php_min,php_max,thhp_min,thhp_max, &
           phm_min,phm_max,thhm_min,thhm_max, &
           has_pp,has_thp, &
           has_php,has_thhp,has_phm,has_thhm, ipass)
      implicit real*8(a-h,o-z)
      real*8 pp(4),pip(4),pim(4)
      logical has_pp,has_thp,has_php,has_thhp,has_phm,has_thhm
      real*8 amag,atheta,deg
      parameter(deg=57.29577951308232d0)

      ipass = 1

!     Proton |p| cut
      if(has_pp) then
        amag = sqrt(pp(2)**2+pp(3)**2+pp(4)**2)
        if(amag.lt.pp_min .or. amag.gt.pp_max) then
          ipass = 0; return
        endif
      endif

!     Proton theta cut
      if(has_thp) then
        amag = sqrt(pp(2)**2+pp(3)**2+pp(4)**2)
        if(amag.gt.0d0) then
          atheta = acos(min(1d0,max(-1d0,pp(4)/amag)))*deg
        else
          atheta = 0d0
        endif
        if(atheta.lt.thp_min .or. atheta.gt.thp_max) then
          ipass = 0; return
        endif
      endif

!     h+ |p| cut
      if(has_php) then
        amag = sqrt(pip(2)**2+pip(3)**2+pip(4)**2)
        if(amag.lt.php_min .or. amag.gt.php_max) then
          ipass = 0; return
        endif
      endif

!     h+ theta cut
      if(has_thhp) then
        amag = sqrt(pip(2)**2+pip(3)**2+pip(4)**2)
        if(amag.gt.0d0) then
          atheta = acos(min(1d0,max(-1d0,pip(4)/amag)))*deg
        else
          atheta = 0d0
        endif
        if(atheta.lt.thhp_min .or. atheta.gt.thhp_max) then
          ipass = 0; return
        endif
      endif

!     h- |p| cut
      if(has_phm) then
        amag = sqrt(pim(2)**2+pim(3)**2+pim(4)**2)
        if(amag.lt.phm_min .or. amag.gt.phm_max) then
          ipass = 0; return
        endif
      endif

!     h- theta cut
      if(has_thhm) then
        amag = sqrt(pim(2)**2+pim(3)**2+pim(4)**2)
        if(amag.gt.0d0) then
          atheta = acos(min(1d0,max(-1d0,pim(4)/amag)))*deg
        else
          atheta = 0d0
        endif
        if(atheta.lt.thhm_min .or. atheta.gt.thhm_max) then
          ipass = 0; return
        endif
      endif

      return
      end


!==============================================================
!     rotz_event  --  Random rotation around beam (z) axis
!     Rotates k2, ph, pp, kgam by uniform random angle in [0,2pi]
!     k1 (beam along z) is invariant under this rotation
!==============================================================
      subroutine rotz_event(k2,ph,pp,kgam,pi,rnd)
      implicit real*8(a-h,o-z)
      real*8 k2(4),ph(4),pp(4),kgam(4)
      real*8 cphi,sphi,px_tmp
      phi_rot = 2d0*pi*rnd
      cphi = cos(phi_rot)
      sphi = sin(phi_rot)
!     Rotate each 4-vector: (E unchanged, px'=px*c-py*s, py'=px*s+py*c)
      px_tmp = k2(2)*cphi - k2(3)*sphi
      k2(3)  = k2(2)*sphi + k2(3)*cphi
      k2(2)  = px_tmp
      px_tmp = ph(2)*cphi - ph(3)*sphi
      ph(3)  = ph(2)*sphi + ph(3)*cphi
      ph(2)  = px_tmp
      px_tmp = pp(2)*cphi - pp(3)*sphi
      pp(3)  = pp(2)*sphi + pp(3)*cphi
      pp(2)  = px_tmp
      px_tmp   = kgam(2)*cphi - kgam(3)*sphi
      kgam(3)  = kgam(2)*sphi + kgam(3)*cphi
      kgam(2)  = px_tmp
      return
      end

!==============================================================
!     write_lund  --  CLAS12 GEMC LUND format
!     Header (10 mandatory + 5 optional fields):
!       npart  1  1  0  0  11  Ebeam  2212  1  sigma  xB  y  W2  Q2  nu
!     Body (14 fields per particle):
!       idx  lifetime  type  pid  parent  daughter  px py pz E mass  vx vy vz
!     type=1: sent to Geant4;  type=0: not sent (beam e-, vector meson)
!     Particle order: 1=beam e-, 2=scattered e-, 3=recoil p,
!                     4=vector meson, 5=h+, 6=h-, [7=photon]
!==============================================================
      subroutine write_lund(lun,iev,k1,ptar,k2,ph,pp,kgam, &
                            has_photon,ivec,iy,pip,pim,ebeam,sigma_ev)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      real*8 k1(4),ptar(4),k2(4),ph(4),pp(4),kgam(4),pip(4),pim(4)
      logical has_photon
      integer*4 iy
      integer pidp, pidm, pidmeson

!     npart: beam e-, scattered e-, recoil p, meson, h+, h-, [photon]
      npart = 6
      if(has_photon) npart = 7

!     Kinematics for optional header fields
      anu_ev = k1(1) - k2(1)
      q2_ev  = 2d0*(k1(1)*k2(1) - k1(2)*k2(2) &
                   -k1(3)*k2(3) - k1(4)*k2(4)) - 2d0*aml2
      w2_ev  = amp2 + 2d0*amp*anu_ev - q2_ev
      yy_ev  = anu_ev / ebeam
      xb_ev  = q2_ev / (2d0*amp*anu_ev)

!     Header: npart 1 1 0 0 11 Ebeam 2212 1 sigma xB y W2 Q2 nu
      write(lun,'(i5,4i3,i5,f8.3,i6,i3,e14.6,5f9.4)') &
            npart, 1, 1, 0, 0, &
            11, ebeam, 2212, 1, sigma_ev, &
            xb_ev, yy_ev, w2_ev, q2_ev, anu_ev

!     Meson PID and daughter masses/PIDs
      if(ivec.eq.3) then
!       phi -> K+ K-
        pidmeson = 333
        amh  = 0.493677d0
        pidp =  321
        pidm = -321
      elseif(ivec.eq.4) then
!       J/psi -> e+ e-
        pidmeson = 443
        amh  = 0.000511d0
        pidp = -11
        pidm =  11
      else
!       rho/omega -> pi+ pi-
        pidmeson = 113
        amh  = 0.13957d0
        pidp =  211
        pidm = -211
      endif

      aml_ev = sqrt(aml2)   ! electron mass

!     Body: idx  lifetime  type  pid  parent  daughter  px py pz E mass  vx vy vz
!     1: beam electron (type=0, not to Geant4)
      write(lun,'(i4,f6.1,i3,i6,2i3,4f12.6,f12.6,3f8.3)') &
            1, 0.0d0, 0, 11, 0, 0, &
            k1(2),k1(3),k1(4),k1(1), aml_ev, 0d0,0d0,0d0
!     2: scattered electron (type=1)
      write(lun,'(i4,f6.1,i3,i6,2i3,4f12.6,f12.6,3f8.3)') &
            2, 0.0d0, 1, 11, 0, 0, &
            k2(2),k2(3),k2(4),k2(1), aml_ev, 0d0,0d0,0d0
!     3: recoil proton (type=1)
      write(lun,'(i4,f6.1,i3,i6,2i3,4f12.6,f12.6,3f8.3)') &
            3, 0.0d0, 1, 2212, 0, 0, &
            pp(2),pp(3),pp(4),pp(1), amp, 0d0,0d0,0d0
!     4: vector meson (type=0, not to Geant4, daughter=5)
      write(lun,'(i4,f6.1,i3,i6,2i3,4f12.6,f12.6,3f8.3)') &
            4, 0.0d0, 0, pidmeson, 0, 5, &
            ph(2),ph(3),ph(4),ph(1), amv, 0d0,0d0,0d0
!     5: positive daughter (type=1, parent=4)
      write(lun,'(i4,f6.1,i3,i6,2i3,4f12.6,f12.6,3f8.3)') &
            5, 0.0d0, 1, pidp, 4, 0, &
            pip(2),pip(3),pip(4),pip(1), amh, 0d0,0d0,0d0
!     6: negative daughter (type=1, parent=4)
      write(lun,'(i4,f6.1,i3,i6,2i3,4f12.6,f12.6,3f8.3)') &
            6, 0.0d0, 1, pidm, 4, 0, &
            pim(2),pim(3),pim(4),pim(1), amh, 0d0,0d0,0d0
!     7: radiated photon (type=1, RC events only)
      if(has_photon) then
        write(lun,'(i4,f6.1,i3,i6,2i3,4f12.6,f12.6,3f8.3)') &
              7, 0.0d0, 1, 22, 0, 0, &
              kgam(2),kgam(3),kgam(4),kgam(1), 0d0, 0d0,0d0,0d0
      endif

      return
      end


!==============================================================
!     conkin
!==============================================================
      subroutine conkin(snuc)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/sxy/s,x,sx,sxp,q2,w2,aly,anu,sqly,an,tamin,tamax,xs,ys
      s=snuc; x=s*(1d0-ys); sx=s-x; sxp=s+x
      w2=amp2+s-q2-x; aly=sx**2+4d0*amp2*q2; sqly=sqrt(aly)
      anu=sx/ap; an=alpha*ys/(8d0*pi**2)*barn
      tamax=(sx+sqly)/ap2; tamin=-q2/amp2/tamax
      return
      end


!==============================================================
!     sample_bw: sample vector meson mass from Breit-Wigner
!       amv0   -- nominal mass (GeV)
!       gamv   -- full width (GeV)
!       iy     -- random seed
!       amwmax -- kinematic upper limit: sqrt(W2) - Mp
!       amv    -- output: sampled mass (GeV)
!
!     Uses: M^2 = M0^2 + M0*Gamma*tan(theta),
!           theta uniform in [theta_min, theta_max]
!           with M_min = 2*M_pi, M_max = min(amwmax, M0+5*Gamma)
!==============================================================
      subroutine sample_bw(amv0,gamv,ivec,iy,amwmax,amv)
      implicit real*8(a-h,o-z)
      real*4 urand
      integer*4 iy, ivec
!     TEMPORARY: delta function (fixed mass, no BW smearing)
!      amv = amv0; return      ! uncomment if you want delta function
!     Sample M from Breit-Wigner with running width (rho only)
!     For rho (ivec=1): Gamma(M) = Gamma0*(p(M)/p(M0))^3*(M0/M)
!     For others: constant width
      if(ivec.eq.3) then
        ampi = 0.493677d0   ! phi -> K+K-
      elseif(ivec.eq.4) then
        ampi = 0.000511d0   ! J/psi -> e+e-
      else
        ampi = 0.13957d0    ! rho/omega -> pi+pi-
      endif
      ammin = 2d0*ampi
      ammax = min(amwmax, amv0+5d0*gamv)
!     Treat as delta function if width < 1 MeV (catches J/psi: 92.9 keV)
      if(ammax.le.ammin .or. gamv.lt.1d-3) then
        amv = amv0
        return
      endif
!     For rho: integrate using constant-width BW as envelope,
!     then apply accept/reject for running width correction
      p0 = sqrt(max(0d0, amv0**2/4d0 - ampi**2))
      thmin = atan((ammin**2-amv0**2)/(amv0*gamv))
      thmax = atan((ammax**2-amv0**2)/(amv0*gamv))
!     Find the maximum of the accept/reject ratio over the sampling range.
!     The constant-width BW is not always above the running-width BW near M0,
!     so we must normalise the ratio by its maximum to get a valid envelope.
      ratio_max = 1d0
      if(ivec.eq.1 .and. p0.gt.0d0) then
        nscan = 500
        do iscan = 0, nscan
          th_s   = thmin + iscan*(thmax-thmin)/dble(nscan)
          amv2_s = amv0**2 + amv0*gamv*tan(th_s)
          amv_s  = sqrt(max(ammin**2, amv2_s))
          pm_s   = sqrt(max(0d0, amv_s**2/4d0 - ampi**2))
          gr_s   = gamv*(pm_s/p0)**3*(amv0/amv_s)
          bwr_s  = 1d0/((amv2_s-amv0**2)**2 + amv0**2*gr_s**2)
          bwc_s  = 1d0/((amv2_s-amv0**2)**2 + amv0**2*gamv**2)
          r_s    = bwr_s/bwc_s*(gr_s/gamv)*(amv_s/amv0)
          if(r_s .gt. ratio_max) ratio_max = r_s
        enddo
        ratio_max = ratio_max * 1.01d0
      endif
      ibw_try = 0
  10  continue
      ibw_try = ibw_try + 1
      if(ibw_try .gt. 10000) then
!       Give up: return nominal mass so the caller can reject the event
        amv = amv0
        return
      endif
      theta = thmin + dble(urand(iy))*(thmax-thmin)
      amv2  = amv0**2 + amv0*gamv*tan(theta)
      amv   = sqrt(max(ammin**2, amv2))
      if(ivec.eq.1) then
!       Running width for rho: Gamma(M) = Gamma0*(p/p0)^3*(M0/M)
        pm = sqrt(max(0d0, amv**2/4d0 - ampi**2))
        if(p0.gt.0d0) then
          gamrun = gamv*(pm/p0)**3*(amv0/amv)
        else
          gamrun = gamv
        endif
!       Accept/reject: ratio normalised by ratio_max ensures valid envelope
        bw_run  = 1d0/((amv2-amv0**2)**2 + amv0**2*gamrun**2)
        bw_const= 1d0/((amv2-amv0**2)**2 + amv0**2*gamv**2)
        ratio   = bw_run/bw_const*(gamrun/gamv)*(amv/amv0)/ratio_max
        if(dble(urand(iy)).gt.ratio) goto 10
      endif
      end

!==============================================================
!     Akushevich physics routines (unchanged)
!==============================================================
      subroutine setcon(ivec,lepton)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/bwpar/amv0,gamv
      dimension amhad(4),gamhad(4)
      data amhad /0.77526d0,  0.78195d0,  1.019412d0, 3.0969d0  /
      data gamhad/0.1502d0,   0.00849d0,  0.004266d0, 0.0000929d0/
      if(lepton.eq.1)aml2=.261112d-6
      if(lepton.eq.2)aml2=.111637d-1
      pi=3.1415926d0; alpha=.729735d-2; barn=.389379d6
      amv0=amhad(ivec); gamv=gamhad(ivec)
      amv=amv0; amp=.938272d0; amp2=amp**2
      ap=2d0*amp; ap2=2d0*amp2; amc2=amp2
      end

      subroutine bornin(sibor)
!
!     Diehl/gagrho virtual photon flux convention:
!       d3sigma/dxB dQ2 dt = (alpha/2pi) * y^2/(1-eps) * (1-xB)/(xB*Q2)
!                            * (sigma_T + eps*sigma_L)
!     where eps = (1-y-1/4*y^2*gamma^2) / (1-y+1/2*y^2+1/4*y^2*gamma^2)
!           gamma^2 = 4*Mp^2*xB^2/Q^2  (Bjorken gamma)
!
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/sxy/s,x,sx,sxp,q2,w2,aly,anu,sqly,an,tamin,tamax,xs,ys
      common/phi/phirad,tdif,phidif,tq,vmax,ivec
      common/sigsig/sigmat,sigmal
      common/pri/ipri
      ipri=0
      call difflt(q2,w2,tdif,sigmal,sigmat)
      gamma2  = 4d0*amp2*xs**2/q2
      eps_num = 1d0 - ys - 0.25d0*ys**2*gamma2
      eps_den = 1d0 - ys + 0.5d0*ys**2 + 0.25d0*ys**2*gamma2
      if(eps_den.le.0d0.or.eps_num.le.0d0) then; sibor=0d0; return; endif
      eps = eps_num/eps_den
      sibor = alpha*barn/(2d0*pi) * ys**2/(1d0-eps) &
              * (1d0-xs)/(xs*q2) * (sigmat + eps*sigmal)
      end

      subroutine bornin_ak(sibor)
!     Akushevich flux convention (kept for rconv computation in RC run)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/sxy/s,x,sx,sxp,q2,w2,aly,anu,sqly,an,tamin,tamax,xs,ys
      common/phi/phirad,tdif,phidif,tq,vmax,ivec
      common/sigsig/sigmat,sigmal
      aga2  = q2/anu**2
      sibor = 2d0*an/(xs*ys**2)*(ys**2*sigmat &
              +2d0*(1d0-ys-.25d0*ys**2*aga2)*(sigmal+sigmat))
      end

      subroutine difflt(q2,w2,t,sigl,sigt)
!
!     Wrapper: converts W2 -> xB, guards thresholds, then dispatches to
!     the particle-specific cross-section functions:
!       ivec=1  rho    -> sigma_T_rho  / sigma_L_rho
!       ivec=2  omega  -> sigma_T_rho  / sigma_L_rho  (same as rho for now)
!       ivec=3  phi    -> sigma_T_phi  / sigma_L_phi
!       ivec=4  jpsi   -> sigma_T_jpsi / sigma_L_jpsi
!     To change a model, edit only the corresponding pair of functions below.
!
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/phi/phirad,tdif,phidif,tq,vmax,ivec
      real*8 sigma_T_rho,  sigma_L_rho
      real*8 sigma_T_phi,  sigma_L_phi
      real*8 sigma_T_jpsi, sigma_L_jpsi

!     threshold check
      if(w2.lt.(amp+amv)**2)then; sigl=0d0; sigt=0d0; return; endif

!     xB from Q2 and W2
      axb = q2/(w2+q2-amp2)
      if(axb.le.0d0.or.axb.ge.1d0)then; sigl=0d0; sigt=0d0; return; endif

!     dispatch by particle type
      if(ivec.eq.4) then
        sigt = sigma_T_jpsi(q2, axb, t)
        sigl = sigma_L_jpsi(q2, axb, t)
      elseif(ivec.eq.3) then
        sigt = sigma_T_phi(q2, axb, t)
        sigl = sigma_L_phi(q2, axb, t)
      else                       ! ivec=1 (rho) and ivec=2 (omega)
        sigt = sigma_T_rho(q2, axb, t)
        sigl = sigma_L_rho(q2, axb, t)
      endif
      end

!======================================================================
!  Particle-specific cross-section functions.
!  Interface: f(q2, xB, t) — all in GeV units.
!  All physics parameters live inside each function.
!  Edit the appropriate pair when fitting to experimental data.
!======================================================================

!----------------------------------------------------------------------
      real*8 function sigma_T_rho(q2,xB,t)
!
!     Transverse rho0 electroproduction cross section.
!     Model: sigma_T = AT * exp(bT*t) * xB/(1-xB) / Q2^3
!
      implicit real*8(a-h,o-z)
!     --- tunable parameters ---
      parameter( AT = 15.d0  )   ! overall normalization
      parameter( bT =  2.00d0 )  ! t-slope (GeV^-2)
!     --------------------------
      if(xB.le.0d0 .or. xB.ge.1d0 .or. q2.le.0d0)then
        sigma_T_rho = 0d0; return
      endif
      sigma_T_rho = AT * exp(bT*t) * xB / (1d0-xB) / q2**3
      end

!----------------------------------------------------------------------
      real*8 function sigma_L_rho(q2,xB,t)
!
!     Longitudinal rho0 electroproduction cross section.
!     Model: sigma_L = AL * exp(bL*t) * xB / Q2^2
!
      implicit real*8(a-h,o-z)
!     --- tunable parameters ---
      parameter( AL = 25.d0  )   ! overall normalization
      parameter( bL =  4.00d0 )  ! t-slope (GeV^-2)
!     --------------------------
      if(xB.le.0d0 .or. xB.ge.1d0 .or. q2.le.0d0)then
        sigma_L_rho = 0d0; return
      endif
      sigma_L_rho = AL * exp(bL*t) * xB / q2**2
      end

!----------------------------------------------------------------------
      real*8 function sigma_T_phi(q2,xB,t)
!
!     Transverse phi electroproduction cross section  dσ_T/dt  [nb/GeV²]
!     Parameters read from input file via common /phimodel/
!
!       σ_T(W,Q²) = c_T(W) · (m²_φ/(m²_φ+Q²))^ν_T          [nb]
!       c_T(W)    = α₁ · (1 − W²_th/W²)^α₂ · W^α₃          [nb]
!
!       dσ_T/dt   = σ_T · bt · exp(bt · t)                     [nb/GeV²]
!       (t < 0, so bt*t < 0 → exponential fall-off)
!
      implicit real*8(a-h,o-z)
      common/phimodel/phi_alf1,phi_alf2,phi_alf3,phi_nuT,phi_bt,phi_cR
!     --- phi mass and threshold ---
      parameter( amph  = 1.019412d0     )  ! m_φ   [GeV]
      parameter( amph2 = amph**2        )  ! m²_φ  [GeV²]
      parameter( amp_  = 0.938272d0     )  ! M_N   [GeV]
      parameter( amp2  = amp_**2        )  ! M²_N  [GeV²]
      parameter( Wth   = 1.96d0         )  ! W_th  [GeV]
      parameter( Wth2  = Wth**2         )  ! W²_th [GeV²]
!     -----------------------------
      if(xB.le.0d0 .or. xB.ge.1d0 .or. q2.le.0d0)then
        sigma_T_phi = 0d0; return
      endif
!     W² from Q² and xB
      w2 = amp2 + q2*(1d0-xB)/xB
      if(w2.le.Wth2)then; sigma_T_phi = 0d0; return; endif
      w  = sqrt(w2)
!     c_T(W) [nb]
      cT = phi_alf1 * (1d0 - Wth2/w2)**phi_alf2 * w**phi_alf3
!     σ_T(W,Q²) [nb]
      sigT = cT / (1d0 + q2/amph2)**phi_nuT
!     dσ_T/dt with exponential t-dependence [nb/GeV²]
      sigma_T_phi = sigT * phi_bt * exp(phi_bt * t)
      end

!----------------------------------------------------------------------
      real*8 function sigma_L_phi(q2,xB,t)
!
!     Longitudinal phi electroproduction cross section  dσ_L/dt  [nb/GeV²]
!     σ_L/σ_T = c_R · Q²/m²_φ,  exponential t-dependence
!     Parameters read from input file via common /phimodel/
!
      implicit real*8(a-h,o-z)
      common/phimodel/phi_alf1,phi_alf2,phi_alf3,phi_nuT,phi_bt,phi_cR
      parameter( amph  = 1.019412d0     )
      parameter( amph2 = amph**2        )
      parameter( amp_  = 0.938272d0     )
      parameter( amp2  = amp_**2        )
      parameter( Wth   = 1.96d0         )
      parameter( Wth2  = Wth**2         )
!     -----------------------------
      if(xB.le.0d0 .or. xB.ge.1d0 .or. q2.le.0d0)then
        sigma_L_phi = 0d0; return
      endif
      w2 = amp2 + q2*(1d0-xB)/xB
      if(w2.le.Wth2)then; sigma_L_phi = 0d0; return; endif
      w  = sqrt(w2)
      cT   = phi_alf1 * (1d0 - Wth2/w2)**phi_alf2 * w**phi_alf3
      sigT = cT / (1d0 + q2/amph2)**phi_nuT
      dsdt = sigT * phi_bt * exp(phi_bt * t)
      sigma_L_phi = phi_cR * (q2/amph2) * dsdt
      end

!----------------------------------------------------------------------
      real*8 function sigma_T_jpsi(q2,xB,t)
!
!     Transverse J/psi electroproduction cross section  dσ_T/dt  [nb/GeV²]
!     Same parametrisation as phi with m_φ → m_J/psi.
!     Matched to Bhawani/lAger formulas.
!
      implicit real*8(a-h,o-z)
      common/jpsimodel/pj_alf1,pj_alf2,pj_alf3,pj_nuT,pj_mg2,pj_cR
!     --- J/psi mass and threshold ---
      parameter( amjp  = 3.0969d0       )  ! m_J/psi [GeV]
      parameter( amjp2 = amjp**2        )  ! m²_J/psi [GeV²]
      parameter( amp_  = 0.938272d0     )  ! M_N [GeV]
      parameter( amp2  = amp_**2        )  ! M²_N
      parameter( Wth2  = (amp_+amjp)**2 )  ! (M_N+m_J/psi)²
!     parameters from input file via common /jpsimodel/
!     -----------------------------
      if(xB.le.0d0 .or. xB.ge.1d0 .or. q2.le.0d0)then
        sigma_T_jpsi = 0d0; return
      endif
      w2 = amp2 + q2*(1d0-xB)/xB
      if(w2.le.Wth2)then; sigma_T_jpsi = 0d0; return; endif
      w  = sqrt(w2)
!     kinematic |t|_min (positive); t argument is negative
      ecm_i  = (w2 + q2 + amp2)/(2d0*w)
      pcm_i  = sqrt(max(0d0, ecm_i**2 - amp2))
      ecm_f  = (w2 + amjp2 - amp2)/(2d0*w)
      pcm_f  = sqrt(max(0d0, ecm_f**2 - amjp2))
      tmin_k = (pcm_i - pcm_f)**2 - ( (q2 + amjp2)/(2d0*w) )**2
      tmax_k = (pcm_i + pcm_f)**2 - ( (q2 + amjp2)/(2d0*w) )**2
      if(-t.lt.tmin_k .or. -t.gt.tmax_k)then
        sigma_T_jpsi = 0d0; return
      endif
      cT   = pj_alf1 * (1d0 - Wth2/w2)**pj_alf2 * w**pj_alf3
      sigT = cT / (1d0 + q2/amjp2)**pj_nuT
      sigma_T_jpsi = sigT * 3d0*(pj_mg2 + tmin_k)**3 / (pj_mg2 - t)**4
      end

!----------------------------------------------------------------------
      real*8 function sigma_L_jpsi(q2,xB,t)
!
!     Longitudinal J/psi electroproduction cross section  dσ_L/dt  [nb/GeV²]
!     Matched to Bhawani/lAger:  σ_L/σ_T = c_R · Q²/m²_J/psi
!
      implicit real*8(a-h,o-z)
      common/jpsimodel/pj_alf1,pj_alf2,pj_alf3,pj_nuT,pj_mg2,pj_cR
      parameter( amjp  = 3.0969d0       )
      parameter( amjp2 = amjp**2        )
      parameter( amp_  = 0.938272d0     )
      parameter( amp2  = amp_**2        )
      parameter( Wth2  = (amp_+amjp)**2 )
!     parameters from input file via common /jpsimodel/
!     (same sigma_T as sigma_T_jpsi, times cR*Q2/m2jp)
!     -----------------------------
      if(xB.le.0d0 .or. xB.ge.1d0 .or. q2.le.0d0)then
        sigma_L_jpsi = 0d0; return
      endif
      w2 = amp2 + q2*(1d0-xB)/xB
      if(w2.le.Wth2)then; sigma_L_jpsi = 0d0; return; endif
      w  = sqrt(w2)
!     kinematic |t|_min (positive); t argument is negative
      ecm_i  = (w2 + q2 + amp2)/(2d0*w)
      pcm_i  = sqrt(max(0d0, ecm_i**2 - amp2))
      ecm_f  = (w2 + amjp2 - amp2)/(2d0*w)
      pcm_f  = sqrt(max(0d0, ecm_f**2 - amjp2))
      tmin_k = (pcm_i - pcm_f)**2 - ( (q2 + amjp2)/(2d0*w) )**2
      tmax_k = (pcm_i + pcm_f)**2 - ( (q2 + amjp2)/(2d0*w) )**2
      if(-t.lt.tmin_k .or. -t.gt.tmax_k)then
        sigma_L_jpsi = 0d0; return
      endif
      cT   = pj_alf1 * (1d0 - Wth2/w2)**pj_alf2 * w**pj_alf3
      sigT = cT / (1d0 + q2/amjp2)**pj_nuT
      dsdt = sigT * 3d0*(pj_mg2 + tmin_k)**3 / (pj_mg2 - t)**4
      sigma_L_jpsi = pj_cR * (q2/amjp2) * dsdt
      end

      double precision function vacpol(t)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      dimension am2(3)
      data am2/.26110d-6,.111637d-1,3.18301d0/
      asuml=0d0
      do 10 i=1,3
        aa2=2d0*am2(i)
        asqlmi=sqrt(t*t+2d0*aa2*t)
        aallmi=log((asqlmi+t)/(asqlmi-t))/asqlmi
  10  asuml=asuml+2d0*(t+aa2)*aallmi/3d0-10d0/9d0 &
                  +4d0*aa2*(1d0-aa2*aallmi)/3d0/t
      if(t.lt.1d0)then
        aaaa=-1.345d-9; abbb=-2.302d-3; accc=4.091d0
      elseif(t.lt.64d0)then
        aaaa=-1.512d-3; abbb=-2.822d-3; accc=1.218d0
      else
        aaaa=-1.1344d-3; abbb=-3.0680d-3; accc=9.9992d-1
      endif
      asumh=-(aaaa+abbb*log(1d0+accc*t))*2d0*pi/alpha
      vacpol=asuml+asumh
      end

      double precision function fspens(x)
      implicit real*8(a-h,o-z)
      af=0d0; aa=1d0; aan=0d0; atch=1d-16
  1   aan=aan+1d0; aa=aa*x; ab=aa/aan**2; af=af+ab
      if(ab-atch)2,2,1
  2   fspens=af
      return
      end

      double precision function fspen(x)
      implicit real*8(a-h,o-z)
      data af1/1.644934d0/
      if(x)8,1,1
  1   if(x-.5d0)2,2,3
  2   fspen=fspens(x); return
  3   if(x-1d0)4,4,5
  4   fspen=af1-log(x)*log(1d0-x+1d-10)-fspens(1d0-x); return
  5   if(x-2d0)6,6,7
  6   fspen=af1-.5d0*log(x)*log((x-1d0)**2/x)+fspens(1d0-1d0/x)
      return
  7   fspen=2d0*af1-.5d0*log(x)**2-fspens(1d0/x); return
  8   if(x+1d0)10,9,9
  9   fspen=-.5d0*log(1d0-x)**2-fspens(x/(x-1d0)); return
 10   fspen=-.5d0*log(1d0-x)*log(x**2/(1d0-x)) &
            -af1+fspens(1d0/(1d0-x))
      return
      end


!==============================================================
!  Exact RC subroutines from Akushevich idiffrad.f
!  qqt -> qqtphi -> rv2ln -> podinl
!  These compute sigma_hard exactly (eq.15 of the paper)
!==============================================================

      subroutine qqt(tai)
      implicit real*8(a-h,o-z)
      external qqtphi
      real*8 qqtphi
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
!     Use 32-point Gauss quadrature for smooth phi integrand
      call dqg32(0d0,2.d0*pi,qqtphi,tai)
      tai=tai/2.d0
      end

      double precision function qqtphi(phi)
      implicit real*8(a-h,o-z)
      external rv2ln
      real*8 rv2ln
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/sxy/s,x,sx,sxp,q2,w2,aly,anu,sqly,an,tamin,tamax,xs,ys
      common/phi/phirad,tdif,phidif,tq,vmax,ivec
      common/ivv/vcurr,cutv
      dimension tlm(4)
      phirad=phi
      tlm(1)=log(xs+tamin)
      tlm(4)=log(xs+tamax)
      tlm(2)=log(xs-q2/s)
      tlm(3)=log(xs+q2/x)
      res=0d0
      do ii=1,3
        ep=1d-10
        call simptx(tlm(ii)+ep,tlm(ii+1)-ep,10,1d-2,rv2ln,re)
        tai=an*alpha/pi*re
        res=res+tai
      enddo
      qqtphi=res
      end

      double precision function rv2ln(taln)
      implicit real*8(a-h,o-z)
      external podinl
      real*8 podinl
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/sxy/s,x,sx,sxp,q2,w2,aly,anu,sqly,an,tamin,tamax,xs,ys
      common/amf2/taa,atm(8,6),sfm0(8)
      common/phi/phirad,tdif,phidif,tq,vmax,ivec
      common/ivv/vcurr,cutv
      ta=exp(taln)-q2/sx
      taa=ta
      sqrtmb=sqrt((ta-tamin)*(tamax-ta)*(s*x*q2-q2**2*amp2-aml2*aly))
      z1=(q2*sxp+ta*(s*sx+ap2*q2)-ap*cos(phirad)*sqrtmb)/aly
      z2=(q2*sxp+ta*(x*sx-ap2*q2)-ap*cos(phirad)*sqrtmb)/aly
      if(abs(z1).lt.1d-10.or.abs(z2).lt.1d-10)then
        rv2ln=0d0; return
      endif
      abb=1.d0/sqly/pi
      bi12=abb/(z1*z2)
      if(abs(bi12).gt.1d20)then; rv2ln=0d0; return; endif
      bi1pi2=abb/z2+abb/z1
      abis=abb/z2**2+abb/z1**2
      abir=abb/z2**2-abb/z1**2
      ahi2=aml2*abis-q2*bi12
      atm(1,1)=4.d0*q2*ahi2
      atm(1,2)=4.d0*ahi2*ta
      atm(1,3)=-2.d0*(2.d0*abb+bi12*ta**2)
      atm(2,1)=2d0*(s*x-amp2*q2)*ahi2/amp2
      atm(2,2)=(2.d0*aml2*abir*sxp-4.d0*amp2*ahi2*ta-bi12*sxp**2*ta+ &
                bi1pi2*sxp*sx+2.d0*ahi2*sx)/(2.d0*amp2)
      atm(2,3)=(2.d0*(2.d0*abb+bi12*ta**2)*amp2-bi12* &
                sx*ta-bi1pi2*sxp)/(2.d0*amp2)
!     Fixed 40-point log-uniform integration: avoids adaptive subdivision
!     singularity at v->0 (podinl ~ 1/v) that causes infinite halving loops.
!     Change of variable u=log(v): integral = int podinl(exp(u))*exp(u) du
      vmin=1d-4
      rvlogmin=log(vmin); rvlogmax=log(vmax)
      du=(rvlogmax-rvlogmin)/39d0
      res=0d0
      do ilogv=1,40
        vv=exp(rvlogmin+(ilogv-1)*du)
        ww=podinl(vv)*vv
        if(ilogv.eq.1.or.ilogv.eq.40) ww=ww*0.5d0
        res=res+ww
      enddo
      res=res*du
      rv2ln=res*(q2/sx+ta)
      end

      double precision function podinl(v)
      implicit real*8(a-h,o-z)
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/sxy/s,x,sx,sxp,q2,w2,aly,anu,sqly,an,tamin,tamax,xs,ys
      common/amf2/ta,atm(8,6),sfm0(8)
      common/phi/phirad,tdif,phidif,tq,vmax,ivec
      common/sigsig/sigmat0,sigmal0
      common/ivv/vcurr,cutv
      dimension sfm(8)
      sxtm=sx+tdif-v
      aak=((2d0*q2+ta*sx)*sxtm-(sx-2d0*ta*amp2)*tq)/aly/2d0
      sqbbk1=sqrt(max(0d0,q2*sxtm**2-sxtm*sx*tq-amp2*tq**2-amv**2*aly))
      sqbbk2=amp*sqrt(max(0d0,(tamax-ta)*(ta-tamin)))
      abbk=sqbbk1*sqbbk2/aly
      d2kvir=2d0*(aak+abbk*cos(phirad-phidif))
      factor=1d0+ta-d2kvir
      if(abs(factor).lt.1d-10)then; podinl=0d0; return; endif
      r=v/factor
      tldq2=q2+r*ta
      tldw2=w2-r*(1.d0+ta)
      tldtd=tdif-r*(ta-d2kvir)
      call difflt(tldq2,tldw2,tldtd,sigmal,sigmat)
      sfm(1)=(sx-r)*sigmat
      sfm(2)=2.d0*ap2/(sx-r)*tldq2*(sigmat+sigmal)
      sfm0(1)=sx*sigmat0
      sfm0(2)=2.d0*ap2/sx*q2*(sigmat0+sigmal0)
      podinl=0.d0
      do isf=1,2
        do irr=1,3
          app=sfm(isf)
          if(irr.eq.1)app=app-sfm0(isf)*(1.d0+r*ta/q2)**2
          pres=app*r**(irr-2)/(q2+r*ta)**2/2.d0
          podinl=podinl-atm(isf,irr)*pres
        enddo
      enddo
      podinl=podinl/factor
      end

!--------------------------------------------------------------
!  sample_v_soft (sig_hard_fix): sample v from the collinear soft
!  kernel dN/dv ~ (1/v) * R(v), R(v) = sigma(W^2-v)/sigma(W^2) the
!  Born threshold suppression at the radiatively shifted W'.
!  Log-uniform proposal in [vmin,vmax]; accept with prob R(v)<=Rmax.
!--------------------------------------------------------------
      subroutine sample_v_soft(vmin_in,vmax_in,q2_in,w2_in,t_in, &
                               iy,vrad_out,iok)
      implicit real*8(a-h,o-z)
      real*4 urand
      integer*4 iy
      call difflt(q2_in,w2_in,t_in,sl0,st0)
      den = st0 + sl0
      if(den.le.0d0)then; iok=0; return; endif
      rvlogmin = log(vmin_in)
      rvlogmax = log(vmax_in)
      do itry=1,1000
        vv = exp(rvlogmin + (rvlogmax-rvlogmin)*dble(urand(iy)))
        call difflt(q2_in,w2_in-vv,t_in,sl,st)
        ratio = (st+sl)/den
        if(dble(urand(iy)).lt.ratio)then
          vrad_out = vv
          iok = 1
          return
        endif
      enddo
      iok = 0
      end

!--------------------------------------------------------------
!  sample_vrad: sample vrad by accept/reject weighted by podinl.
!  Bug #2 fix: replaces log-flat sampling with physical distribution.
!  podinl is declared external here, NOT in the main program,
!  to avoid interfering with simpux/rv2ln function-argument passing.
!--------------------------------------------------------------
      subroutine sample_vrad(vmin_in,vmax_in,tamin_in,tamax_in, &
                             phidif_in,iy,vrad_out,iok)
      implicit real*8(a-h,o-z)
      external podinl
      real*8 podinl
      real*4 urand
      integer*4 iy
      common/cmp/pi,alpha,amp,amp2,ap,ap2,aml2,amc2,amv,barn
      common/sxy/s,x,sx,sxp,q2,w2,aly,anu,sqly,an,tamin,tamax,xs,ys
      common/amf2/taa,atm(8,6),sfm0(8)
      common/phi/phirad,tdif,phidif,tq,vmax,ivec

!     Set up atm at near-collinear ta with phirad = phidif
      phirad = phidif_in
      ta_c = tamin_in + 1d-3*(tamax_in-tamin_in)
      if(ta_c.gt.tamax_in) ta_c = 0.5d0*(tamin_in+tamax_in)

!     Populate atm for ta_c (mirrors rv2ln logic)
      ta   = ta_c
      taa  = ta
      sqrtmb = sqrt(max(0d0,(ta-tamin_in)*(tamax_in-ta) &
                    *(s*x*q2-q2**2*amp2-aml2*aly)))
      z1 = (q2*sxp+ta*(s*sx+ap2*q2)-ap*cos(phirad)*sqrtmb)/aly
      z2 = (q2*sxp+ta*(x*sx-ap2*q2)-ap*cos(phirad)*sqrtmb)/aly
      if(abs(z1).lt.1d-30.or.abs(z2).lt.1d-30)then
        iok = 0
        return
      endif
      abb    = 1.d0/sqly/pi
      bi12   = abb/(z1*z2)
      bi1pi2 = abb/z2+abb/z1
      abis   = abb/z2**2+abb/z1**2
      abir   = abb/z2**2-abb/z1**2
      ahi2   = aml2*abis-q2*bi12
      atm(1,1) = 4.d0*q2*ahi2
      atm(1,2) = 4.d0*ahi2*ta
      atm(1,3) = -2.d0*(2.d0*abb+bi12*ta**2)
      atm(2,1) = 2d0*(s*x-amp2*q2)*ahi2/amp2
      atm(2,2) = (2.d0*aml2*abir*sxp-4.d0*amp2*ahi2*ta &
                  -bi12*sxp**2*ta+bi1pi2*sxp*sx &
                  +2.d0*ahi2*sx)/(2.d0*amp2)
      atm(2,3) = (2.d0*(2.d0*abb+bi12*ta**2)*amp2 &
                  -bi12*sx*ta-bi1pi2*sxp)/(2.d0*amp2)

!     Scan podinl over v to find local maximum
      rvlogmin = log(vmin_in)
      rvlogmax = log(vmax_in)
      pmax = 0d0
      do iscan = 1, 50
        frac  = dble(iscan-1)/49d0
        vscan = exp(rvlogmin + frac*(rvlogmax-rvlogmin))
        pval  = podinl(vscan)
        if(pval.gt.pmax) pmax = pval
      enddo
      pmax = pmax*2d0   ! safety margin

      if(pmax.le.0d0)then
        iok = 0
        return
      endif

!     Accept/reject loop
      vrad_out = exp(rvlogmin+dble(urand(iy))*(rvlogmax-rvlogmin))
      do itry = 1, 10000
        vtry = exp(rvlogmin+dble(urand(iy))*(rvlogmax-rvlogmin))
        pval = max(0d0, podinl(vtry))
        if(dble(urand(iy)).lt.pval/pmax)then
          vrad_out = vtry
          exit
        endif
      enddo
      iok = 1
      end

      subroutine simpsx(a,b,np,ep,func,res)
      implicit real*8(a-h,o-z)
      external func
      step=(b-a)/np
      call simps(a,b,step,ep,1d-18,func,ra,res,r2,r3)
      end

      subroutine simptx(a,b,np,ep,func,res)
      implicit real*8(a-h,o-z)
      external func
      step=(b-a)/np
      call simpt(a,b,step,ep,1d-18,func,ra,res,r2,r3)
      end

      subroutine simpux(a,b,np,ep,func,res)
      implicit real*8(a-h,o-z)
      external func
      step=(b-a)/np
      call simpu(a,b,step,ep,1d-18,func,ra,res,r2,r3)
      end

      subroutine dqg32(xl,xu,fct,y)
      double precision xl,xu,y,a,b,c,fct
      a=.5d0*(xu+xl); b=xu-xl
      c=.49863193092474078d0*b
      y=.35093050047350483d-2*(fct(a+c)+fct(a-c))
      c=.49280575577263417d0*b
      y=y+.8137197365452835d-2*(fct(a+c)+fct(a-c))
      c=.48238112779375322d0*b
      y=y+.12696032654631030d-1*(fct(a+c)+fct(a-c))
      c=.46745303796886984d0*b
      y=y+.17136931456510717d-1*(fct(a+c)+fct(a-c))
      c=.44816057788302606d0*b
      y=y+.21417949011113340d-1*(fct(a+c)+fct(a-c))
      c=.42468380686628499d0*b
      y=y+.25499029631188088d-1*(fct(a+c)+fct(a-c))
      c=.39724189798397120d0*b
      y=y+.29342046739267774d-1*(fct(a+c)+fct(a-c))
      c=.36609105937014484d0*b
      y=y+.32911111388180923d-1*(fct(a+c)+fct(a-c))
      c=.33152213346510760d0*b
      y=y+.36172897054424253d-1*(fct(a+c)+fct(a-c))
      c=.29385787862038116d0*b
      y=y+.39096947893535153d-1*(fct(a+c)+fct(a-c))
      c=.25344995446611470d0*b
      y=y+.41655962113473378d-1*(fct(a+c)+fct(a-c))
      c=.21067563806531767d0*b
      y=y+.43826046502201906d-1*(fct(a+c)+fct(a-c))
      c=.16593430114106382d0*b
      y=y+.45586939347881942d-1*(fct(a+c)+fct(a-c))
      c=.11964368112606854d0*b
      y=y+.46922199540402283d-1*(fct(a+c)+fct(a-c))
      c=.7223598079139825d-1*b
      y=y+.47819360039637430d-1*(fct(a+c)+fct(a-c))
      c=.24153832843869158d-1*b
      y=b*(y+.48270044257363900d-1*(fct(a+c)+fct(a-c)))
      end

      subroutine simps(a1,b1,h1,reps1,aeps1,funct,x,ai,aih,aiabs)
      implicit real*8(a-h,o-z)
      dimension f(7),p(5)
      integer nhalv
      h=dsign(h1,b1-a1); s=dsign(1.d0,h)
      a=a1; b=b1; ai=0.d0; aih=0.d0; aiabs=0.d0
      p(2)=4.d0; p(4)=4.d0; p(3)=2.d0; p(5)=1.d0
      if(b-a) 1,2,1
    1 reps=dabs(reps1); aeps=dabs(aeps1)
      do 3 k=1,7
  3   f(k)=10.d16
      x=a; c=0.d0; f(1)=funct(x)/3.d0; nhalv=0
    4 x0=x
      if((x0+4.*h-b)*s) 5,5,6
    6 h=(b-x0)/4.; if(h) 7,2,7
    7 do 8 k=2,7
  8   f(k)=10.d16
      c=1.d0
    5 di2=f(1); di3=dabs(f(1))
      do 9 k=2,5
      x=x+h
      if((x-b)*s) 23,24,24
   24 x=b
   23 if(f(k)-10.d16) 10,11,10
   11 f(k)=funct(x)/3.
   10 di2=di2+p(k)*f(k)
    9 di3=di3+p(k)*abs(f(k))
      di1=(f(1)+4.*f(3)+f(5))*2.*h
      di2=di2*h; di3=di3*h
      if(reps) 12,13,12
   13 if(aeps) 12,14,12
   12 eps=dabs((aiabs+di3)*reps)
      if(eps-aeps) 15,16,16
   15 eps=aeps
   16 delta=dabs(di2-di1)
      if(delta-eps) 20,21,21
   20 if(delta-eps/8.) 17,14,14
   17 h=2.*h; f(1)=f(5); f(2)=f(6); f(3)=f(7)
      do 19 k=4,7
  19  f(k)=10.d16
      go to 18
   14 f(1)=f(5); f(3)=f(6); f(5)=f(7)
      f(2)=10.d16; f(4)=10.d16; f(6)=10.d16; f(7)=10.d16
   18 di1=di2+(di2-di1)/15.
      ai=ai+di1; aih=aih+di2; aiabs=aiabs+di3
      go to 22
   21 nhalv=nhalv+1
      if(nhalv.gt.10) go to 14
      h=h/2.; f(7)=f(5); f(6)=f(4); f(5)=f(3)
      f(3)=f(2); f(2)=10.d16; f(4)=10.d16
      x=x0; c=0.d0
      go to 5
   22 if(c) 2,4,2
    2 return
      end

      subroutine simpt(a1,b1,h1,reps1,aeps1,funct,x,ai,aih,aiabs)
      implicit real*8(a-h,o-z)
      dimension f(7),p(5)
      integer nhalv
      h=dsign(h1,b1-a1); s=dsign(1.d0,h)
      a=a1; b=b1; ai=0.d0; aih=0.d0; aiabs=0.d0
      p(2)=4.d0; p(4)=4.d0; p(3)=2.d0; p(5)=1.d0
      if(b-a) 1,2,1
    1 reps=dabs(reps1); aeps=dabs(aeps1)
      do 3 k=1,7
  3   f(k)=10.d16
      x=a; c=0.d0; f(1)=funct(x)/3.; nhalv=0
    4 x0=x
      if((x0+4.*h-b)*s) 5,5,6
    6 h=(b-x0)/4.; if(h) 7,2,7
    7 do 8 k=2,7
  8   f(k)=10.d16
      c=1.d0
    5 di2=f(1); di3=dabs(f(1))
      do 9 k=2,5
      x=x+h
      if((x-b)*s) 23,24,24
   24 x=b
   23 if(f(k)-10.d16) 10,11,10
   11 f(k)=funct(x)/3.
   10 di2=di2+p(k)*f(k)
    9 di3=di3+p(k)*abs(f(k))
      di1=(f(1)+4.*f(3)+f(5))*2.*h
      di2=di2*h; di3=di3*h
      if(reps) 12,13,12
   13 if(aeps) 12,14,12
   12 eps=dabs((aiabs+di3)*reps)
      if(eps-aeps) 15,16,16
   15 eps=aeps
   16 delta=dabs(di2-di1)
      if(delta-eps) 20,21,21
   20 if(delta-eps/8.) 17,14,14
   17 h=2.*h; f(1)=f(5); f(2)=f(6); f(3)=f(7)
      do 19 k=4,7
  19  f(k)=10.d16
      go to 18
   14 f(1)=f(5); f(3)=f(6); f(5)=f(7)
      f(2)=10.d16; f(4)=10.d16; f(6)=10.d16; f(7)=10.d16
   18 di1=di2+(di2-di1)/15.
      ai=ai+di1; aih=aih+di2; aiabs=aiabs+di3
      go to 22
   21 nhalv=nhalv+1
      if(nhalv.gt.10) go to 14
      h=h/2.; f(7)=f(5); f(6)=f(4); f(5)=f(3)
      f(3)=f(2); f(2)=10.d16; f(4)=10.d16
      x=x0; c=0.d0
      go to 5
   22 if(c) 2,4,2
    2 return
      end

      subroutine simpu(a1,b1,h1,reps1,aeps1,funct,x,ai,aih,aiabs)
      implicit real*8(a-h,o-z)
      dimension f(7),p(5)
      integer nhalv
      h=dsign(h1,b1-a1); s=dsign(1.d0,h)
      a=a1; b=b1; ai=0.d0; aih=0.d0; aiabs=0.d0
      p(2)=4.d0; p(4)=4.d0; p(3)=2.d0; p(5)=1.d0
      if(b-a) 1,2,1
    1 reps=dabs(reps1); aeps=dabs(aeps1)
      do 3 k=1,7
  3   f(k)=10.d16
      x=a; c=0.d0; f(1)=funct(x)/3.; nhalv=0
    4 x0=x
      if((x0+4.*h-b)*s) 5,5,6
    6 h=(b-x0)/4.; if(h) 7,2,7
    7 do 8 k=2,7
  8   f(k)=10.d16
      c=1.d0
    5 di2=f(1); di3=dabs(f(1))
      do 9 k=2,5
      x=x+h
      if((x-b)*s) 23,24,24
   24 x=b
   23 if(f(k)-10.d16) 10,11,10
   11 f(k)=funct(x)/3.
   10 di2=di2+p(k)*f(k)
    9 di3=di3+p(k)*abs(f(k))
      di1=(f(1)+4.*f(3)+f(5))*2.*h
      di2=di2*h; di3=di3*h
      if(reps) 12,13,12
   13 if(aeps) 12,14,12
   12 eps=dabs((aiabs+di3)*reps)
      if(eps-aeps) 15,16,16
   15 eps=aeps
   16 delta=dabs(di2-di1)
      if(delta-eps) 20,21,21
   20 if(delta-eps/8.) 17,14,14
   17 h=2.*h; f(1)=f(5); f(2)=f(6); f(3)=f(7)
      do 19 k=4,7
  19  f(k)=10.d16
      go to 18
   14 f(1)=f(5); f(3)=f(6); f(5)=f(7)
      f(2)=10.d16; f(4)=10.d16; f(6)=10.d16; f(7)=10.d16
   18 di1=di2+(di2-di1)/15.
      ai=ai+di1; aih=aih+di2; aiabs=aiabs+di3
      go to 22
   21 nhalv=nhalv+1
      if(nhalv.gt.10) go to 14
      h=h/2.; f(7)=f(5); f(6)=f(4); f(5)=f(3)
      f(3)=f(2); f(2)=10.d16; f(4)=10.d16
      x=x0; c=0.d0
      go to 5
   22 if(c) 2,4,2
    2 return
      end

!==============================================================
!     SPIN DENSITY MATRIX ELEMENT (SDME) ANGULAR DISTRIBUTION
!     Schilling-Wolf formalism for electroproduction
!     Reference: K. Schilling, P. Seyboth, G. Wolf,
!       Nucl. Phys. B 15, 397 (1970)
!     HERMES: Eur. Phys. J. C 62, 659 (2009)
!
!     15-element SDME array r(15) ordering:
!       r(1)  = r04_00        unpolarized (sigma_T + eps*sigma_L)
!       r(2)  = Re r04_10
!       r(3)  = r04_1-1
!       r(4)  = r1_11         transverse linear pol. (cos 2Phi)
!       r(5)  = r1_00
!       r(6)  = Re r1_10
!       r(7)  = r1_1-1
!       r(8)  = Im r2_10      transverse linear pol. (sin 2Phi)
!       r(9)  = Im r2_1-1
!       r(10) = r5_11         T-L interference (cos Phi)
!       r(11) = r5_00
!       r(12) = Re r5_10
!       r(13) = r5_1-1
!       r(14) = Im r6_10      T-L interference (sin Phi)
!       r(15) = Im r6_1-1
!
!     Angles (helicity frame):
!       costh  = cos(Theta) of h+ in meson rest frame
!       phid   = phi of h+ decay plane vs production plane
!       Phi    = azimuth of production plane vs lepton plane
!       eps    = virtual photon polarization parameter
!==============================================================


!==============================================================
!     sdme_W: full angular distribution W(costh, phid, Phi)
!     Returns the normalized W value (should be >= 0)
!==============================================================
      function sdme_W(costh, phid, aPhi, eps, r)
      implicit real*8(a-h,o-z)
      real*8 r(15), sdme_W

      sinth2 = 1d0 - costh*costh
      sinth  = sqrt(max(0d0, sinth2))
      sin2th = 2d0*sinth*costh       ! sin(2*Theta)
      cos2ph = cos(2d0*phid)
      sin2ph = sin(2d0*phid)
      cosph  = cos(phid)
      sinph  = sin(phid)
      cos2Phi = cos(2d0*aPhi)
      sin2Phi = sin(2d0*aPhi)
      cosPhi  = cos(aPhi)
      sinPhi  = sin(aPhi)

      sq2 = sqrt(2d0)

!     ── W^0: unpolarized part (alpha=04) ──────────────────
      w0 = 0.5d0*(1d0 - r(1)) &
         + 0.5d0*(3d0*r(1) - 1d0)*costh*costh &
         - sq2*r(2)*sin2th*cosph &
         - r(3)*sinth2*cos2ph

!     ── W^1: transverse lin. pol. (alpha=1, couples to cos2Phi)
      w1 = r(4)*sinth2 + r(5)*costh*costh &
         - sq2*r(6)*sin2th*cosph &
         - r(7)*sinth2*cos2ph

!     ── W^2: transverse lin. pol. (alpha=2, couples to sin2Phi)
      w2 = sq2*r(8)*sin2th*sinph &
         + r(9)*sinth2*sin2ph

!     ── W^5: T-L interference (alpha=5, couples to cosPhi)
      w5 = r(10)*sinth2 + r(11)*costh*costh &
         - sq2*r(12)*sin2th*cosph &
         - r(13)*sinth2*cos2ph

!     ── W^6: T-L interference (alpha=6, couples to sinPhi)
      w6 = sq2*r(14)*sin2th*sinph &
         + r(15)*sinth2*sin2ph

!     ── Full angular distribution ─────────────────────────
!     W = (3/8pi^2) * { W0 - eps*cos2Phi*W1 - eps*sin2Phi*W2
!                       + sqrt(2*eps*(1+eps))*cosPhi*W5
!                       + sqrt(2*eps*(1+eps))*sinPhi*W6 }

      eTL = sqrt(max(0d0, 2d0*eps*(1d0+eps)))

      sdme_W = (3d0/(8d0*acos(-1d0)**2)) * ( &
               w0 &
             - eps*cos2Phi*w1 &
             - eps*sin2Phi*w2 &
             + eTL*cosPhi*w5 &
             + eTL*sinPhi*w6 )

!     Protect against small negative from numerics
      if(sdme_W .lt. 0d0) sdme_W = 0d0

      return
      end


!==============================================================
!     sdme_Wmax: conservative upper bound for accept/reject
!     Scans a coarse grid to find max of W, then adds 20% margin
!==============================================================
      function sdme_Wmax(eps, r)
      implicit real*8(a-h,o-z)
      real*8 r(15), sdme_Wmax, sdme_W
      external sdme_W

      aPI = acos(-1d0)
      wmax_val = 0d0
      nct  = 40
      nphi = 40
      nPHI = 40

      do i = 0, nct
        acosth = -1d0 + 2d0*dble(i)/dble(nct)
        do j = 0, nphi
          aphid = 2d0*aPI*dble(j)/dble(nphi)
          do k = 0, nPHI
            aPhiV = 2d0*aPI*dble(k)/dble(nPHI)
            wval = sdme_W(acosth, aphid, aPhiV, eps, r)
            if(wval .gt. wmax_val) wmax_val = wval
          enddo
        enddo
      enddo

!     Add 20% safety margin
      sdme_Wmax = 1.2d0 * wmax_val
      if(sdme_Wmax .lt. 1d-20) sdme_Wmax = 1d0

      return
      end


!==============================================================
!     calc_epsilon: virtual photon polarization parameter
!       eps = (1 - y - y^2*Q^2/(4*nu^2))
!           / (1 - y + y^2/2 + y^2*Q^2/(4*nu^2))
!     Input: ebeam = beam energy, y = nu/ebeam, Q2
!==============================================================
      function calc_epsilon(ebeam, y, Q2)
      implicit real*8(a-h,o-z)
      real*8 calc_epsilon

      anu = y*ebeam
      if(anu.gt.1d-10) then
        gamma2 = Q2/anu**2
      else
        calc_epsilon = 1d0
        return
      endif

      anum = 1d0 - y - 0.25d0*y*y*gamma2
      aden = 1d0 - y + 0.5d0*y*y + 0.25d0*y*y*gamma2

      if(aden.gt.1d-10) then
        calc_epsilon = anum/aden
      else
        calc_epsilon = 1d0
      endif

!     Physical bounds
      calc_epsilon = max(0d0, min(1d0, calc_epsilon))

      return
      end


!==============================================================
!     sample_sdme_angles: Accept/reject MC sampling of (cosΘ, φ, Φ)
!     from the full Schilling-Wolf W distribution.
!
!     Input:  eps    = virtual photon polarization
!             r(15)  = SDME array
!             iy     = random seed
!             wmax   = upper bound on W (precomputed by sdme_Wmax)
!     Output: costh  = cos(Theta) of h+ in helicity frame
!             phid   = phi of h+ (decay plane vs production plane)
!             aPhiOut= Phi (production plane vs lepton plane)
!==============================================================
      subroutine sample_sdme_angles(costh, phid, aPhiOut, &
                                    eps, r, wmax, iy)
      implicit real*8(a-h,o-z)
      real*4 urand
      real*8 r(15), sdme_W
      external sdme_W
      integer*4 iy

      aPI = acos(-1d0)

!     Accept/reject loop
 10   continue
        costh   = 2d0*dble(urand(iy)) - 1d0
        phid    = 2d0*aPI*dble(urand(iy))
        aPhiOut = 2d0*aPI*dble(urand(iy))
        rtest   = dble(urand(iy))

        wval = sdme_W(costh, phid, aPhiOut, eps, r)

        if(rtest*wmax .gt. wval) goto 10

      return
      end


!==============================================================
!     init_sdme: initialize SDME array to default values
!     mode = 0: all zero (isotropic: W = 3/(8 pi^2) * 1/2)
!     mode = 1: SCHC + NPE (photoproduction-like transverse)
!     mode = 2: user-specified (no change, assumes r already set)
!==============================================================
      subroutine init_sdme(r, mode)
      implicit real*8(a-h,o-z)
      real*8 r(15)
      integer mode

!     Zero everything first
      do i = 1, 15
        r(i) = 0d0
      enddo

      if(mode.eq.1) then
!       SCHC + NPE for photoproduction:
!       r04_00 = 0  (no longitudinal photon)
!       r1_11  = 1/2
!       r1_1-1 = 1/2
!       Im r2_1-1 = -1/2
!       All others zero.
!       This gives W0 ~ sin^2(Theta), with cos2Phi and sin2Phi
!       modulations from the virtual photon transverse polarization
        r(1)  = 0d0       ! r04_00
        r(4)  = 0.5d0     ! r1_11
        r(7)  = 0.5d0     ! r1_1-1
        r(9)  = -0.5d0    ! Im r2_1-1
      endif

      return
      end


      FUNCTION URAND(IY)
!     LCG random number generator.
!     The original 32-bit overflow relied on implicit INTEGER*4 wrap-around,
!     which gfortran -O2 computes in 64-bit, breaking the period.
!     Fix: use explicit 64-bit arithmetic and mask to 31 bits.
      INTEGER*4 IY
      INTEGER*8 IY8, IA8, IC8, M8
      REAL S
      DATA S,IA8,IC8,M8/.46566128E-9,843314861_8,453816693_8,2147483648_8/
      IY8 = INT(IY,8)
      IY8 = MOD(IY8*IA8 + IC8, M8)
      IF(IY8.LT.0) IY8 = IY8 + M8
      IY = INT(IY8,4)
      URAND=FLOAT(IY)*S
      END
