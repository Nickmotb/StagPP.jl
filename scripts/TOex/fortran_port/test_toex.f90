! Unit tests of the unified TOex port, linked against the scratch watercarb objects.
!  parity: toex_coh_state_u and toex_solid_curve on states.txt (StagPP reference in j_*.txt)
!  closure: both engines on a P-T-H grid, element closure and the solid sink's consistency
!  cell: synthetic cell through gather -> solve -> write-back, conservation per tracer
!  speed: per-cell cost of legacy vs unified on representative mantle cells
program test_toex
  use mpi
  use control,      only: me0
  use refstat,      only: Stag_fO2, Stag_metal, Stag_fO2step, Stag_sH2O, Stag_sH2Oa, Stag_H2OSdens, &
                          Stag_H2OMdens, EM_Xhz, EM_Xbs
  use PerplexStuff, only: load_s_fO2_data
  use volatilestuff
  implicit none
  character(len=256) :: tdir, mode
  integer :: ierr
  call MPI_Init(ierr)
  me0 = .true.
  call get_command_argument(1, mode)
  call get_command_argument(2, tdir)
  Stag_fO2      = trim(tdir)//'/StagfO2.dat'
  Stag_metal    = trim(tdir)//'/StagMetal.dat'
  Stag_fO2step  = trim(tdir)//'/StagfO2step.dat'
  Stag_sH2O     = trim(tdir)//'/StagH2O.dat'
  Stag_sH2Oa    = trim(tdir)//'/StagH2Oa.dat'
  Stag_H2OSdens = trim(tdir)//'/StagSdens.dat'
  Stag_H2OMdens = trim(tdir)//'/StagMdens.dat'
  ! watercarb defaults (defaults_include.f90): the fO2 blend weights are built from these
  EM_Xhz = [43.43, 45.93, 8.34, 0.9, 1.0, 0.01, 0.3]
  EM_Xbs = [50.42, 9.77, 7.1, 12.54, 16.8, 2.23, 0.07]
  call load_s_fO2_data(.true., .true.)
  call toex_bind_solid_fO2
  call toex_bind_sH2Oa
  select case (trim(mode))
  case ('parity');  call t_parity
  case ('closure'); call t_closure
  case ('cell');    call t_cell
  case ('speed');   call t_speed
  case ('debug');   call t_debug
  end select
  call MPI_Finalize(ierr)

contains

  ! Inner solve at the Julia states, written for compare.jl
  subroutine t_parity
    type(toex_ctx_t) :: ctx
    type(toex_coh_t) :: cs
    real(8) :: f, P, T, pb, phi, TOex, TC, H, mx, C(TOEX_NAMAX)
    integer :: k, u1, u2, u3, ios
    open(newunit=u1, file='states.txt', status='old')
    open(newunit=u2, file='f_coh.txt', status='replace')
    open(newunit=u3, file='f_curve.txt', status='replace')
    do
       read(u1,*,iostat=ios) k, f, P, T, pb, phi, TOex, TC, H
       if (ios /= 0) exit
       toex_unified = .true.
       call toex_build_ctx(P, T, pb, phi, TOex, TC, H, ctx, mx)
       toex_phi_seed_ok = .false.;  toex_u_seed = VQNAN;  zd_Vwarm = -1.d0
       call toex_coh_state_u(f, ctx, cs)
       write(u2,'(i6,22es25.16e3)') k, cs%aC, cs%n_fluid, cs%XH2O, cs%XCO2, cs%XCH4, cs%XCO, cs%XH2, &
            cs%nC_fluid, cs%nC_melt, cs%nC_solid, cs%nC_carb, cs%molH_fl, cs%nH_solid, cs%nH2_melt, &
            cs%nH_H2O_melt, cs%a_solid, cs%fH2, cs%fex, cs%S_fluid, dble(cs%it_ideality), cs%vol_fluid, &
            cs%res_ideality
       call toex_solid_curve(P, T, pb, C)
       write(u3,'(i6,4es25.16e3)') k, C(1), C(6), C(11), C(16)
    end do
    close(u1);  close(u2);  close(u3)
  end subroutine t_parity

  ! Both engines over a grid; unified gets the total water, legacy the free water above C(1)
  subroutine t_closure
    type(toex_result_t) :: r, rl
    real(8) :: Ps(6), Hs(5), ph(2), C(TOEX_NAMAX), P, T, phi, H, Hf, TOex, TC, pb
    real(8) :: eH, eC, eO, eS, wH, wC, wO, wS, Ms
    integer :: i, j, l, m, ntot, nconv, nlconv, nfl, kk
    Ps = [0.5d0, 2.d0, 5.d0, 10.d0, 20.d0, 60.d0]
    Hs = [1.d-5, 1.d-4, 5.d-4, 2.d-3, 1.d-2]
    ph = [0.d0, 0.02d0]
    TOex = 1.d-3;  TC = 1.d-4
    ntot = 0;  nconv = 0;  nlconv = 0;  nfl = 0
    wH = 0.d0;  wC = 0.d0;  wO = 0.d0;  wS = 0.d0
    print '(a)', '   P GPa      T K   pbas   phi       H      | conv  dFMQ*   a_solid  solid/H   melt/H  fluid/H |' // &
                 ' legacy conv  f_uni-f_leg'
    do i = 1, 6
       do kk = 1, 2
          P = Ps(i)
          T = merge(1600.d0 + 12.d0*min(P,23.5d0) + 5.d0*max(P-23.5d0,0.d0), min(700.d0+60.d0*P, 1300.d0+12.d0*P), kk==1)
          do l = 1, 2
             phi = ph(l);  if (kk == 2) phi = 0.d0
             do j = 1, 5
                do m = 1, 2
                   pb = merge(0.d0, 1.d0, m==1);  H = Hs(j)
                   toex_unified = .true.
                   call toex_partition(P, T, pb, phi, TOex, TC, H, r)
                   ntot = ntot + 1
                   if (r%converged .and. r%ok) nconv = nconv + 1
                   if (r%mH_fluid > 0.d0) nfl = nfl + 1
                   ! closures in grams on the 1 kg reference cell, relative to each budget
                   eH = abs(r%mH_fluid + r%mH_melt + r%mH_solid - 1.d3*H)/(1.d3*H)
                   eC = abs(r%mC_gph + r%mC_carb + r%mC_melt + r%mC_fluid - 1.d3*TC)/(1.d3*TC)
                   eO = abs(r%sOex + r%mOex + r%carbonO + r%carbO + r%cohsink + r%meltH2ex - 1.d3*TOex)/(1.d3*TOex)
                   Ms = 1.d3*(1.d0 - phi)
                   eS = abs(r%mH_solid - 1.d-2*Ms*r%C_solid)/max(1.d-2*Ms*r%C_solid_sat, 1.d-30)
                   wH = max(wH, eH);  wC = max(wC, eC);  wO = max(wO, eO);  wS = max(wS, eS)
                   ! legacy on the same cell: only the water above the solid's a = 1 capacity is free
                   call toex_solid_curve(P, T, pb, C)
                   Hf = max(H - (1.d0 - phi)*1.d-2*C(1), 0.d0)
                   toex_unified = .false.
                   call toex_partition(P, T, pb, phi, TOex, TC, Hf, rl)
                   if (rl%converged .and. rl%ok) nlconv = nlconv + 1
                   if (m == 1 .and. (j == 1 .or. j == 3 .or. j == 5)) &
                      print '(f7.1,f9.0,f6.1,f6.2,es10.2,a,l5,f8.3,f10.4,3f9.4,a,l5,f12.4)', P, T, pb, phi, H, ' |', &
                         r%converged, r%fO2 - toex_logfO2_EMOG(T, 1.d4*P), r%aH2O_solid, &
                         r%mH_solid/(1.d3*H), r%mH_melt/(1.d3*H), r%mH_fluid/(1.d3*H), ' |', rl%converged, r%fO2 - rl%fO2
                end do
             end do
          end do
       end do
    end do
    print '(a)', '   (* dFMQ column is fO2 - EMOG here)'
    print '(a,i4,a,i4,a,i4,a,i4)', ' unified converged ', nconv, ' of ', ntot, '   legacy converged ', nlconv, &
         '   cells with a fluid ', nfl
    print '(a,es9.2,a,es9.2,a,es9.2,a,es9.2)', ' worst closure  H ', wH, '  C ', wC, '  O ', wO, &
         '   solid water vs Ms*C(a) ', wS
  end subroutine t_closure

  ! One cell of 12 tracers: harzburgite, basalt and melt, through toex_react_cell
  subroutine t_cell
    use tracerstuff, only: tra, celltra, TYPEPOS, MASSPOS, TEMPPOS, ntracervar, traceelement_water, &
         traceelement_carbon, harz_bas_TTG_cc, tt_solid_harzburgite, tt_solid_basalt, &
         tt_molten_harzburgite, tt_molten_basalt, tt_solid_TTG, tt_solid_cc, tt_molten_TTG, tt_molten_cc
    use refstat,     only: Pref_adiabat
    use grid,        only: gl
    use volatileprocs, only: toex_react_cell
    integer, parameter :: NT = 12
    real(8) :: w0, c0, o0, w1, c1, o1, Ci(NT), Cv(TOEX_NAMAX), lC(TOEX_NAMAX), bs, rmin, rmax
    integer :: i, ieng, cas
    real(8) :: Pc, Tc, Hb, Hf
    character(len=40) :: lbl
    traceelement_water = .true.;  traceelement_carbon = .true.;  harz_bas_TTG_cc = .false.
    tt_solid_harzburgite = 0.;  tt_solid_basalt = 1.;  tt_solid_TTG = 1.1;  tt_solid_cc = 1.2
    tt_molten_harzburgite = 2.;  tt_molten_basalt = 3.;  tt_molten_TTG = 3.1;  tt_molten_cc = 3.2
    TYPEPOS = 1;  MASSPOS = 2;  TEMPPOS = 3;  H2OPOS = 4;  H2OFREEPOS = 5;  WSOLPOS = 6
    OxCPOS = 7;  ReCPOS = 8;  fO2POS = 9;  F3FTPOS = 10;  OEXPOS = 11;  CFLUIDPOS = 12;  OFLUIDPOS = 13
    ntracervar = 13
    allocate(tra(ntracervar, NT, 1));  allocate(celltra(0:0,0:0,0:0,1))
    celltra(0,0,0,1)%n = NT;  celltra(0,0,0,1)%nother = 0
    allocate(celltra(0,0,0,1)%traNum(NT))
    do i = 1, NT;  celltra(0,0,0,1)%traNum(i) = i;  end do
    gl%iz0 = 0;  allocate(Pref_adiabat(0:4))
    do cas = 1, 3
       select case (cas)
       case (1);  Pc = 3.d0;  Tc = 1636.d0;  Hb = 0.01d0;  Hf = 0.d0;   lbl = 'adiabat 3 GPa, 100 ppm, 2 melt'
       case (2);  Pc = 3.d0;  Tc = 880.d0;   Hb = 1.0d0;   Hf = 0.5d0;  lbl = 'slab 3 GPa, wet, fluid'
       case (3);  Pc = 15.d0; Tc = 1780.d0;  Hb = 0.05d0;  Hf = 0.d0;   lbl = 'transition zone 15 GPa'
       end select
       Pref_adiabat = Pc*1.d9
       do ieng = 1, 2
          toex_unified = (ieng == 1)
          do i = 1, NT
             tra(:,i,1) = 0.
             tra(MASSPOS,i,1) = 1.0 + 0.1*i                    ! kg
             tra(TEMPPOS,i,1) = Tc
             if (i <= 5) then;  tra(TYPEPOS,i,1) = 0.              ! harzburgite
             else if (i <= 10) then;  tra(TYPEPOS,i,1) = 1.        ! basalt
             else;  tra(TYPEPOS,i,1) = 3.;  end if                 ! molten basalt
             if (cas == 2 .and. i > 10) tra(TYPEPOS,i,1) = 1.      ! no melt in the slab
             tra(H2OPOS,i,1)     = Hb*(0.5 + 0.1*i)                ! wt%, uneven on purpose
             tra(H2OFREEPOS,i,1) = Hf
             tra(WSOLPOS,i,1)    = 1.0
             tra(OxCPOS,i,1)     = 0.005;  tra(ReCPOS,i,1) = 0.005
             tra(OEXPOS,i,1)     = 0.1 + 0.01*i
             tra(fO2POS,i,1)     = 0.
          end do
          call totals(NT, w0, c0, o0)
          call toex_react_cell(0,0,0,1)
          call totals(NT, w1, c1, o1)
          ! solid tracers: H2OPOS over their own C_i(a) must be one number (unified)
          rmin = huge(1.d0);  rmax = 0.d0
          do i = 1, NT
             if (tra(TYPEPOS,i,1) >= 2.) cycle
             bs = tra(TYPEPOS,i,1)
             call toex_solid_curve(Pc, Tc, bs, Cv)
             lC = 0.d0;  lC(1:toex_sa_nA) = log(max(Cv(1:toex_sa_nA), 1.d-12))
             Ci(i) = toex_solid_h2o_wt(toex_last_react_a, toex_sa_la, lC, toex_sa_nA)
             if (Ci(i) > 0.d0) then
                rmin = min(rmin, tra(H2OPOS,i,1)/Ci(i));  rmax = max(rmax, tra(H2OPOS,i,1)/Ci(i))
             end if
          end do
          print '(a,a,a,l2,a,3es10.2,a,f8.4,a,f8.5,a,f8.5)', ' ', lbl(1:30), merge(' unified', ' legacy ', ieng==1), &
               toex_last_react_conv, '  |d(H,C,O)|/total', abs(w1-w0)/w0, abs(c1-c0)/c0, abs(o1-o0)/o0, &
               '  a_solid', toex_last_react_a, '  bound/C_i(a) min', rmin, ' max', rmax
       end do
    end do
  end subroutine t_cell

  subroutine totals(NT, w, c, o)
    use tracerstuff, only: tra, MASSPOS
    integer, intent(in)  :: NT
    real(8), intent(out) :: w, c, o
    integer :: j
    w = 0.d0;  c = 0.d0;  o = 0.d0
    do j = 1, NT
       w = w + (tra(H2OPOS,j,1) + tra(H2OFREEPOS,j,1))*tra(MASSPOS,j,1)
       c = c + (tra(OxCPOS,j,1) + tra(ReCPOS,j,1) + tra(CFLUIDPOS,j,1))*tra(MASSPOS,j,1)
       o = o + (tra(OEXPOS,j,1) + tra(OFLUIDPOS,j,1))*tra(MASSPOS,j,1)
    end do
  end subroutine totals

  ! Cost per cell, legacy vs unified, on the same physical cells: cold start and a warm start
  !  0.05 log units off the previous answer (what a time step hands the solver)
  subroutine t_speed
    ! engine 1 legacy (free water), 2 unified with the solid (total water), 3 unified with the solid's
    !  capacity zeroed and the free water only: the new algorithm on the legacy physics
    integer, parameter :: NP = 24, NREP = 10, NE = 3
    type(toex_result_t) :: r
    real(8) :: P, T, phi, H, Hf, C(TOEX_NAMAX), fsol(NE,NP,4), t0, t1, cost(NE,2,4), C1(NP,4)
    integer(8) :: ncoh(NE,2,4), neos(NE,2,4), nres(3,NE,2,4)
    integer :: ieng, iw, k, rep, reg, nfail(NE,4)
    real(4), allocatable :: tab_keep(:,:,:,:,:)
    character(len=34) :: rname(4)
    rname = [character(len=34) :: 'adiabat, 300 ppm water, no melt', 'adiabat <= 6 GPa, 300 ppm, 1% melt', &
             'slab, 1 wt% water (fluid)', 'adiabat, 300 ppm, no carbon']
    cost = 0.d0;  ncoh = 0;  neos = 0;  nfail = 0;  nres = 0
    tab_keep = toex_sa_tab
    do reg = 1, 4
       do k = 1, NP
          call cellstate(reg, k, P, T, phi, H);  call toex_solid_curve(P, T, 0.d0, C);  C1(k,reg) = C(1)
       end do
    end do
    ! warm-up: the IDV table and the ZD09 constants are built on first use
    do k = 1, NP
       call cellstate(1, k, P, T, phi, H)
       toex_unified = .false.;  call toex_partition(P, T, 0.d0, phi, 1.d-3, 1.d-4, 0.d0, r)
       toex_unified = .true.;   call toex_partition(P, T, 0.d0, phi, 1.d-3, 1.d-4, H, r)
    end do
    do reg = 1, 4
       do ieng = 1, NE
          toex_unified = (ieng >= 2)
          if (ieng == 3) toex_sa_tab = 0.          ! outside the timing: the copy is 19 MB
          do iw = 1, 2
             do rep = 1, NREP
                toex_n_coh = 0;  toex_n_eos = 0;  toex_u_nres = 0
                call cpu_time(t0)
                do k = 1, NP
                   call cellstate(reg, k, P, T, phi, H)
                   Hf = max(H - (1.d0 - phi)*1.d-2*C1(k,reg), 0.d0)
                   if (ieng == 2) Hf = H
                   if (iw == 1) then
                      call toex_partition(P, T, 0.d0, phi, 1.d-3, merge(0.d0, 1.d-4, reg==4), Hf, r)
                      fsol(ieng,k,reg) = r%fO2
                   else
                      call toex_partition(P, T, 0.d0, phi, 1.d-3, merge(0.d0, 1.d-4, reg==4), Hf, r, &
                                          f_guess=fsol(ieng,k,reg) + 0.05d0)
                   end if
                   if (rep == 1 .and. iw == 1 .and. .not.(r%converged .and. r%ok)) nfail(ieng,reg) = nfail(ieng,reg) + 1
                end do
                call cpu_time(t1)
                cost(ieng,iw,reg) = cost(ieng,iw,reg) + (t1 - t0)/(NP*NREP)
                ncoh(ieng,iw,reg) = ncoh(ieng,iw,reg) + toex_n_coh
                neos(ieng,iw,reg) = neos(ieng,iw,reg) + toex_n_eos
                nres(:,ieng,iw,reg) = nres(:,ieng,iw,reg) + toex_u_nres
             end do
          end do
          if (ieng == 3) toex_sa_tab = tab_keep
       end do
    end do
    toex_sa_tab = tab_keep
    print '(a)', ' ms per cell solve (single core)        start   legacy   unified  unified,no-solid' // &
                 '   unified/legacy   solid share of unified'
    do reg = 1, 4
       do iw = 1, 2
          print '(1x,a34,a7,f9.3,f10.3,f13.3,f17.1,f20.0,a)', rname(reg), merge(' cold ', ' warm ', iw==1), &
               1.d3*cost(1,iw,reg), 1.d3*cost(2,iw,reg), 1.d3*cost(3,iw,reg), cost(2,iw,reg)/cost(1,iw,reg), &
               1.d2*(cost(2,iw,reg) - cost(3,iw,reg))/cost(2,iw,reg), ' %'
       end do
    end do
    print '(a,4i3,a,4i3,a,4i3)', ' unconverged cells: legacy', nfail(1,:), '   unified', nfail(2,:), '   no-solid', nfail(3,:)
    print '(a)', ' per cell (warm start):         coh_state calls  leg/uni     EOS calls leg/uni' // &
                 '     unified residuals per coh_state call: carbon  hydrogen  fluid-n'
    do reg = 1, 4
       print '(1x,a34,2f8.1,2f10.1,3f12.1)', rname(reg), &
            dble(ncoh(1,2,reg))/(NP*NREP), dble(ncoh(2,2,reg))/(NP*NREP), &
            dble(neos(1,2,reg))/(NP*NREP), dble(neos(2,2,reg))/(NP*NREP), &
            dble(nres(:,2,2,reg))/max(dble(ncoh(2,2,reg)),1.d0)
    end do
  end subroutine t_speed

  ! Representative cells: column from 0.3 to 40 GPa (adiabat), or a slab column
  subroutine cellstate(reg, k, P, T, phi, H)
    integer, intent(in) :: reg, k
    real(8), intent(out) :: P, T, phi, H
    real(8) :: Tad
    P = 0.3d0*(40.d0/0.3d0)**(dble(k-1)/23.d0)
    if (reg == 2) P = 0.3d0*(6.d0/0.3d0)**(dble(k-1)/23.d0)
    Tad = 1600.d0 + 12.d0*min(P,23.5d0) + 5.d0*max(P-23.5d0,0.d0)
    T = Tad;  phi = 0.d0;  H = 3.d-4
    if (reg == 2) phi = 0.01d0
    if (reg == 3) then;  T = min(700.d0 + 60.d0*P, Tad - 300.d0);  H = 1.d-2;  end if
  end subroutine cellstate

  subroutine t_debug
    real(8) :: Cv(TOEX_NAMAX), lC(TOEX_NAMAX), pp(5), cc(5), aa(2)
    integer :: q, ja
    type(toex_ctx_t) :: ctx
    type(toex_state_t) :: st
    real(8) :: mx, f, pb
    integer :: i, m
    pp = [0.d0, 0.2d0, 0.58d0, 0.8d0, 1.d0];  aa = [1.d0, 0.3867d0]
    do ja = 1, 2
       do q = 1, 5
          call toex_solid_curve(3.d0, 880.d0, pp(q), Cv)
          lC = 0.d0;  lC(1:toex_sa_nA) = log(max(Cv(1:toex_sa_nA), 1.d-12))
          cc(q) = toex_solid_h2o_wt(aa(ja), toex_sa_la, lC, toex_sa_nA)
       end do
       print '(a,f7.4,a,5f9.4,a,f8.4)', ' 3 GPa 880 K  a =', aa(ja), '  C(p = 0, .2, .58, .8, 1) wt%', cc, &
            '   mass mix 0.42*C(0)+0.58*C(1)', 0.42d0*cc(1) + 0.58d0*cc(5)
    end do
    do m = 1, 2
       pb = merge(0.d0, 1.d0, m==1)
       call toex_build_ctx(2.d0, 1624.d0, pb, 0.d0, 1.d-3, 1.d-4, 1.d-4, ctx, mx)
       print '(a,f4.1,a,es10.3,a,es10.3,a,l2,a,2f8.4)', ' pbas', pb, '  maxsOex', ctx%maxsOex, '  minsOex', ctx%minsOex, &
            '  sfO2_ready', ctx%sfO2_ready, '  wcol', ctx%wcol
       do i = 0, 8
          f = -120.d0 + 20.d0*i
          call toex_derive(f, ctx, st)
          print '(a,f7.1,a,es11.3,a,es11.3,a,es11.3,a,es11.3)', '   f', f, '  total', st%total, '  sOex', st%sOex, &
               '  fs', st%fs, '  fex', st%fex
       end do
    end do
  end subroutine t_debug

end program test_toex
