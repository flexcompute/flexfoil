program setbl_driver
  include 'XFOIL.INC'
  real :: setbl_samples(5)

  call bench_setbl(setbl_samples)

  call prepare_fixture()
  call setbl

  write(*,'(A)') '{'
  call emit_state('setbl', .true.)
  write(*,'(A)') '  "perf": {'
  call emit_perf_case('setbl', 20, setbl_samples, .false.)
  write(*,'(A)') '  }'
  write(*,'(A)') '}'

contains

  subroutine prepare_fixture()
    include 'XFOIL.INC'
    integer :: ibl

    call init_xfoil()
    call gen_naca0012(160)
    call pangen(.false.)
    lgamu = .false.
    call ggcalc
    lgamu = .true.
    alfa = 4.0d0 * atan(1.0d0) / 45.0d0
    call set_alpha_state()
    call xywake
    call qwcalc
    call qiset
    call stfind
    call iblpan
    call xicalc
    call iblsys
    call uicalc
    call qdcalc
    call blpini

    do ibl = 1, nbl(1)
      uedg(ibl,1) = uinv(ibl,1)
    end do
    do ibl = 1, nbl(2)
      uedg(ibl,2) = uinv(ibl,2)
    end do
  end subroutine prepare_fixture

  subroutine bench_setbl(samples)
    include 'XFOIL.INC'
    real, intent(out) :: samples(5)
    real :: t0, t1
    real :: thet_save(ivx,2), dstr_save(ivx,2), uedg_save(ivx,2), mass_save(ivx,2), ctau_save(ivx,2)
    real :: xssitr_save(2)
    logical :: lblini_save
    integer :: itran_save(2)
    integer :: sample_idx, iter

    call prepare_fixture()
    call save_bl_state(thet_save, dstr_save, uedg_save, mass_save, ctau_save, itran_save, xssitr_save, lblini_save)

    do sample_idx = 1, 5
      call cpu_time(t0)
      do iter = 1, 20
        call restore_bl_state(thet_save, dstr_save, uedg_save, mass_save, ctau_save, itran_save, xssitr_save, lblini_save)
        call setbl
      end do
      call cpu_time(t1)
      samples(sample_idx) = t1 - t0
    end do
  end subroutine bench_setbl

  subroutine save_bl_state(thet_out, dstr_out, uedg_out, mass_out, ctau_out, itran_out, xssitr_out, lblini_out)
    include 'XFOIL.INC'
    real, intent(out) :: thet_out(ivx,2), dstr_out(ivx,2), uedg_out(ivx,2), mass_out(ivx,2), ctau_out(ivx,2)
    real, intent(out) :: xssitr_out(2)
    logical, intent(out) :: lblini_out
    integer, intent(out) :: itran_out(2)

    thet_out = thet
    dstr_out = dstr
    uedg_out = uedg
    mass_out = mass
    ctau_out = ctau
    itran_out = itran
    xssitr_out = xssitr
    lblini_out = lblini
  end subroutine save_bl_state

  subroutine restore_bl_state(thet_in, dstr_in, uedg_in, mass_in, ctau_in, itran_in, xssitr_in, lblini_in)
    include 'XFOIL.INC'
    real, intent(in) :: thet_in(ivx,2), dstr_in(ivx,2), uedg_in(ivx,2), mass_in(ivx,2), ctau_in(ivx,2)
    real, intent(in) :: xssitr_in(2)
    logical, intent(in) :: lblini_in
    integer, intent(in) :: itran_in(2)

    thet = thet_in
    dstr = dstr_in
    uedg = uedg_in
    mass = mass_in
    ctau = ctau_in
    itran = itran_in
    xssitr = xssitr_in
    lblini = lblini_in
  end subroutine restore_bl_state

  subroutine emit_state(key, trailing)
    include 'XFOIL.INC'
    character(len=*), intent(in) :: key
    logical, intent(in) :: trailing

    write(*,'(A,A,A)') '  "', trim(key), '": {'
    write(*,'(A,I0,A)') '    "nbl_upper": ', nbl(1), ','
    write(*,'(A,I0,A)') '    "nbl_lower": ', nbl(2), ','
    write(*,'(A,I0,A)') '    "iblte_upper": ', iblte(1), ','
    write(*,'(A,I0,A)') '    "iblte_lower": ', iblte(2), ','
    write(*,'(A,I0,A)') '    "itran_upper": ', itran(1), ','
    write(*,'(A,I0,A)') '    "itran_lower": ', itran(2), ','
    write(*,'(A,ES24.16E3,A)') '    "xssitr_upper": ', xssitr(1), ','
    write(*,'(A,ES24.16E3,A)') '    "xssitr_lower": ', xssitr(2), ','
    if (lblini) then
      write(*,'(A)') '    "lblini": true,'
    else
      write(*,'(A)') '    "lblini": false,'
    end if
    call emit_real_array('upper_x', xssi(1:nbl(1),1), nbl(1), .true., 4)
    call emit_real_array('lower_x', xssi(1:nbl(2),2), nbl(2), .true., 4)
    call emit_real_array('upper_theta', thet(1:nbl(1),1), nbl(1), .true., 4)
    call emit_real_array('lower_theta', thet(1:nbl(2),2), nbl(2), .true., 4)
    call emit_real_array('upper_dstr', dstr(1:nbl(1),1), nbl(1), .true., 4)
    call emit_real_array('lower_dstr', dstr(1:nbl(2),2), nbl(2), .true., 4)
    call emit_real_array('upper_uedg', uedg(1:nbl(1),1), nbl(1), .true., 4)
    call emit_real_array('lower_uedg', uedg(1:nbl(2),2), nbl(2), .true., 4)
    call emit_real_array('upper_ctau', ctau(1:nbl(1),1), nbl(1), .true., 4)
    call emit_real_array('lower_ctau', ctau(1:nbl(2),2), nbl(2), .true., 4)
    call emit_real_array('upper_mass', mass(1:nbl(1),1), nbl(1), .true., 4)
    call emit_real_array('lower_mass', mass(1:nbl(2),2), nbl(2), .false., 4)
    if (trailing) then
      write(*,'(A)') '  },'
    else
      write(*,'(A)') '  }'
    end if
  end subroutine emit_state

  subroutine init_xfoil()
    include 'XFOIL.INC'
    logical :: ldbg
    integer :: ludbg, idbgcall, idbgiter
    common /XDEBUG/ ldbg, ludbg, idbgcall, idbgiter
    real :: pi_local

    ldbg = .false.
    ludbg = 6
    idbgcall = 0
    idbgiter = 0

    pi_local = 4.0d0 * atan(1.0d0)
    pi = pi_local
    hopi = 0.5d0 / pi_local
    qopi = 0.25d0 / pi_local
    dtor = pi_local / 180.0d0
    gamma = 1.4d0
    gamm1 = gamma - 1.0d0

    lvisc = .true.
    lalfa = .true.
    lwake = .false.
    lpacc = .false.
    lblini = .false.
    lipan = .false.
    lqaij = .false.
    ladij = .false.
    lwdij = .false.
    lvconv = .false.
    lgamu = .false.
    lqinu = .false.
    sharp = .false.
    n = 0
    nb = 0
    nw = 0
    ist = 0
    qinf = 1.0d0
    matyp = 1
    minf = 0.0d0
    minf1 = 0.0d0
    reinf = 1.0d6
    reinf1 = 1.0d6
    reinf_cl = 1.0d6
    alfa = 0.0d0
    cosa = 1.0d0
    sina = 0.0d0
    xcmref = 0.25d0
    ycmref = 0.0d0
    waklen = 1.0d0
    cvpar = 1.0d0
    cterat = 0.15d0
    ctrrat = 0.2d0
    xsref1 = 1.0d0
    xsref2 = 1.0d0
    xpref1 = 1.0d0
    xpref2 = 1.0d0
    acrit(1) = 9.0d0
    acrit(2) = 9.0d0
    xstrip(1) = 1.0d0
    xstrip(2) = 1.0d0
    itran(1:2) = 0
    xssitr(1:2) = 0.0d0
    xoctr(1:2) = 0.0d0
    gam(1:iqx) = 0.0d0
    gam_a(1:iqx) = 0.0d0
    qinv(1:izx) = 0.0d0
    qinv_a(1:izx) = 0.0d0
    qinvu(1:izx,1) = 0.0d0
    qinvu(1:izx,2) = 0.0d0
    sig(1:izx) = 0.0d0
    xssi(1:ivx,1:2) = 0.0d0
    uedg(1:ivx,1:2) = 0.0d0
    thet(1:ivx,1:2) = 0.0d0
    dstr(1:ivx,1:2) = 0.0d0
    mass(1:ivx,1:2) = 0.0d0
    ctau(1:ivx,1:2) = 0.0d0
    tau(1:ivx,1:2) = 0.0d0
    dis(1:ivx,1:2) = 0.0d0
    ctq(1:ivx,1:2) = 0.0d0
    delt(1:ivx,1:2) = 0.0d0
    tstr(1:ivx,1:2) = 0.0d0
  end subroutine init_xfoil

  subroutine gen_naca0012(npanel)
    include 'XFOIL.INC'
    integer, intent(in) :: npanel
    integer :: i, nhalf
    real :: beta, xx, yt, pi_local, t

    pi_local = 4.0d0 * atan(1.0d0)
    t = 0.12d0
    nhalf = npanel / 2
    n = 0
    do i = 0, nhalf
      beta = pi_local * dble(i) / dble(nhalf)
      xx = 0.5d0 * (1.0d0 - cos(beta))
      yt = 5.0d0 * t * (0.2969d0 * sqrt(xx) - 0.126d0 * xx - 0.3516d0 * xx**2 &
        + 0.2843d0 * xx**3 - 0.1036d0 * xx**4)
      n = n + 1
      x(n) = xx
      y(n) = yt
    end do
    do i = 1, nhalf
      beta = pi_local * dble(nhalf - i) / dble(nhalf)
      xx = 0.5d0 * (1.0d0 - cos(beta))
      yt = 5.0d0 * t * (0.2969d0 * sqrt(xx) - 0.126d0 * xx - 0.3516d0 * xx**2 &
        + 0.2843d0 * xx**3 - 0.1036d0 * xx**4)
      n = n + 1
      x(n) = xx
      y(n) = -yt
    end do
    name = 'NACA 0012'
  end subroutine gen_naca0012

  subroutine set_alpha_state()
    include 'XFOIL.INC'
    integer :: i

    cosa = cos(alfa)
    sina = sin(alfa)
    do i = 1, n
      gam(i) = cosa * gamu(i,1) + sina * gamu(i,2)
      gam_a(i) = -sina * gamu(i,1) + cosa * gamu(i,2)
    end do
    psio = cosa * gamu(n+1,1) + sina * gamu(n+1,2)
    call set_te_geometry()
    call qiset
  end subroutine set_alpha_state

  subroutine set_te_geometry()
    include 'XFOIL.INC'
    real :: dxte, dyte, dxs, dys, scs, sds

    dxte = x(1) - x(n)
    dyte = y(1) - y(n)
    dxs = 0.5d0 * (-xp(1) + xp(n))
    dys = 0.5d0 * (-yp(1) + yp(n))

    ante = dxs * dyte - dys * dxte
    aste = dxs * dxte + dys * dyte
    dste = sqrt(dxte**2 + dyte**2)
    sharp = dste .lt. 0.0001d0 * chord

    if (sharp) then
      scs = 1.0d0
      sds = 0.0d0
    else
      scs = ante / dste
      sds = aste / dste
    end if

    sigte = 0.5d0 * (gam(1) - gam(n)) * scs
    gamte = -0.5d0 * (gam(1) - gam(n)) * sds
    sigte_a = 0.5d0 * (gam_a(1) - gam_a(n)) * scs
    gamte_a = -0.5d0 * (gam_a(1) - gam_a(n)) * sds
  end subroutine set_te_geometry

  subroutine emit_perf_case(name, inner_loops, samples, trailing)
    character(len=*), intent(in) :: name
    integer, intent(in) :: inner_loops
    real, intent(in) :: samples(5)
    logical, intent(in) :: trailing
    real :: sorted(5)

    sorted = samples
    call sort_real(sorted, 5)

    write(*,'(A,A,A)') '    "', trim(name), '": {'
    write(*,'(A,I0,A)') '      "inner_loops": ', inner_loops, ','
    call emit_real_array('samples_seconds', sorted, 5, .true., 6)
    write(*,'(A,ES24.16E3)') '      "median_seconds": ', sorted(3)
    if (trailing) then
      write(*,'(A)') '    },'
    else
      write(*,'(A)') '    }'
    end if
  end subroutine emit_perf_case

  subroutine emit_real_array(name, arr, nvals, trailing, indent)
    character(len=*), intent(in) :: name
    integer, intent(in) :: nvals, indent
    real, intent(in) :: arr(*)
    logical, intent(in) :: trailing
    integer :: i

    write(*,'(A)',advance='no') repeat(' ', indent) // '"' // trim(name) // '": ['
    do i = 1, nvals
      if (i > 1) write(*,'(A)',advance='no') ', '
      write(*,'(ES24.16E3)',advance='no') arr(i)
    end do
    if (trailing) then
      write(*,'(A)') '],'
    else
      write(*,'(A)') ']'
    end if
  end subroutine emit_real_array

  subroutine sort_real(arr, nvals)
    integer, intent(in) :: nvals
    real, intent(inout) :: arr(nvals)
    integer :: i, j
    real :: tmp

    do i = 1, nvals - 1
      do j = i + 1, nvals
        if (arr(j) < arr(i)) then
          tmp = arr(i)
          arr(i) = arr(j)
          arr(j) = tmp
        end if
      end do
    end do
  end subroutine sort_real

end program setbl_driver
