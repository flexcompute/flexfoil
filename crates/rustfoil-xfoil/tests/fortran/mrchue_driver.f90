program mrchue_driver
  include 'XFOIL.INC'

  call prepare_fixture()
  call blpini
  call mrchue

  write(*,'(A)') '{'
  call emit_state('mrchue')
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
    alfa = 15.0d0 * atan(1.0d0) / 45.0d0
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

    do ibl = 1, nbl(1)
      uedg(ibl,1) = uinv(ibl,1)
    end do
    do ibl = 1, nbl(2)
      uedg(ibl,2) = uinv(ibl,2)
    end do
  end subroutine prepare_fixture

  subroutine emit_state(key)
    include 'XFOIL.INC'
    character(len=*), intent(in) :: key

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
    write(*,'(A)') '  }'
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
    scs = (dxte*dxs + dyte*dys) / max(1.0d-12, sqrt((dxte*dxte + dyte*dyte) * (dxs*dxs + dys*dys)))
    sds = (dxte*dys - dyte*dxs) / max(1.0d-12, sqrt((dxte*dxte + dyte*dyte) * (dxs*dxs + dys*dys)))
    sharp = abs(sds) .lt. 1.0d-4
    if (sharp) then
      dste = 0.0d0
      ante = 0.0d0
      aste = 0.0d0
    else
      dste = sqrt(dxte*dxte + dyte*dyte)
      aste = atan2(dyte, dxte)
      ante = max(0.0d0, -0.5d0*scs)
    end if
  end subroutine set_te_geometry

  subroutine emit_real_array(name, values, count, trailing, indent)
    character(len=*), intent(in) :: name
    real, intent(in) :: values(*)
    integer, intent(in) :: count
    logical, intent(in) :: trailing
    integer, intent(in) :: indent
    integer :: i
    character(len=16) :: pad

    pad = '                '
    write(*,'(A,A,A)', advance='no') pad(1:indent), '"', trim(name)
    write(*,'(A)', advance='no') '": ['
    do i = 1, count
      if (i .gt. 1) write(*,'(A)', advance='no') ', '
      write(*,'(ES24.16E3)', advance='no') values(i)
    end do
    if (trailing) then
      write(*,'(A)') '],'
    else
      write(*,'(A)') ']'
    end if
  end subroutine emit_real_array

end program mrchue_driver
