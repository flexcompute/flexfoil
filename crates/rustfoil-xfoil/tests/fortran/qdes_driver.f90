program qdes_driver
  include 'XFOIL.INC'
  character(len=32) :: arg
  real :: alpha_deg, perturb, pi_local
  integer :: niterq
  real :: base_x(iqx), base_y(iqx)
  integer :: i, ile

  alpha_deg = 4.0d0
  perturb = 0.02d0
  niterq = 2

  call getarg(1, arg)
  if (len_trim(arg) > 0) read(arg, *) alpha_deg
  call getarg(2, arg)
  if (len_trim(arg) > 0) read(arg, *) perturb
  call getarg(3, arg)
  if (len_trim(arg) > 0) read(arg, *) niterq

  call init_xfoil()
  call gen_naca0012(80)
  base_x(1:n) = x(1:n)
  base_y(1:n) = y(1:n)

  lgamu = .false.
  call ggcalc
  lgamu = .true.

  pi_local = 4.0d0 * atan(1.0d0)
  alfa = alpha_deg * pi_local / 180.0d0
  call set_alpha_state()

  if (nsp .ne. n) then
    lqspec = .false.
    liqset = .false.
  end if

  algam = alfa
  clgam = cl
  cmgam = cm
  chx = xte - xle
  chy = yte - yle
  chsq = chx**2 + chy**2
  nsp = n
  do i = 1, nsp
    qgamm(i) = gam(i)
    sspec(i) = s(i) / s(n)
    xspoc(i) = ((x(i) - xle) * chx + (y(i) - yle) * chy) / chsq
    yspoc(i) = ((y(i) - yle) * chx - (x(i) - xle) * chy) / chsq
  end do
  ssple = sle / s(n)
  nqsp = 1
  kqtarg = 1

  call gamqsp(1)

  ile = 1
  do i = 2, n
    if (x(i) .lt. x(ile)) ile = i
  end do
  iq1 = 1
  iq2 = ile
  liqset = .true.

  call perturb_upper_qspec(1, perturb, iq1, iq2)
  call splqsp(1)
  call mixed(1, niterq)

  write(*,'(A)') '{'
  write(*,'(A,F10.4,A)') '  "alpha_deg": ', alpha_deg, ','
  write(*,'(A,F10.4,A)') '  "perturb": ', perturb, ','
  write(*,'(A,I0,A)') '  "iterations": ', niterq, ','
  call emit_real_array('base_x', base_x, n, .true., 2)
  call emit_real_array('base_y', base_y, n, .true., 2)
  call emit_targets(iq1, iq2, ile, n)
  call emit_real_array('output_x', x, n, .true., 2)
  call emit_real_array('output_y', y, n, .true., 2)
  write(*,'(A,ES24.16E3,A)') '  "cl": ', cl, ','
  write(*,'(A,ES24.16E3,A)') '  "cm": ', cm, ','
  write(*,'(A,ES24.16E3)') '  "cdp": ', cdp
  write(*,'(A)') '}'

contains

  subroutine perturb_upper_qspec(kqsp, amplitude, start_idx, end_idx)
    include 'XFOIL.INC'
    integer, intent(in) :: kqsp, start_idx, end_idx
    real, intent(in) :: amplitude
    integer :: i
    real :: fs, blend

    do i = start_idx, end_idx
      fs = (s(i) - s(start_idx)) / max(s(end_idx) - s(start_idx), 1.0d-12)
      blend = sin(4.0d0 * atan(1.0d0) * fs)**2
      qspec(i, kqsp) = qspec(i, kqsp) * (1.0d0 + amplitude * blend)
    end do
  end subroutine perturb_upper_qspec

  subroutine emit_targets(start_idx, end_idx, le_idx, npts)
    include 'XFOIL.INC'
    integer, intent(in) :: start_idx, end_idx, le_idx, npts
    integer :: i, j, count
    integer :: lower_count
    real :: out_x(iqx), out_q(iqx), lower_x(iqx), lower_q(iqx)

    count = end_idx - start_idx + 1
    do i = 1, count
      j = end_idx - i + 1
      out_x(i) = x(j)
      out_q(i) = qspec(j, 1)
    end do

    call emit_real_array('target_upper_x', out_x, count, .true., 2)
    call emit_real_array('target_upper_q', out_q, count, .true., 2)

    lower_count = npts - le_idx + 1
    do i = 1, lower_count
      j = le_idx + i - 1
      lower_x(i) = x(j)
      lower_q(i) = qspec(j, 1)
    end do

    call emit_real_array('target_lower_x', lower_x, lower_count, .true., 2)
    call emit_real_array('target_lower_q', lower_q, lower_count, .true., 2)
  end subroutine emit_targets

  subroutine emit_real_array(name, arr, nvals, trailing, indent)
    character(len=*), intent(in) :: name
    real, intent(in) :: arr(*)
    integer, intent(in) :: nvals, indent
    logical, intent(in) :: trailing
    integer :: i
    character(len=16) :: spacer

    spacer = '                '
    write(*,'(A)',advance='no') spacer(1:indent)//'"'//trim(name)//'": ['
    do i = 1, nvals
      if (i .gt. 1) write(*,'(A)',advance='no') ', '
      write(*,'(ES24.16E3)',advance='no') arr(i)
    end do
    if (trailing) then
      write(*,'(A)') '],'
    else
      write(*,'(A)') ']'
    end if
  end subroutine emit_real_array

  subroutine init_xfoil()
    include 'XFOIL.INC'
    logical :: ldbg
    integer :: ludbg, idbgcall, idbgiter
    common /XDEBUG/ ldbg, ludbg, idbgcall, idbgiter
    real :: pi_init

    ldbg = .false.
    ludbg = 6
    idbgcall = 0
    idbgiter = 0

    pi_init = 4.0d0 * atan(1.0d0)
    pi = pi_init
    hopi = 0.5d0 / pi_init
    qopi = 0.25d0 / pi_init
    dtor = pi_init / 180.0d0
    gamma = 1.4d0
    gamm1 = gamma - 1.0d0

    lvisc = .false.
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
    liqset = .false.
    sharp = .false.

    qinf = 1.0d0
    minf = 0.0d0
    minf1 = 0.0d0
    reinf = 1.0d6
    reinf1 = 1.0d6
    reinf_cl = 1.0d6
    xcmref = 0.25d0
    ycmref = 0.0d0
    waklen = 1.0d0
    acrit(1) = 9.0d0
    acrit(2) = 9.0d0
    xstrip(1) = 1.0d0
    xstrip(2) = 1.0d0
  end subroutine init_xfoil

  subroutine gen_naca0012(npanel)
    include 'XFOIL.INC'
    integer, intent(in) :: npanel
    integer :: i, nhalf
    real :: beta, xx, yt, t

    t = 0.12d0
    nhalf = npanel/2
    n = 0
    do i = nhalf, 0, -1
      beta = pi * dble(i) / dble(nhalf)
      xx = 0.5d0 * (1.0d0 - cos(beta))
      yt = 5.0d0*t*(0.2969d0*sqrt(xx) - 0.126d0*xx - 0.3516d0*xx**2 &
        + 0.2843d0*xx**3 - 0.1036d0*xx**4)
      n = n + 1
      x(n) = xx
      y(n) = yt
    end do
    do i = 1, nhalf
      beta = pi * dble(i) / dble(nhalf)
      xx = 0.5d0 * (1.0d0 - cos(beta))
      yt = 5.0d0*t*(0.2969d0*sqrt(xx) - 0.126d0*xx - 0.3516d0*xx**2 &
        + 0.2843d0*xx**3 - 0.1036d0*xx**4)
      n = n + 1
      x(n) = xx
      y(n) = -yt
    end do
    name = 'NACA 0012'
    call finalize_geometry()
  end subroutine gen_naca0012

  subroutine finalize_geometry()
    include 'XFOIL.INC'

    call scalc(x, y, s, n)
    call segspl(x, xp, s, n)
    call segspl(y, yp, s, n)
    call ncalc(x, y, s, n, nx, ny)
    call lefind(sle, x, xp, y, yp, s, n)
    xle = seval(sle, x, xp, s, n)
    yle = seval(sle, y, yp, s, n)
    xte = 0.5d0 * (x(1) + x(n))
    yte = 0.5d0 * (y(1) + y(n))
    chord = sqrt((xte - xle)**2 + (yte - yle)**2)
    call set_te_geometry()
    call apcalc
    lqinu = .false.
    lwake = .false.
    lqaij = .false.
    ladij = .false.
    lwdij = .false.
  end subroutine finalize_geometry

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
  end subroutine set_te_geometry

  subroutine splqsp(kqsp)
    include 'XFOIL.INC'
    integer, intent(in) :: kqsp
    integer :: i

    call splind(qspec(2, kqsp), qspecp(2, kqsp), sspec(2), nsp - 2, -999.0, -999.0)

    i = 1
    call splind(qspec(i, kqsp), qspecp(i, kqsp), sspec(i), 2, -999.0, qspecp(i + 1, kqsp))

    i = nsp - 1
    call splind(qspec(i, kqsp), qspecp(i, kqsp), sspec(i), 2, qspecp(i, kqsp), -999.0)
  end subroutine splqsp

  subroutine gamqsp(kqsp)
    include 'XFOIL.INC'
    integer, intent(in) :: kqsp
    integer :: i

    alqsp(kqsp) = algam
    clqsp(kqsp) = clgam
    cmqsp(kqsp) = cmgam

    do i = 1, nsp
      qspec(i, kqsp) = qgamm(i)
    end do

    qdof0 = 0.0d0
    qdof1 = 0.0d0
    qdof2 = 0.0d0
    qdof3 = 0.0d0

    call splqsp(kqsp)

    if (.not. liqset) then
      iq1 = 1
      iq2 = nsp
    end if
  end subroutine gamqsp

  subroutine mixed(kqsp, niterq)
    include 'XFOIL.INC'
    integer, intent(in) :: kqsp, niterq
    integer :: i, j, iter, inmax, igmax
    real :: bwt, cosa_loc, sina_loc, fs, res, xbis, ybis, qbis
    real :: ds1, ds2, dsmin, ag1, ag2, abis, cbis, sbis
    real :: dnmax, dgmax

    bwt = 0.1d0
    cosa_loc = cos(alfa)
    sina_loc = sin(alfa)
    call scalc(x, y, s, n)

    do i = 1, n
      qf0(i) = 0.0d0
      qf1(i) = 0.0d0
      qf2(i) = 0.0d0
      qf3(i) = 0.0d0
    end do

    do i = iq1, iq2
      fs = (s(i) - s(iq1)) / max(s(iq2) - s(iq1), 1.0d-12)
      qf0(i) = 1.0d0 - fs
      qf1(i) = fs
      qf2(i) = 0.0d0
      qf3(i) = 0.0d0
      gam(i) = qspec(i, kqsp) + qdof0 * qf0(i) + qdof1 * qf1(i)
    end do

    do iter = 1, niterq
      do i = 1, n + 5
        do j = 1, n + 5
          q(i, j) = 0.0d0
        end do
      end do

      call ncalc(x, y, s, n, nx, ny)

      do i = 1, n
        call psilin(i, x(i), y(i), nx(i), ny(i), psi, psi_ni, .true., .false.)
        dzdn(i) = dzdn(i) + psi_ni

        do j = 1, iq1 - 1
          q(i, j) = q(i, j) + dzdg(j)
        end do
        do j = iq1, iq2
          q(i, j) = q(i, j) + dzdn(j)
        end do
        do j = iq2 + 1, n
          q(i, j) = q(i, j) + dzdg(j)
        end do

        dq(i) = psio - psi
        q(i, n + 1) = q(i, n + 1) - 1.0d0
        q(i, n + 2) = q(i, n + 2) + z_qdof0
        q(i, n + 3) = q(i, n + 3) + z_qdof1
        q(i, n + 4) = q(i, n + 4) + z_qdof2
        q(i, n + 5) = q(i, n + 5) + z_qdof3
      end do

      dq(n + 1) = -(gam(1) + gam(n))
      call gamlin(n + 1, 1, 1.0)
      call gamlin(n + 1, n, 1.0)

      if (sharp) then
        ag1 = atan2(-yp(1), -xp(1))
        ag2 = atanc(yp(n), xp(n), ag1)
        abis = 0.5d0 * (ag1 + ag2)
        cbis = cos(abis)
        sbis = sin(abis)

        ds1 = sqrt((x(1) - x(2))**2 + (y(1) - y(2))**2)
        ds2 = sqrt((x(n) - x(n - 1))**2 + (y(n) - y(n - 1))**2)
        dsmin = min(ds1, ds2)

        xbis = xte - bwt * dsmin * cbis
        ybis = yte - bwt * dsmin * sbis
        call psilin(0, xbis, ybis, -sbis, cbis, psi, qbis, .false., .true.)
        res = qbis

        do j = 1, n + 5
          q(n, j) = 0.0d0
        end do
        do j = 1, n
          call gamlin(n, j, dqdg(j))
          q(n, j) = dqdg(j)
        end do
        dq(n) = -res
      end if

      q(n + 2, iq1) = 1.0d0
      dq(n + 2) = 0.0d0
      q(n + 3, iq2) = 1.0d0
      dq(n + 3) = 0.0d0
      q(n + 4, n + 4) = 1.0d0
      dq(n + 4) = -qdof2
      q(n + 5, n + 5) = 1.0d0
      dq(n + 5) = -qdof3

      call gauss(iqx, n + 5, q, dq, 1)

      inmax = 0
      igmax = 0
      dnmax = 0.0d0
      dgmax = 0.0d0

      do i = 1, iq1 - 1
        gam(i) = gam(i) + dq(i)
        if (abs(dq(i)) .gt. abs(dgmax)) then
          dgmax = dq(i)
          igmax = i
        end if
      end do

      do i = iq1, iq2
        x(i) = x(i) + nx(i) * dq(i)
        y(i) = y(i) + ny(i) * dq(i)
        if (abs(dq(i)) .gt. abs(dnmax)) then
          dnmax = dq(i)
          inmax = i
        end if
      end do

      do i = iq2 + 1, n
        gam(i) = gam(i) + dq(i)
        if (abs(dq(i)) .gt. abs(dgmax)) then
          dgmax = dq(i)
          igmax = i
        end if
      end do

      psio = psio + dq(n + 1)
      qdof0 = qdof0 + dq(n + 2)
      qdof1 = qdof1 + dq(n + 3)
      qdof2 = qdof2 + dq(n + 4)
      qdof3 = qdof3 + dq(n + 5)

      cosa_loc = cos(alfa)
      sina_loc = sin(alfa)
      call scalc(x, y, s, n)
      call segspl(x, xp, s, n)
      call segspl(y, yp, s, n)

      do i = iq1, iq2
        gam(i) = qspec(i, kqsp) + qdof0 * qf0(i) + qdof1 * qf1(i) + qdof2 * qf2(i) + qdof3 * qf3(i)
      end do

      call set_te_geometry()

      if (abs(dnmax) .lt. 5.0d-5 .and. abs(dgmax) .lt. 5.0d-4) then
        return
      end if
    end do
  end subroutine mixed

  subroutine gamlin(i, j, coef)
    include 'XFOIL.INC'
    integer, intent(in) :: i, j
    real, intent(in) :: coef

    if (j .ge. iq1 .and. j .le. iq2) then
      q(i, n + 2) = q(i, n + 2) + coef * qf0(j)
      q(i, n + 3) = q(i, n + 3) + coef * qf1(j)
      q(i, n + 4) = q(i, n + 4) + coef * qf2(j)
      q(i, n + 5) = q(i, n + 5) + coef * qf3(j)
    else
      q(i, j) = q(i, j) + coef
    end if
  end subroutine gamlin

end program qdes_driver
