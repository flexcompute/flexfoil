program qdcalc_driver
  include 'XFOIL.INC'
  real :: samples(5)

  call bench_qdcalc(samples)

  call init_xfoil()
  call gen_naca0012(40)
  lgamu = .false.
  call ggcalc
  lgamu = .true.
  alfa = 4.0d0 * atan(1.0d0) / 45.0d0
  call set_alpha_state()
  call xywake
  call qwcalc
  call qiset
  call qdcalc

  write(*,'(A)') '{'
  write(*,'(A,ES24.16E3,A)') '  "xle": ', xle, ','
  write(*,'(A,ES24.16E3,A)') '  "yle": ', yle, ','
  write(*,'(A,ES24.16E3,A)') '  "xte": ', xte, ','
  write(*,'(A,ES24.16E3,A)') '  "yte": ', yte, ','
  write(*,'(A,ES24.16E3,A)') '  "chord": ', chord, ','
  call emit_real_array('wake_x', x(n+1:n+nw), nw, .true., 2)
  call emit_real_array('wake_y', y(n+1:n+nw), nw, .true., 2)
  call emit_wake_state(.true.)
  call emit_dij_flat(.true.)
  write(*,'(A,I0,A)') '  "matrix_size": ', n + nw, ','
  write(*,'(A)') '  "perf": {'
  write(*,'(A,I0,A)') '    "inner_loops": ', 100, ','
  call emit_real_array('samples_seconds', samples, 5, .true., 4)
  call sort_real(samples, 5)
  write(*,'(A,ES24.16E3)') '    "median_seconds": ', samples(3)
  write(*,'(A)') '  }'
  write(*,'(A)') '}'

contains

  subroutine emit_wake_state(trailing)
    include 'XFOIL.INC'
    logical, intent(in) :: trailing
    real :: wake_s(iwx), wake_nx_out(iwx), wake_ny_out(iwx), wake_apanel(iwx)
    real :: wake_qinvu_0_out(iwx), wake_qinvu_90_out(iwx), wake_qinv_out(iwx), wake_qinv_a_out(iwx)

    wake_s(1:nw) = s(n+1:n+nw) - s(n)
    wake_nx_out(1:nw) = nx(n+1:n+nw)
    wake_ny_out(1:nw) = ny(n+1:n+nw)
    wake_apanel(1:nw) = 0.0
    if (nw .gt. 1) wake_apanel(1:nw-1) = apanel(n+1:n+nw-1)
    wake_qinvu_0_out(1:nw) = qinvu(n+1:n+nw,1)
    wake_qinvu_90_out(1:nw) = qinvu(n+1:n+nw,2)
    wake_qinv_out(1:nw) = qinv(n+1:n+nw)
    wake_qinv_a_out(1:nw) = qinv_a(n+1:n+nw)

    call emit_real_array('wake_s', wake_s, nw, .true., 2)
    call emit_real_array('wake_nx', wake_nx_out, nw, .true., 2)
    call emit_real_array('wake_ny', wake_ny_out, nw, .true., 2)
    call emit_real_array('wake_apanel', wake_apanel, nw, .true., 2)
    call emit_real_array('wake_qinvu_0', wake_qinvu_0_out, nw, .true., 2)
    call emit_real_array('wake_qinvu_90', wake_qinvu_90_out, nw, .true., 2)
    call emit_real_array('wake_qinv', wake_qinv_out, nw, .true., 2)
    call emit_real_array('wake_qinv_a', wake_qinv_a_out, nw, trailing, 2)
  end subroutine emit_wake_state

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

    pi_local = 4.0d0*atan(1.0d0)
    pi = pi_local
    hopi = 0.5d0 / pi_local
    qopi = 0.25d0 / pi_local
    dtor = pi_local / 180.0d0
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
    gam(1:iqx) = 0.0d0
    gam_a(1:iqx) = 0.0d0
    qinv(1:izx) = 0.0d0
    qinv_a(1:izx) = 0.0d0
    qinvu(1:izx,1) = 0.0d0
    qinvu(1:izx,2) = 0.0d0
    sig(1:izx) = 0.0d0
  end subroutine init_xfoil

  subroutine gen_naca0012(npanel)
    include 'XFOIL.INC'
    integer, intent(in) :: npanel
    integer :: i, nhalf
    real :: beta, xx, yt, pi, t

    pi = 4.0d0*atan(1.0d0)
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
    sigte_a = 0.5d0 * (gam_a(1) - gam_a(n)) * scs
    gamte_a = -0.5d0 * (gam_a(1) - gam_a(n)) * sds
  end subroutine set_te_geometry

  subroutine bench_qdcalc(samples)
    real, intent(out) :: samples(5)
    integer :: sample_idx, iter
    real :: t0, t1

    call init_xfoil()
    call gen_naca0012(40)
    lgamu = .false.
    call ggcalc
    lgamu = .true.
    alfa = 4.0d0 * atan(1.0d0) / 45.0d0
    call set_alpha_state()
    call xywake
    call qwcalc
    call qiset

    do sample_idx = 1, 5
      call cpu_time(t0)
      do iter = 1, 100
        ladij = .false.
        lwdij = .false.
        call qdcalc
      end do
      call cpu_time(t1)
      samples(sample_idx) = t1 - t0
    end do
  end subroutine bench_qdcalc

  subroutine emit_dij_flat(trailing)
    include 'XFOIL.INC'
    logical, intent(in) :: trailing
    real :: dij_flat((n+nw)*(n+nw))
    integer :: i, j, idx, nvals

    nvals = n + nw
    idx = 0
    do i = 1, nvals
      do j = 1, nvals
        idx = idx + 1
        dij_flat(idx) = dij(i,j)
      end do
    end do
    call emit_real_array('dij_flat', dij_flat, idx, trailing, 2)
  end subroutine emit_dij_flat

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

end program qdcalc_driver
