program state_topology_driver

  write(*,'(A)') '{'
  call emit_stfind()
  write(*,'(A)') ','
  call emit_iblpan()
  write(*,'(A)') ','
  call emit_xicalc()
  write(*,'(A)') ','
  call emit_uicalc()
  write(*,'(A)') ','
  call emit_qvfue()
  write(*,'(A)') ','
  call emit_gamqv()
  write(*,'(A)') ','
  call emit_ueset()
  write(*,'(A)') ','
  call emit_dsset()
  write(*,'(A)') ','
  call emit_stmove()
  write(*,'(A)') ','
  call emit_perf()
  write(*,'(A)') '}'

contains

  subroutine emit_stfind()
    include 'XFOIL.INC'

    call init_raw_fixture()
    call stfind

    write(*,'(A)') '  "stfind": {'
    write(*,'(A,I0,A)') '    "ist": ', ist - 1, ','
    write(*,'(A,ES24.16E3,A)') '    "sst": ', sst, ','
    write(*,'(A,ES24.16E3,A)') '    "sst_go": ', sst_go, ','
    write(*,'(A,ES24.16E3)') '    "sst_gp": ', sst_gp
    write(*,'(A)') '  }'
  end subroutine emit_stfind

  subroutine emit_iblpan()
    include 'XFOIL.INC'
    integer :: upper_ipan(IVX)
    integer :: upper_len

    call init_raw_fixture()
    call stfind
    call iblpan
    call build_upper_plot_int(ipan(1,1), ipan(1,2), iblte(1), iblte(2), nw, upper_ipan, upper_len)

    write(*,'(A)') '  "iblpan": {'
    write(*,'(A,I0,A)') '    "nbl_upper": ', upper_len, ','
    write(*,'(A,I0,A)') '    "nbl_lower": ', nbl(2), ','
    write(*,'(A,I0,A)') '    "iblte_upper": ', iblte(1), ','
    write(*,'(A,I0,A)') '    "iblte_lower": ', iblte(2), ','
    call emit_int_array('ipan_upper', upper_ipan(2:upper_len) - 1, upper_len - 1, .true., 4)
    call emit_int_array('ipan_lower', ipan(2:nbl(2),2) - 1, nbl(2) - 1, .false., 4)
    write(*,'(A)') '  }'
  end subroutine emit_iblpan

  subroutine emit_xicalc()
    include 'XFOIL.INC'
    real :: upper_x(IVX)
    integer :: upper_len

    call init_raw_fixture()
    call stfind
    call iblpan
    call xicalc
    call build_upper_plot_real(xssi(1,1), xssi(1,2), iblte(1), iblte(2), nw, upper_x, upper_len)

    write(*,'(A)') '  "xicalc": {'
    call emit_real_array('upper_x', upper_x, upper_len, .true., 4)
    call emit_real_array('lower_x', xssi(1:nbl(2),2), nbl(2), .true., 4)
    call emit_real_array('wgap', wgap(1:nw), nw, .false., 4)
    write(*,'(A)') '  }'
  end subroutine emit_xicalc

  subroutine emit_uicalc()
    include 'XFOIL.INC'
    real :: upper_uinv(IVX)
    real :: upper_uinv_a(IVX)
    integer :: upper_len

    call init_raw_fixture()
    call stfind
    call iblpan
    call xicalc
    call uicalc
    call build_upper_plot_real(uinv(1,1), uinv(1,2), iblte(1), iblte(2), nw, upper_uinv, upper_len)
    call build_upper_plot_real(uinv_a(1,1), uinv_a(1,2), iblte(1), iblte(2), nw, upper_uinv_a, upper_len)

    write(*,'(A)') '  "uicalc": {'
    call emit_real_array('upper_uinv', upper_uinv, upper_len, .true., 4)
    call emit_real_array('lower_uinv', uinv(1:nbl(2),2), nbl(2), .true., 4)
    call emit_real_array('upper_uinv_a', upper_uinv_a, upper_len, .true., 4)
    call emit_real_array('lower_uinv_a', uinv_a(1:nbl(2),2), nbl(2), .false., 4)
    write(*,'(A)') '  }'
  end subroutine emit_uicalc

  subroutine emit_qvfue()
    include 'XFOIL.INC'

    call prepare_fixture()
    call qvfue

    write(*,'(A)') '  "qvfue": {'
    call emit_real_array('qvis', qvis(1:n+nw), n+nw, .false., 4)
    write(*,'(A)') '  }'
  end subroutine emit_qvfue

  subroutine emit_gamqv()
    include 'XFOIL.INC'

    call prepare_fixture()
    call qvfue
    call gamqv

    write(*,'(A)') '  "gamqv": {'
    call emit_real_array('gam', gam(1:n), n, .true., 4)
    call emit_real_array('gam_a', gam_a(1:n), n, .false., 4)
    write(*,'(A)') '  }'
  end subroutine emit_gamqv

  subroutine emit_ueset()
    include 'XFOIL.INC'
    real :: upper_uedg(IVX)
    integer :: upper_len

    call prepare_fixture()
    call ueset
    call build_upper_plot_real(uedg(1,1), uedg(1,2), iblte(1), iblte(2), nw, upper_uedg, upper_len)

    write(*,'(A)') '  "ueset": {'
    call emit_real_array('upper_uedg', upper_uedg, upper_len, .true., 4)
    call emit_real_array('lower_uedg', uedg(1:nbl(2),2), nbl(2), .false., 4)
    write(*,'(A)') '  }'
  end subroutine emit_ueset

  subroutine emit_dsset()
    include 'XFOIL.INC'
    real :: upper_dstr(IVX)
    integer :: upper_len

    call prepare_fixture()
    call ueset
    call dsset
    call build_upper_plot_real(dstr(1,1), dstr(1,2), iblte(1), iblte(2), nw, upper_dstr, upper_len)

    write(*,'(A)') '  "dsset": {'
    call emit_real_array('upper_dstr', upper_dstr, upper_len, .true., 4)
    call emit_real_array('lower_dstr', dstr(1:nbl(2),2), nbl(2), .false., 4)
    write(*,'(A)') '  }'
  end subroutine emit_dsset

  subroutine emit_stmove()
    include 'XFOIL.INC'
    integer :: istold, idif, is, ibl, i
    real :: ueps, dudx

    call prepare_fixture()
    gam(1) = 0.72d0
    gam(2) = 0.28d0
    gam(3) = -0.05d0
    gam(4) = -0.25d0
    gam(5) = -0.40d0
    gam(6) = -0.20d0

    istold = ist
    call stfind

    if (istold .eq. ist) then
      call xicalc
    else
      call iblpan
      call uicalc
      call xicalc

      if (ist .gt. istold) then
        idif = ist - istold

        itran(1) = itran(1) + idif
        itran(2) = itran(2) - idif

        do ibl = nbl(1), idif + 2, -1
          ctau(ibl,1) = ctau(ibl-idif,1)
          thet(ibl,1) = thet(ibl-idif,1)
          dstr(ibl,1) = dstr(ibl-idif,1)
          uedg(ibl,1) = uedg(ibl-idif,1)
        end do

        dudx = uedg(idif+2,1) / xssi(idif+2,1)
        do ibl = idif + 1, 2, -1
          ctau(ibl,1) = ctau(idif+2,1)
          thet(ibl,1) = thet(idif+2,1)
          dstr(ibl,1) = dstr(idif+2,1)
          uedg(ibl,1) = dudx * xssi(ibl,1)
        end do

        do ibl = 2, nbl(2)
          ctau(ibl,2) = ctau(ibl+idif,2)
          thet(ibl,2) = thet(ibl+idif,2)
          dstr(ibl,2) = dstr(ibl+idif,2)
          uedg(ibl,2) = uedg(ibl+idif,2)
        end do
      else
        idif = istold - ist

        itran(1) = itran(1) - idif
        itran(2) = itran(2) + idif

        do ibl = nbl(2), idif + 2, -1
          ctau(ibl,2) = ctau(ibl-idif,2)
          thet(ibl,2) = thet(ibl-idif,2)
          dstr(ibl,2) = dstr(ibl-idif,2)
          uedg(ibl,2) = uedg(ibl-idif,2)
        end do

        dudx = uedg(idif+2,2) / xssi(idif+2,2)
        do ibl = idif + 1, 2, -1
          ctau(ibl,2) = ctau(idif+2,2)
          thet(ibl,2) = thet(idif+2,2)
          dstr(ibl,2) = dstr(idif+2,2)
          uedg(ibl,2) = dudx * xssi(ibl,2)
        end do

        do ibl = 2, nbl(1)
          ctau(ibl,1) = ctau(ibl+idif,1)
          thet(ibl,1) = thet(ibl+idif,1)
          dstr(ibl,1) = dstr(ibl+idif,1)
          uedg(ibl,1) = uedg(ibl+idif,1)
        end do
      end if

      ueps = 1.0e-7
      do is = 1, 2
        do ibl = 2, nbl(is)
          i = ipan(ibl,is)
          if (uedg(ibl,is) .le. ueps) then
            uedg(ibl,is) = ueps
            qvis(i) = vti(ibl,is) * ueps
            gam(i) = vti(ibl,is) * ueps
          end if
        end do
      end do
    end if

    do is = 1, 2
      do ibl = 2, nbl(is)
        mass(ibl,is) = dstr(ibl,is) * uedg(ibl,is)
      end do
    end do

    write(*,'(A)') '  "stmove": {'
    write(*,'(A,I0,A)') '    "ist": ', ist - 1, ','
    write(*,'(A,I0,A)') '    "nbl_upper": ', nbl(1), ','
    write(*,'(A,I0,A)') '    "nbl_lower": ', nbl(2), ','
    call emit_int_array('ipan_upper', ipan(2:nbl(1),1) - 1, nbl(1) - 1, .true., 4)
    call emit_int_array('ipan_lower', ipan(2:nbl(2),2) - 1, nbl(2) - 1, .true., 4)
    call emit_real_array('upper_x', xssi(1:nbl(1),1), nbl(1), .true., 4)
    call emit_real_array('lower_x', xssi(1:nbl(2),2), nbl(2), .true., 4)
    call emit_real_array('upper_theta', thet(1:nbl(1),1), nbl(1), .true., 4)
    call emit_real_array('lower_theta', thet(1:nbl(2),2), nbl(2), .true., 4)
    call emit_real_array('upper_dstr', dstr(1:nbl(1),1), nbl(1), .true., 4)
    call emit_real_array('lower_dstr', dstr(1:nbl(2),2), nbl(2), .true., 4)
    call emit_real_array('upper_uedg', uedg(1:nbl(1),1), nbl(1), .true., 4)
    call emit_real_array('lower_uedg', uedg(1:nbl(2),2), nbl(2), .true., 4)
    call emit_real_array('upper_mass', mass(1:nbl(1),1), nbl(1), .true., 4)
    call emit_real_array('lower_mass', mass(1:nbl(2),2), nbl(2), .false., 4)
    write(*,'(A)') '  }'
  end subroutine emit_stmove

  subroutine emit_perf()
    real :: samples(5)

    write(*,'(A)') '  "perf": {'

    call bench_stfind(samples)
    call emit_perf_case('stfind', 5000, samples, .true.)
    call bench_iblpan(samples)
    call emit_perf_case('iblpan', 5000, samples, .true.)
    call bench_xicalc(samples)
    call emit_perf_case('xicalc', 5000, samples, .true.)
    call bench_uicalc(samples)
    call emit_perf_case('uicalc', 5000, samples, .true.)
    call bench_qvfue(samples)
    call emit_perf_case('qvfue', 5000, samples, .true.)
    call bench_gamqv(samples)
    call emit_perf_case('gamqv', 5000, samples, .true.)
    call bench_ueset(samples)
    call emit_perf_case('ueset', 5000, samples, .true.)
    call bench_dsset(samples)
    call emit_perf_case('dsset', 5000, samples, .false.)

    write(*,'(A)') '  }'
  end subroutine emit_perf

  subroutine bench_stfind(samples)
    real, intent(out) :: samples(5)
    integer :: sample_idx, iter
    real :: t0, t1

    do sample_idx = 1, 5
      call cpu_time(t0)
      do iter = 1, 5000
        call init_raw_fixture()
        call stfind
      end do
      call cpu_time(t1)
      samples(sample_idx) = t1 - t0
    end do
  end subroutine bench_stfind

  subroutine bench_iblpan(samples)
    real, intent(out) :: samples(5)
    integer :: sample_idx, iter
    real :: t0, t1

    do sample_idx = 1, 5
      call cpu_time(t0)
      do iter = 1, 5000
        call init_raw_fixture()
        call stfind
        call iblpan
      end do
      call cpu_time(t1)
      samples(sample_idx) = t1 - t0
    end do
  end subroutine bench_iblpan

  subroutine bench_xicalc(samples)
    real, intent(out) :: samples(5)
    integer :: sample_idx, iter
    real :: t0, t1

    do sample_idx = 1, 5
      call cpu_time(t0)
      do iter = 1, 5000
        call init_raw_fixture()
        call stfind
        call iblpan
        call xicalc
      end do
      call cpu_time(t1)
      samples(sample_idx) = t1 - t0
    end do
  end subroutine bench_xicalc

  subroutine bench_uicalc(samples)
    real, intent(out) :: samples(5)
    integer :: sample_idx, iter
    real :: t0, t1

    do sample_idx = 1, 5
      call cpu_time(t0)
      do iter = 1, 5000
        call init_raw_fixture()
        call stfind
        call iblpan
        call xicalc
        call uicalc
      end do
      call cpu_time(t1)
      samples(sample_idx) = t1 - t0
    end do
  end subroutine bench_uicalc

  subroutine bench_qvfue(samples)
    real, intent(out) :: samples(5)
    integer :: sample_idx, iter
    real :: t0, t1

    do sample_idx = 1, 5
      call cpu_time(t0)
      do iter = 1, 5000
        call prepare_fixture()
        call qvfue
      end do
      call cpu_time(t1)
      samples(sample_idx) = t1 - t0
    end do
  end subroutine bench_qvfue

  subroutine bench_gamqv(samples)
    real, intent(out) :: samples(5)
    integer :: sample_idx, iter
    real :: t0, t1

    do sample_idx = 1, 5
      call cpu_time(t0)
      do iter = 1, 5000
        call prepare_fixture()
        call qvfue
        call gamqv
      end do
      call cpu_time(t1)
      samples(sample_idx) = t1 - t0
    end do
  end subroutine bench_gamqv

  subroutine bench_ueset(samples)
    real, intent(out) :: samples(5)
    integer :: sample_idx, iter
    real :: t0, t1

    do sample_idx = 1, 5
      call cpu_time(t0)
      do iter = 1, 5000
        call prepare_fixture()
        call ueset
      end do
      call cpu_time(t1)
      samples(sample_idx) = t1 - t0
    end do
  end subroutine bench_ueset

  subroutine bench_dsset(samples)
    real, intent(out) :: samples(5)
    integer :: sample_idx, iter
    real :: t0, t1

    do sample_idx = 1, 5
      call cpu_time(t0)
      do iter = 1, 5000
        call prepare_fixture()
        call ueset
        call dsset
      end do
      call cpu_time(t1)
      samples(sample_idx) = t1 - t0
    end do
  end subroutine bench_dsset

  subroutine emit_perf_case(name, inner_loops, samples, trailing)
    character(len=*), intent(in) :: name
    integer, intent(in) :: inner_loops
    real, intent(in) :: samples(5)
    logical, intent(in) :: trailing
    real :: sorted(5)
    integer :: i

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

  subroutine init_raw_fixture()
    include 'XFOIL.INC'
    logical :: ldbg
    integer :: ludbg, idbgcall, idbgiter
    common /XDEBUG/ ldbg, ludbg, idbgcall, idbgiter
    integer :: i, j

    ldbg = .false.
    ludbg = 6
    idbgcall = 0
    idbgiter = 0

    sharp = .true.
    lipan = .false.
    lblini = .false.
    lwake = .true.
    lwdij = .true.
    qinf = 1.0d0
    alfa = 0.0d0
    cosa = 1.0d0
    sina = 0.0d0
    ante = 0.02d0
    nw = 2
    n = 6

    x(1) = 1.00d0
    y(1) = 0.00d0
    x(2) = 0.75d0
    y(2) = 0.08d0
    x(3) = 0.40d0
    y(3) = 0.12d0
    x(4) = 0.00d0
    y(4) = 0.00d0
    x(5) = 0.35d0
    y(5) = -0.07d0
    x(6) = 0.85d0
    y(6) = -0.02d0
    x(7) = 1.05d0
    y(7) = 0.00d0
    x(8) = 1.25d0
    y(8) = 0.00d0

    s(1) = 0.0d0
    do i = 2, n
      s(i) = s(i-1) + sqrt((x(i) - x(i-1))**2 + (y(i) - y(i-1))**2)
    end do
    s(n+1) = s(n)
    s(n+2) = s(n+1) + sqrt((x(n+2) - x(n+1))**2 + (y(n+2) - y(n+1))**2)

    xp(1) = 1.0d0
    yp(1) = 0.1d0
    xp(n) = 1.0d0
    yp(n) = -0.1d0

    gam(1) = 0.72d0
    gam(2) = 0.38d0
    gam(3) = 0.12d0
    gam(4) = -0.20d0
    gam(5) = -0.55d0
    gam(6) = -0.25d0

    qinv(1) = 0.72d0
    qinv(2) = 0.38d0
    qinv(3) = 0.12d0
    qinv(4) = -0.20d0
    qinv(5) = -0.55d0
    qinv(6) = -0.25d0
    qinv(7) = 0.22d0
    qinv(8) = 0.18d0

    qinv_a(1) = 0.10d0
    qinv_a(2) = 0.08d0
    qinv_a(3) = 0.03d0
    qinv_a(4) = -0.04d0
    qinv_a(5) = -0.07d0
    qinv_a(6) = -0.02d0
    qinv_a(7) = 0.01d0
    qinv_a(8) = 0.01d0

    qvis(1:n+nw) = qinv(1:n+nw)
    gam_a(1:n) = qinv_a(1:n)

    do i = 1, n + nw
      do j = 1, n + nw
        dij(i,j) = 5.0d-4 * dble(i) - 3.0d-4 * dble(j)
      end do
    end do
  end subroutine init_raw_fixture

  subroutine prepare_fixture()
    include 'XFOIL.INC'
    integer :: is, ibl

    call init_raw_fixture()
    call stfind
    call iblpan
    call xicalc
    call uicalc

    itran(1) = 2
    itran(2) = 2

    do is = 1, 2
      ctau(1,is) = 0.0d0
      thet(1,is) = 0.0d0
      dstr(1,is) = 0.0d0
      mass(1,is) = 0.0d0
      uedg(1,is) = 0.0d0
      do ibl = 2, nbl(is)
        ctau(ibl,is) = 0.03d0
        thet(ibl,is) = 0.008d0 + 0.002d0*dble(ibl) + 0.0005d0*dble(is)
        dstr(ibl,is) = 0.012d0 + 0.003d0*dble(ibl) + 0.0008d0*dble(is)
        mass(ibl,is) = 0.025d0 + 0.004d0*dble(ibl) + 0.001d0*dble(is)
        uedg(ibl,is) = uinv(ibl,is) + 0.03d0*dble(ibl)
      end do
    end do
  end subroutine prepare_fixture

  subroutine build_upper_plot_real(upper_airfoil, lower_full, iblte_upper, iblte_lower, nwake, output, out_len)
    real, intent(in) :: upper_airfoil(*), lower_full(*)
    integer, intent(in) :: iblte_upper, iblte_lower, nwake
    real, intent(out) :: output(*)
    integer, intent(out) :: out_len
    integer :: i

    out_len = iblte_upper + nwake
    do i = 1, iblte_upper
      output(i) = upper_airfoil(i)
    end do
    do i = 1, nwake
      output(iblte_upper + i) = lower_full(iblte_lower + i)
    end do
  end subroutine build_upper_plot_real

  subroutine build_upper_plot_int(upper_airfoil, lower_full, iblte_upper, iblte_lower, nwake, output, out_len)
    integer, intent(in) :: upper_airfoil(*), lower_full(*)
    integer, intent(in) :: iblte_upper, iblte_lower, nwake
    integer, intent(out) :: output(*)
    integer, intent(out) :: out_len
    integer :: i

    out_len = iblte_upper + nwake
    do i = 1, iblte_upper
      output(i) = upper_airfoil(i)
    end do
    do i = 1, nwake
      output(iblte_upper + i) = lower_full(iblte_lower + i)
    end do
  end subroutine build_upper_plot_int

  subroutine emit_real_array(name, arr, nvals, trailing, indent)
    character(len=*), intent(in) :: name
    integer, intent(in) :: nvals, indent
    real, intent(in) :: arr(nvals)
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

  subroutine emit_int_array(name, arr, nvals, trailing, indent)
    character(len=*), intent(in) :: name
    integer, intent(in) :: nvals, indent
    integer, intent(in) :: arr(nvals)
    logical, intent(in) :: trailing
    integer :: i

    write(*,'(A)',advance='no') repeat(' ', indent) // '"' // trim(name) // '": ['
    do i = 1, nvals
      if (i > 1) write(*,'(A)',advance='no') ', '
      write(*,'(I0)',advance='no') arr(i)
    end do
    if (trailing) then
      write(*,'(A)') '],'
    else
      write(*,'(A)') ']'
    end if
  end subroutine emit_int_array

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

end program state_topology_driver
