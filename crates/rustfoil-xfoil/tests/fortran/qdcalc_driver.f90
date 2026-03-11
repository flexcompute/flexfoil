program qdcalc_driver
  include 'XFOIL.INC'
  real :: samples(5)

  call init_xfoil()
  call gen_naca0012(40)
  call pangen(.false.)
  lgamu = .false.
  call ggcalc
  lgamu = .true.
  alfa = 4.0d0 * atan(1.0d0) / 45.0d0
  call xywake
  call qwcalc
  call qiset
  call qdcalc

  write(*,'(A)') '{'
  call emit_real_array('wake_x', x(n+1:n+nw), nw, .true., 2)
  call emit_real_array('wake_y', y(n+1:n+nw), nw, .true., 2)
  call emit_diag_sample(.true.)
  call emit_row0_sample(.true.)
  write(*,'(A,I0,A)') '  "matrix_size": ', n + nw, ','
  call bench_qdcalc(samples)
  write(*,'(A)') '  "perf": {'
  write(*,'(A,I0,A)') '    "inner_loops": ', 100, ','
  call emit_real_array('samples_seconds', samples, 5, .true., 4)
  call sort_real(samples, 5)
  write(*,'(A,ES24.16E3)') '    "median_seconds": ', samples(3)
  write(*,'(A)') '  }'
  write(*,'(A)') '}'

contains

  subroutine init_xfoil()
    include 'XFOIL.INC'
    logical :: ldbg
    integer :: ludbg, idbgcall, idbgiter
    common /XDEBUG/ ldbg, ludbg, idbgcall, idbgiter

    ldbg = .false.
    ludbg = 6
    idbgcall = 0
    idbgiter = 0

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
    qinf = 1.0d0
    minf = 0.0d0
    minf1 = 0.0d0
    alfa = 0.0d0
    cosa = 1.0d0
    sina = 0.0d0
    xcmref = 0.25d0
    ycmref = 0.0d0
    waklen = 1.0d0
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
    do i = 0, nhalf
      beta = pi * dble(i) / dble(nhalf)
      xx = 0.5d0 * (1.0d0 - cos(beta))
      yt = 5.0d0*t*(0.2969d0*sqrt(xx) - 0.126d0*xx - 0.3516d0*xx**2 &
        + 0.2843d0*xx**3 - 0.1036d0*xx**4)
      n = n + 1
      x(n) = xx
      y(n) = yt
    end do
    do i = 1, nhalf
      beta = pi * dble(nhalf - i) / dble(nhalf)
      xx = 0.5d0 * (1.0d0 - cos(beta))
      yt = 5.0d0*t*(0.2969d0*sqrt(xx) - 0.126d0*xx - 0.3516d0*xx**2 &
        + 0.2843d0*xx**3 - 0.1036d0*xx**4)
      n = n + 1
      x(n) = xx
      y(n) = -yt
    end do
    name = 'NACA 0012'
  end subroutine gen_naca0012

  subroutine bench_qdcalc(samples)
    real, intent(out) :: samples(5)
    integer :: sample_idx, iter
    real :: t0, t1

    do sample_idx = 1, 5
      call cpu_time(t0)
      do iter = 1, 100
        call init_xfoil()
        call gen_naca0012(40)
        call pangen(.false.)
        lgamu = .false.
        call ggcalc
        lgamu = .true.
        alfa = 4.0d0 * atan(1.0d0) / 45.0d0
        call xywake
        call qwcalc
        call qiset
        call qdcalc
      end do
      call cpu_time(t1)
      samples(sample_idx) = t1 - t0
    end do
  end subroutine bench_qdcalc

  subroutine emit_diag_sample(trailing)
    include 'XFOIL.INC'
    logical, intent(in) :: trailing
    real :: diag_sample(20)
    integer :: i, nvals

    nvals = min(n + nw, 20)
    do i = 1, nvals
      diag_sample(i) = dij(i,i)
    end do
    call emit_real_array('diag_sample', diag_sample, nvals, trailing, 2)
  end subroutine emit_diag_sample

  subroutine emit_row0_sample(trailing)
    include 'XFOIL.INC'
    logical, intent(in) :: trailing
    real :: row_sample(20)
    integer :: i, nvals

    nvals = min(n + nw, 20)
    do i = 1, nvals
      row_sample(i) = dij(1,i)
    end do
    call emit_real_array('row0_sample', row_sample, nvals, trailing, 2)
  end subroutine emit_row0_sample

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
