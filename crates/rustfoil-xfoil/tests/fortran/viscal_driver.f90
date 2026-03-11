program viscal_driver
  include 'XFOIL.INC'
  character(len=32) :: arg
  real :: alpha_deg, pi

  pi = 4.0d0 * atan(1.0d0)
  alpha_deg = 4.0d0
  call getarg(1, arg)
  if (len_trim(arg) > 0) read(arg, *) alpha_deg

  call init_xfoil()
  call gen_naca0012(160)
  call pangen(.false.)
  lgamu = .false.
  call ggcalc
  lgamu = .true.

  alfa = alpha_deg * pi / 180.0d0
  cosa = cos(alfa)
  sina = sin(alfa)
  call viscal(10)

  write(*,'(A)') '{'
  write(*,'(A)') '  "cases": ['
  write(*,'(A)') '    {'
  write(*,'(A,F10.4,A)') '      "alpha_deg": ', alpha_deg, ','
  if (lvconv) then
    write(*,'(A)') '      "converged": true,'
  else
    write(*,'(A)') '      "converged": false,'
  end if
  write(*,'(A,I0,A)') '      "iterations": ', 10, ','
  write(*,'(A,ES24.16E3,A)') '      "cl": ', cl, ','
  write(*,'(A,ES24.16E3,A)') '      "cd": ', cd, ','
  write(*,'(A,ES24.16E3,A)') '      "cm": ', cm, ','
  write(*,'(A,ES24.16E3,A)') '      "x_tr_upper": ', xoctr(1), ','
  write(*,'(A,ES24.16E3)') '      "x_tr_lower": ', xoctr(2)
  write(*,'(A)') '    }'
  write(*,'(A)') '  ]'
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

    qinf = 1.0d0
    minf = 0.0d0
    minf1 = 0.0d0
    reinf = 1.0d6
    reinf1 = 1.0d6
    reinf_cl = 1.0d6
    xcmref = 0.25d0
    ycmref = 0.0d0
    acrit(1) = 9.0d0
    acrit(2) = 9.0d0
    xstrip(1) = 1.0d0
    xstrip(2) = 1.0d0
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

end program viscal_driver
