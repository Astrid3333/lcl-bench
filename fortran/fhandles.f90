module fhandles
  use kinds, only: dp
  use fjson, only: jval, jnum_req, jget
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_quiet_nan, ieee_value
  implicit none

  real(dp), save :: s_ep, s_sg, s_Dd, s_ma, s_re, s_dwa, s_dwb

contains

  subroutine load_params(name, e, p, np, ok, msg)
    character(len=*), intent(in) :: name
    type(jval), intent(in) :: e
    real(dp), intent(out) :: p(:)
    integer, intent(out) :: np
    logical, intent(out) :: ok
    character(len=*), intent(out) :: msg
    type(jval) :: pr
    integer :: k
    real(dp) :: nan

    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    do k = 1, size(p)
      p(k) = nan
    end do
    np = 0
    ok = .true.
    msg = ''
    pr = jget(e, 'params')

    select case (trim(name))
    case ('logistic')
      np = 1
      p(1) = jnum_req(pr, 'r')
    case ('henon')
      np = 2
      p(1) = jnum_req(pr, 'a')
      p(2) = jnum_req(pr, 'b')
    case ('clifford', 'tinkerbell')
      np = 4
      p(1) = jnum_req(pr, 'a')
      p(2) = jnum_req(pr, 'b')
      p(3) = jnum_req(pr, 'c')
      p(4) = jnum_req(pr, 'd')
    case ('lorenz')
      np = 3
      p(1) = jnum_req(pr, 'sigma')
      p(2) = jnum_req(pr, 'rho')
      p(3) = jnum_req(pr, 'beta')
    case ('rossler')
      np = 3
      p(1) = jnum_req(pr, 'a')
      p(2) = jnum_req(pr, 'b')
      p(3) = jnum_req(pr, 'c')
    case ('thomas')
      np = 1
      p(1) = jnum_req(pr, 'b')
    case ('aizawa')
      np = 6
      p(1) = jnum_req(pr, 'a')
      p(2) = jnum_req(pr, 'b')
      p(3) = jnum_req(pr, 'c')
      p(4) = jnum_req(pr, 'd')
      p(5) = jnum_req(pr, 'e')
      p(6) = jnum_req(pr, 'f')
    case default
      ok = .false.
      msg = 'no params mapping for ' // trim(name)
      np = 0
      return
    end select

    do k = 1, np
      if (ieee_is_nan(p(k))) then
        ok = .false.
        msg = 'missing or non-numeric param for ' // trim(name)
        np = 0
        return
      end if
    end do
  end subroutine load_params

  subroutine map_step(name, p, x, y)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: p(:)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:)
    real(dp) :: t

    select case (trim(name))
    case ('logistic')
      y(1) = (p(1) * x(1)) * (1.0_dp - x(1))
    case ('henon')
      y(1) = (1.0_dp - (p(1) * (x(1) * x(1)))) + x(2)
      y(2) = p(2) * x(1)
    case ('clifford')
      t = p(1) * x(2)
      y(1) = sin(t) + p(3) * cos(p(1) * x(1))
      t = p(2) * x(1)
      y(2) = sin(t) + p(4) * cos(p(2) * x(2))
    case ('tinkerbell')
      y(1) = ((x(1) * x(1) - x(2) * x(2)) + p(1) * x(1)) + p(2) * x(2)
      y(2) = ((2.0_dp * x(1)) * x(2) + p(3) * x(1)) + p(4) * x(2)
    case default
      y(1) = 0.0_dp
    end select
  end subroutine map_step

  subroutine rhs(name, p, x, y)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: p(:)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: y(:)

    select case (trim(name))
    case ('lorenz')
      y(1) = p(1) * (x(2) - x(1))
      y(2) = x(1) * (p(2) - x(3)) - x(2)
      y(3) = x(1) * x(2) - p(3) * x(3)
    case ('rossler')
      y(1) = (-x(2)) - x(3)
      y(2) = x(1) + p(1) * x(2)
      y(3) = p(2) + x(3) * (x(1) - p(3))
    case ('thomas')
      y(1) = sin(x(2)) - p(1) * x(1)
      y(2) = sin(x(3)) - p(1) * x(2)
      y(3) = sin(x(1)) - p(1) * x(3)
    case ('aizawa')
      y(1) = ((x(3) - p(2)) * x(1)) - p(4) * x(2)
      y(2) = p(4) * x(1) + ((x(3) - p(2)) * x(2))
      y(3) = (((p(3) + p(1) * x(3)) - ((x(3) * x(3) * x(3)) / 3.0_dp)) &
             - ((x(1) * x(1) + x(2) * x(2)) * (1.0_dp + p(5) * x(3)))) &
             + ((p(6) * x(3)) * (x(1) * x(1) * x(1)))
    case default
      y(1) = 0.0_dp
      y(2) = 0.0_dp
      y(3) = 0.0_dp
    end select
  end subroutine rhs

  subroutine load_pot(name, e, ok, msg)
    character(len=*), intent(in) :: name
    type(jval), intent(in) :: e
    logical, intent(out) :: ok
    character(len=*), intent(out) :: msg
    type(jval) :: pr
    real(dp) :: nan
    integer :: k
    real(dp) :: vals(3)

    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    s_ep = nan
    s_sg = nan
    s_Dd = nan
    s_ma = nan
    s_re = nan
    s_dwa = nan
    s_dwb = nan
    ok = .true.
    msg = ''
    pr = jget(e, 'params')
    select case (trim(name))
    case ('lennard_jones')
      s_ep = jnum_req(pr, 'eps')
      s_sg = jnum_req(pr, 'sigma')
      vals(1) = s_ep
      vals(2) = s_sg
      if (any(ieee_is_nan(vals(1:2)))) then
        ok = .false.
        msg = 'missing or non-numeric param for ' // trim(name)
      end if
    case ('morse')
      s_Dd = jnum_req(pr, 'De')
      s_ma = jnum_req(pr, 'a')
      s_re = jnum_req(pr, 're')
      vals(1) = s_Dd
      vals(2) = s_ma
      vals(3) = s_re
      if (any(ieee_is_nan(vals))) then
        ok = .false.
        msg = 'missing or non-numeric param for ' // trim(name)
      end if
    case ('double_well')
      s_dwa = jnum_req(pr, 'a')
      s_dwb = jnum_req(pr, 'b')
      vals(1) = s_dwa
      vals(2) = s_dwb
      if (any(ieee_is_nan(vals(1:2)))) then
        ok = .false.
        msg = 'missing or non-numeric param for ' // trim(name)
      end if
    case default
      ok = .false.
      msg = 'load_pot: unknown potential ' // trim(name)
    end select
    k = 0
  end subroutine load_pot

  function lj_v(r) result(v)
    real(dp), intent(in) :: r
    real(dp) :: v
    real(dp) :: y

    y = s_sg / r
    v = (4.0_dp * s_ep) * (y**12 - y**6)
  end function lj_v

  function morse_v(r) result(v)
    real(dp), intent(in) :: r
    real(dp) :: v
    real(dp) :: t

    t = 1.0_dp - exp(-s_ma * (r - s_re))
    v = (s_Dd * (t * t)) - s_Dd
  end function morse_v

  function dw_v(x) result(v)
    real(dp), intent(in) :: x
    real(dp) :: v

    v = s_dwa * (x**4) - s_dwb * (x**2)
  end function dw_v

end module fhandles
