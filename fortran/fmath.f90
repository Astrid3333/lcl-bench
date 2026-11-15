module fmath
  use kinds, only: dp, i8
  implicit none
  private

  real(dp), parameter :: PI = 3.141592653589793238462643383279502884_dp
  real(dp), parameter :: LN2 = 0.6931471805599453_dp

  public :: pi_value, cordic_gain, cordic_rot
  public :: hastings_sin, hastings_atan, hastings_exp, hastings_log1p
  public :: lcg_next, lcg_normals, lcg_period
  public :: egyptian_mult, egyptian_fraction, igcd
  public :: babylonian_sqrt, aryabhata_sin
  public :: ext_gcd, chinese_remainder, chinese_elimination, vigesimal
  public :: ms_step, wh_step, wh_step16

contains

  pure function pi_value() result(p)
    real(dp) :: p
    p = PI
  end function

  pure function cordic_gain(n) result(k)
    integer, intent(in) :: n
    real(dp) :: k
    integer :: i
    k = 1.0_dp
    do i = 0, n - 1
      k = k * sqrt(1.0_dp + 2.0_dp**(-2 * i))
    end do
  end function

  subroutine cordic_rot(theta, atan_tab, x, y, zres)
    real(dp), intent(in) :: theta
    real(dp), intent(in) :: atan_tab(:)
    real(dp), intent(out) :: x, y, zres
    real(dp) :: t, z, sigma, f, xn, yn, sign_cos
    integer :: i
    x = 1.0_dp
    y = 0.0_dp
    t = theta
    do while (t > PI)
      t = t - 2.0_dp * PI
    end do
    do while (t < -PI)
      t = t + 2.0_dp * PI
    end do
    sign_cos = 1.0_dp
    if (t > PI / 2.0_dp) then
      t = PI - t
      sign_cos = -1.0_dp
    else if (t < -PI / 2.0_dp) then
      t = -PI - t
      sign_cos = -1.0_dp
    end if
    z = t
    do i = 1, size(atan_tab)
      if (z >= 0.0_dp) then
        sigma = 1.0_dp
      else
        sigma = -1.0_dp
      end if
      f = 2.0_dp**(-(i - 1))
      xn = x - sigma * y * f
      yn = y + sigma * x * f
      z = z - sigma * atan_tab(i)
      x = xn
      y = yn
    end do
    x = x * sign_cos
    zres = z
  end subroutine

  pure function hpoly(c, xx) result(r)
    real(dp), intent(in) :: c(:), xx
    real(dp) :: r
    integer :: i
    r = c(size(c))
    do i = size(c) - 1, 1, -1
      r = r * xx + c(i)
    end do
  end function

  pure function hastings_sin(x) result(v)
    real(dp), intent(in) :: x
    real(dp) :: v, s, xx, kk
    real(dp), parameter :: c(6) = [1.0_dp, -1.0_dp/6.0_dp, 1.0_dp/120.0_dp, &
         -1.0_dp/5040.0_dp, 1.0_dp/362880.0_dp, -1.0_dp/39916800.0_dp]
    s = 1.0_dp
    xx = x
    if (xx < 0.0_dp) then
      s = -1.0_dp
      xx = -xx
    end if
    kk = floor(xx / (2.0_dp * PI))
    xx = xx - kk * 2.0_dp * PI
    if (xx > PI) then
      xx = PI - xx
    end if
    if (xx > PI / 2.0_dp) then
      xx = PI - xx
    end if
    v = s * xx * hpoly(c, xx * xx)
  end function

  pure function hastings_atan(x) result(v)
    real(dp), intent(in) :: x
    real(dp) :: v, s, xx, t, inv
    real(dp), parameter :: c(6) = [1.0_dp, -1.0_dp/3.0_dp, 1.0_dp/5.0_dp, &
         -1.0_dp/7.0_dp, 1.0_dp/9.0_dp, -1.0_dp/11.0_dp]
    s = 1.0_dp
    xx = x
    if (xx < 0.0_dp) then
      s = -1.0_dp
      xx = -xx
    end if
    if (xx >= 0.5_dp .and. xx <= 2.0_dp) then
      t = (xx - 1.0_dp) / (xx + 1.0_dp)
      v = s * (PI / 4.0_dp + t * hpoly(c, t * t))
    else if (xx > 1.0_dp) then
      inv = 1.0_dp / xx
      v = s * (PI / 2.0_dp - inv * hpoly(c, inv * inv))
    else
      v = s * xx * hpoly(c, xx * xx)
    end if
  end function

  pure function hastings_exp(x) result(v)
    real(dp), intent(in) :: x
    real(dp) :: v
    if (x < 0.0_dp) then
      v = 1.0_dp / exp_pos(-x)
    else
      v = exp_pos(x)
    end if
  end function

  pure function exp_pos(x) result(v)
    real(dp), intent(in) :: x
    real(dp) :: v, r, y
    integer :: k
    real(dp), parameter :: c(6) = [1.0_dp, 1.0_dp/2.0_dp, 1.0_dp/6.0_dp, &
         1.0_dp/24.0_dp, 1.0_dp/120.0_dp, 1.0_dp/720.0_dp]
    if (x > 1.0_dp) then
      k = int(floor(x / LN2))
      r = x - k * LN2
      v = 2.0_dp**k * (1.0_dp + r * hpoly(c, r))
    else
      k = 0
      r = x
      do while (r > 0.25_dp)
        r = r / 2.0_dp
        k = k + 1
      end do
      y = 1.0_dp + r * hpoly(c, r)
      do while (k > 0)
        y = y * y
        k = k - 1
      end do
      v = y
    end if
  end function

  pure function hastings_log1p(x) result(v)
    real(dp), intent(in) :: x
    real(dp) :: v, t
    real(dp), parameter :: c(10) = [1.0_dp, -1.0_dp/2.0_dp, 1.0_dp/3.0_dp, &
         -1.0_dp/4.0_dp, 1.0_dp/5.0_dp, -1.0_dp/6.0_dp, 1.0_dp/7.0_dp, &
         -1.0_dp/8.0_dp, 1.0_dp/9.0_dp, -1.0_dp/10.0_dp]
    if (x > 1.0_dp) then
      t = (x - 1.0_dp) / (x + 1.0_dp)
      v = LN2 + t * hpoly(c, t)
    else if (x >= 0.25_dp) then
      t = (x - 1.0_dp) / 2.0_dp
      v = LN2 + t * hpoly(c, t)
    else
      v = x * hpoly(c, x)
    end if
  end function

  pure function lcg_next(a, c, m, x) result(y)
    integer(i8), intent(in) :: a, c, m, x
    integer(i8) :: y
    y = modulo(a * x + c, m)
  end function

  subroutine lcg_normals(a, m, box, count, minnorm)
    integer(i8), intent(in) :: a, m
    integer, intent(in) :: box
    integer, intent(out) :: count
    real(dp), intent(out) :: minnorm
    integer(i8) :: a2, h1, h2, h3
    integer :: i2, i3
    real(dp) :: nrm
    a2 = modulo(a * a, m)
    count = 0
    minnorm = huge(1.0_dp)
    do i3 = -box, box
      h3 = int(i3, i8)
      do i2 = -box, box
        h2 = int(i2, i8)
        h1 = modulo(-(h2 * a + h3 * a2), m)
        if (2 * h1 > m) h1 = h1 - m
        if (abs(h1) <= box) then
          if (.not. (h1 == 0 .and. h2 == 0 .and. h3 == 0)) then
            count = count + 1
            nrm = sqrt(real(h1 * h1 + h2 * h2 + h3 * h3, dp))
            if (nrm < minnorm) minnorm = nrm
          end if
        end if
      end do
    end do
  end subroutine

  function lcg_period(a, c, m, x0) result(p)
    integer(i8), intent(in) :: a, c, m, x0
    integer :: p
    integer(i8) :: x, xstart
    logical, allocatable :: seen(:)
    xstart = modulo(x0, m)
    x = xstart
    p = 0
    allocate(seen(int(m)))
    seen = .false.
    do while (p < int(m))
      if (seen(int(x) + 1)) exit
      seen(int(x) + 1) = .true.
      x = modulo(a * x + c, m)
      p = p + 1
    end do
  end function

  pure function egyptian_mult(a, b) result(res)
    integer(i8), intent(in) :: a, b
    integer(i8) :: res, aa, bb
    logical :: neg
    aa = a
    bb = b
    neg = .false.
    if (aa < 0_i8) then
      aa = -aa
      neg = .not. neg
    end if
    bb = b
    if (bb < 0_i8) then
      bb = -bb
      neg = .not. neg
    end if
    res = 0_i8
    do while (bb > 0_i8)
      if (modulo(bb, 2_i8) == 1_i8) res = res + aa
      aa = 2_i8 * aa
      bb = bb / 2_i8
    end do
    if (neg) res = -res
  end function

  pure function igcd(a, b) result(g)
    integer(i8), intent(in) :: a, b
    integer(i8) :: g, x, y, tmp
    x = abs(a)
    y = abs(b)
    do while (y /= 0_i8)
      tmp = x
      x = y
      y = modulo(tmp, y)
    end do
    g = x
  end function

  subroutine egyptian_fraction(num, den, max_terms, dens, nd)
    integer(i8), intent(in) :: num, den
    integer, intent(in) :: max_terms
    integer(i8), intent(out) :: dens(:)
    integer, intent(out) :: nd
    integer(i8) :: n, d, m, nn, nd2, g
    n = num
    d = den
    nd = 0
    do while (n > 0_i8 .and. nd < max_terms)
      m = (d + n - 1_i8) / n
      nd = nd + 1
      dens(nd) = m
      nn = n * m - d
      nd2 = d * m
      g = igcd(nn, nd2)
      if (g > 0_i8) then
        n = nn / g
        d = nd2 / g
      else
        n = 0_i8
      end if
    end do
  end subroutine

  pure function babylonian_sqrt(x, iterations) result(v)
    real(dp), intent(in) :: x
    integer, intent(in) :: iterations
    real(dp) :: v
    integer :: i
    if (x < 0.0_dp) then
      v = -1.0_dp
      return
    end if
    if (x == 0.0_dp) then
      v = 0.0_dp
      return
    end if
    v = x / 2.0_dp
    if (v == 0.0_dp) v = 1.0_dp
    do i = 1, iterations
      v = 0.5_dp * (v + x / v)
    end do
  end function

  pure function aryabhata_sin(theta, pasos) result(v)
    real(dp), intent(in) :: theta
    integer, intent(in) :: pasos
    real(dp) :: v, th, paso, kk, s_prev, s_cur, s_next
    integer :: n, idx, i
    th = theta
    do while (th < 0.0_dp)
      th = th + 2.0_dp * PI
    end do
    do while (th > 2.0_dp * PI)
      th = th - 2.0_dp * PI
    end do
    n = pasos
    if (n <= 0) n = 360
    paso = 2.0_dp * PI / real(n, dp)
    idx = int(th / paso)
    if (idx < 0) idx = 0
    if (idx > n) idx = n
    kk = 2.0_dp * (1.0_dp - cos(paso))
    s_prev = 0.0_dp
    s_cur = sin(paso)
    if (idx == 0) then
      v = 0.0_dp
      return
    end if
    if (idx == 1) then
      v = s_cur
      return
    end if
    do i = 1, idx - 1
      s_next = 2.0_dp * s_cur - s_prev - kk * s_cur
      s_prev = s_cur
      s_cur = s_next
    end do
    v = s_cur
  end function

  subroutine ext_gcd(a, b, g, x, y)
    integer(i8), intent(in) :: a, b
    integer(i8), intent(out) :: g, x, y
    integer(i8) :: old_r, r, old_s, s, old_t, t, q, tmp
    old_r = a
    r = b
    old_s = 1_i8
    s = 0_i8
    old_t = 0_i8
    t = 1_i8
    do while (r /= 0_i8)
      q = old_r / r
      if (modulo(old_r, r) /= 0_i8) then
        if ((old_r < 0_i8) .neqv. (r < 0_i8)) q = q - 1_i8
      end if
      tmp = old_r - q * r
      old_r = r
      r = tmp
      tmp = old_s - q * s
      old_s = s
      s = tmp
      tmp = old_t - q * t
      old_t = t
      t = tmp
    end do
    g = old_r
    x = old_s
    y = old_t
  end subroutine

  subroutine chinese_remainder(residues, moduli, x, mm, ok)
    integer(i8), intent(in) :: residues(:), moduli(:)
    integer(i8), intent(out) :: x, mm
    logical, intent(out) :: ok
    integer :: i, n
    integer(i8) :: mi, g, inv, yy
    n = size(moduli)
    mm = 1_i8
    do i = 1, n
      mm = mm * moduli(i)
    end do
    x = 0_i8
    ok = .true.
    do i = 1, n
      mi = mm / moduli(i)
      call ext_gcd(modulo(mi, moduli(i)), moduli(i), g, inv, yy)
      if (g /= 1_i8) then
        ok = .false.
        return
      end if
      x = x + residues(i) * mi * modulo(inv, moduli(i))
    end do
    x = modulo(x, mm)
  end subroutine

  subroutine chinese_elimination(a, b, x)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), intent(out) :: x(:)
    integer :: n, col, fila, piv, k, i, j
    real(dp) :: tmp, factor, acc
    real(dp) :: aw(size(a, 1), size(a, 2)), bw(size(b))
    n = size(b)
    aw = a
    bw = b
    x = 0.0_dp
    do col = 1, n
      piv = col
      do fila = col + 1, n
        if (abs(aw(fila, col)) > abs(aw(piv, col))) piv = fila
      end do
      if (abs(aw(piv, col)) < 1.0e-14_dp) return
      if (piv /= col) then
        do k = 1, n
          tmp = aw(col, k)
          aw(col, k) = aw(piv, k)
          aw(piv, k) = tmp
        end do
        tmp = bw(col)
        bw(col) = bw(piv)
        bw(piv) = tmp
      end if
      do fila = col + 1, n
        factor = aw(fila, col) / aw(col, col)
        do k = col, n
          aw(fila, k) = aw(fila, k) - factor * aw(col, k)
        end do
        bw(fila) = bw(fila) - factor * bw(col)
      end do
    end do
    do i = n, 1, -1
      acc = bw(i)
      do j = i + 1, n
        acc = acc - aw(i, j) * x(j)
      end do
      x(i) = acc / aw(i, i)
    end do
  end subroutine

  subroutine vigesimal(valor, max_digitos, digs, nd)
    integer(i8), intent(in) :: valor
    integer, intent(in) :: max_digitos
    integer(i8), intent(out) :: digs(:)
    integer, intent(out) :: nd
    integer(i8) :: v
    v = valor
    if (v == 0_i8) then
      digs(1) = 0_i8
      nd = 1
      return
    end if
    nd = 0
    do while (v > 0_i8 .and. nd < max_digitos)
      nd = nd + 1
      digs(nd) = modulo(v, 20_i8)
      v = v / 20_i8
    end do
  end subroutine

  pure function ms_step(x, digits) result(y)
    integer(i8), intent(in) :: x
    integer, intent(in) :: digits
    integer(i8) :: y, p
    p = 10_i8 ** (digits / 2)
    y = modulo((x * x) / p, 10_i8 ** digits)
  end function

  subroutine wh_step(s, mults, mods, r)
    integer(i8), intent(inout) :: s(3)
    integer(i8), intent(in) :: mults(3), mods(3)
    real(dp), intent(out) :: r
    integer :: i
    do i = 1, 3
      s(i) = modulo(mults(i) * s(i), mods(i))
    end do
    r = modulo(s(1) / real(mods(1), dp) + s(2) / real(mods(2), dp) + &
               s(3) / real(mods(3), dp), 1.0_dp)
  end subroutine

  pure function floor_div(a, b) result(q)
    integer(i8), intent(in) :: a, b
    integer(i8) :: q
    q = a / b
    if (modulo(a, b) /= 0_i8) then
      if ((a < 0_i8) .neqv. (b < 0_i8)) q = q - 1_i8
    end if
  end function

  subroutine wh_step16(s, mults, mods, r)
    integer(i8), intent(inout) :: s(3)
    integer(i8), intent(in) :: mults(3), mods(3)
    real(dp), intent(out) :: r
    s(1) = mults(1) * modulo(s(1), 177_i8) - 2_i8 * floor_div(s(1), 177_i8)
    s(2) = mults(2) * modulo(s(2), 176_i8) - 35_i8 * floor_div(s(2), 176_i8)
    s(3) = mults(3) * modulo(s(3), 178_i8) - 63_i8 * floor_div(s(3), 178_i8)
    r = modulo(s(1) / real(mods(1), dp) + s(2) / real(mods(2), dp) + &
               s(3) / real(mods(3), dp), 1.0_dp)
  end subroutine

end module fmath
