module fkernels
  use kinds, only: dp, i8
  use fhandles, only: map_step, rhs
  use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
  implicit none

  interface
    function pot_vfun(x) result(v)
      import :: dp
      real(dp), intent(in) :: x
      real(dp) :: v
    end function
  end interface

contains

  function lyap_map(name, p, x0, n, transient) result(lam)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: p(:)
    real(dp), intent(in) :: x0(:)
    integer, intent(in) :: n, transient
    real(dp) :: lam
    real(dp), allocatable :: x(:), y(:), xt(:), yt(:), d(:)
    real(dp) :: acc, d0, t, nd, s
    integer :: i, j, dim

    d0 = 1.0e-8_dp
    dim = size(x0)
    allocate(x(dim), y(dim), xt(dim), yt(dim), d(dim))
    x = x0
    do i = 1, transient
      call map_step(name, p, x, xt)
      x = xt
    end do
    t = d0 / sqrt(real(dim, dp))
    y = x + t
    acc = 0.0_dp
    do i = 1, n
      call map_step(name, p, x, xt)
      call map_step(name, p, y, yt)
      x = xt
      y = yt
      d = y - x
      s = 0.0_dp
      do j = 1, dim
        s = s + d(j) * d(j)
      end do
      nd = sqrt(s)
      acc = acc + log(nd / d0)
      y = x + d0 * (d / nd)
    end do
    lam = acc / real(n, dp)
  end function lyap_map

  function lyap_ode(name, p, x0, dt, t_transient, t_lyap) result(lam)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: p(:)
    real(dp), intent(in) :: x0(:)
    real(dp), intent(in) :: dt, t_transient, t_lyap
    real(dp) :: lam
    real(dp), allocatable :: x(:), y(:), d(:)
    real(dp) :: h, acc, d0, t, nd, s
    integer :: i, j, n, n_tr, dim

    d0 = 1.0e-8_dp
    h = dt
    dim = size(x0)
    allocate(x(dim), y(dim), d(dim))
    x = x0
    n_tr = nint(t_transient / h)
    do i = 1, n_tr
      call rk4_step(x)
    end do
    t = d0 / sqrt(real(dim, dp))
    y = x + t
    n = nint(t_lyap / h)
    acc = 0.0_dp
    do i = 1, n
      call rk4_step(x)
      call rk4_step(y)
      d = y - x
      s = 0.0_dp
      do j = 1, dim
        s = s + d(j) * d(j)
      end do
      nd = sqrt(s)
      acc = acc + log(nd / d0)
      y = x + d0 * (d / nd)
    end do
    lam = acc / (real(n, dp) * h)

  contains

    subroutine rk4_step(z)
      real(dp), intent(inout) :: z(:)
      real(dp) :: k1(size(z)), k2(size(z)), k3(size(z)), k4(size(z))
      real(dp) :: tmp(size(z))

      call rhs(name, p, z, k1)
      tmp = z + (h / 2.0_dp) * k1
      call rhs(name, p, tmp, k2)
      tmp = z + (h / 2.0_dp) * k2
      call rhs(name, p, tmp, k3)
      tmp = z + h * k3
      call rhs(name, p, tmp, k4)
      z = z + (h / 6.0_dp) * (k1 + 2.0_dp * k2 + 2.0_dp * k3 + k4)
    end subroutine rk4_step

  end function lyap_ode

  function box_counting_dim(pts, sizes) result(dimv)
    real(dp), intent(in) :: pts(:,:), sizes(:)
    real(dp) :: dimv
    integer(i8), allocatable :: keys(:)
    real(dp), allocatable :: cnts(:), logx(:), logy(:)
    real(dp) :: lo(size(pts, 2))
    real(dp) :: sx, sy, mx, my, sxy, sxx
    integer(i8) :: k1, k2, kcnt
    integer :: i, m, np, ns, ngood
    integer :: iq1, iq2

    np = size(pts, 1)
    ns = size(sizes)
    allocate(keys(np), cnts(ns), logx(ns), logy(ns))
    do m = 1, size(pts, 2)
      lo(m) = minval(pts(:, m))
    end do
    do i = 1, ns
      do m = 1, np
        iq1 = int(floor((pts(m, 1) - lo(1)) / sizes(i)))
        iq2 = int(floor((pts(m, 2) - lo(2)) / sizes(i)))
        k1 = int(iq1, i8)
        k2 = int(iq2, i8)
        keys(m) = k1 * 1000000000_i8 + k2
      end do
      call hsort_i8(keys)
      kcnt = 1_i8
      do m = 2, np
        if (keys(m) /= keys(m - 1)) kcnt = kcnt + 1_i8
      end do
      cnts(i) = real(kcnt, dp)
    end do

    ngood = 0
    do i = 1, ns
      if (cnts(i) > 1.0_dp) ngood = ngood + 1
    end do
    if (ngood < 2) then
      dimv = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if

    sx = 0.0_dp
    sy = 0.0_dp
    do i = 1, ns
      if (cnts(i) > 1.0_dp) then
        logx(i) = log(1.0_dp / sizes(i))
        logy(i) = log(cnts(i))
        sx = sx + logx(i)
        sy = sy + logy(i)
      end if
    end do
    mx = sx / real(ngood, dp)
    my = sy / real(ngood, dp)
    sxy = 0.0_dp
    sxx = 0.0_dp
    do i = 1, ns
      if (cnts(i) > 1.0_dp) then
        sxy = sxy + (logx(i) - mx) * (logy(i) - my)
        sxx = sxx + (logx(i) - mx) ** 2
      end if
    end do
    dimv = sxy / sxx
  end function box_counting_dim

  subroutine corr_dim_grassberger(pts, n_radii, c_lo, c_hi, dimv, ok, msg)
    real(dp), intent(in) :: pts(:,:)
    integer, intent(in) :: n_radii
    real(dp), intent(in) :: c_lo, c_hi
    real(dp), intent(out) :: dimv
    logical, intent(out) :: ok
    character(len=*), intent(out) :: msg
    real(dp), allocatable :: d(:), rr(:), lx(:), ly(:)
    logical, allocatable :: mk(:)
    real(dp) :: dmin, dmax, lmin, span, cnt, cval, x
    real(dp) :: sx, sy, mx, my, sxy, sxx
    integer :: i, j, t, m, np, npairs, k, lo_i, hi_i, mid, nmask

    ok = .true.
    msg = ''
    np = size(pts, 1)
    npairs = np * (np - 1) / 2
    allocate(d(npairs))
    t = 0
    do i = 1, np - 1
      do j = i + 1, np
        t = t + 1
        d(t) = sqrt((pts(j, 1) - pts(i, 1)) ** 2 + (pts(j, 2) - pts(i, 2)) ** 2)
      end do
    end do

    m = 0
    do t = 1, npairs
      if (d(t) > 0.0_dp) then
        m = m + 1
        d(m) = d(t)
      end if
    end do
    if (m < 10) then
      ok = .false.
      dimv = 0.0_dp
      msg = 'corr_dim: too few nonzero pairs'
      return
    end if

    dmin = d(1)
    dmax = d(1)
    do t = 2, m
      if (d(t) < dmin) dmin = d(t)
      if (d(t) > dmax) dmax = d(t)
    end do

    allocate(rr(n_radii), lx(n_radii), ly(n_radii), mk(n_radii))
    lmin = log(dmin)
    span = (log(dmax) - lmin) / real(n_radii - 1, dp)
    do k = 1, n_radii
      rr(k) = exp(lmin + real(k - 1, dp) * span)
    end do

    call hsort_dp(d(1:m))

    nmask = 0
    sx = 0.0_dp
    sy = 0.0_dp
    do k = 1, n_radii
      lo_i = 1
      hi_i = m
      cnt = 0.0_dp
      do while (lo_i <= hi_i)
        mid = (lo_i + hi_i) / 2
        if (d(mid) <= rr(k)) then
          cnt = real(mid, dp)
          lo_i = mid + 1
        else
          hi_i = mid - 1
        end if
      end do
      cval = cnt / real(m, dp)
      if (cval >= c_lo .and. cval <= c_hi) then
        mk(k) = .true.
        lx(k) = log(rr(k))
        ly(k) = log(cval)
        sx = sx + lx(k)
        sy = sy + ly(k)
        nmask = nmask + 1
      else
        mk(k) = .false.
        lx(k) = 0.0_dp
        ly(k) = 0.0_dp
      end if
    end do

    if (nmask < 3) then
      ok = .false.
      dimv = 0.0_dp
      msg = 'corr_dim: empty scaling window'
      return
    end if

    mx = sx / real(nmask, dp)
    my = sy / real(nmask, dp)
    sxy = 0.0_dp
    sxx = 0.0_dp
    do k = 1, n_radii
      if (mk(k)) then
        x = lx(k) - mx
        sxy = sxy + x * (ly(k) - my)
        sxx = sxx + x * x
      end if
    end do
    dimv = sxy / sxx
  end subroutine corr_dim_grassberger

  subroutine fd_eigenvalues(g, V, k, ev, ok, msg)
    real(dp), intent(in) :: g(:), V(:)
    integer, intent(in) :: k
    real(dp), intent(out) :: ev(k)
    logical, intent(out) :: ok
    character(len=*), intent(out) :: msg
    real(dp), allocatable :: dd(:), ee(:), z(:,:), work(:)
    real(dp) :: h, c
    integer :: n, info, j

    ok = .true.
    msg = ''
    n = size(g)
    h = g(2) - g(1)
    c = 0.5_dp / (h * h)
    allocate(dd(n), ee(n - 1), z(n, k), work(4 * n))
    do j = 1, n
      dd(j) = 2.0_dp * c + V(j)
    end do
    do j = 1, n - 1
      ee(j) = -c
    end do
    call dstev('N', n, dd, ee, z, n, work, info)
    if (info /= 0) then
      ok = .false.
      ev = 0.0_dp
      write(msg, '(a,i0)') 'dstev failed with info=', info
      return
    end if
    do j = 1, k
      ev(j) = dd(j)
    end do
  end subroutine fd_eigenvalues

  subroutine golden_min(fun, lo, hi, xmin, fmin)
    procedure(pot_vfun) :: fun
    real(dp), intent(in) :: lo, hi
    real(dp), intent(out) :: xmin, fmin
    real(dp) :: a, b, c, d, fc, fd, gm
    integer :: it

    gm = 0.5_dp * (3.0_dp - sqrt(5.0_dp))
    a = lo
    b = hi
    c = a + gm * (b - a)
    d = b - gm * (b - a)
    fc = fun(c)
    fd = fun(d)
    it = 0
    do while ((b - a) > 1.0e-14_dp * (abs(a) + abs(b) + 1.0_dp) .and. it < 500)
      if (fc < fd) then
        b = d
        d = c
        fd = fc
        c = a + gm * (b - a)
        fc = fun(c)
      else
        a = c
        c = d
        fc = fd
        d = b - gm * (b - a)
        fd = fun(d)
      end if
      it = it + 1
    end do
    if (fc <= fd) then
      xmin = c
      fmin = fc
    else
      xmin = d
      fmin = fd
    end if
  end subroutine golden_min

  subroutine hsort_dp(a)
    real(dp), intent(inout) :: a(:)
    integer :: n, i, ir
    real(dp) :: tmp

    n = size(a)
    if (n < 2) return
    do i = n / 2, 1, -1
      call sift_down(i, n)
    end do
    do ir = n, 2, -1
      tmp = a(1)
      a(1) = a(ir)
      a(ir) = tmp
      call sift_down(1, ir - 1)
    end do

  contains

    subroutine sift_down(root0, bottom0)
      integer, intent(in) :: root0, bottom0
      integer :: root, child, bottom
      real(dp) :: t

      root = root0
      bottom = bottom0
      do
        child = 2 * root
        if (child > bottom) exit
        if (child < bottom) then
          if (a(child) < a(child + 1)) child = child + 1
        end if
        if (a(root) >= a(child)) exit
        t = a(root)
        a(root) = a(child)
        a(child) = t
        root = child
      end do
    end subroutine sift_down

  end subroutine hsort_dp

  subroutine hsort_i8(a)
    integer(i8), intent(inout) :: a(:)
    integer :: n, i, ir
    integer(i8) :: tmp

    n = size(a)
    if (n < 2) return
    do i = n / 2, 1, -1
      call sift_down(i, n)
    end do
    do ir = n, 2, -1
      tmp = a(1)
      a(1) = a(ir)
      a(ir) = tmp
      call sift_down(1, ir - 1)
    end do

  contains

    subroutine sift_down(root0, bottom0)
      integer, intent(in) :: root0, bottom0
      integer :: root, child, bottom
      integer(i8) :: t

      root = root0
      bottom = bottom0
      do
        child = 2 * root
        if (child > bottom) exit
        if (child < bottom) then
          if (a(child) < a(child + 1)) child = child + 1
        end if
        if (a(root) >= a(child)) exit
        t = a(root)
        a(root) = a(child)
        a(child) = t
        root = child
      end do
    end subroutine sift_down

  end subroutine hsort_i8

end module fkernels
