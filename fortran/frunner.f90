program f_runner
  use kinds, only: dp
  use fjson
  use fhandles
  use fkernels
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite, &
       ieee_quiet_nan, ieee_value
  implicit none

  integer, parameter :: MAXCHK = 64
  integer, parameter :: MAXRUN = 64
  integer, parameter :: NVPMAX = 8

  type :: check_t
    character(len=64) :: quantity = ''
    real(dp) :: computed = 0.0_dp
    real(dp) :: reference = 0.0_dp
    real(dp) :: tol = 0.0_dp
    real(dp) :: rel_err = 0.0_dp
    logical :: pass = .true.
  end type check_t

  type :: run_t
    character(len=64) :: name = ''
    character(len=64) :: family = ''
    character(len=512) :: errmsg = ''
    character(len=512) :: orbit_path = ''
    logical :: pass = .true.
    logical :: is_err = .false.
    logical :: has_orbit = .false.
    logical :: has_lyap1 = .false.
    logical :: has_dim = .false.
    real(dp) :: lyap1_obs = 0.0_dp
    real(dp) :: dim_obs = 0.0_dp
    type(check_t) :: chk(MAXCHK)
    integer :: nchk = 0
  end type run_t

  character(len=64) :: sel
  character(len=512) :: outdir, emsg, jpath
  type(jval) :: reg1, reg2, ents, ent, top, ra
  type(run_t) :: runs(MAXRUN)
  logical :: ok, ok2
  integer :: nrun, npass, ei, ri, nent, u, ios
  character(len=8) :: cdate
  character(len=10) :: ctime
  character(len=15) :: ts

  sel = 'all'
  outdir = 'results'
  if (command_argument_count() >= 1) call get_command_argument(1, sel)
  if (command_argument_count() >= 2) call get_command_argument(2, outdir)

  call jload('registries/attractors.json', reg1, ok, emsg)
  if (.not. ok) call jload('../registries/attractors.json', reg1, ok, emsg)
  if (.not. ok) then
    write(*, '(a)') 'ERROR: ' // trim(emsg)
    error stop 1
  end if
  call jload('registries/potentials.json', reg2, ok2, emsg)
  if (.not. ok2) call jload('../registries/potentials.json', reg2, ok2, emsg)
  if (.not. ok2) then
    write(*, '(a)') 'ERROR: ' // trim(emsg)
    error stop 1
  end if

  call execute_command_line('mkdir -p "' // trim(outdir) // '"', exitstat=ios)

  nrun = 0
  call jcopy(ents, jget(reg1, 'entries'))
  nent = jalen(ents)
  do ei = 1, nent
    call jcopy(ent, jat(ents, ei))
    if (selected(ent, sel)) call run_one(ent, outdir, runs, nrun)
  end do
  call jcopy(ents, jget(reg2, 'entries'))
  nent = jalen(ents)
  do ei = 1, nent
    call jcopy(ent, jat(ents, ei))
    if (selected(ent, sel)) call run_one(ent, outdir, runs, nrun)
  end do

  npass = 0
  do ri = 1, nrun
    if (runs(ri)%pass) npass = npass + 1
  end do
  write(*, '(/,a,i0,a,i0,a)') 'Summary: ', npass, '/', nrun, ' entries PASS'

  call date_and_time(date=cdate, time=ctime)
  ts = cdate(1:8) // 'T' // ctime(1:6)
  jpath = trim(outdir) // '/registry_run_f_' // ts // '.json'
  call jcopy(top, jmk_obj())
  call jset(top, 'timestamp', jmk_str(trim(ts)))
  call jset(top, 'selector', jmk_str(trim(sel)))
  call jset(top, 'language', jmk_str('fortran'))
  call jset(top, 'n_pass', jmk_num(real(npass, dp)))
  call jset(top, 'n_total', jmk_num(real(nrun, dp)))
  call jcopy(ra, jmk_arr())
  do ri = 1, nrun
    call jpush(ra, run_json(runs(ri)))
  end do
  call jset(top, 'runs', ra)
  open(newunit=u, file=trim(jpath), status='replace', action='write')
  call jdump(u, top)
  write(u, '(a)') ''
  close(u)
  write(*, '(a)') 'Saved: ' // trim(jpath)

contains

  function nanv() result(x)
    real(dp) :: x
    x = ieee_value(0.0_dp, ieee_quiet_nan)
  end function

  function atnum(o, i) result(x)
    type(jval), intent(in) :: o
    integer, intent(in) :: i
    real(dp) :: x
    type(jval) :: c
    call jcopy(c, jat(o, i))
    if (c%t == JV_NUM) then
      x = c%v
    else
      x = nanv()
    end if
  end function

  function lower(s) result(o)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: o
    integer :: i, c
    o = s
    do i = 1, len(s)
      c = iachar(o(i:i))
      if (c >= iachar('A') .and. c <= iachar('Z')) o(i:i) = achar(c + 32)
    end do
  end function

  function linspace(a, b, n) result(g)
    real(dp), intent(in) :: a, b
    integer, intent(in) :: n
    real(dp), allocatable :: g(:)
    real(dp) :: stp
    integer :: i
    allocate(g(max(n, 0)))
    if (n < 2) return
    stp = (b - a) / real(n - 1, dp)
    do i = 1, n
      g(i) = a + stp * real(i - 1, dp)
    end do
  end function

  function sgn3(x) result(s)
    real(dp), intent(in) :: x
    real(dp) :: s
    if (x > 0.0_dp) then
      s = 1.0_dp
    else if (x < 0.0_dp) then
      s = -1.0_dp
    else
      s = 0.0_dp
    end if
  end function

  function fg(x, p) result(s)
    real(dp), intent(in) :: x
    integer, intent(in) :: p
    character(len=:), allocatable :: s
    character(len=97) :: c, f, mant, frac, edig
    character(len=80) :: digits
    integer :: i, e10, ddec, ios, ia, nd
    real(dp) :: a
    logical :: neg

    if (ieee_is_nan(x)) then
      s = 'nan'
      return
    end if
    if (.not. ieee_is_finite(x)) then
      if (x < 0.0_dp) then
        s = '-inf'
      else
        s = 'inf'
      end if
      return
    end if
    neg = x < 0.0_dp
    a = abs(x)
    if (a == 0.0_dp) then
      s = '0'
      return
    end if
    write(f, '("(es0.",i0,"e2))")') p - 1
    write(c, f) a
    i = index(c, 'E')
    if (i == 0) i = index(c, 'e')
    mant = c(:i-1)
    digits = ''
    nd = 0
    i = 1
    do while (i <= len_trim(mant))
      if (mant(i:i) /= '.') then
        nd = nd + 1
        digits(nd:nd) = mant(i:i)
      end if
      i = i + 1
    end do
    i = index(c, 'E')
    if (i == 0) i = index(c, 'e')
    read(c(i + 1:), *, iostat=ios) e10
    if (e10 < -4 .or. e10 >= p) then
      frac = digits(2:p)
      do while (len_trim(frac) > 0)
        if (frac(len_trim(frac):len_trim(frac)) /= '0') exit
        frac(len_trim(frac):len_trim(frac)) = ' '
      end do
      if (len_trim(frac) > 0) then
        c = digits(1:1) // '.' // trim(frac)
      else
        c = digits(1:1)
      end if
      if (e10 < 0) then
        ia = -e10
        if (ia < 100) then
          write(edig, '(i2.2)') ia
        else
          write(edig, '(i0)') ia
        end if
        if (neg) then
          s = '-' // c(1:len_trim(c)) // 'e-' // trim(edig)
        else
          s = c(1:len_trim(c)) // 'e-' // trim(edig)
        end if
      else
        ia = e10
        if (ia < 100) then
          write(edig, '(i2.2)') ia
        else
          write(edig, '(i0)') ia
        end if
        if (neg) then
          s = '-' // c(1:len_trim(c)) // 'e+' // trim(edig)
        else
          s = c(1:len_trim(c)) // 'e+' // trim(edig)
        end if
      end if
    else
      ddec = p - 1 - e10
      write(f, '("(f0.",i0,")")') ddec
      write(c, f) a
      c = adjustl(c)
      if (index(c, '.') > 0) then
        do while (len_trim(c) > 0)
          if (c(len_trim(c):len_trim(c)) /= '0') exit
          c(len_trim(c):len_trim(c)) = ' '
        end do
        if (len_trim(c) > 0) then
          if (c(len_trim(c):len_trim(c)) == '.') c(len_trim(c):len_trim(c)) = ' '
        end if
      end if
      if (c(1:1) == '.') c = '0' // c(1:len(c)-1)
      if (neg) then
        s = '-' // trim(c)
      else
        s = trim(c)
      end if
    end if
  end function

  subroutine init_run(r, name, family)
    type(run_t), intent(inout) :: r
    character(len=*), intent(in) :: name, family
    r%name = name
    r%family = family
    r%errmsg = ''
    r%orbit_path = ''
    r%pass = .true.
    r%is_err = .false.
    r%has_orbit = .false.
    r%has_lyap1 = .false.
    r%has_dim = .false.
    r%lyap1_obs = 0.0_dp
    r%dim_obs = 0.0_dp
    r%nchk = 0
  end subroutine

  subroutine add_check(r, quantity, computed, reference, tol, rel_err, ok)
    type(run_t), intent(inout) :: r
    character(len=*), intent(in) :: quantity
    real(dp), intent(in) :: computed, reference, tol, rel_err
    logical, intent(in) :: ok
    if (r%nchk >= MAXCHK) return
    r%nchk = r%nchk + 1
    r%chk(r%nchk)%quantity = quantity
    r%chk(r%nchk)%computed = computed
    r%chk(r%nchk)%reference = reference
    r%chk(r%nchk)%tol = tol
    r%chk(r%nchk)%rel_err = rel_err
    r%chk(r%nchk)%pass = ok
  end subroutine

  subroutine add_obs(r, quantity, computed)
    type(run_t), intent(inout) :: r
    character(len=*), intent(in) :: quantity
    real(dp), intent(in) :: computed
    call add_check(r, quantity, computed, nanv(), nanv(), nanv(), .true.)
  end subroutine

  subroutine set_error(r, msg)
    type(run_t), intent(inout) :: r
    character(len=*), intent(in) :: msg
    r%is_err = .true.
    r%errmsg = msg
    r%pass = .false.
    r%nchk = 1
    r%has_orbit = .false.
    r%has_lyap1 = .false.
    r%has_dim = .false.
    r%chk(1)%quantity = 'error'
    r%chk(1)%computed = nanv()
    r%chk(1)%reference = nanv()
    r%chk(1)%tol = nanv()
    r%chk(1)%rel_err = nanv()
    r%chk(1)%pass = .false.
  end subroutine

  subroutine finish(r)
    type(run_t), intent(inout) :: r
    logical :: p
    integer :: i
    p = .true.
    do i = 1, r%nchk
      if (.not. r%chk(i)%pass) p = .false.
    end do
    r%pass = p
  end subroutine

  function pass_str(b) result(s)
    logical, intent(in) :: b
    character(len=4) :: s
    if (b) then
      s = 'PASS'
    else
      s = 'FAIL'
    end if
  end function

  subroutine printrun(r)
    type(run_t), intent(in) :: r
    character(len=1024) :: line
    integer :: i
    if (r%is_err) then
      write(*, '(a)') '  ERROR: ' // trim(r%errmsg)
      write(*, '(a)') '  => FAIL'
      return
    end if
    do i = 1, r%nchk
      if (ieee_is_nan(r%chk(i)%reference)) then
        line = '  ' // pass_str(r%chk(i)%pass) // ' ' // &
             trim(r%chk(i)%quantity) // ' (observed) = ' // &
             trim(fg(r%chk(i)%computed, 12))
      else
        line = '  ' // pass_str(r%chk(i)%pass) // ' ' // &
             trim(r%chk(i)%quantity) // ': computed=' // &
             trim(fg(r%chk(i)%computed, 12)) // ' ref=' // &
             trim(fg(r%chk(i)%reference, 12)) // ' rel_err=' // &
             trim(fg(r%chk(i)%rel_err, 3)) // ' tol=' // &
             trim(fg(r%chk(i)%tol, 3))
      end if
      write(*, '(a)') trim(line)
    end do
    if (r%pass) then
      write(*, '(a)') '  => PASS'
    else
      write(*, '(a)') '  => FAIL'
    end if
  end subroutine

  function check_json(c) result(o)
    type(check_t), intent(in) :: c
    type(jval) :: o
    call jcopy(o, jmk_obj())
    call jset(o, 'quantity', jmk_str(trim(c%quantity)))
    call jset(o, 'computed', jmk_num(c%computed))
    call jset(o, 'reference', jmk_num(c%reference))
    call jset(o, 'tol', jmk_num(c%tol))
    call jset(o, 'rel_err', jmk_num(c%rel_err))
    call jset(o, 'pass', jmk_bool(c%pass))
  end function

  function run_json(r) result(o)
    type(run_t), intent(in) :: r
    type(jval) :: o, ca
    integer :: i
    call jcopy(o, jmk_obj())
    call jset(o, 'name', jmk_str(trim(r%name)))
    call jset(o, 'family', jmk_str(trim(r%family)))
    call jset(o, 'pass', jmk_bool(r%pass))
    call jcopy(ca, jmk_arr())
    do i = 1, r%nchk
      call jpush(ca, check_json(r%chk(i)))
    end do
    call jset(o, 'checks', ca)
    if (r%has_orbit) call jset(o, 'orbit_path', jmk_str(trim(r%orbit_path)))
    if (r%has_lyap1) call jset(o, 'lyap1_observed', jmk_num(r%lyap1_obs))
    if (r%has_dim) call jset(o, 'dim_observed', jmk_num(r%dim_obs))
    if (r%is_err) call jset(o, 'error', jmk_str(trim(r%errmsg)))
  end function

  function selected(e, s) result(y)
    type(jval), intent(in) :: e
    character(len=*), intent(in) :: s
    logical :: y
    character(len=64) :: fam, nmx
    if (s == 'all') then
      y = .true.
      return
    end if
    fam = jstr(e, 'family', '')
    if (s == 'attractors') then
      y = fam /= 'potential'
      return
    end if
    if (s == 'potentials') then
      y = fam == 'potential'
      return
    end if
    nmx = jstr(e, 'name', '')
    y = (s == nmx) .or. (s == fam)
  end function

  subroutine seed_rng()
    integer, allocatable :: s(:)
    integer :: n, j
    call random_seed(size=n)
    allocate(s(n))
    do j = 1, n
      s(j) = 20261023 + 37 * j
    end do
    call random_seed(put=s)
    deallocate(s)
  end subroutine

  subroutine run_one(ent, outdir, runs, nrun)
    type(jval), intent(in) :: ent
    character(len=*), intent(in) :: outdir
    type(run_t), intent(inout) :: runs(MAXRUN)
    integer, intent(inout) :: nrun
    type(run_t) :: r
    logical :: ok
    character(len=512) :: emsg
    character(len=64) :: fam, nmx

    if (nrun >= MAXRUN) return
    nmx = jstr(ent, 'name', '')
    fam = jstr(ent, 'family', '')
    call init_run(r, nmx, fam)
    write(*, '(a)') '[' // trim(fam) // '] ' // trim(nmx)
    emsg = ''
    select case (trim(fam))
    case ('map')
      call eval_map(ent, outdir, r, ok, emsg)
    case ('ifs')
      call eval_ifs(ent, r, ok, emsg)
    case ('ode')
      call eval_ode(ent, r, ok, emsg)
    case ('neural')
      call eval_neural(ent, r, ok, emsg)
    case ('potential')
      call eval_potential(ent, r, ok, emsg)
    case default
      ok = .false.
      emsg = 'unknown family: ' // trim(fam)
    end select
    if (.not. ok) then
      call set_error(r, emsg)
    else
      call finish(r)
    end if
    call printrun(r)
    nrun = nrun + 1
    runs(nrun) = r
  end subroutine

  subroutine eval_map(e, outdir, r, ok, emsg)
    type(jval), intent(in) :: e
    character(len=*), intent(in) :: outdir
    type(run_t), intent(inout) :: r
    logical, intent(out) :: ok
    character(len=*), intent(out) :: emsg
    type(jval) :: ref, oh, ly, cd, proto, valsv, row, cell
    real(dp) :: p(NVPMAX)
    real(dp) :: lam, md, val, tol, rel, dimv, rlt, abst, clo, chi
    real(dp), allocatable :: x0(:), orbit(:, :), pts(:, :), x(:), xt(:)
    integer :: norb, nlyap, ntr, d, i, j, k, u, ios, npts2, nrad, ptr2, np, need
    logical :: cok
    character(len=512) :: cmsg, path
    character(len=64) :: nmx

    nmx = jstr(e, 'name', '')
    call load_params(nmx, e, p, np, ok, emsg)
    if (.not. ok) return
    x0 = jtovec(e, 'x0')
    d = size(x0)
    if (nmx == 'logistic') then
      need = 1
    else
      need = 2
    end if
    if (d < need) then
      ok = .false.
      emsg = 'x0 too short for ' // trim(nmx)
      return
    end if
    norb = jint(e, 'n_orbit', 0)
    allocate(orbit(norb, d), x(d), xt(d))
    x = x0
    do i = 1, norb
      call map_step(nmx, p, x, xt)
      x = xt
      orbit(i, :) = x
    end do
    do i = 1, norb
      do j = 1, d
        if (.not. ieee_is_finite(orbit(i, j))) then
          ok = .false.
          emsg = 'Orbit of ' // trim(nmx) // ' is non-finite'
          return
        end if
      end do
    end do
    call execute_command_line('mkdir -p "' // trim(outdir) // '/orbits"', &
         exitstat=ios)
    path = trim(outdir) // '/orbits/' // trim(nmx) // '_orbit.txt'
    open(newunit=u, file=trim(path), status='replace', action='write')
    do i = 1, norb
      write(u, '(*(g0.17,1x))') orbit(i, 1:d)
    end do
    close(u)
    r%orbit_path = path
    r%has_orbit = .true.

    nlyap = jint(e, 'n_lyap', 0)
    ntr = jint(e, 'transient', 0)
    lam = lyap_map(nmx, p, x0, nlyap, ntr)
    r%lyap1_obs = lam
    r%has_lyap1 = .true.

    call jcopy(ref, jget(e, 'ref'))
    call jcopy(oh, jget(ref, 'orbit_head'))
    if (oh%t == JV_OBJ) then
      k = jint(oh, 'n', 0)
      call jcopy(valsv, jget(oh, 'values'))
      md = 0.0_dp
      do i = 1, min(k, jalen(valsv), norb)
        call jcopy(row, jat(valsv, i))
        do j = 1, min(d, jalen(row))
          call jcopy(cell, jat(row, j))
          if (cell%t /= JV_NUM) then
            ok = .false.
            emsg = 'orbit_head values not numeric'
            return
          end if
          md = max(md, abs(orbit(i, j) - cell%v))
        end do
      end do
      abst = jnum_req(oh, 'abs_tol')
      call add_check(r, 'orbit_head_maxabs', md, 0.0_dp, abst, md, md <= abst)
    end if
    call jcopy(ly, jget(ref, 'lyap1'))
    if (ly%t == JV_OBJ) then
      val = jnum_req(ly, 'value')
      abst = jnum_req(ly, 'abs_tol')
      rel = abs(lam - val) / max(abs(val), epsilon(1.0_dp))
      call add_check(r, 'lyap1', lam, val, abst, rel, abs(lam - val) <= abst)
    else
      call add_obs(r, 'lyap1_observed', lam)
    end if
    call jcopy(cd, jget(ref, 'corr_dim'))
    if (cd%t == JV_OBJ) then
      call jcopy(proto, jget(cd, 'protocol'))
      ptr2 = jint(proto, 'transient', 0)
      npts2 = jint(proto, 'n_points', 0)
      nrad = jint(proto, 'n_radii', 0)
      clo = jnum_req(proto, 'c_lo')
      chi = jnum_req(proto, 'c_hi')
      x = x0
      do i = 1, ptr2
        call map_step(nmx, p, x, xt)
        x = xt
      end do
      allocate(pts(npts2, d))
      do i = 1, npts2
        call map_step(nmx, p, x, xt)
        x = xt
        pts(i, :) = x
      end do
      call corr_dim_grassberger(pts, nrad, clo, chi, dimv, cok, cmsg)
      if (.not. cok) then
        ok = .false.
        emsg = cmsg
        return
      end if
      val = jnum_req(cd, 'value')
      rlt = jnum_req(cd, 'rel_tol')
      tol = rlt * abs(val)
      rel = abs(dimv - val) / abs(val)
      call add_check(r, 'corr_dim', dimv, val, tol, rel, rel <= rlt)
    end if
    ok = .true.
    emsg = ''
  end subroutine

  subroutine eval_ifs(e, r, ok, emsg)
    type(jval), intent(in) :: e
    type(run_t), intent(inout) :: r
    logical, intent(out) :: ok
    character(len=*), intent(out) :: emsg
    type(jval) :: pr, ref, dimr, vert, maps, mrow, row
    real(dp), allocatable :: x0(:), cp(:), pv(:), sizes(:), pts(:, :), vv(:, :)
    real(dp) :: x(2), u, dimv, val, tol, rel, rlt, ps, acc
    real(dp) :: ma, mb, mc, mdd, me, mf, x1, x2
    integer :: n, bi, i, k, nc, nv
    character(len=64) :: nmx

    nmx = jstr(e, 'name', '')
    n = jint(e, 'n_points', 0)
    bi = jint(e, 'burn_in', 0)
    x0 = jtovec(e, 'x0')
    if (size(x0) /= 2) then
      ok = .false.
      emsg = 'bad x0 for ' // trim(nmx)
      return
    end if
    x = x0(1:2)
    allocate(pts(n, 2))
    call jcopy(pr, jget(e, 'params'))
    call seed_rng()
    if (jhas(pr, 'vertices')) then
      call jcopy(vert, jget(pr, 'vertices'))
      nv = jalen(vert)
      allocate(vv(nv, 2), pv(nv), cp(nv))
      do i = 1, nv
        call jcopy(row, jat(vert, i))
        vv(i, 1) = atnum(row, 1)
        vv(i, 2) = atnum(row, 2)
      end do
      pv = jtovec(pr, 'p')
      if (size(pv) /= nv) then
        ok = .false.
        emsg = 'bad p for ' // trim(nmx)
        return
      end if
      ps = 0.0_dp
      do i = 1, nv
        ps = ps + pv(i)
      end do
      acc = 0.0_dp
      do i = 1, nv
        acc = acc + pv(i) / ps
        cp(i) = acc
      end do
      do i = 0, n + bi - 1
        call random_number(u)
        k = 1
        do while (k < nv .and. cp(k) < u)
          k = k + 1
        end do
        x = (x + vv(k, :)) / 2.0_dp
        if (i >= bi) then
          pts(i - bi + 1, 1) = x(1)
          pts(i - bi + 1, 2) = x(2)
        end if
      end do
    else if (jhas(pr, 'maps')) then
      call jcopy(maps, jget(pr, 'maps'))
      nc = jalen(maps)
      allocate(pv(nc), cp(nc))
      do i = 1, nc
        pv(i) = jnum(jat(maps, i), 'p', nanv())
      end do
      ps = 0.0_dp
      do i = 1, nc
        ps = ps + pv(i)
      end do
      acc = 0.0_dp
      do i = 1, nc
        acc = acc + pv(i) / ps
        cp(i) = acc
      end do
      do i = 0, n + bi - 1
        call random_number(u)
        k = 1
        do while (k < nc .and. cp(k) < u)
          k = k + 1
        end do
        call jcopy(mrow, jat(maps, k))
        ma = jnum(mrow, 'a', nanv())
        mb = jnum(mrow, 'b', nanv())
        mc = jnum(mrow, 'c', nanv())
        mdd = jnum(mrow, 'd', nanv())
        me = jnum(mrow, 'e', nanv())
        mf = jnum(mrow, 'f', nanv())
        x1 = (ma * x(1) + mb * x(2)) + me
        x2 = (mc * x(1) + mdd * x(2)) + mf
        x(1) = x1
        x(2) = x2
        if (i >= bi) then
          pts(i - bi + 1, 1) = x(1)
          pts(i - bi + 1, 2) = x(2)
        end if
      end do
    else
      ok = .false.
      emsg = 'Unknown IFS param shape for ' // trim(nmx)
      return
    end if

    sizes = jtovec(e, 'box_sizes')
    dimv = box_counting_dim(pts, sizes)
    r%has_dim = .true.
    r%dim_obs = dimv
    call jcopy(ref, jget(e, 'ref'))
    call jcopy(dimr, jget(ref, 'dim'))
    if (dimr%t == JV_OBJ) then
      val = jnum_req(dimr, 'value')
      rlt = jnum_req(dimr, 'rel_tol')
      tol = rlt * abs(val)
      rel = abs(dimv - val) / abs(val)
      call add_check(r, 'box_dim', dimv, val, tol, rel, rel <= rlt)
    else
      call add_obs(r, 'box_dim_observed', dimv)
    end if
    ok = .true.
    emsg = ''
  end subroutine

  subroutine eval_ode(e, r, ok, emsg)
    type(jval), intent(in) :: e
    type(run_t), intent(inout) :: r
    logical, intent(out) :: ok
    character(len=*), intent(out) :: emsg
    type(jval) :: ref, ly
    real(dp) :: p(NVPMAX)
    real(dp), allocatable :: x0(:)
    real(dp) :: lam, val, tol, rel, rlt, dt, tt, tl
    integer :: np, d
    character(len=64) :: nmx

    nmx = jstr(e, 'name', '')
    call load_params(nmx, e, p, np, ok, emsg)
    if (.not. ok) return
    x0 = jtovec(e, 'x0')
    d = size(x0)
    if (d /= 3) then
      ok = .false.
      emsg = 'bad x0 for ' // trim(nmx)
      return
    end if
    dt = jnum_req(e, 'dt_lyap')
    tt = jnum_req(e, 't_transient')
    tl = jnum_req(e, 't_lyap')
    lam = lyap_ode(nmx, p, x0, dt, tt, tl)
    r%lyap1_obs = lam
    r%has_lyap1 = .true.
    call jcopy(ref, jget(e, 'ref'))
    call jcopy(ly, jget(ref, 'lyap1'))
    if (ly%t == JV_OBJ) then
      val = jnum_req(ly, 'value')
      if (jhas(ly, 'rel_tol')) then
        rlt = jnum_req(ly, 'rel_tol')
        tol = rlt * abs(val)
        rel = abs(lam - val) / abs(val)
        call add_check(r, 'lyap1', lam, val, tol, rel, rel <= rlt)
      else
        tol = jnum_req(ly, 'abs_tol')
        rel = abs(lam - val) / max(abs(val), epsilon(1.0_dp))
        call add_check(r, 'lyap1', lam, val, tol, rel, abs(lam - val) <= tol)
      end if
    else
      call add_obs(r, 'lyap1_observed', lam)
    end if
    ok = .true.
    emsg = ''
  end subroutine

  subroutine eval_neural(e, r, ok, emsg)
    type(jval), intent(in) :: e
    type(run_t), intent(inout) :: r
    logical, intent(out) :: ok
    character(len=*), intent(out) :: emsg
    type(jval) :: pr, pats, row
    real(dp), allocatable :: pp(:, :), w(:, :), s(:)
    real(dp) :: acc
    integer :: n, npat, i, j, k
    logical :: okp
    character(len=64) :: nmx, qty

    nmx = jstr(e, 'name', '')
    call jcopy(pr, jget(e, 'params'))
    n = jint(pr, 'n', 0)
    call jcopy(pats, jget(pr, 'patterns'))
    npat = jalen(pats)
    if (n < 1 .or. npat < 1) then
      ok = .false.
      emsg = 'missing params for ' // trim(nmx)
      return
    end if
    allocate(pp(npat, n), w(n, n), s(n))
    do i = 1, npat
      call jcopy(row, jat(pats, i))
      do j = 1, n
        pp(i, j) = atnum(row, j)
      end do
    end do
    do i = 1, n
      do j = 1, n
        acc = 0.0_dp
        do k = 1, npat
          acc = acc + pp(k, i) * pp(k, j)
        end do
        w(i, j) = acc / real(n, dp)
      end do
      w(i, i) = 0.0_dp
    end do
    do i = 1, npat
      do j = 1, n
        acc = 0.0_dp
        do k = 1, n
          acc = acc + w(j, k) * pp(i, k)
        end do
        s(j) = acc
      end do
      okp = .true.
      do j = 1, n
        if (sgn3(s(j)) /= pp(i, j)) okp = .false.
      end do
      write(qty, '(a,i0)') 'fixed_point_pattern_', i
      if (okp) then
        call add_check(r, qty, 1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, .true.)
      else
        call add_check(r, qty, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, .false.)
      end if
    end do
    ok = .true.
    emsg = ''
  end subroutine

  subroutine eval_potential(e, r, ok, emsg)
    type(jval), intent(in) :: e
    type(run_t), intent(inout) :: r
    logical, intent(out) :: ok
    character(len=*), intent(out) :: emsg
    type(jval) :: ref, pr, xg, eigref, omega, force, rmin, rmina, vmin
    real(dp), allocatable :: rr(:), gv(:), vg(:), en(:), ev(:), evals(:), rc(:)
    real(dp) :: lo, hi, xmn, vmn, kpar, mpar, lam, hh, vp, vm, fn, fa, f1
    real(dp) :: val, tol, rel, rlt, w, hstep, x1, x2, avec, kap, r1
    integer :: ng, k, j
    character(len=64) :: nmx, nraw, qty

    nraw = jstr(e, 'name', '')
    nmx = lower(nraw)
    call jcopy(ref, jget(e, 'ref'))
    call jcopy(pr, jget(e, 'params'))
    select case (trim(nmx))
    case ('lennard_jones', 'morse', 'double_well')
      call load_pot(nmx, e, ok, emsg)
      if (.not. ok) return
      rr = jtovec(e, 'r_range')
      if (size(rr) < 2) then
        ok = .false.
        emsg = 'bad r_range for ' // trim(nraw)
        return
      end if
      lo = rr(1)
      hi = rr(2)
      if (trim(nmx) == 'lennard_jones') then
        call golden_min(lj_v, lo, hi, xmn, vmn)
      else if (trim(nmx) == 'morse') then
        call golden_min(morse_v, lo, hi, xmn, vmn)
      else
        call golden_min(dw_v, 0.0_dp, hi, xmn, vmn)
      end if
      call jcopy(rmin, jget(ref, 'r_min'))
      if (rmin%t == JV_OBJ) then
        val = jnum_req(rmin, 'value')
        tol = jnum_req(rmin, 'abs_tol')
        rel = abs(xmn - val) / abs(val)
        call add_check(r, 'r_min', xmn, val, tol, rel, abs(xmn - val) <= tol)
      end if
      call jcopy(rmina, jget(ref, 'r_min_abs'))
      if (rmina%t == JV_OBJ) then
        val = jnum_req(rmina, 'value')
        tol = jnum_req(rmina, 'abs_tol')
        rel = abs(abs(xmn) - val) / abs(val)
        call add_check(r, 'r_min_abs', abs(xmn), val, tol, rel, &
             abs(abs(xmn) - val) <= tol)
      end if
      call jcopy(vmin, jget(ref, 'V_min'))
      if (vmin%t == JV_OBJ) then
        val = jnum_req(vmin, 'value')
        tol = jnum_req(vmin, 'abs_tol')
        rel = abs(vmn - val) / max(abs(val), epsilon(1.0_dp))
        call add_check(r, 'V_min', vmn, val, tol, rel, abs(vmn - val) <= tol)
      end if

    case ('harmonic', 'hydrogen', 'poschl_teller')
      call jcopy(xg, jget(e, 'x_grid'))
      x1 = jnum(xg, 'xmin', nanv())
      x2 = jnum(xg, 'xmax', nanv())
      ng = jint(xg, 'n', 0)
      if (ng < 2) then
        ok = .false.
        emsg = 'bad x_grid for ' // trim(nraw)
        return
      end if
      allocate(gv(ng), vg(ng))
      select case (trim(nmx))
      case ('harmonic')
        gv = linspace(x1, x2, ng)
        kpar = jnum_req(pr, 'k')
        do j = 1, ng
          vg(j) = (0.5_dp * kpar) * (gv(j) * gv(j))
        end do
      case ('hydrogen')
        hstep = x2 / real(ng + 1, dp)
        do j = 1, ng
          gv(j) = real(j, dp) * hstep
          vg(j) = -1.0_dp / gv(j)
        end do
      case ('poschl_teller')
        gv = linspace(x1, x2, ng)
        lam = jnum_req(pr, 'lambda')
        do j = 1, ng
          vg(j) = (-0.5_dp * lam * (lam + 1.0_dp)) / &
               (cosh(gv(j)) * cosh(gv(j)))
        end do
      end select
      en = jtovec(e, 'eigen_n')
      k = size(en)
      if (k < 1) then
        ok = .false.
        emsg = 'bad eigen_n for ' // trim(nraw)
        return
      end if
      allocate(ev(k))
      call fd_eigenvalues(gv, vg, k, ev, ok, emsg)
      if (.not. ok) return
      call jcopy(eigref, jget(ref, 'eigenvalues'))
      evals = jtovec(eigref, 'value')
      if (size(evals) < k) then
        ok = .false.
        emsg = 'eigenvalues value list too short for ' // trim(nraw)
        return
      end if
      rlt = jnum_req(eigref, 'rel_tol')
      do j = 1, k
        val = evals(j)
        tol = rlt * abs(val)
        rel = abs(ev(j) - val) / abs(val)
        qty = 'E_' // trim(fg(en(j), 12))
        call add_check(r, qty, ev(j), val, tol, rel, rel <= rlt)
      end do
      if (trim(nmx) == 'harmonic') then
        call jcopy(omega, jget(ref, 'omega'))
        if (omega%t == JV_OBJ) then
          kpar = jnum_req(pr, 'k')
          mpar = jnum_req(pr, 'm')
          w = sqrt(kpar / mpar)
          val = jnum_req(omega, 'value')
          tol = jnum_req(omega, 'abs_tol')
          rel = abs(w - val) / abs(val)
          call add_check(r, 'omega', w, val, tol, rel, abs(w - val) <= tol)
        end if
      end if

    case ('yukawa')
      avec = jnum_req(pr, 'A')
      kap = jnum_req(pr, 'kappa')
      hh = 1.0e-5_dp
      rc = jtovec(e, 'r_check')
      call jcopy(force, jget(ref, 'force'))
      do j = 1, size(rc)
        r1 = rc(j)
        vp = (avec * exp(-kap * (r1 + hh))) / (r1 + hh)
        vm = (avec * exp(-kap * (r1 - hh))) / (r1 - hh)
        fn = -(vp - vm) / (2.0_dp * hh)
        fa = (avec * exp(-kap * r1)) * ((1.0_dp / (r1 * r1)) + (kap / r1))
        val = fa
        tol = jnum_req(force, 'abs_tol')
        rel = abs(fn - val) / abs(val)
        qty = 'force_r_' // trim(fg(r1, 6))
        call add_check(r, qty, fn, val, tol, rel, abs(fn - val) <= tol)
      end do
      vp = (avec * exp(-kap * (1.0_dp + hh))) / (1.0_dp + hh)
      vm = (avec * exp(-kap * (1.0_dp - hh))) / (1.0_dp - hh)
      f1 = -(vp - vm) / (2.0_dp * hh)
      val = jnum_req(force, 'value_at_r1')
      tol = jnum_req(force, 'abs_tol')
      rel = abs(f1 - val) / abs(val)
      call add_check(r, 'force_r_1_ref', f1, val, tol, rel, abs(f1 - val) <= tol)

    case default
      ok = .false.
      emsg = 'Unknown potential: ' // trim(nraw)
      return
    end select

    if (r%nchk == 0) then
      ok = .false.
      emsg = 'No checks produced for ' // trim(nraw)
      return
    end if
    ok = .true.
    emsg = ''
  end subroutine

end program
