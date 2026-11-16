program f_runner
  use kinds, only: dp, i8
  use fjson
  use fhandles
  use fkernels
  use fmath
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
  type(jval) :: reg1, reg2, reg3, reg4, reg5, reg6, reg7, reg8
  type(jval) :: reg9, reg10, reg11, reg12, ents, ent, top, ra
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
  call jload('registries/energy_drift.json', reg3, ok2, emsg)
  if (.not. ok2) call jload('../registries/energy_drift.json', reg3, ok2, emsg)
  if (.not. ok2) then
    write(*, '(a)') 'ERROR: ' // trim(emsg)
    error stop 1
  end if
  call jload('registries/cordic.json', reg4, ok2, emsg)
  if (.not. ok2) call jload('../registries/cordic.json', reg4, ok2, emsg)
  if (.not. ok2) then
    write(*, '(a)') 'ERROR: ' // trim(emsg)
    error stop 1
  end if
  call jload('registries/hastings.json', reg5, ok2, emsg)
  if (.not. ok2) call jload('../registries/hastings.json', reg5, ok2, emsg)
  if (.not. ok2) then
    write(*, '(a)') 'ERROR: ' // trim(emsg)
    error stop 1
  end if
  call jload('registries/lcg_lehmer.json', reg6, ok2, emsg)
  if (.not. ok2) call jload('../registries/lcg_lehmer.json', reg6, ok2, emsg)
  if (.not. ok2) then
    write(*, '(a)') 'ERROR: ' // trim(emsg)
    error stop 1
  end if
  call jload('registries/matematica_ancestral.json', reg7, ok2, emsg)
  if (.not. ok2) call jload('../registries/matematica_ancestral.json', reg7, ok2, emsg)
  if (.not. ok2) then
    write(*, '(a)') 'ERROR: ' // trim(emsg)
    error stop 1
  end if
  call jload('registries/prng_legacy.json', reg8, ok2, emsg)
  if (.not. ok2) call jload('../registries/prng_legacy.json', reg8, ok2, emsg)
  if (.not. ok2) then
    write(*, '(a)') 'ERROR: ' // trim(emsg)
    error stop 1
  end if

  call jload('registries/goertzel_fft.json', reg9, ok2, emsg)
  if (.not. ok2) call jload('../registries/goertzel_fft.json', reg9, ok2, emsg)
  if (.not. ok2) then
    write(*, '(a)') 'ERROR: ' // trim(emsg)
    error stop 1
  end if

  call jload('registries/archimedes_pi.json', reg10, ok2, emsg)
  if (.not. ok2) call jload('../registries/archimedes_pi.json', reg10, ok2, emsg)
  if (.not. ok2) then
    write(*, '(a)') 'ERROR: ' // trim(emsg)
    error stop 1
  end if

  call jload('registries/middle_square_weyl.json', reg11, ok2, emsg)
  if (.not. ok2) call jload('../registries/middle_square_weyl.json', reg11, ok2, emsg)
  if (.not. ok2) then
    write(*, '(a)') 'ERROR: ' // trim(emsg)
    error stop 1
  end if

  call jload('registries/hybrid_prng.json', reg12, ok2, emsg)
  if (.not. ok2) call jload('../registries/hybrid_prng.json', reg12, ok2, emsg)
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
  call jcopy(ents, jget(reg3, 'entries'))
  nent = jalen(ents)
  do ei = 1, nent
    call jcopy(ent, jat(ents, ei))
    if (selected(ent, sel)) call run_one(ent, outdir, runs, nrun)
  end do
  call jcopy(ents, jget(reg4, 'entries'))
  nent = jalen(ents)
  do ei = 1, nent
    call jcopy(ent, jat(ents, ei))
    if (selected(ent, sel)) call run_one(ent, outdir, runs, nrun)
  end do
  call jcopy(ents, jget(reg5, 'entries'))
  nent = jalen(ents)
  do ei = 1, nent
    call jcopy(ent, jat(ents, ei))
    if (selected(ent, sel)) call run_one(ent, outdir, runs, nrun)
  end do
  call jcopy(ents, jget(reg6, 'entries'))
  nent = jalen(ents)
  do ei = 1, nent
    call jcopy(ent, jat(ents, ei))
    if (selected(ent, sel)) call run_one(ent, outdir, runs, nrun)
  end do
  call jcopy(ents, jget(reg7, 'entries'))
  nent = jalen(ents)
  do ei = 1, nent
    call jcopy(ent, jat(ents, ei))
    if (selected(ent, sel)) call run_one(ent, outdir, runs, nrun)
  end do
  call jcopy(ents, jget(reg8, 'entries'))
  nent = jalen(ents)
  do ei = 1, nent
    call jcopy(ent, jat(ents, ei))
    if (selected(ent, sel)) call run_one(ent, outdir, runs, nrun)
  end do
  call jcopy(ents, jget(reg9, 'entries'))
  nent = jalen(ents)
  do ei = 1, nent
    call jcopy(ent, jat(ents, ei))
    if (selected(ent, sel)) call run_one(ent, outdir, runs, nrun)
  end do
  call jcopy(ents, jget(reg10, 'entries'))
  nent = jalen(ents)
  do ei = 1, nent
    call jcopy(ent, jat(ents, ei))
    if (selected(ent, sel)) call run_one(ent, outdir, runs, nrun)
  end do
  call jcopy(ents, jget(reg11, 'entries'))
  nent = jalen(ents)
  do ei = 1, nent
    call jcopy(ent, jat(ents, ei))
    if (selected(ent, sel)) call run_one(ent, outdir, runs, nrun)
  end do
  call jcopy(ents, jget(reg12, 'entries'))
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
    case ('energy_drift')
      call eval_energy_drift(ent, r, ok, emsg)
    case ('cordic')
      call eval_cordic(ent, r, ok, emsg)
    case ('hastings')
      call eval_hastings(ent, r, ok, emsg)
    case ('lcg')
      call eval_lcg(ent, r, ok, emsg)
    case ('ancestral')
      call eval_ancestral(ent, r, ok, emsg)
    case ('prng_legacy')
      call eval_prng_legacy(ent, r, ok, emsg)
    case ('dsp_bin')
      call eval_dsp_bin(ent, r, ok, emsg)
    case ('archimedes_pi')
      call eval_archimedes_pi(ent, r, ok, emsg)
    case ('prng_bin')
      call eval_prng_bin(ent, r, ok, emsg)
    case ('prng')
      call eval_hybrid_prng(ent, r, ok, emsg)
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

  subroutine eval_energy_drift(e, r, ok, emsg)
    ! Double-well trajectory: velocity Verlet vs RK4, relative energy drift.
    type(jval), intent(in) :: e
    type(run_t), intent(inout) :: r
    logical, intent(out) :: ok
    character(len=*), intent(out) :: emsg
    type(jval) :: pr, ref, ck
    real(dp) :: a, b, mass, h, T, x, v, E0, acc, accn, en, rel, maxr
    real(dp) :: xv, vv, val, tol, rlt
    real(dp) :: y1, y2, k11, k12, k21, k22, k31, k32, k41, k42
    integer :: n, k
    character(len=64) :: keys(5)
    real(dp) :: vals(5)
    logical :: has_abs

    call jcopy(pr, jget(e, 'params'))
    a = jnum_req(pr, 'a')
    b = jnum_req(pr, 'b')
    mass = jnum_req(e, 'mass')
    h = jnum_req(e, 'h')
    T = jnum_req(e, 'T')
    ! x0/v0 are arrays; extract element 1
    block
      type(jval) :: arr, cell0
      call jcopy(arr, jget(e, 'x0'))
      call jcopy(cell0, jat(arr, 1))
      x = cell0%v
      call jcopy(arr, jget(e, 'v0'))
      call jcopy(cell0, jat(arr, 1))
      v = cell0%v
    end block
    n = nint(T / h)

    E0 = 0.5_dp * mass * v * v + (a * x**4 - b * x**2)

    ! velocity Verlet kick-drift-kick
    xv = x
    vv = v
    maxr = 0.0_dp
    do k = 1, n
      acc = (-4.0_dp * a * xv**3 + 2.0_dp * b * xv) / mass
      xv = xv + vv * h + 0.5_dp * acc * h * h
      accn = (-4.0_dp * a * xv**3 + 2.0_dp * b * xv) / mass
      vv = vv + 0.5_dp * (acc + accn) * h
      en = 0.5_dp * mass * vv * vv + (a * xv**4 - b * xv**2)
      rel = abs(en - E0) / abs(E0)
      if (rel > maxr) maxr = rel
    end do
    vals(2) = rel
    vals(3) = maxr

    ! RK4 on (x, v)
    y1 = x
    y2 = v
    maxr = 0.0_dp
    do k = 1, n
      k11 = y2
      k12 = (-4.0_dp * a * y1**3 + 2.0_dp * b * y1) / mass
      k21 = y2 + 0.5_dp * h * k12
      k22 = (-4.0_dp * a * (y1 + 0.5_dp * h * k11)**3 &
             + 2.0_dp * b * (y1 + 0.5_dp * h * k11)) / mass
      k31 = y2 + 0.5_dp * h * k22
      k32 = (-4.0_dp * a * (y1 + 0.5_dp * h * k21)**3 &
             + 2.0_dp * b * (y1 + 0.5_dp * h * k21)) / mass
      k41 = y2 + h * k32
      k42 = (-4.0_dp * a * (y1 + h * k31)**3 + 2.0_dp * b * (y1 + h * k31)) / mass
      y1 = y1 + (h / 6.0_dp) * (k11 + 2.0_dp * k21 + 2.0_dp * k31 + k41)
      y2 = y2 + (h / 6.0_dp) * (k12 + 2.0_dp * k22 + 2.0_dp * k32 + k42)
      en = 0.5_dp * mass * y2 * y2 + (a * y1**4 - b * y1**2)
      rel = abs(en - E0) / abs(E0)
      if (rel > maxr) maxr = rel
    end do
    vals(4) = rel
    vals(5) = maxr
    vals(1) = E0

    keys(1) = 'E0'
    keys(2) = 'verlet_final_rel_drift'
    keys(3) = 'verlet_max_rel_drift'
    keys(4) = 'rk4_final_rel_drift'
    keys(5) = 'rk4_max_rel_drift'

    call jcopy(ref, jget(e, 'ref'))
    do k = 1, 5
      call jcopy(ck, jget(ref, trim(keys(k))))
      if (ck%t /= JV_OBJ) then
        call add_obs(r, trim(keys(k)) // '_observed', vals(k))
        cycle
      end if
      val = jnum_req(ck, 'value')
      has_abs = jhas(ck, 'abs_tol')
      if (has_abs) then
        tol = jnum_req(ck, 'abs_tol')
      else
        rlt = jnum_req(ck, 'rel_tol')
        tol = rlt * abs(val)
      end if
      rel = abs(vals(k) - val) / max(abs(val), epsilon(1.0_dp))
      call add_check(r, trim(keys(k)), vals(k), val, tol, rel, &
                     abs(vals(k) - val) <= tol)
    end do
    ok = .true.
    emsg = ''
  end subroutine

  function mtag(m) result(s)
    real(dp), intent(in) :: m
    character(len=32) :: s
    character(len=64) :: buf
    integer :: n
    write(buf, '(f20.6)') m
    buf = adjustl(buf)
    n = len_trim(buf)
    if (index(buf, '.') > 0) then
      do while (n > 1)
        if (buf(n:n) == '0') then
          n = n - 1
        else
          exit
        end if
      end do
      if (n >= 1) then
        if (buf(n:n) == '.') n = n - 1
      end if
    end if
    s = buf(1:n)
  end function

  function ref_tol(ref, key) result(tol)
    type(jval), intent(in) :: ref
    character(len=*), intent(in) :: key
    real(dp) :: tol
    type(jval) :: ck
    call jcopy(ck, jget(ref, key))
    if (jhas(ck, 'abs_tol')) then
      tol = jnum_req(ck, 'abs_tol')
    else
      tol = jnum_req(ck, 'rel_tol') * abs(jnum_req(ck, 'value'))
    end if
  end function

  subroutine check_ref(r, ref, key, computed)
    type(run_t), intent(inout) :: r
    type(jval), intent(in) :: ref
    character(len=*), intent(in) :: key
    real(dp), intent(in) :: computed
    type(jval) :: ck, cv
    real(dp) :: val, tol, rel
    call jcopy(ck, jget(ref, key))
    if (ck%t /= JV_OBJ) then
      call add_obs(r, trim(key) // '_observed', computed)
      return
    end if
    if (.not. jhas(ck, 'value')) then
      call add_obs(r, trim(key) // '_observed', computed)
      return
    end if
    call jcopy(cv, jget(ck, 'value'))
    if (cv%t /= JV_NUM) then
      ! value ausente o null: comparar solo registra el valor (pass=true)
      call add_obs(r, trim(key) // '_observed', computed)
      return
    end if
    tol = ref_tol(ref, key)
    val = jnum_req(ck, 'value')
    rel = abs(computed - val) / max(abs(val), 1.0_dp)
    call add_check(r, trim(key), computed, val, tol, rel, &
                   abs(computed - val) <= tol)
  end subroutine

  subroutine vec_from(o, v)
    type(jval), intent(in) :: o
    real(dp), allocatable, intent(out) :: v(:)
    type(jval) :: cell
    integer :: j, n
    n = jalen(o)
    allocate(v(max(n, 1)))
    do j = 1, n
      call jcopy(cell, jat(o, j))
      v(j) = cell%v
    end do
  end subroutine

  subroutine mat_from(o, m)
    type(jval), intent(in) :: o
    real(dp), allocatable, intent(out) :: m(:, :)
    type(jval) :: row, cell
    integer :: i, j, nr, nc
    nr = jalen(o)
    nc = 0
    do i = 1, nr
      call jcopy(row, jat(o, i))
      nc = max(nc, jalen(row))
    end do
    allocate(m(max(nr, 1), max(nc, 1)))
    do i = 1, nr
      call jcopy(row, jat(o, i))
      do j = 1, jalen(row)
        call jcopy(cell, jat(row, j))
        m(i, j) = cell%v
      end do
    end do
  end subroutine

  subroutine eval_cordic(e, r, ok, emsg)
    type(jval), intent(in) :: e
    type(run_t), intent(inout) :: r
    logical, intent(out) :: ok
    character(len=*), intent(out) :: emsg
    type(jval) :: pr, ref, mv, cell
    real(dp), allocatable :: atan_tab(:)
    real(dp) :: kg, xr, yr, zr, theta, theta_r, xq, yq, eq, ec, es
    real(dp) :: max_cos, max_sin, max_res, max_q15, cos_tol, sin_tol, tol
    integer :: n_iter, n_samples, nmult, i, k, theta_f
    character(len=64) :: tag

    call jcopy(pr, jget(e, 'params'))
    call jcopy(ref, jget(e, 'ref'))
    n_iter = jint(pr, 'iterations', 0)
    n_samples = jint(pr, 'samples', 0)
    call jcopy(mv, jget(pr, 'angle_multipliers'))
    nmult = jalen(mv)
    allocate(atan_tab(max(n_iter, 1)))
    do i = 1, n_iter
      atan_tab(i) = atan(2.0_dp**(-(i - 1)))
    end do
    kg = cordic_gain(n_iter)
    call check_ref(r, ref, 'gain_K', kg)
    cos_tol = ref_tol(ref, 'cos_err_max')
    sin_tol = ref_tol(ref, 'sin_err_max')
    max_cos = 0.0_dp
    max_sin = 0.0_dp
    max_res = 0.0_dp
    max_q15 = 0.0_dp
    do i = 1, nmult
      call jcopy(cell, jat(mv, i))
      theta = cell%v * pi_value()
      call cordic_rot(theta, atan_tab, xr, yr, zr)
      ec = abs(xr / kg - cos(theta))
      es = abs(yr / kg - sin(theta))
      tag = 'm' // trim(mtag(cell%v))
      call add_check(r, 'cos_err_' // trim(tag), ec, 0.0_dp, cos_tol, ec, ec <= cos_tol)
      call add_check(r, 'sin_err_' // trim(tag), es, 0.0_dp, sin_tol, es, es <= sin_tol)
      theta_f = int(theta / pi_value() * 32768.0_dp)
      theta_r = (real(theta_f, dp) / 32768.0_dp) * pi_value()
      call cordic_rot(theta_r, atan_tab, xq, yq, zr)
      xq = real(int(xq * 32768.0_dp), dp)
      yq = real(int(yq * 32768.0_dp), dp)
      eq = max(abs(xq / (32768.0_dp * kg) - cos(theta_r)), &
               abs(yq / (32768.0_dp * kg) - sin(theta_r)))
      if (eq > max_q15) max_q15 = eq
    end do
    do k = 0, n_samples - 1
      theta = -pi_value() + 2.0_dp * pi_value() * real(k, dp) / &
              real(n_samples - 1, dp)
      call cordic_rot(theta, atan_tab, xr, yr, zr)
      max_cos = max(max_cos, abs(xr / kg - cos(theta)))
      max_sin = max(max_sin, abs(yr / kg - sin(theta)))
      max_res = max(max_res, abs(zr))
    end do
    call add_check(r, 'cos_err_max', max_cos, 0.0_dp, cos_tol, max_cos, max_cos <= cos_tol)
    call add_check(r, 'sin_err_max', max_sin, 0.0_dp, sin_tol, max_sin, max_sin <= sin_tol)
    tol = ref_tol(ref, 'residual_max')
    call add_check(r, 'residual_max', max_res, 0.0_dp, tol, max_res, max_res <= tol)
    tol = ref_tol(ref, 'q15_err_max')
    call add_check(r, 'q15_err_max', max_q15, 0.0_dp, tol, max_q15, max_q15 <= tol)
    ok = .true.
    emsg = ''
  end subroutine

  subroutine eval_hastings(e, r, ok, emsg)
    type(jval), intent(in) :: e
    type(run_t), intent(inout) :: r
    logical, intent(out) :: ok
    character(len=*), intent(out) :: emsg
    type(jval) :: pr, ref, rr
    real(dp), allocatable :: rng(:)
    real(dp) :: err, tol, v, a, b, xx, ex, pi2, pi6
    integer :: n_samples, k

    call jcopy(pr, jget(e, 'params'))
    call jcopy(ref, jget(e, 'ref'))
    n_samples = jint(pr, 'samples', 0)
    pi2 = pi_value() / 2.0_dp
    pi6 = pi_value() / 6.0_dp

    call jcopy(rr, jget(pr, 'sin_range'))
    call vec_from(rr, rng)
    a = rng(1)
    b = rng(2)
    err = 0.0_dp
    do k = 0, n_samples - 1
      xx = a + (b - a) * real(k, dp) / real(n_samples - 1, dp)
      err = max(err, abs(hastings_sin(xx) - sin(xx)))
    end do
    tol = ref_tol(ref, 'sin_err_max')
    call add_check(r, 'sin_err_max', err, 0.0_dp, tol, err, err <= tol)

    err = abs(hastings_sin(pi6) - sin(pi6))
    tol = ref_tol(ref, 'sin_err_pi6')
    call add_check(r, 'sin_err_pi6', err, 0.0_dp, tol, err, err <= tol)

    err = abs(hastings_sin(pi2) - sin(pi2))
    tol = ref_tol(ref, 'sin_err_pi2')
    call add_check(r, 'sin_err_pi2', err, 0.0_dp, tol, err, err <= tol)

    v = hastings_sin(pi2)
    call check_ref(r, ref, 'sin_val_pi2', v)

    call jcopy(rr, jget(pr, 'atan_range'))
    call vec_from(rr, rng)
    a = rng(1)
    b = rng(2)
    err = 0.0_dp
    do k = 0, n_samples - 1
      xx = a + (b - a) * real(k, dp) / real(n_samples - 1, dp)
      err = max(err, abs(hastings_atan(xx) - atan(xx)))
    end do
    tol = ref_tol(ref, 'atan_err_max')
    call add_check(r, 'atan_err_max', err, 0.0_dp, tol, err, err <= tol)

    err = abs(hastings_atan(1.0_dp) - atan(1.0_dp))
    tol = ref_tol(ref, 'atan_err_one')
    call add_check(r, 'atan_err_one', err, 0.0_dp, tol, err, err <= tol)

    call jcopy(rr, jget(pr, 'exp_range'))
    call vec_from(rr, rng)
    a = rng(1)
    b = rng(2)
    err = 0.0_dp
    do k = 0, n_samples - 1
      xx = a + (b - a) * real(k, dp) / real(n_samples - 1, dp)
      ex = exp(xx)
      err = max(err, abs(hastings_exp(xx) - ex) / ex)
    end do
    tol = ref_tol(ref, 'exp_rel_err_max')
    call add_check(r, 'exp_rel_err_max', err, 0.0_dp, tol, err, err <= tol)

    err = abs(hastings_exp(1.0_dp) - exp(1.0_dp))
    tol = ref_tol(ref, 'exp_err_one')
    call add_check(r, 'exp_err_one', err, 0.0_dp, tol, err, err <= tol)

    v = hastings_exp(1.0_dp)
    call check_ref(r, ref, 'exp_val_one', v)

    call jcopy(rr, jget(pr, 'log1p_range'))
    call vec_from(rr, rng)
    a = rng(1)
    b = rng(2)
    err = 0.0_dp
    do k = 0, n_samples - 1
      xx = a + (b - a) * real(k, dp) / real(n_samples - 1, dp)
      err = max(err, abs(hastings_log1p(xx) - log(1.0_dp + xx)))
    end do
    tol = ref_tol(ref, 'log1p_err_max')
    call add_check(r, 'log1p_err_max', err, 0.0_dp, tol, err, err <= tol)

    err = abs(hastings_log1p(1.0_dp) - log(1.0_dp + 1.0_dp))
    tol = ref_tol(ref, 'log1p_err_one')
    call add_check(r, 'log1p_err_one', err, 0.0_dp, tol, err, err <= tol)

    v = hastings_log1p(1.0_dp)
    call check_ref(r, ref, 'log1p_val_one', v)
    ok = .true.
    emsg = ''
  end subroutine

  subroutine eval_lcg(e, r, ok, emsg)
    type(jval), intent(in) :: e
    type(run_t), intent(inout) :: r
    logical, intent(out) :: ok
    character(len=*), intent(out) :: emsg
    type(jval) :: pr, ref, rc, mc, pco, pcmo
    integer(i8) :: a_r, c_r, m_r, a_m, c_m, m_m
    integer(i8) :: x1, x2, x3, seq_r(4), seq_m(4), uniq(256)
    integer(i8) :: a_p, c_p, m_p, x_p, a_pm, c_pm, m_pm, x_pm
    integer :: box, n_samples, i, j, k, nuniq, count_r, count_m
    integer :: per_r, per_m, found
    real(dp) :: cc0, cc1, cc2, kk, ki, max_dev, min_norm_r, min_norm_m, tol
    character(len=64) :: nmx

    call jcopy(pr, jget(e, 'params'))
    call jcopy(ref, jget(e, 'ref'))
    call jcopy(rc, jget(pr, 'randu'))
    a_r = int(jnum_req(rc, 'a'), i8)
    c_r = int(jnum_req(rc, 'c'), i8)
    m_r = int(jnum_req(rc, 'm'), i8)
    call jcopy(mc, jget(pr, 'minstd'))
    a_m = int(jnum_req(mc, 'a'), i8)
    c_m = int(jnum_req(mc, 'c'), i8)
    m_m = int(jnum_req(mc, 'm'), i8)
    box = jint(pr, 'normal_box_h', 0)
    n_samples = jint(pr, 'marsaglia_samples', 0)
    call jcopy(pco, jget(pr, 'marsaglia_coeffs'))
    cc0 = atnum(pco, 1)
    cc1 = atnum(pco, 2)
    cc2 = atnum(pco, 3)

    seq_r(1) = lcg_next(a_r, c_r, m_r, 1_i8)
    do i = 2, 4
      seq_r(i) = lcg_next(a_r, c_r, m_r, seq_r(i - 1))
    end do
    seq_m(1) = lcg_next(a_m, c_m, m_m, 1_i8)
    do i = 2, 4
      seq_m(i) = lcg_next(a_m, c_m, m_m, seq_m(i - 1))
    end do
    do i = 1, 4
      write(nmx, '("randu_x",i0)') i
      call check_ref(r, ref, trim(nmx), real(seq_r(i), dp))
      write(nmx, '("minstd_x",i0)') i
      call check_ref(r, ref, trim(nmx), real(seq_m(i), dp))
    end do

    x1 = lcg_next(a_r, c_r, m_r, 1_i8)
    x2 = lcg_next(a_r, c_r, m_r, x1)
    x3 = lcg_next(a_r, c_r, m_r, x2)
    nuniq = 0
    max_dev = 0.0_dp
    do k = 1, n_samples
      kk = (cc0 * real(x1, dp) + cc1 * real(x2, dp) + cc2 * real(x3, dp)) / &
           real(m_r, dp)
      ki = real(nint(kk), dp)
      max_dev = max(max_dev, abs(kk - ki))
      found = 0
      do j = 1, nuniq
        if (real(uniq(j), dp) == ki) then
          found = 1
          exit
        end if
      end do
      if (found == 0 .and. nuniq < 256) then
        nuniq = nuniq + 1
        uniq(nuniq) = int(ki, i8)
      end if
      x1 = x2
      x2 = x3
      x3 = lcg_next(a_r, c_r, m_r, x3)
    end do
    call check_ref(r, ref, 'randu_plane_count', real(nuniq, dp))
    tol = ref_tol(ref, 'marsaglia_integrality')
    call add_check(r, 'marsaglia_integrality', max_dev, 0.0_dp, tol, max_dev, max_dev <= tol)

    call lcg_normals(a_r, m_r, box, count_r, min_norm_r)
    call lcg_normals(a_m, m_m, box, count_m, min_norm_m)
    call check_ref(r, ref, 'randu_normal_count', real(count_r, dp))
    call check_ref(r, ref, 'minstd_normal_count', real(count_m, dp))
    call check_ref(r, ref, 'randu_min_normal_norm', min_norm_r)

    call jcopy(pco, jget(pr, 'period_case'))
    a_p = int(jnum_req(pco, 'a'), i8)
    c_p = int(jnum_req(pco, 'c'), i8)
    m_p = int(jnum_req(pco, 'm'), i8)
    x_p = int(jnum_req(pco, 'x0'), i8)
    call jcopy(pcmo, jget(pr, 'period_case_minstd'))
    a_pm = int(jnum_req(pcmo, 'a'), i8)
    c_pm = int(jnum_req(pcmo, 'c'), i8)
    m_pm = int(jnum_req(pcmo, 'm'), i8)
    x_pm = int(jnum_req(pcmo, 'x0'), i8)
    per_r = lcg_period(a_p, c_p, m_p, x_p)
    per_m = lcg_period(a_pm, c_pm, m_pm, x_pm)
    call check_ref(r, ref, 'randu_period', real(per_r, dp))
    call check_ref(r, ref, 'minstd_period', real(per_m, dp))
    ok = .true.
    emsg = ''
  end subroutine

  subroutine eval_ancestral(e, r, ok, emsg)
    type(jval), intent(in) :: e
    type(run_t), intent(inout) :: r
    logical, intent(out) :: ok
    character(len=*), intent(out) :: emsg
    type(jval) :: pr, ref, sqp, arp, tcrp, lins, arr
    real(dp), allocatable :: amat(:, :), bvec(:), dens(:)
    real(dp), allocatable :: xvec(:)
    integer(i8), allocatable :: resid(:), modu(:), densi(:), digs(:)
    integer(i8) :: mult(2), frac(2), x_crt, mm_crt, vval
    integer :: nterms, ndig, i, j, n, ndig_max, nres, nmod
    real(dp) :: xval, sq, err, err_tab, err_sin, theta, paso, s_ary, idx_f
    real(dp) :: acc, resid_max, sumv, tol
    logical :: cok
    character(len=64) :: nmx

    call jcopy(pr, jget(e, 'params'))
    call jcopy(ref, jget(e, 'ref'))
    call jcopy(arr, jget(pr, 'egyptian_multiply'))
    mult(1) = int(atnum(arr, 1), i8)
    mult(2) = int(atnum(arr, 2), i8)
    call jcopy(arr, jget(pr, 'egyptian_fraction'))
    frac(1) = int(atnum(arr, 1), i8)
    frac(2) = int(atnum(arr, 2), i8)

    call check_ref(r, ref, 'egyptian_product', real(egyptian_mult(mult(1), mult(2)), dp))
    call check_ref(r, ref, 'egyptian_product_signed', &
                   real(egyptian_mult(-mult(1), mult(2)), dp))

    allocate(densi(32))
    call egyptian_fraction(frac(1), frac(2), 32, densi, nterms)
    allocate(dens(max(nterms, 1)))
    do i = 1, nterms
      dens(i) = real(densi(i), dp)
    end do
    do i = 1, min(nterms, 2)
      write(nmx, '("egyptian_den_",i0)') i
      call check_ref(r, ref, trim(nmx), real(densi(i), dp))
    end do
    call check_ref(r, ref, 'egyptian_terms', real(nterms, dp))
    sumv = 0.0_dp
    do i = 1, nterms
      sumv = sumv + 1.0_dp / dens(i)
    end do
    err = abs(sumv - real(frac(1), dp) / real(frac(2), dp))
    tol = ref_tol(ref, 'egyptian_sum_err')
    call add_check(r, 'egyptian_sum_err', err, 0.0_dp, tol, err, err <= tol)

    call jcopy(sqp, jget(pr, 'babylonian_sqrt'))
    xval = jnum_req(sqp, 'x')
    sq = babylonian_sqrt(xval, jint(sqp, 'iterations', 0))
    err = abs(sq - sqrt(xval))
    tol = ref_tol(ref, 'babylonian_sqrt_err')
    call add_check(r, 'babylonian_sqrt_err', err, 0.0_dp, tol, err, err <= tol)

    call jcopy(arp, jget(pr, 'aryabhata'))
    theta = jnum_req(arp, 'theta')
    ndig_max = jint(arp, 'steps', 360)
    paso = 2.0_dp * pi_value() / real(ndig_max, dp)
    s_ary = aryabhata_sin(theta, ndig_max)
    idx_f = real(int(theta / paso), dp)
    err_tab = abs(s_ary - sin(idx_f * paso))
    tol = ref_tol(ref, 'aryabhata_table_err')
    call add_check(r, 'aryabhata_table_err', err_tab, 0.0_dp, tol, err_tab, err_tab <= tol)
    err_sin = abs(s_ary - sin(theta))
    tol = ref_tol(ref, 'aryabhata_sin_err')
    call add_check(r, 'aryabhata_sin_err', err_sin, 0.0_dp, tol, err_sin, err_sin <= tol)

    call jcopy(tcrp, jget(pr, 'tcr'))
    call jcopy(arr, jget(tcrp, 'residues'))
    nres = jalen(arr)
    allocate(resid(max(nres, 1)))
    do i = 1, nres
      resid(i) = int(atnum(arr, i), i8)
    end do
    call jcopy(arr, jget(tcrp, 'moduli'))
    nmod = jalen(arr)
    allocate(modu(max(nmod, 1)))
    do i = 1, nmod
      modu(i) = int(atnum(arr, i), i8)
    end do
    call chinese_remainder(resid(1:nres), modu(1:nmod), x_crt, mm_crt, cok)
    call check_ref(r, ref, 'tcr_x', real(x_crt, dp))
    call check_ref(r, ref, 'tcr_M', real(mm_crt, dp))

    vval = int(jnum_req(pr, 'vigesimal'), i8)
    allocate(digs(16))
    call vigesimal(vval, 16, digs, ndig)
    do i = 1, min(ndig, 3)
      write(nmx, '("vigesimal_d",i0)') i - 1
      call check_ref(r, ref, trim(nmx), real(digs(i), dp))
    end do
    call check_ref(r, ref, 'vigesimal_digits', real(ndig, dp))
    sumv = 0.0_dp
    do i = 1, ndig
      sumv = sumv + real(digs(i), dp) * 20.0_dp**real(i - 1, dp)
    end do
    call check_ref(r, ref, 'vigesimal_value', sumv)

    call jcopy(lins, jget(pr, 'linear_system'))
    call jcopy(arr, jget(lins, 'A'))
    call mat_from(arr, amat)
    call jcopy(arr, jget(lins, 'b'))
    call vec_from(arr, bvec)
    n = size(bvec)
    allocate(xvec(max(n, 1)))
    call chinese_elimination(amat, bvec, xvec)
    do i = 1, n
      write(nmx, '("linear_x",i0)') i
      call check_ref(r, ref, trim(nmx), xvec(i))
    end do
    resid_max = 0.0_dp
    do i = 1, n
      acc = 0.0_dp
      do j = 1, n
        acc = acc + amat(i, j) * xvec(j)
      end do
      resid_max = max(resid_max, abs(acc - bvec(i)))
    end do
    tol = ref_tol(ref, 'linear_residual_max')
    call add_check(r, 'linear_residual_max', resid_max, 0.0_dp, tol, resid_max, &
                   resid_max <= tol)
    ok = .true.
    emsg = ''
  end subroutine

  subroutine eval_prng_legacy(e, r, ok, emsg)
    type(jval), intent(in) :: e
    type(run_t), intent(inout) :: r
    logical, intent(out) :: ok
    character(len=*), intent(out) :: emsg
    type(jval) :: pr, ref, arr
    integer(i8) :: x, y, seed, mults(3), mods(3), sstate(3), pstate(3)
    integer(i8) :: bad(16), periods(3), lcmv, g
    integer :: digits, nseq, nstates, nseeds, nsamples, i, k, s, nbad, n
    integer :: nfixed_nz, nout, best, best_seed
    logical :: isbad, zero_fixed
    logical, allocatable :: seen(:), seen2(:), seenout(:)
    real(dp) :: rval, r16, rmin, rmax, maxdiff, viol
    character(len=64) :: eno, nmx

    call jcopy(pr, jget(e, 'params'))
    call jcopy(ref, jget(e, 'ref'))
    eno = jstr(e, 'name', '')

    if (trim(eno) == 'middle_square_vn') then
      digits = jint(pr, 'digits', 4)
      seed = int(jnum_req(pr, 'seed'), i8)
      nseq = jint(pr, 'seq_values', 0)
      nstates = 10 ** digits
      x = seed
      do k = 1, nseq
        x = ms_step(x, digits)
        write(nmx, '("seq_",i0)') k
        call check_ref(r, ref, trim(nmx), real(x, dp))
      end do

      nfixed_nz = 0
      zero_fixed = .false.
      do s = 0, nstates - 1
        y = ms_step(int(s, i8), digits)
        if (y == int(s, i8) .and. s /= 0) nfixed_nz = nfixed_nz + 1
      end do
      if (ms_step(0_i8, digits) == 0_i8) zero_fixed = .true.
      call check_ref(r, ref, 'fixed_points_nonzero_count', real(nfixed_nz, dp))
      call check_ref(r, ref, 'zero_is_fixed', merge(1.0_dp, 0.0_dp, zero_fixed))

      allocate(seenout(nstates))
      seenout = .false.
      nout = 0
      do s = 0, nstates - 1
        y = ms_step(int(s, i8), digits)
        if (.not. seenout(int(y) + 1)) then
          seenout(int(y) + 1) = .true.
          nout = nout + 1
        end if
      end do
      write(nmx, '("distinct_outputs_n",i0)') digits
      call check_ref(r, ref, trim(nmx), real(nout, dp))

      allocate(seen(nstates))
      seen = .false.
      x = seed
      k = 0
      do while (.not. seen(int(x) + 1))
        seen(int(x) + 1) = .true.
        x = ms_step(x, digits)
        k = k + 1
      end do
      call check_ref(r, ref, 'period_from_seed', real(k, dp))

      call jcopy(arr, jget(pr, 'bad_two_digit_targets'))
      nbad = jalen(arr)
      do i = 1, nbad
        bad(i) = int(atnum(arr, i), i8)
      end do
      nseeds = jint(pr, 'two_digit_seeds', 0)
      allocate(seen2(100))
      best = 0
      best_seed = -1
      do s = 0, nseeds - 1
        x = int(s, i8)
        seen2 = .false.
        n = 0
        do while (.true.)
          isbad = .false.
          do i = 1, nbad
            if (x == bad(i)) isbad = .true.
          end do
          if (isbad) exit
          if (seen2(int(x) + 1)) exit
          seen2(int(x) + 1) = .true.
          x = ms_step(x, 2)
          n = n + 1
        end do
        if (n > best) then
          best = n
          best_seed = s
        end if
      end do
      call check_ref(r, ref, 'two_digit_max_length', real(best, dp))
      call check_ref(r, ref, 'two_digit_max_length_seed', real(best_seed, dp))
    else if (trim(eno) == 'wichmann_hill') then
      call jcopy(arr, jget(pr, 'seeds'))
      do i = 1, 3
        sstate(i) = int(atnum(arr, i), i8)
      end do
      call jcopy(arr, jget(pr, 'multipliers'))
      do i = 1, 3
        mults(i) = int(atnum(arr, i), i8)
      end do
      call jcopy(arr, jget(pr, 'moduli'))
      do i = 1, 3
        mods(i) = int(atnum(arr, i), i8)
      end do
      nsamples = jint(pr, 'samples', 0)

      do i = 1, 3
        allocate(seen(int(mods(i))))
        seen = .false.
        x = 1_i8
        k = 0
        do while (.not. seen(int(x) + 1))
          seen(int(x) + 1) = .true.
          x = modulo(mults(i) * x, mods(i))
          k = k + 1
        end do
        periods(i) = int(k, i8)
        write(nmx, '("sub_period_",i0)') i
        call check_ref(r, ref, trim(nmx), real(periods(i), dp))
        deallocate(seen)
      end do
      lcmv = periods(1)
      do i = 2, 3
        g = igcd(lcmv, periods(i))
        lcmv = (lcmv / g) * periods(i)
      end do
      call check_ref(r, ref, 'combined_period_lcm', real(lcmv, dp))

      pstate = sstate
      rmin = huge(1.0_dp)
      rmax = -huge(1.0_dp)
      maxdiff = 0.0_dp
      viol = 0.0_dp
      do k = 1, nsamples
        call wh_step(sstate, mults, mods, rval)
        call wh_step16(pstate, mults, mods, r16)
        maxdiff = max(maxdiff, abs(rval - r16))
        rmin = min(rmin, rval)
        rmax = max(rmax, rval)
        viol = max(viol, max(0.0_dp, rval - 1.0_dp) + max(0.0_dp, -rval))
        if (k <= 4) then
          write(nmx, '("r_",i0)') k - 1
          call check_ref(r, ref, trim(nmx), rval)
        end if
      end do
      call check_ref(r, ref, 'r_min', rmin)
      call check_ref(r, ref, 'r_max', rmax)
      call check_ref(r, ref, 'r_range_violation', viol)
      call check_ref(r, ref, 'variant16_max_diff', maxdiff)
    else
      ok = .false.
      emsg = 'unknown prng_legacy entry: ' // trim(eno)
      return
    end if
    ok = .true.
    emsg = ''
  end subroutine

  ! ---------------------------------------------------------------- dsp_bin
  ! Python: registry_run.py eval_dsp_bin (396-442). Firma sin outdir.
  ! Checks: {K,alpha,beta}_k<k> con ref N*A/2, tol rt*ref; variant_spread.
  ! Todos null-ref = add_check con ok=.true. (patron mkcheck(None)).
  subroutine eval_dsp_bin(e, r, ok, emsg)
    type(jval), intent(in) :: e
    type(run_t), intent(inout) :: r
    logical, intent(out) :: ok
    character(len=*), intent(out) :: emsg
    integer, parameter :: MAXT = 64
    integer :: n, nt, i, k, kk, ntmp, ntmp2, ntmp3
    integer :: tone_k(MAXT)
    real(dp) :: tone_a(MAXT), tone_phi(MAXT)
    real(dp), allocatable :: x(:)
    real(dp) :: re_s, im_s, mag2, df, phase, arg, dphi
    real(dp) :: rt, st
    real(dp), parameter :: PI2 = 6.283185307179586476925286766559_dp
    real(dp), allocatable :: tmpk(:), tmpa(:), tmpp(:)
    real(dp) :: w, c, refv, mk, ma, mb
    real(dp) :: spread
    n = jint(e, "N", 8192)
    allocate(x(0:n))
    nt = min(size(tone_k), 3)
    tmpk = jtovec(e, "tone_k")
    tmpa = jtovec(e, "tone_A")
    tmpp = jtovec(e, "tone_phi")
    ntmp = min(size(tmpk), MAXT)
    ntmp2 = min(size(tmpa), MAXT)
    ntmp3 = min(size(tmpp), MAXT)
    do i = 1, ntmp
       tone_k(i) = nint(tmpk(i))
    end do
    do i = 1, ntmp2
       tone_a(i) = tmpa(i)
    end do
    do i = 1, ntmp3
       tone_phi(i) = tmpp(i)
    end do
    x = 0.0_dp
    do i = 0, n-1
      do kk = 1, nt
        x(i) = x(i) + tone_a(kk)*cos(PI2*tone_k(kk)*real(i, dp)/real(n, dp) + tone_phi(kk))
      end do
    end do
    rt = jnum(e, "ref.magnitude.rel_tol", 1.0e-8_dp)
    st = jnum(e, "ref.variant_spread.abs_tol", 1.0e-8_dp)
    spread = 0.0_dp
    do kk = 1, nt
      k = tone_k(kk)
      w = PI2*real(k, dp)/real(n, dp)
      c = 2.0_dp*cos(w)
      refv = real(n, dp)*tone_a(kk)/2.0_dp
      mk = goertzel_k(x, c)
      ma = goertzel_alpha(x, c)
      mb = goertzel_beta(x, c)
      call add_check(r, "K_k"//trim(itostr(k)), mk, refv, rt*refv, abs(mk-refv)/refv, abs(mk-refv)/refv <= rt)
      call add_check(r, "alpha_k"//trim(itostr(k)), ma, refv, rt*refv, abs(ma-refv)/refv, abs(ma-refv)/refv <= rt)
      call add_check(r, "beta_k"//trim(itostr(k)), mb, refv, rt*refv, abs(mb-refv)/refv, abs(mb-refv)/refv <= rt)
      spread = max(spread, abs(ma-mk)/refv, abs(mb-mk)/refv)
    end do
    call add_check(r, "variant_spread", spread, 0.0_dp, st, spread, spread <= st)
    ok = .true.
    emsg = ''
  end subroutine

  ! ------------------------------------------------------- archimedes_pi
  ! Python: registry_run.py eval_archimedes_pi (730-817). Firma sin outdir.
  ! pi_est = r(n - levels) tras romberg in-place. Checks pi_stable/pi_naive
  ! con pi_ref/tol; variant_spread ref 0.0.
  subroutine eval_archimedes_pi(e, r, ok, emsg)
    type(jval), intent(in) :: e
    type(run_t), intent(inout) :: r
    logical, intent(out) :: ok
    character(len=*), intent(out) :: emsg
    integer :: iterations, levels, i, j, m, n, na
    real(dp) :: t, pi_ref, tol, pi_s, pi_n, diff
    real(dp), allocatable :: a(:), rr(:)
    real(dp) :: denom, st
    iterations = jint(e, "params.iterations", 10)
    levels = jint(e, "params.levels", 3)
    pi_ref = jnum(e, "ref.pi.value", 3.141592653589793_dp)
    tol = jnum(e, "ref.pi.abs_tol", 1.0e-8_dp)
    st = jnum(e, "ref.variant_spread.abs_tol", 1.0e-8_dp)
    na = iterations + 1
    allocate(a(na), rr(na))
    ! stable
    t = 1.0_dp/sqrt(3.0_dp)
    n = 6
    a(1) = real(n, dp)*t
    do i = 1, iterations
      t = t/(sqrt(1.0_dp+t*t)+1.0_dp)
      n = 2*n
      a(i+1) = real(n, dp)*t
    end do
    rr(1:na) = a(1:na)
    do j = 1, levels
      denom = real(4**j, dp) - 1.0_dp
      m = na - j + 1
      do i = 1, m-1
        rr(i) = rr(i+1) + (rr(i+1)-rr(i))/denom
      end do
    end do
    pi_s = rr(na - levels)
    ! naive
    t = 1.0_dp/sqrt(3.0_dp)
    n = 6
    a(1) = real(n, dp)*t
    do i = 1, iterations
      t = (sqrt(1.0_dp+t*t)-1.0_dp)/t
      n = 2*n
      a(i+1) = real(n, dp)*t
    end do
    rr(1:na) = a(1:na)
    do j = 1, levels
      denom = real(4**j, dp) - 1.0_dp
      m = na - j + 1
      do i = 1, m-1
        rr(i) = rr(i+1) + (rr(i+1)-rr(i))/denom
      end do
    end do
    pi_n = rr(na - levels)
    diff = abs(pi_s - pi_n)/abs(pi_ref)
    call add_check(r, "pi_stable", pi_s, pi_ref, tol, abs(pi_s-pi_ref)/pi_ref, abs(pi_s-pi_ref) <= tol)
    call add_check(r, "pi_naive", pi_n, pi_ref, tol, abs(pi_n-pi_ref)/pi_ref, abs(pi_n-pi_ref) <= tol)
    call add_check(r, "variant_spread", diff, 0.0_dp, st, diff, diff <= st)
    ok = .true.
    emsg = ''
  end subroutine

  ! ---------------------------------------------------------- prng_bin
  ! Python: registry_run.py eval_middle_square_weyl (452-582). Sin outdir.
  ! weyl/tcr: x = mid_ext((x*x+w)**2) (doble cuadrado, como _mid4(x*x+w)).
  ! seen 0-based global; x0 en 0..modulus-1; hit-flag para ciclo.
  ! weyl_period/tcr_period/variant_spread: null-ref condicional.
  ! mid_ext(q, digits) extrae digitos medios de q (NO cuadra); ext_gcd SUBRUTINA.
  subroutine eval_prng_bin(e, r, ok, emsg)
    type(jval), intent(in) :: e
    type(run_t), intent(inout) :: r
    logical, intent(out) :: ok
    character(len=*), intent(out) :: emsg
    integer :: digits, modulus, weyl_s, tcr_p1, tcr_p2, max_iter
    integer :: i, j, x, w, w1, w2, k, idx, tr, trw, trcrt
    integer :: n_cycles, zero_sink, weyl_period, tcr_period
    integer(i8) :: inv, g, yy, q
    integer :: m
    integer, allocatable :: seen(:), seenw(:), seen_crt(:)
    integer, allocatable :: tr_arr(:)
    real(dp) :: zs_tol, zs_rel, zs_err
    real(dp) :: ref_c, tol_c, ref_z, st, spread, wp_ref, wp_tol, tp_ref, tp_tol
    logical :: hit
    digits = jint(e, "params.digits", 4)
    modulus = jint(e, "params.modulus", 10000)
    weyl_s = jint(e, "params.weyl_s", 6329)
    tcr_p1 = jint(e, "params.tcr_p1", 9973)
    tcr_p2 = jint(e, "params.tcr_p2", 9967)
    max_iter = jint(e, "params.max_iter", 1000000)
    ref_c = jnum(e, "ref.middle_square_cycles.value", 0.0_dp)
    tol_c = jnum(e, "ref.middle_square_cycles.abs_tol", 0.0_dp)
    ref_z = jnum(e, "ref.middle_square_zero_sink.value", 0.0_dp)
    zs_tol = jnum(e, "ref.middle_square_zero_sink.abs_tol", 0.0_dp)
    st = jnum(e, "ref.variant_spread.abs_tol", 1.0e-6_dp)
    ! --- middle-square cycles (exhaustivo, seen global tipo Python)
    allocate(seen(0:modulus-1), tr_arr(modulus))
    seen = 0
    n_cycles = 0
    do i = 0, modulus-1
      x = i
      tr = 0
      do
        if (seen(x) /= 0) exit
        seen(x) = 1
        tr = tr + 1
        tr_arr(tr) = x
        x = int(mid_ext(int(x,8)*int(x,8), digits), 4)
      end do
      hit = .false.
      do idx = 1, tr
        if (tr_arr(idx) == x) then
          hit = .true.
          exit
        end if
      end do
      if (hit) n_cycles = n_cycles + 1
    end do
    call add_check(r, "middle_square_cycles", real(n_cycles, dp), ref_c, tol_c, &
      abs(real(n_cycles, dp)-ref_c)/max(abs(ref_c), 1.0_dp), &
      abs(real(n_cycles, dp)-ref_c) <= tol_c)
    ! --- middle-square zero sink (exhaustivo)
    zero_sink = 0
    do i = 0, modulus-1
      x = i
      do j = 1, 50
        x = int(mid_ext(int(x,8)*int(x,8), digits), 4)
      end do
      if (x == 0) zero_sink = zero_sink + 1
    end do
    zs_rel = real(zero_sink, dp)/real(modulus, dp)
    zs_err = abs(zs_rel - ref_z)
    call add_check(r, "middle_square_zero_sink", zs_rel, ref_z, zs_tol, zs_err, zs_err <= zs_tol)
    ! --- weyl_period (null-ref condicional, como Python)
    x = 1234
    w = 0
    allocate(seenw(0:modulus-1))
    seenw = 0
    trw = 0
    weyl_period = 0
    do i = 1, max_iter
      w = modulo(w + weyl_s, modulus)
      q = int(x,8)*int(x,8) + int(w,8)
      x = int(mid_ext(q*q, digits), 4)
      if (seenw(x) /= 0) then
        do idx = 1, trw
          if (tr_arr(idx) == x) then
            weyl_period = i - idx + 1
            exit
          end if
        end do
        exit
      end if
      seenw(x) = 1
      trw = trw + 1
      tr_arr(trw) = x
    end do
    if (.not. ieee_is_nan(jnum(e, "ref.weyl_period.value", nanv()))) then
      wp_ref = jnum(e, "ref.weyl_period.value", 0.0_dp)
      wp_tol = jnum(e, "ref.weyl_period.abs_tol", 0.0_dp)
      call add_check(r, "weyl_period", real(weyl_period, dp), wp_ref, wp_tol, &
        abs(real(weyl_period, dp)-wp_ref)/max(abs(wp_ref), 1.0_dp), &
        abs(real(weyl_period, dp)-wp_ref) <= wp_tol)
    else
      call add_check(r, "weyl_period", real(weyl_period, dp), nanv(), nanv(), nanv(), .true.)
    end if
    ! --- tcr_period (null-ref condicional, como Python)
    x = 1234
    w1 = 0
    w2 = 0
    call ext_gcd(modulo(int(tcr_p1,8), int(tcr_p2,8)), int(tcr_p2,8), g, inv, yy)
    inv = modulo(inv, int(tcr_p2,8))
    m = tcr_p1*tcr_p2
    allocate(seen_crt(0:modulus-1))
    seen_crt = 0
    trcrt = 0
    tcr_period = 0
    do i = 1, max_iter
      w1 = modulo(w1 + weyl_s, tcr_p1)
      w2 = modulo(w2 + weyl_s, tcr_p2)
      k = int(modulo(int(w2 - w1,8)*inv, int(tcr_p2,8)))
      w = modulo(w1 + tcr_p1*k, m)
      q = int(x,8)*int(x,8) + int(w,8)
      x = int(mid_ext(q*q, digits), 4)
      if (seen_crt(x) /= 0) then
        do idx = 1, trcrt
          if (tr_arr(idx) == x) then
            tcr_period = i - idx + 1
            exit
          end if
        end do
        exit
      end if
      seen_crt(x) = 1
      trcrt = trcrt + 1
      tr_arr(trcrt) = x
    end do
    if (.not. ieee_is_nan(jnum(e, "ref.tcr_period.value", nanv()))) then
      tp_ref = jnum(e, "ref.tcr_period.value", 0.0_dp)
      tp_tol = jnum(e, "ref.tcr_period.abs_tol", 0.0_dp)
      call add_check(r, "tcr_period", real(tcr_period, dp), tp_ref, tp_tol, &
        abs(real(tcr_period, dp)-tp_ref)/max(abs(tp_ref), 1.0_dp), &
        abs(real(tcr_period, dp)-tp_ref) <= tp_tol)
    else
      call add_check(r, "tcr_period", real(tcr_period, dp), nanv(), nanv(), nanv(), .true.)
    end if
    ! --- variant_spread (null-ref si value es null, como Python)
    spread = abs(real(weyl_period, dp) - real(tcr_period, dp))/max(abs(real(tcr_period, dp)), 1.0_dp)
    if (.not. ieee_is_nan(jnum(e, "ref.variant_spread.value", nanv()))) then
      call add_check(r, "variant_spread", spread, &
        jnum(e, "ref.variant_spread.value", 0.0_dp), st, spread, spread <= st)
    else
      call add_check(r, "variant_spread", spread, nanv(), nanv(), nanv(), .true.)
    end if
    ok = .true.
    emsg = ''
  end subroutine

  ! ---------------------------------------------------------- hybrid_prng
  ! Python: registry_run.py eval_hybrid_prng (585-727). Sin outdir.
  ! ms_next_4 = medios 4 de x*x (ciclos/zero_sink); chi2 usa (x*x) %% 10000.
  ! r64 = u64/2**64 (Python: float(r)/float(2**64-1), ambos -> 2**64).
  ! chi2/chi2_bins/chi2_samples: null-ref condicional (value null -> ok).
  subroutine eval_hybrid_prng(e, r, ok, emsg)
    type(jval), intent(in) :: e
    type(run_t), intent(inout) :: r
    logical, intent(out) :: ok
    character(len=*), intent(out) :: emsg
    integer :: seed, n_samples, nbins
    integer :: i, j, x, b, idx
    integer(i8) :: mw, ms1, ms2, ms_comb, hi, lo
    integer :: n_cycles, zero_sink, tr
    integer, allocatable :: seen(:), hist(:), tr_arr(:)
    real(dp) :: zs_tol, zs_rel, zs_err, r64, expected, chi2
    real(dp) :: ref_c, tol_c, ref_z, ref_chi2, ref_bins, ref_samples
    real(dp) :: tol_chi2, tol_bins, tol_samples
    logical :: hit
    seed = jint(e, "params.seed", 1234)
    n_samples = jint(e, "params.n_samples", 10000)
    nbins = jint(e, "params.nbins", 10)
    ref_c = jnum(e, "ref.middle_square_cycles.value", 0.0_dp)
    tol_c = jnum(e, "ref.middle_square_cycles.abs_tol", 0.0_dp)
    ref_z = jnum(e, "ref.middle_square_zero_sink.value", 0.0_dp)
    zs_tol = jnum(e, "ref.middle_square_zero_sink.abs_tol", 0.0_dp)
    ! --- middle-square cycles (10000, seen global)
    allocate(seen(0:9999), tr_arr(10000))
    seen = 0
    n_cycles = 0
    do i = 0, 9999
      x = i
      tr = 0
      do
        if (seen(x) /= 0) exit
        seen(x) = 1
        tr = tr + 1
        tr_arr(tr) = x
        x = ms_next_4(x)
      end do
      hit = .false.
      do idx = 1, tr
        if (tr_arr(idx) == x) then
          hit = .true.
          exit
        end if
      end do
      if (hit) n_cycles = n_cycles + 1
    end do
    call add_check(r, "middle_square_cycles", real(n_cycles, dp), ref_c, tol_c, &
      abs(real(n_cycles, dp)-ref_c)/max(abs(ref_c), 1.0_dp), &
      abs(real(n_cycles, dp)-ref_c) <= tol_c)
    ! --- zero sink (10000, 50 pasos, 4 digitos)
    zero_sink = 0
    do i = 0, 9999
      x = i
      do j = 1, 50
        x = ms_next_4(x)
      end do
      if (x == 0) zero_sink = zero_sink + 1
    end do
    zs_rel = real(zero_sink, dp)/10000.0_dp
    zs_err = abs(zs_rel - ref_z)
    call add_check(r, "middle_square_zero_sink", zs_rel, ref_z, zs_tol, zs_err, zs_err <= zs_tol)
    ! --- chi2 (middle-square ultimos-4 + msws64 xor), orden Python
    allocate(hist(0:nbins-1))
    hist = 0
    ms1 = int(modulo(seed, 9999) + 1, i8)
    ms2 = int(modulo(ishft(seed, -16), 9999) + 1, i8)
    mw = int(ior(seed, 1), i8)
    do i = 1, n_samples
      ms1 = modulo(ms1*ms1, 10000_i8)
      ms2 = modulo(ms2*ms2, 10000_i8)
      call msws64(mw)
      ms_comb = ior(ishft(ms1, 48), ishft(ms2, 32))
      ms_comb = ieor(mw, ms_comb)
      hi = ishft(ms_comb, -32)
      lo = iand(ms_comb, 4294967295_i8)
      r64 = (real(hi, dp)*4294967296.0_dp + real(lo, dp)) / 18446744073709551616.0_dp
      b = int(r64*real(nbins, dp))
      if (b >= nbins) b = nbins - 1
      if (b < 0) b = 0
      hist(b) = hist(b) + 1
    end do
    expected = real(n_samples, dp)/real(nbins, dp)
    chi2 = 0.0_dp
    do b = 0, nbins-1
      chi2 = chi2 + (real(hist(b), dp) - expected)**2 / expected
    end do
    ! chi2 null-ref condicional (value null -> check pasada)
    if (.not. ieee_is_nan(jnum(e, "ref.chi2.value", nanv()))) then
      ref_chi2 = jnum(e, "ref.chi2.value", 0.0_dp)
      tol_chi2 = jnum(e, "ref.chi2.abs_tol", 1.0e-6_dp)
      call add_check(r, "chi2", chi2, ref_chi2, tol_chi2, &
        abs(chi2-ref_chi2)/max(abs(ref_chi2),1.0_dp), &
        abs(chi2-ref_chi2) <= tol_chi2)
    else
      call add_check(r, "chi2", chi2, nanv(), nanv(), nanv(), .true.)
    end if
    ! chi2_bins null-ref condicional
    if (.not. ieee_is_nan(jnum(e, "ref.chi2_bins.value", nanv()))) then
      ref_bins = jnum(e, "ref.chi2_bins.value", 0.0_dp)
      tol_bins = jnum(e, "ref.chi2_bins.abs_tol", 0.0_dp)
      call add_check(r, "chi2_bins", real(nbins, dp), ref_bins, tol_bins, &
        abs(real(nbins,dp)-ref_bins)/max(abs(ref_bins),1.0_dp), &
        abs(real(nbins,dp)-ref_bins) <= tol_bins)
    else
      call add_check(r, "chi2_bins", real(nbins, dp), nanv(), nanv(), nanv(), .true.)
    end if
    ! chi2_samples null-ref condicional
    if (.not. ieee_is_nan(jnum(e, "ref.chi2_samples.value", nanv()))) then
      ref_samples = jnum(e, "ref.chi2_samples.value", 0.0_dp)
      tol_samples = jnum(e, "ref.chi2_samples.abs_tol", 0.0_dp)
      call add_check(r, "chi2_samples", real(n_samples, dp), ref_samples, tol_samples, &
        abs(real(n_samples,dp)-ref_samples)/max(abs(ref_samples),1.0_dp), &
        abs(real(n_samples,dp)-ref_samples) <= tol_samples)
    else
      call add_check(r, "chi2_samples", real(n_samples, dp), nanv(), nanv(), nanv(), .true.)
    end if
    ok = .true.
    emsg = ''
  end subroutine

  ! ------------------------------------------------------- ms_next_4
  ! 4-digitos fijos, sin digits del JSON (hardcodeado Python 600).
  pure integer function ms_next_4(x) result(y)
    integer, intent(in) :: x
    character(len=8) :: s
    write(s, '(I0.8)') x*x
    read(s(3:6), *) y
  end function

end program
