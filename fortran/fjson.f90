module fjson
  use kinds, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite, &
       ieee_quiet_nan, ieee_value
  implicit none

  integer, parameter :: JV_NULL = 0
  integer, parameter :: JV_BOOL = 1
  integer, parameter :: JV_NUM = 2
  integer, parameter :: JV_STR = 3
  integer, parameter :: JV_ARR = 4
  integer, parameter :: JV_OBJ = 5

  type :: jval
    integer :: t = JV_NULL
    logical :: b = .false.
    real(dp) :: v = 0.0_dp
    character(len=:), allocatable :: str
    character(len=64), allocatable :: keys(:)
    type(jval), allocatable :: kids(:)
    integer :: n = 0
  end type jval

  type :: pstate
    character(len=:), allocatable :: s
    integer :: i = 1
    logical :: err = .false.
    character(len=160) :: msg = ''
  end type pstate

  private :: pstate, set_err, skip_ws, is_ws, expect_lit, parse_value, &
       parse_number, parse_string, is_numchar, hex4, utf8, jfind, grow, &
       app_val, wnum, wstr, hex2

contains

  subroutine jload(path, root, ok, msg)
    character(len=*), intent(in) :: path
    type(jval), intent(out) :: root
    logical, intent(out) :: ok
    character(len=*), intent(out) :: msg
    integer :: u, ios, sz
    character(len=:), allocatable :: buf
    type(pstate) :: st

    ok = .false.
    msg = ''
    root%t = JV_NULL
    open(newunit=u, file=path, status='old', action='read', &
         access='stream', form='unformatted', iostat=ios)
    if (ios /= 0) then
      msg = 'cannot open ' // trim(path)
      return
    end if
    inquire(u, size=sz)
    if (sz <= 0) then
      close(u)
      msg = 'empty file ' // trim(path)
      return
    end if
    allocate(character(len=sz) :: buf)
    read(u, iostat=ios) buf
    close(u)
    if (ios /= 0) then
      msg = 'cannot read ' // trim(path)
      return
    end if
    st = pstate(s=buf, i=1, err=.false., msg='')
    call skip_ws(st)
    call parse_value(st, root)
    if (.not. st%err) then
      call skip_ws(st)
      if (st%i <= len(st%s)) call set_err(st, 'trailing junk')
    end if
    if (st%err) then
      msg = st%msg
      root%t = JV_NULL
      return
    end if
    ok = .true.
  end subroutine jload

  ! Walk a dotted path "a.b.c" level-by-level.
  ! Intermediate miss or non-object → JV_NULL.
  ! Single-segment keys behave exactly as before.
  function jget(o, key) result(v)
    type(jval), intent(in) :: o
    character(len=*), intent(in) :: key
    type(jval) :: v, cur, nxt
    integer :: k, p, n
    character(len=len(key)) :: seg, rest

    call jcopy(cur, o)
    rest = key
    do
      n = len_trim(rest)
      if (n <= 0) then
        v%t = JV_NULL
        return
      end if
      p = index(rest, '.')
      if (p <= 0) exit          ! last segment
      if (p == 1) then          ! empty segment
        v%t = JV_NULL
        return
      end if
      if (cur%t /= JV_OBJ) then
        v%t = JV_NULL
        return
      end if
      seg = rest(1:p-1)
      k = jfind(cur, seg(1:len_trim(seg)))
      if (k <= 0) then
        v%t = JV_NULL
        return
      end if
      call jcopy(nxt, cur%kids(k))  ! temp: avoid aliasing cur with cur%kids
      call jcopy(cur, nxt)
      rest = rest(p+1:)
    end do
    if (cur%t /= JV_OBJ) then
      v%t = JV_NULL
      return
    end if
    k = jfind(cur, rest(1:n))
    if (k > 0) then
      call jcopy(v, cur%kids(k))
    else
      v%t = JV_NULL
    end if
  end function jget

  ! True if the final key exists on the dotted path (value may be null).
  logical function jhas(o, key) result(r)
    type(jval), intent(in) :: o
    character(len=*), intent(in) :: key
    type(jval) :: cur, nxt
    integer :: k, p, n
    character(len=len(key)) :: seg, rest
    call jcopy(cur, o)
    rest = key
    do
      n = len_trim(rest)
      if (n <= 0) then
        r = .false.
        return
      end if
      p = index(rest, '.')
      if (p <= 0) exit
      if (p == 1) then
        r = .false.
        return
      end if
      if (cur%t /= JV_OBJ) then
        r = .false.
        return
      end if
      seg = rest(1:p-1)
      k = jfind(cur, seg(1:len_trim(seg)))
      if (k <= 0) then
        r = .false.
        return
      end if
      call jcopy(nxt, cur%kids(k))  ! temp: avoid aliasing
      call jcopy(cur, nxt)
      rest = rest(p+1:)
    end do
    if (cur%t /= JV_OBJ) then
      r = .false.
      return
    end if
    r = jfind(cur, rest(1:n)) > 0
  end function jhas

  real(dp) function jnum(o, key, dflt) result(r)
    type(jval), intent(in) :: o
    character(len=*), intent(in) :: key
    real(dp), intent(in) :: dflt
    type(jval) :: v
    call jcopy(v, jget(o, key))
    if (v%t == JV_NUM) then
      r = v%v
    else
      r = dflt
    end if
  end function jnum

  real(dp) function jnum_req(o, key) result(r)
    type(jval), intent(in) :: o
    character(len=*), intent(in) :: key
    type(jval) :: v
    real(dp) :: q
    q = ieee_value(0.0_dp, ieee_quiet_nan)
    call jcopy(v, jget(o, key))
    if (v%t == JV_NUM) then
      r = v%v
    else
      r = q
    end if
  end function jnum_req

  integer function jint(o, key, dflt) result(r)
    type(jval), intent(in) :: o
    character(len=*), intent(in) :: key
    integer, intent(in) :: dflt
    real(dp) :: x
    x = jnum(o, key, real(dflt, dp))
    if (ieee_is_nan(x)) then
      r = dflt
    else
      r = nint(x)
    end if
  end function jint

  function jstr(o, key, dflt) result(s)
    type(jval), intent(in) :: o
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: dflt
    character(len=:), allocatable :: s
    type(jval) :: v
    call jcopy(v, jget(o, key))
    if (v%t == JV_STR) then
      s = v%str
    else
      s = dflt
    end if
  end function jstr

  function jtovec(o, key) result(x)
    type(jval), intent(in) :: o
    character(len=*), intent(in) :: key
    real(dp), allocatable :: x(:)
    type(jval) :: v
    integer :: k
    call jcopy(v, jget(o, key))
    if (v%t /= JV_ARR) then
      allocate(x(0))
      return
    end if
    allocate(x(v%n))
    do k = 1, v%n
      if (v%kids(k)%t == JV_NUM) then
        x(k) = v%kids(k)%v
      else
        x(k) = ieee_value(0.0_dp, ieee_quiet_nan)
      end if
    end do
  end function jtovec

  function jat(o, i) result(v)
    type(jval), intent(in) :: o
    integer, intent(in) :: i
    type(jval) :: v
    if (o%t == JV_ARR .and. i >= 1 .and. i <= o%n) then
      call jcopy(v, o%kids(i))
    else
      v%t = JV_NULL
    end if
  end function jat

  integer function jalen(o) result(r)
    type(jval), intent(in) :: o
    if (o%t == JV_ARR) then
      r = o%n
    else
      r = 0
    end if
  end function jalen

  type(jval) function jmk_null() result(r)
    r%t = JV_NULL
  end function jmk_null

  type(jval) function jmk_bool(b) result(r)
    logical, intent(in) :: b
    r%t = JV_BOOL
    r%b = b
  end function jmk_bool

  type(jval) function jmk_num(x) result(r)
    real(dp), intent(in) :: x
    r%t = JV_NUM
    r%v = x
  end function jmk_num

  type(jval) function jmk_str(s) result(r)
    character(len=*), intent(in) :: s
    r%t = JV_STR
    r%str = s
  end function jmk_str

  type(jval) function jmk_arr() result(r)
    r%t = JV_ARR
    r%n = 0
  end function jmk_arr

  type(jval) function jmk_obj() result(r)
    r%t = JV_OBJ
    r%n = 0
  end function jmk_obj

  subroutine jpush(parent, child)
    type(jval), intent(inout) :: parent
    type(jval), intent(in) :: child
    call grow(parent)
    parent%n = parent%n + 1
    call jcopy(parent%kids(parent%n), child)
  end subroutine jpush

  subroutine jset(parent, key, val)
    type(jval), intent(inout) :: parent
    character(len=*), intent(in) :: key
    type(jval), intent(in) :: val
    integer :: k
    k = jfind(parent, key)
    if (k > 0) then
      call jcopy(parent%kids(k), val)
      return
    end if
    call grow(parent)
    parent%n = parent%n + 1
    parent%keys(parent%n) = key
    call jcopy(parent%kids(parent%n), val)
  end subroutine jset

  recursive subroutine jcopy(d, s)
    type(jval), intent(out) :: d
    type(jval), intent(in) :: s
    integer :: i
    d%t = s%t
    d%b = s%b
    d%v = s%v
    d%n = s%n
    if (allocated(s%str)) d%str = s%str
    if (allocated(s%keys)) then
      allocate(d%keys(size(s%keys)))
      d%keys = s%keys
    end if
    if (allocated(s%kids)) then
      allocate(d%kids(size(s%kids)))
      do i = 1, size(s%kids)
        call jcopy(d%kids(i), s%kids(i))
      end do
    end if
  end subroutine jcopy

  subroutine jdump(u, root)
    integer, intent(in) :: u
    type(jval), intent(in) :: root
    character(len=:), allocatable :: buf
    buf = ''
    call app_val(buf, root, 0)
    write(u, '(a)', advance='no') buf
  end subroutine jdump

  subroutine set_err(st, m)
    type(pstate), intent(inout) :: st
    character(len=*), intent(in) :: m
    if (st%err) return
    st%err = .true.
    write(st%msg, '(a,i0)') trim(m) // ' at ', st%i
  end subroutine set_err

  logical function is_ws(c) result(r)
    character(len=1), intent(in) :: c
    r = c == ' ' .or. c == achar(9) .or. c == achar(10) .or. c == achar(13)
  end function is_ws

  subroutine skip_ws(st)
    type(pstate), intent(inout) :: st
    do while (st%i <= len(st%s))
      if (.not. is_ws(st%s(st%i:st%i))) exit
      st%i = st%i + 1
    end do
  end subroutine skip_ws

  subroutine expect_lit(st, lit)
    type(pstate), intent(inout) :: st
    character(len=*), intent(in) :: lit
    integer :: L
    L = len(lit)
    if (st%i + L - 1 > len(st%s)) then
      call set_err(st, 'unexpected end')
      return
    end if
    if (st%s(st%i:st%i + L - 1) /= lit) then
      call set_err(st, 'bad literal')
      return
    end if
    st%i = st%i + L
  end subroutine expect_lit

  recursive subroutine parse_value(st, v)
    type(pstate), intent(inout) :: st
    type(jval), intent(out) :: v
    character(len=1) :: c
    character(len=:), allocatable :: key
    type(jval) :: kid

    v%t = JV_NULL
    v%n = 0
    call skip_ws(st)
    if (st%err) return
    if (st%i > len(st%s)) then
      call set_err(st, 'unexpected end')
      return
    end if
    c = st%s(st%i:st%i)
    if (c == 'n') then
      call expect_lit(st, 'null')
      v%t = JV_NULL
    else if (c == 't') then
      call expect_lit(st, 'true')
      v%t = JV_BOOL
      v%b = .true.
    else if (c == 'f') then
      call expect_lit(st, 'false')
      v%t = JV_BOOL
      v%b = .false.
    else if (c == '"') then
      st%i = st%i + 1
      call parse_string(st, v%str)
      if (st%err) return
      v%t = JV_STR
    else if (c == '[') then
      st%i = st%i + 1
      v%t = JV_ARR
      v%n = 0
      call skip_ws(st)
      if (st%err) return
      if (st%i <= len(st%s)) then
        if (st%s(st%i:st%i) == ']') then
          st%i = st%i + 1
          return
        end if
      end if
      do
        call parse_value(st, kid)
        if (st%err) return
        call jpush(v, kid)
        call skip_ws(st)
        if (st%err) return
        if (st%i > len(st%s)) then
          call set_err(st, 'unterminated array')
          return
        end if
        c = st%s(st%i:st%i)
        if (c == ',') then
          st%i = st%i + 1
        else if (c == ']') then
          st%i = st%i + 1
          exit
        else
          call set_err(st, 'expected , or ]')
          return
        end if
      end do
    else if (c == '{') then
      st%i = st%i + 1
      v%t = JV_OBJ
      v%n = 0
      call skip_ws(st)
      if (st%err) return
      if (st%i <= len(st%s)) then
        if (st%s(st%i:st%i) == '}') then
          st%i = st%i + 1
          return
        end if
      end if
      do
        call skip_ws(st)
        if (st%err) return
        if (st%i > len(st%s)) then
          call set_err(st, 'unterminated object')
          return
        end if
        if (st%s(st%i:st%i) /= '"') then
          call set_err(st, 'expected string key')
          return
        end if
        st%i = st%i + 1
        call parse_string(st, key)
        if (st%err) return
        call skip_ws(st)
        if (st%err) return
        if (st%i > len(st%s)) then
          call set_err(st, 'unterminated object')
          return
        end if
        if (st%s(st%i:st%i) /= ':') then
          call set_err(st, 'expected :')
          return
        end if
        st%i = st%i + 1
        call parse_value(st, kid)
        if (st%err) return
        call jset(v, key, kid)
        call skip_ws(st)
        if (st%err) return
        if (st%i > len(st%s)) then
          call set_err(st, 'unterminated object')
          return
        end if
        c = st%s(st%i:st%i)
        if (c == ',') then
          st%i = st%i + 1
        else if (c == '}') then
          st%i = st%i + 1
          exit
        else
          call set_err(st, 'expected , or }')
          return
        end if
      end do
    else if (c == '-' .or. (c >= '0' .and. c <= '9')) then
      call parse_number(st, v%v)
      if (st%err) return
      v%t = JV_NUM
    else
      call set_err(st, 'unexpected character')
    end if
  end subroutine parse_value

  logical function is_numchar(c) result(r)
    character(len=1), intent(in) :: c
    r = (c >= '0' .and. c <= '9') .or. c == '-' .or. c == '+' .or. &
         c == '.' .or. c == 'e' .or. c == 'E'
  end function is_numchar

  subroutine parse_number(st, val)
    type(pstate), intent(inout) :: st
    real(dp), intent(out) :: val
    integer :: start, ios
    character(len=:), allocatable :: tok
    start = st%i
    do while (st%i <= len(st%s))
      if (.not. is_numchar(st%s(st%i:st%i))) exit
      st%i = st%i + 1
    end do
    tok = st%s(start:st%i - 1)
    read(tok, *, iostat=ios) val
    if (ios /= 0) then
      call set_err(st, 'bad number')
      val = 0.0_dp
    end if
  end subroutine parse_number

  subroutine parse_string(st, out)
    type(pstate), intent(inout) :: st
    character(len=:), allocatable, intent(out) :: out
    character(len=:), allocatable :: buf
    character(len=1) :: c, c2
    integer :: cp, cp2
    logical :: ok
    buf = ''
    do
      if (st%i > len(st%s)) then
        call set_err(st, 'unterminated string')
        return
      end if
      c = st%s(st%i:st%i)
      st%i = st%i + 1
      if (c == '"') exit
      if (c == '\') then
        if (st%i > len(st%s)) then
          call set_err(st, 'unterminated escape')
          return
        end if
        c2 = st%s(st%i:st%i)
        st%i = st%i + 1
        if (c2 == '"') then
          buf = buf // '"'
        else if (c2 == '\') then
          buf = buf // '\'
        else if (c2 == '/') then
          buf = buf // '/'
        else if (c2 == 'b') then
          buf = buf // achar(8)
        else if (c2 == 'f') then
          buf = buf // achar(12)
        else if (c2 == 'n') then
          buf = buf // achar(10)
        else if (c2 == 'r') then
          buf = buf // achar(13)
        else if (c2 == 't') then
          buf = buf // achar(9)
        else if (c2 == 'u') then
          call hex4(st, cp, ok)
          if (.not. ok) return
          if (cp >= 55296 .and. cp <= 56319) then
            if (st%i + 5 <= len(st%s)) then
              if (st%s(st%i:st%i + 1) == '\u') then
                st%i = st%i + 2
                call hex4(st, cp2, ok)
                if (.not. ok) return
                if (cp2 >= 56320 .and. cp2 <= 57343) then
                  cp = 65536 + (cp - 55296) * 1024 + (cp2 - 56320)
                else
                  call set_err(st, 'bad surrogate pair')
                  return
                end if
              else
                call set_err(st, 'lone surrogate')
                return
              end if
            else
              call set_err(st, 'lone surrogate')
              return
            end if
          else if (cp >= 56320 .and. cp <= 57343) then
            call set_err(st, 'lone low surrogate')
            return
          end if
          buf = buf // utf8(cp)
        else
          call set_err(st, 'bad escape')
          return
        end if
      else
        buf = buf // c
      end if
    end do
    out = buf
  end subroutine parse_string

  subroutine hex4(st, cp, ok)
    type(pstate), intent(inout) :: st
    integer, intent(out) :: cp
    logical, intent(out) :: ok
    integer :: k, d, v
    ok = .false.
    cp = 0
    if (st%i + 3 > len(st%s)) then
      call set_err(st, 'truncated \\u')
      return
    end if
    do k = 0, 3
      d = iachar(st%s(st%i + k:st%i + k))
      if (d >= 48 .and. d <= 57) then
        v = d - 48
      else if (d >= 97 .and. d <= 102) then
        v = d - 87
      else if (d >= 65 .and. d <= 70) then
        v = d - 55
      else
        call set_err(st, 'bad hex digit')
        return
      end if
      cp = cp * 16 + v
    end do
    st%i = st%i + 4
    ok = .true.
  end subroutine hex4

  function utf8(cp) result(s)
    integer, intent(in) :: cp
    character(len=:), allocatable :: s
    integer :: b1, b2, b3, b4
    if (cp < 128) then
      s = char(cp)
    else if (cp < 2048) then
      b1 = 192 + cp / 64
      b2 = 128 + mod(cp, 64)
      s = char(b1) // char(b2)
    else if (cp < 65536) then
      b1 = 224 + cp / 4096
      b2 = 128 + mod(cp / 64, 64)
      b3 = 128 + mod(cp, 64)
      s = char(b1) // char(b2) // char(b3)
    else
      b1 = 240 + cp / 262144
      b2 = 128 + mod(cp / 4096, 64)
      b3 = 128 + mod(cp / 64, 64)
      b4 = 128 + mod(cp, 64)
      s = char(b1) // char(b2) // char(b3) // char(b4)
    end if
  end function utf8

  integer function jfind(o, key) result(k)
    type(jval), intent(in) :: o
    character(len=*), intent(in) :: key
    k = 0
    if (o%t /= JV_OBJ) return
    do k = 1, o%n
      if (o%keys(k) == key) return
    end do
    k = 0
  end function jfind

  subroutine grow(parent)
    type(jval), intent(inout) :: parent
    type(jval), allocatable :: tmp(:)
    character(len=64), allocatable :: ktmp(:)
    integer :: cap, i
    if (.not. allocated(parent%kids)) then
      allocate(parent%kids(8))
      if (parent%t == JV_OBJ) allocate(parent%keys(8))
      return
    end if
    if (parent%n < size(parent%kids)) return
    cap = size(parent%kids)
    allocate(tmp(2 * cap))
    do i = 1, parent%n
      call jcopy(tmp(i), parent%kids(i))
    end do
    call move_alloc(tmp, parent%kids)
    if (allocated(parent%keys)) then
      allocate(ktmp(2 * cap))
      ktmp(1:parent%n) = parent%keys(1:parent%n)
      call move_alloc(ktmp, parent%keys)
    end if
  end subroutine grow

  recursive subroutine app_val(buf, v, ind)
    character(len=:), allocatable, intent(inout) :: buf
    type(jval), intent(in) :: v
    integer, intent(in) :: ind
    integer :: k

    select case (v%t)
    case (JV_NULL)
      buf = buf // 'null'
    case (JV_BOOL)
      if (v%b) then
        buf = buf // 'true'
      else
        buf = buf // 'false'
      end if
    case (JV_NUM)
      buf = buf // wnum(v%v)
    case (JV_STR)
      if (allocated(v%str)) then
        buf = buf // wstr(v%str)
      else
        buf = buf // '""'
      end if
    case (JV_ARR)
      if (v%n == 0) then
        buf = buf // '[]'
        return
      end if
      buf = buf // '['
      do k = 1, v%n
        buf = buf // achar(10) // repeat(' ', 2 * (ind + 1))
        call app_val(buf, v%kids(k), ind + 1)
        if (k < v%n) buf = buf // ','
      end do
      buf = buf // achar(10) // repeat(' ', 2 * ind) // ']'
    case (JV_OBJ)
      if (v%n == 0) then
        buf = buf // '{}'
        return
      end if
      buf = buf // '{'
      do k = 1, v%n
        buf = buf // achar(10) // repeat(' ', 2 * (ind + 1)) // &
             wstr(trim(v%keys(k))) // ': '
        call app_val(buf, v%kids(k), ind + 1)
        if (k < v%n) buf = buf // ','
      end do
      buf = buf // achar(10) // repeat(' ', 2 * ind) // '}'
    case default
      buf = buf // 'null'
    end select
  end subroutine app_val

  function wnum(x) result(s)
    real(dp), intent(in) :: x
    character(len=:), allocatable :: s
    character(len=40) :: c
    if (.not. ieee_is_finite(x)) then
      s = 'null'
      return
    end if
    if (ieee_is_nan(x)) then
      s = 'null'
      return
    end if
    if (abs(x) < 1.0e9_dp .and. x == real(int(x), dp)) then
      write(c, '(i0)') int(x)
      s = trim(c)
    else
      write(c, '(es26.17e3)') x
      s = trim(adjustl(c))
    end if
  end function wnum

  function wstr(s) result(o)
    character(len=*), intent(in) :: s
    character(len=:), allocatable :: o
    integer :: i, a
    character(len=1) :: c
    o = '"'
    do i = 1, len(s)
      c = s(i:i)
      a = iachar(c)
      if (c == '"') then
        o = o // '\"'
      else if (c == '\') then
        o = o // '\\'
      else if (a == 8) then
        o = o // '\b'
      else if (a == 12) then
        o = o // '\f'
      else if (a == 10) then
        o = o // '\n'
      else if (a == 13) then
        o = o // '\r'
      else if (a == 9) then
        o = o // '\t'
      else if (a < 32) then
        o = o // '\u00' // hex2(a)
      else
        o = o // c
      end if
    end do
    o = o // '"'
  end function wstr

  function hex2(k) result(s)
    integer, intent(in) :: k
    character(len=2) :: s
    character(len=16), parameter :: HC = '0123456789ABCDEF'
    integer :: hi, lo
    hi = k / 16
    lo = mod(k, 16)
    s = HC(hi + 1:hi + 1) // HC(lo + 1:lo + 1)
  end function hex2

end module fjson
