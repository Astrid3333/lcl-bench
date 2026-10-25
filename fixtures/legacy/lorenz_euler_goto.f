c     Lorenz attractor, forward Euler, double precision, with a classic
c     divergence guard written as a block IF + GOTO to the trailing
c     continue (label 200). The guard never trips at this dt.
      program lorenz_euler_goto
      implicit none
      integer i, n
      parameter (n = 100000)
      double precision x, y, z, t, dt, sigma, rho, beta
      double precision dx, dy, dz
      sigma = 10.0d0
      rho = 28.0d0
      beta = 2.666666666666667d0
      x = 1.0d0
      y = 1.0d0
      z = 1.0d0
      dt = 0.01d0
      t = 0.0d0
      do i = 1, n
        dx = sigma * (y - x)
        dy = x * (rho - z) - y
        dz = x * y - beta * z
        x = x + dt * dx
        y = y + dt * dy
        z = z + dt * dz
        t = t + dt
        if (abs(x) .gt. 1.0d30 .or. abs(y) .gt. 1.0d30) then
          goto 200
        end if
      end do
  200 continue
      write(*,*) x, y, z
      end
