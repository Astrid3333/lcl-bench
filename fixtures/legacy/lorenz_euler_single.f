c     Lorenz attractor, forward Euler, bare REAL = single precision.
c     Classic legacy defects: no REAL*8, first-order step.
      program lorenz_euler
      implicit none
      integer i, n
      parameter (n = 50000)
      real x, y, z, t, dt, sigma, rho, beta
      real dx, dy, dz
      sigma = 10.0
      rho = 28.0
      beta = 2.6666667
      x = 1.0
      y = 1.0
      z = 1.0
      dt = 0.01
      t = 0.0
      do i = 1, n
        dx = sigma * (y - x)
        dy = x * (rho - z) - y
        dz = x * y - beta * z
        x = x + dt * dx
        y = y + dt * dy
        z = z + dt * dz
        t = t + dt
      end do
      write(*,*) x, y, z
      end
