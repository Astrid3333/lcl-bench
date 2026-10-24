c     Lorenz attractor, forward Euler, single precision, dt far too large.
      program lorenz_euler_huge
      implicit none
      integer i, n
      parameter (n = 20000)
      real x, y, z, dt, sigma, rho, beta
      real dx, dy, dz
      sigma = 10.0
      rho = 28.0
      beta = 2.6666667
      x = 1.0
      y = 1.0
      z = 1.0
      dt = 0.1
      do i = 1, n
        dx = sigma * (y - x)
        dy = x * (rho - z) - y
        dz = x * y - beta * z
        x = x + dt * dx
        y = y + dt * dy
        z = z + dt * dz
      end do
      write(*,*) x, y, z
      end
