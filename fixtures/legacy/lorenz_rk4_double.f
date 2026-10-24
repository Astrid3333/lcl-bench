c     Lorenz attractor, classic Runge-Kutta 4, double precision.
c     Correct reference-era code: REAL*8 + RK4 stages.
      program lorenz_rk4
      implicit none
      integer i, n
      parameter (n = 100000)
      double precision x, y, z, t, dt, sigma, rho, beta
      double precision k1x, k1y, k1z, k2x, k2y, k2z
      double precision k3x, k3y, k3z, k4x, k4y, k4z
      sigma = 10.0d0
      rho = 28.0d0
      beta = 2.666666666666667d0
      x = 1.0d0
      y = 1.0d0
      z = 1.0d0
      dt = 0.01d0
      t = 0.0d0
      do i = 1, n
        k1x = dt * sigma * (y - x)
        k1y = dt * (x * (rho - z) - y)
        k1z = dt * (x * y - beta * z)
        k2x = dt * sigma * ((y + 0.5d0*k1y) - (x + 0.5d0*k1x))
        k2y = dt * ((x + 0.5d0*k1x) * (rho - (z + 0.5d0*k1z))
     &            - (y + 0.5d0*k1y))
        k2z = dt * ((x + 0.5d0*k1x) * (y + 0.5d0*k1y)
     &            - beta * (z + 0.5d0*k1z))
        k3x = dt * sigma * ((y + 0.5d0*k2y) - (x + 0.5d0*k2x))
        k3y = dt * ((x + 0.5d0*k2x) * (rho - (z + 0.5d0*k2z))
     &            - (y + 0.5d0*k2y))
        k3z = dt * ((x + 0.5d0*k2x) * (y + 0.5d0*k2y)
     &            - beta * (z + 0.5d0*k2z))
        k4x = dt * sigma * ((y + k3y) - (x + k3x))
        k4y = dt * ((x + k3x) * (rho - (z + k3z)) - (y + k3y))
        k4z = dt * ((x + k3x) * (y + k3y) - beta * (z + k3z))
        x = x + (k1x + 2.0d0*k2x + 2.0d0*k3x + k4x) / 6.0d0
        y = y + (k1y + 2.0d0*k2y + 2.0d0*k3y + k4y) / 6.0d0
        z = z + (k1z + 2.0d0*k2z + 2.0d0*k3z + k4z) / 6.0d0
        t = t + dt
      end do
      write(*,*) x, y, z
      end
