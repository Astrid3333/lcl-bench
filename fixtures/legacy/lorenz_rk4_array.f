c     Lorenz attractor, classic Runge-Kutta 4, state held in a 3-element
c     array. Exercises array reads/writes and array-element assignments;
c     no IF and no GOTO (structured translation path).
      program lorenz_rk4_array
      implicit none
      integer i, n
      parameter (n = 100000)
      double precision y(3), k1(3), k2(3), k3(3), k4(3)
      double precision dt, sigma, rho, beta
      sigma = 10.0d0
      rho = 28.0d0
      beta = 2.666666666666667d0
      y(1) = 1.0d0
      y(2) = 1.0d0
      y(3) = 1.0d0
      dt = 0.01d0
      do i = 1, n
        k1(1) = dt * sigma * (y(2) - y(1))
        k1(2) = dt * (y(1) * (rho - y(3)) - y(2))
        k1(3) = dt * (y(1) * y(2) - beta * y(3))
        k2(1) = dt * sigma * ((y(2) + 0.5d0*k1(2)) - (y(1) + 0.5d0*k1(1)))
        k2(2) = dt * ((y(1) + 0.5d0*k1(1)) * (rho - (y(3) + 0.5d0*k1(3)))
     &            - (y(2) + 0.5d0*k1(2)))
        k2(3) = dt * ((y(1) + 0.5d0*k1(1)) * (y(2) + 0.5d0*k1(2))
     &            - beta * (y(3) + 0.5d0*k1(3)))
        k3(1) = dt * sigma * ((y(2) + 0.5d0*k2(2)) - (y(1) + 0.5d0*k2(1)))
        k3(2) = dt * ((y(1) + 0.5d0*k2(1)) * (rho - (y(3) + 0.5d0*k2(3)))
     &            - (y(2) + 0.5d0*k2(2)))
        k3(3) = dt * ((y(1) + 0.5d0*k2(1)) * (y(2) + 0.5d0*k2(2))
     &            - beta * (y(3) + 0.5d0*k2(3)))
        k4(1) = dt * sigma * ((y(2) + k3(2)) - (y(1) + k3(1)))
        k4(2) = dt * ((y(1) + k3(1)) * (rho - (y(3) + k3(3)))
     &            - (y(2) + k3(2)))
        k4(3) = dt * ((y(1) + k3(1)) * (y(2) + k3(2))
     &            - beta * (y(3) + k3(3)))
        y(1) = y(1) + (k1(1) + 2.0d0*k2(1) + 2.0d0*k3(1) + k4(1))/6.0d0
        y(2) = y(2) + (k1(2) + 2.0d0*k2(2) + 2.0d0*k3(2) + k4(2))/6.0d0
        y(3) = y(3) + (k1(3) + 2.0d0*k2(3) + 2.0d0*k3(3) + k4(3))/6.0d0
      end do
      write(*,*) y(1), y(2), y(3)
      end
