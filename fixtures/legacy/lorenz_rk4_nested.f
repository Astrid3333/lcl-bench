c     Lorenz attractor, classic Runge-Kutta 4 with a nested do-loop:
c     the outer do steps time, the inner do combines the four stages
c     over the state vector with an array index taken from the inner
c     loop variable. Exercises do-inside-do, array reads/writes indexed
c     by a loop variable, and a double end do; no IF and no GOTO.
      program lorenz_rk4_nested
      implicit none
      integer i, j, n
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
        do j = 1, 3
          y(j) = y(j) + (k1(j) + 2.0d0*k2(j) + 2.0d0*k3(j)
     &            + k4(j))/6.0d0
        end do
      end do
      write(*,*) y(1), y(2), y(3)
      end
