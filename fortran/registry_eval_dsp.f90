! registry_eval_dsp.f90 - Evaluador de Goertzel para dsp_bin (Fortran)
!
! Calcula magnitudes de Goertzel (K, alpha, beta) para tonos dados
! e imprime los resultados como números separados por espacios.
!
! Uso:
!   ./registry_eval_dsp <N> <num_tones> <k1> <A1> <phi1> ...
!
! Salida: para cada tono, imprime: mag_K mag_alpha mag_beta refv
!         luego imprime: spread

program registry_eval_dsp
    implicit none
    integer, parameter :: dp = selected_real_kind(15, 30)
    integer :: N, num_tones, i, t
    real(dp), allocatable :: k(:), A(:), phi(:), x(:)
    real(dp) :: w, c, s1, s2, s0, s, d, ee, b, sN, sP
    real(dp) :: mag_K, mag_alpha, mag_beta, refv
    real(dp) :: max_spread, spread
    character(len=32) :: arg_str

    ! Leer argumentos de línea de comandos (como strings, luego convertir)
    if (command_argument_count() < 3) then
        print *, "Uso: registry_eval_dsp <N> <num_tones> <k1> <A1> <phi1> ..."
        stop 1
    end if

    call get_command_argument(1, arg_str)
    read(arg_str, *) N
    call get_command_argument(2, arg_str)
    read(arg_str, *) num_tones

    allocate(k(num_tones), A(num_tones), phi(num_tones), x(N))

    do i = 1, num_tones
        call get_command_argument(2 + 3*(i-1) + 1, arg_str)
        read(arg_str, *) k(i)
        call get_command_argument(2 + 3*(i-1) + 2, arg_str)
        read(arg_str, *) A(i)
        call get_command_argument(2 + 3*(i-1) + 3, arg_str)
        read(arg_str, *) phi(i)
    end do

    ! Generar señal x[n] = sum A_i*cos(2*pi*k_i*n/N + phi_i)
    x = 0.0_dp
    do t = 1, num_tones
        do i = 1, N
            x(i) = x(i) + A(t) * cos(2.0_dp * 3.141592653589793238462643383279502884_dp * k(t) * (i-1) / N + phi(t))
        end do
    end do

    ! Para cada tono, calcular las 3 formulaciones de Goertzel
    max_spread = 0.0_dp
    do t = 1, num_tones
        w = 2.0_dp * 3.141592653589793238462643383279502884_dp * k(t) / N
        c = 2.0_dp * cos(w)

        ! Goertzel K
        s1 = 0.0_dp; s2 = 0.0_dp
        do i = 1, N
            s0 = x(i) + c * s1 - s2
            s2 = s1
            s1 = s0
        end do
        mag_K = sqrt(s1**2 + s2**2 - c * s1 * s2)

        ! Goertzel alpha
        s = 0.0_dp; d = 0.0_dp
        do i = 1, N
            d = d + x(i) + (c - 2.0_dp) * s
            s = s + d
        end do
        sN = s; sP = s - d
        mag_alpha = sqrt(sN**2 + sP**2 - c * sN * sP)

        ! Goertzel beta
        s = 0.0_dp; ee = 0.0_dp; b = c + 2.0_dp
        do i = 1, N
            ee = x(i) + b * s - ee
            s = ee - s
        end do
        sN = s; sP = ee - s
        mag_beta = sqrt(sN**2 + sP**2 - c * sN * sP)

        ! Referencia
        refv = N * A(t) / 2.0_dp

        ! Imprimir: mag_K mag_alpha mag_beta refv
        print '(4(f0.15, 1x))', mag_K, mag_alpha, mag_beta, refv

        ! Calcular spread para este tono
        spread = max(abs(mag_alpha - mag_K) / refv, abs(mag_beta - mag_K) / refv)
        max_spread = max(max_spread, spread)
    end do

    ! Imprimir spread máximo
    print '(f0.15)', max_spread
end program registry_eval_dsp