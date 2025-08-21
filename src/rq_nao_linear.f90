program residuos_quadraticos_nao_lineares
    !use lapack, only: dgesv
    implicit NONE(external)

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            integer, intent(in) :: n, nrhs, lda, ldb
            double precision, intent(in out) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgesv
    end interface

    integer, parameter :: dp = 8
    integer :: io, io_status, i, j, num_colunas, num_linhas, nrhs, lda, ldb, info, n
    integer :: max_iteracoes
    integer, allocatable :: ipiv(:)
    character(len=512) :: io_msg, filename, arg, fmt_string
    logical :: file_exists, done
    real :: work(50)
    real(dp), allocatable :: t(:), y(:), x(:), novo_x(:), d(:)
    real(dp), allocatable :: A(:, :), AT(:, :), ATA(:, :), b(:), ATb(:)
    real(dp) :: alpha, beta, delta, gama
    real(dp) :: norma_atual, norma_anterior
    real(dp), parameter :: tolerancia = 1.0d-8

    num_linhas = 100
    num_colunas = 4
    alpha = 1.0d0
    beta = 1.0d0
    delta = 1.0d0
    gama = 1.0d0

    if (command_argument_count() < 2) then
        filename = "input_nao_linear.txt"
        arg = "4"
    else
        call get_command_argument(1, filename)
        call get_command_argument(2, arg)
    end if
    read (arg, *, iostat=io_status) n

    inquire (file=trim(filename), exist=file_exists)
    if (.not. file_exists) then
        print *, 'Arquivo não encontrado:', trim(filename)
        stop
    end if

    allocate (t(num_linhas), y(num_linhas))
    allocate (x(num_colunas), novo_x(num_colunas), d(num_colunas))
    allocate (b(num_linhas))
    allocate (A(num_linhas, num_colunas))
    allocate (AT(num_colunas, num_linhas))
    allocate (ATA(num_colunas, num_colunas))
    allocate (ATb(num_colunas))
    allocate (ipiv(num_colunas))

    open (newunit=io, file=trim(filename), status="old", action="read")
    do i = 1, num_linhas
        read (io, *, iostat=io_status, iomsg=io_msg) t(i), y(i)
        if (io_status /= 0) then
            print *, "Erro ao ler linha", i, ":", trim(io_msg)
            stop
        end if
    end do
    close (io)

    x = [alpha, beta, delta, gama]
    max_iteracoes = 100
    done = .false.
    norma_anterior = huge(1.0d0)

!do while error > tolerancia

    do i = 1, max_iteracoes
        call calcula(A, num_linhas, num_colunas, t, y, x, b)
        norma_atual = sqrt(sum(b*b))

        if (norma_atual < tolerancia) then
            done = .true.
            print *, "Convergência alcançada!"
            print *, "Número de iterações: ", i
            exit
        end if

        norma_anterior = norma_atual

        AT = transpose(A)
        ATA = matmul(AT, A)
        ATb = matmul(AT, b)

        print *, i
        call dgesv(num_colunas, 1, ATA, num_colunas, ipiv, ATb, num_colunas, info)

        if (info /= 0) then
            print *, "Erro ao resolver o sistema linear: código", info
            exit
        end if

        d = ATb*-(1.0d0)

        write (*, *) x
        write (*, *) d
        x = x + d
        read (*, *)
    end do

    if (.not. done) then
        print *, "Atenção: Máximo de iterações atingido sem convergência"
    end if

    print *, "Número de iterações: ", i
    print *, "Resultado:"
    print *, "  alpha =", x(1)
    print *, "   beta =", x(2)
    print *, "  delta =", x(3)
    print *, "   gama =", x(4)

    deallocate (t, y, x, novo_x, d, A, b, AT, ATA, ATb, ipiv)

CONTAINS
    subroutine calcula(A, num_linhas, num_colunas, t, y, x, b)
        implicit none(external)
        integer, intent(in) :: num_linhas, num_colunas
        real(dp), intent(in) :: t(num_linhas), y(num_linhas), x(num_colunas)
        real(dp), intent(out) :: A(num_linhas, num_colunas), b(num_linhas)
        integer :: i
        real(dp) :: alpha, beta, delta, gama
        real(dp) :: exp_bt, sin_gt, cos_gt

        alpha = x(1)
        beta = x(2)
        delta = x(3)
        gama = x(4)

        do i = 1, num_linhas
            exp_bt = exp(beta*t(i))  !
            sin_gt = sin(gama*t(i))
            cos_gt = cos(gama*t(i))

            ! Calcula o valor de b (resíduo): b_i = y_i - f(t_i)
            b(i) = y(i) - (alpha*exp_bt + delta*sin_gt)

            ! Jacobiana negativa (para o método de Newton)
            A(i, 1) = exp_bt
            A(i, 2) = alpha*t(i)*exp_bt
            A(i, 3) = sin_gt
            A(i, 4) = delta*t(i)*cos_gt

        end do

    end subroutine calcula

end program residuos_quadraticos_nao_lineares
