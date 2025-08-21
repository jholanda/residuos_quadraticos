program residuos_quadraticos_linear
    implicit none(external)

    integer :: io, io_status, i, j, n, m, nrhs, lda, ldb, info
    integer, allocatable :: ipiv(:)
    character(len=512) :: io_msg, filename, arg, fmt_string
    logical :: file_exists
    real :: work(50)
    real(kind=8), allocatable :: A(:, :), AT(:, :), ATA(:, :), b(:), ATb(:), x(:)

    external dgesv
    if (command_argument_count() < 2) then
        filename = "input_linear.txt"
        arg = "2"
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

    m = 0
    open (newunit=io, file=trim(filename), status="old", action="read", iostat=io_status)
    do
        read (io, *, iostat=io_status)
        if (io_status /= 0) exit
        m = m + 1
    end do
    close (io)

    allocate (A(m, n))
    allocate (AT(n, m))
    allocate (ATA(n, n))
    allocate (b(m))
    allocate (ATb(n))
    allocate (x(n))
    allocate (ipiv(n))

    open (newunit=io, file=trim(filename), status="old", action="read")
    do i = 1, m
        A(i, 1) = 1.0
        read (io, *, iostat=io_status, iomsg=io_msg) A(i, 2), b(i)

        do j = 3, n
            A(i, j) = A(i, j - 1)*A(i, 2)
        end do
    end do
    close (io)

    AT = transpose(A)
    ATA = matmul(AT, A)
    ATb = matmul(AT, b)

    !call dgels('N', n, n, 1, ATA, n, ATb, n, work, 50, info)

    call dgesv(n, 1, ATA, n, ipiv, ATb, n, info)

    if (info /= 0) then
        print *, "Erro ao resolver o sistema: código", info
        stop
    end if

    x = ATb(1:n)

    open (unit=20, file=trim(filename)//"_grau"//trim(arg)//".txt", status="replace", action="write")
    write (20, "(15F18.3)") (x(j), j=1, n)
    flush (20)
    close (20)

    deallocate (A, AT, ATA, b, ATb, x, ipiv)

    print *, "Coeficientes calculados com sucesso e salvos em: ", trim(filename)//"_grau"//trim(arg)//".txt"

end program residuos_quadraticos_linear
