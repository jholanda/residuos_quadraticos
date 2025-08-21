program gera_input_nao_linear
    implicit none
    integer :: i, out
    real :: t, r, noise, y
    real :: alpha, beta, delta, gama, k

    open (newunit=out, file="input_nao_linear.txt", status="replace", action="write")

    noise = 0.6
    alpha = 1.0
    beta = -0.2
    delta = 2.5
    gama = 0.3

    do i = 1, 100
        t = real(i)
        call random_number(noise)

        y = alpha*exp(beta*t) + delta*sin(gama*t) + noise
        write (out, "(F5.1,1X,F12.2)") t, y

    end do

    close (out)

end program gera_input_nao_linear
