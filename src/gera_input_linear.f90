program gera_input_linear
    implicit none

    integer :: i, out, a, b, k
    real :: x, r, fator_aleatorio, noise, polinomio

    open (unit=out, file="input_linear.txt", status="replace", action="write")

    a = 1
    b = 2
    k = 5

    do i = 1, 100, +2
        x = x + 1
        !call random_number(r)
        !noise = k * r * (a*x + b)
        call random_number(r)
        call random_number(fator_aleatorio)
        noise = (k*2.0*r)*(fator_aleatorio*a*x + b**2)
        polinomio = a*x + b + noise
        write (out, "(F8.2,1X,F8.2)") x, polinomio
    end do

    close (out)

end program gera_input_linear
