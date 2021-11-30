program algoritmoQD
implicit none

integer n, max_iter
real(8) tol
real(8), allocatable, dimension(:) :: P
complex(8), allocatable, dimension(:) :: raices

n = 3
max_iter = 100
tol = 0.0001
allocate(P(0:n), raices(1:n))

P(0) = -9.
P(1) = -9.
P(2) = 1.
P(3) = 1.

call QD(P, n, tol, max_iter, raices)
call mostrar_raices(raices, n, tol)

deallocate(P, raices)

contains
!###############################################################################
! P: coeficientes
! n: grado del polinomio
subroutine QD(P, n, tol, max_iter, raices)
    integer, intent(in) :: n, max_iter
    real(8), intent(in) :: tol
    real(8), intent(in), dimension(0:n) :: P
    complex(8), intent(out), dimension(1:n) :: raices

    real(8), dimension(1:n) :: q, q_ant
    real(8), dimension(0:n) :: e
    integer i, iter
    real(8) ei, u, v

    ! Inicializacion
    q(:) = 0.
    q(1) = -1. * P(n-1) / P(n)

    e(0) = 0.
    e(n) = 0.

    do i=1, n-1
        e(i) = P(n-i-1) / P(n-i)
    end do

    iter = 1
    ei = tol * 2.

    open(2, file='tabla.dat')
    write(*, '(A)') 'Generando tabla en "tabla.dat"'

    ! Cabecera
    write(2,'(a)', advance='no') 'iteración        e0         '
    do i=1, n
        write(2,'(a, i1, a, i1, a)', advance='no') 'q', i,'         e', i,'         '
    end do
    write(2,*)

    do while((iter <= max_iter) .and. (ei >= tol))
        ! Fila de tabla
        write(2,'(i5, a)', advance='no') iter, '      '
        do i=1, n
            write(2,'(f22.6)', advance='no')  q(i)
        end do
        write(2,*)
        do i=0, n
            write(2,'(f22.6)', advance='no') e(i)
        end do
        write(2,*)

        ! Calculo de los nuevos valores
        q_ant = q
        do i=1, n
            q(i) = e(i) - e(i-1) + q(i)
        end do

        do i=1, n-1
            e(i) = e(i) * q(i+1) / q(i)
        end do

        ei = maxval(abs(e))
        iter = iter + 1
    end do

    ! Comprobamos la existencia de raices complejas
    if(iter > max_iter) then
        do i=0, n-1
            if (abs(e(i)) >= tol) then
                u = q(i) + q(i+1)
                v = -(q_ant(i) * q(i+1))
                CALL bairstow (P, N, u, v, tol)
                CALL resolvente(-u, -v, raices(i), raices(i+1))
            else
                raices(i+1) = cmplx(q(i+1), 0., kind=8)
            end if
        end do
    else
        raices(:) = q(:)
    end if

    close(2)
end subroutine QD

!###############################################################################

subroutine bairstow(Pol, n, u, v, tol)
    integer, intent(in) :: n
    real(8), dimension(0:n), intent(in) :: Pol
    real(8), intent(inout) :: u, v
    real(8), intent(in) :: tol

    real(8), dimension(0:2) :: q
    real(8), dimension(0:3) :: p
    real(8) h, k, e
    integer i

    e = tol * 2.

    do while(e >= tol)
        q(:) = 0.
        p(:) = 0.

        do i=n, 1, -1
            q(0) = Pol(i) + u * q(1) + v * q(2)
            q(2) = q(1)
            q(1) = q(0)

            p(0) = q(0) + u * p(1) + v * p(2)
            p(3) = p(2)
            p(2) = p(1)
            p(1) = p(0)
        end do
        q(0) = Pol(0) + u * q(1) + v * q(2)

        h = (q(0) * p(3) - q(1) * p(2)) / (p(2)**2. - p(1) * p(3))
        k = (q(1) * p(1) - q(0) * p(2)) / (p(2)**2. - p(1) * p(3))

        u = u + h
        v = v + k

        e = maxval(abs(q(0:1)))
    end do
end subroutine bairstow

!###############################################################################

subroutine resolvente(u, v, c1, c2)
    real(8), intent(in) :: u, v
    complex(8), intent(out) :: c1, c2

    c1 = -u + sqrt(u**2 - 4. * v) / 2.
    c2 = -u - sqrt(u**2 - 4. * v) / 2.
end subroutine resolvente

!###############################################################################

subroutine mostrar_raices(raices, n, tol)
    integer, intent(in) :: n
    complex(8), dimension(1:n), intent(in) :: raices
    real(8), intent(in) :: tol

    integer i

    open(2, file='raices.dat')

    write(*, '(A)') 'Las raices obtenidas se encuentran en "raices.dat"'
    do i=1, n
        ! Raiz real
        if(imag(raices(i)) == 0.) then
            write(2, '(A, I2, A, F20.10)') 'X', i, '= ', real(raices(i))

        ! Raiz  solo imaginaria
        else if(real(raices(i)) <= tol) then
            write(2, '(A, I2, A, F20.10, A)') 'X', i, '= ', imag(raices(i)), 'i'
        ! Raiz compleja
        else
            write(2, '(A, I2, A, 2F20.10, A)') 'X', i, '= ', real(raices(i)), imag(raices(i)), 'i'
        end if
    end do

    close(2)

end subroutine
end program

