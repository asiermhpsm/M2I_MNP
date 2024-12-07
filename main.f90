program main
    ! Resuelve un sistema de ecuaciones lineales Ax=b usando el método de QR
    ! Las matrices deben estar definidas en un archivo de texto con el siguiente formato:
    ! n
    ! a11 a12 ... a1n
    ! a21 a22 ... a2n
    ! ...
    ! an1 an2 ... ann
    ! b1
    ! b2
    ! ...
    ! bn
    ! Para compilar: gfortran module_QR.f90 main.f90
    ! Para ejecutar: .\a < matriz_flexibilidad.dat
    use module_QR, only: QR, real64
    implicit none

    real(real64), allocatable :: A(:,:), b(:,:), x(:)
    integer :: i, n

    ! Definicion del sistema de ecuaciones
    read*, n
    allocate(A(n,n), b(n,1), x(n))
    print *, "Matriz A:"
    do i = 1, n
        read*, A(i,:)
        print*,A(i,:)
    end do
    print *, "Vector b:"
    do i = 1, n
        read*, b(i,:)
        print*,b(i,:)
    end do

    ! Resolucion del sistema Ax=b
    call QR(A,b)
    call solSist(A,b,x)

    ! Impresion de la solucion
    print *, "La solucion del sistema es:"
    do i=1,size(x)
        print *, x(i)
    end do

    contains
        
        subroutine solSist(A, b, x)
            ! Resuelve un sistema triangular superior (cuadrado)
            ! A: matriz de coeficientes
            ! b: vector (columna, i.e matriz nx1) de términos independientes
            ! x: vector solucion
            real(real64), intent(in) :: A(:,:), b(:,:)
            real(real64), intent(out) :: x(:)
            real(real64) :: aux
            integer :: i, j, n

            if (size(A,1) /= size(A,2) .or. size(A,1) /= size(b,1) .or. size(b,2) /= 1) then
                ERROR STOP "El sistema triangular a resolver no está bien planteado ya que las dimensiones no coinciden"
            end if

            n = size(A,1)
            x(n) = b(n,1)/A(n,n)
            do i=n-1,1,-1
                aux = 0
                do j=i+1,n
                    aux = aux + A(i,j)*x(j)
                end do
                x(i) = (b(i,1)-aux)/A(i,i)
            end do
            
        end subroutine solSist

end program main
