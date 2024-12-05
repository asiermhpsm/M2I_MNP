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
    ! Para compilar: gfortran main.f90
    ! Para ejecutar: .\a < matriz_flexibilidad.dat
    use iso_Fortran_env, only: real64
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
        subroutine QR(A, b)
            ! Transforma la matriz A y el vector b para resolver un sistema triangular superior usando QR
            ! A: matriz de coeficientes
            ! b: vector (columna, i.e matriz nx1) de términos independientes
            real(real64) :: A(:,:), b(:,:), Hamp(size(A,1),size(A,2))
            real(real64), allocatable :: v(:), H(:,:)
            integer :: k, n

            if (size(A,1) /= size(A,2) .or. size(A,1) /= size(b,1) .or. size(b,2) /= 1) then
                ERROR STOP "El sistema no está bien planteado ya que las dimensiones de los coeficientes no coinciden"
            end if

            n = size(A,1)

            do k = 1, n-1
                if (esVectTriang(A, k)) then
                    cycle
                end if
                ! Calcula el vector vk = ak - ||ak||e1
                v = A(k:,k)
                v(1) = v(1) - norm2(v)
                
                ! Calcula la matriz de Householder (ampliada) Hamp
                H = -2/dot_product(v, v)*matmul(reshape(v, [size(v), 1]), reshape(v, [1, size(v)]))
                do i=1,n
                    H(i,i) = 1 + H(i,i)
                end do
                call HouseholderAmpliada(H, n, Hamp)
                
                ! Actualiza las matrices A y b
                A = matmul(Hamp, A)
                b = matmul(Hamp, b)

            end do

        end subroutine QR

        subroutine HouseholderAmpliada(H, n, Hamp)
            ! Dada una matriz de Householder H, devuelve la matriz de Householder ampliada de tamaño nxn Hamp
            ! H: matriz de Householder
            ! n: tamaño de la matriz ampliada
            ! Hamp: matriz de Householder ampliada
            real(real64), intent(in), allocatable :: H(:,:)
            integer, intent(in) :: n
            real(real64), intent(out) :: Hamp(:,:)
            integer :: i, j, dif

            dif = n - size(H,1)
            do i=1,n
                do j=1,n
                    if (i <= dif .or. j <= dif) then
                        if (i==j) then
                            Hamp(i,j) = 1
                        else
                            Hamp(i,j) = 0
                        end if
                    else
                        Hamp(i,j) = H(i-dif,j-dif)
                    end if
                end do
            end do
            
        end subroutine HouseholderAmpliada
        
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

        logical function esVectTriang(matriz, k)
            ! Devuelve true si la fila k de la matriz tiene todo 0s a partir de la fila k+1
            real(real64), intent(in) :: matriz(:,:)
            integer, intent(in) :: k
            integer :: i

            esVectTriang = .true.
            do i = k+1, size(matriz,1)
                if (abs(matriz(i, k)) > 1E-20) then
                    esVectTriang = .false.
                    return
                end if
            end do
        end function esVectTriang

end program main
