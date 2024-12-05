program main
    use iso_Fortran_env, only: real64
    implicit none

    real(real64), allocatable :: A(:,:), b(:,:), x(:)
    integer :: i, j
    
    A = reshape(real([5,7,3,2,5,12,18,4,3], real64), [3,3], order=[2,1])
    b = reshape(real([8,9,2], real64), [3,1])
    call QR(A,b)

    call solSist(A,b,x)

    contains
        subroutine QR(A, b)
            ! Transforma la matriz A y el vector b para resolver un sistema triangular superior usando QR
            ! A: matriz de coeficientes
            ! b: vector (columna, i.e matriz nx1) de términos independientes
            real(real64) :: A(:,:), b(:,:), Hamp(size(A,1),size(A,2))
            real(real64), allocatable :: v(:), H(:,:)
            integer :: k

            if (size(A,1) /= size(A,2) .or. size(A,1) /= size(b,1) .or. size(b,2) /= 1) then
                ERROR STOP "El sistema no está bien planteado ya que las dimensiones de los coeficientes no coinciden"
            end if

            do k = 1, size(A,2)-1
                ! Calcula el vector vk = ak - ||ak||e1
                v = A(k:,k)
                v(1) = v(1) - norm2(v)
                

                ! Calcula la matriz de Householder (ampliada) Hamp
                H = -2/dot_product(v, v)*matmul(reshape(v, [size(v), 1]), reshape(v, [1, size(v)]))
                do i=1,size(H,1)
                    H(i,i) = 1 + H(i,i)
                end do
                call HouseholderAmpliada(H, size(A,1), Hamp)
                

                ! Actualiza las matrices A y b
                A = matmul(Hamp, A)
                b = matmul(Hamp, b)

            end do

        end subroutine QR

        subroutine HouseholderAmpliada(H, n, Hamp)
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
            real(real64), allocatable :: A(:,:), b(:,:), x(:)
            real(real64) :: aux
            integer :: i, j, n

            if (size(A,1) /= size(A,2) .or. size(A,1) /= size(b,1) .or. size(b,2) /= 1) then
                ERROR STOP "El sistema triangular a resolver no está bien planteado ya que las dimensiones no coinciden"
            end if

            n = size(A,1)
            print *, "X:"
            allocate(x(n))
            x(n) = b(n,1)/A(n,n)
            print *, n, x(n)
            do i=n-1,1,-1
                aux = 0
                do j=i+1,n
                    aux = aux + A(i,j)*x(j)
                end do
                x(i) = (b(i,1)-aux)/A(i,i)
                print *, i, x(i)
            end do
            
        end subroutine solSist

end program main
