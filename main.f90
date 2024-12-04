program main
    use iso_Fortran_env, only: real64
    implicit none

    real(real64), allocatable :: A(:,:), b(:,:), H(:,:), Hamp(:,:)
    integer :: i, j
    
    A = reshape(real([1,2,3,4,5,6,7,8,9], real64), [3,3], order=[2,1])
    b = reshape(real([5,6], real64), [2,1])
    call QR(A)



    contains
        subroutine QR(A)
            ! Resuleve el sistema Ax = b usando la factorización QR
            ! A: matriz de coeficientes
            real(real64), allocatable :: A(:,:), v(:), H(:,:), Hamp(:,:)
            integer :: k

            do k = 1, size(A,2)-1
                ! Calcula el vector vk = ak - ||ak||e1
                v = A(k:,k)
                v(1) = v(1) - norm2(v)

                ! Calcula la matriz de Householder H
                call Householder(v, H)
                allocate(Hamp(size(A,1),size(A,1)))
                call HouseholderAmpliada(H, size(A,1), Hamp)
                

                ! Actualiza la matriz A
                A = matmul(Hamp, A)
                print *, "A:-----------------"
                do i=1, size(A,1)
                    do j=1, size(A,2)
                        print *, i, j, A(i,j)
                    end do
                end do

                deallocate(Hamp)
            end do

        end subroutine QR

        subroutine Householder(v, H)
            ! Calcula la matriz de Householder H = I - 2/(vT*v) * v*vT
            ! v: vector sobre el que se calcula la matriz
            ! H: matriz de Householder
            real(real64), intent(in), allocatable :: v(:)
            real(real64), intent(inout), allocatable :: H(:,:)
            real(real64), allocatable :: aux
            integer :: i,j

            H = matmul(reshape(v, [size(v), 1]), reshape(v, [1, size(v)]))
            aux = 2/dot_product(v, v)
            do i=1, size(H,1)
                do j=1, size(H,2)
                    if (i == j) then
                        H(i,j) = 1 - aux*H(i,j)
                    else
                        H(i,j) = -aux*H(i,j)
                    end if
                end do
            end do
            
        end subroutine Householder


        subroutine HouseholderAmpliada(H, n, Hamp)
            real(real64), intent(in), allocatable :: H(:,:)
            integer, intent(in) :: n
            real(real64), intent(inout), allocatable :: Hamp(:,:)
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
        

end program main
