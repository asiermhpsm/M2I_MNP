program leer_matriz
    implicit none
    integer :: filas, columnas, i, j
    real, allocatable :: matriz(:,:)
    character(len=100) :: nombre_fichero

    ! Especificar el nombre del archivo
    nombre_fichero = 'matriz.txt'

    ! Abrir el archivo
    open(unit=10, file=nombre_fichero, status='old', action='read')

    ! Leer las dimensiones de la matriz
    read(10, *) filas, columnas

    ! Asignar memoria para la matriz
    allocate(matriz(filas, columnas))

    ! Leer los datos de la matriz
    do i = 1, filas
        read(10, *) (matriz(i, j), j = 1, columnas)
    end do

    ! Cerrar el archivo
    close(10)

    ! Imprimir la matriz leída
    print *, "Matriz leída:"
    do i = 1, filas
        print *, (matriz(i, j), j = 1, columnas)
    end do

    ! Liberar memoria
    deallocate(matriz)
end program leer_matriz
