module integral

    contains

    subroutine Riemann(a,b,h,f,Area)

        real(8), intent(in) ::a,b,h
        interface
            function f(x)
                real(8) :: x
                real(8) :: f
            end function
        end interface
        real(8), intent(inout)  :: Area


        real(8), allocatable    :: X(:)
        real(8)     :: dx
        integer     :: i,n

        ! Calculo los puntos en el intervalo
        n = (b-a)/h  !Numero de intervalos que puedo hacer (si la division no es entera me quedo con el entero)
        !   a = 0 |---|---|---|-| b=1  h = 0.3
        !   n = 1/0.3 = 3
        dx = (b-a)/n
        !   dx = 0.3333333333333
        allocate(X(0:n))
        do i = 0,n
            X(i) = a + i*dx
        enddo
        ! Método de Riemann con punto medio
        Area = 0
        do i = 1, n !bucle por intervalos
            Area = Area + dx*f((x(i)+x(i-1))/2)
        enddo
    end subroutine

    subroutine trapecio(a,b,h,f,area)
        real(8), intent(in) ::a,b,h
        interface
            function f(x)
                real(8) :: x
                real(8) :: f
            end function
        end interface
        real(8), intent(inout)  :: Area

        real(8), allocatable:: X(:)
        real(8) :: dx, sum  
        integer::i,n 
        n = (b-a)/h
        dx = (b-a)/n 
        
        
        allocate(X(0:n))
        do i = 0,n
            X(i) = a + i*dx
        enddo
        sum=0
        area=0
        do i=1, n-1
            sum=sum +2*f(x(i))
        enddo
        area= (dx/2d0)*(f(a)+sum+f(b))
    end subroutine 

    subroutine simpson(a, b, h, f, area)
        implicit none
        real(8), intent(in) ::a,b,h
        interface
            function f(x)
                real(8) :: x
                real(8) :: f
            end function
        end interface
        real(8), intent(inout)  :: Area

        real(8), allocatable:: X(:)
        real(8) :: dx, sum1, sum2 
        integer :: i,k,p, n 
        n = (b-a)/h 
        dx= (b-a)/n
        allocate(x(0:n))
            
        do i=1,n 
            x(i)=a+i*dx
        enddo

    sum1=0
    sum2=0
    area=0
    if ((n/2)*2==n) then
        k=n-1
        p=n-2
    else
        k=n-2
        p=n-1
    endif
    do i=1, k, 2
        sum1=sum1+4*f(x(i))
    enddo
    do i=2,p,2
        sum2=sum2+2*f(x(i))
    enddo
    area=(dx/3d0)*(f(a)+sum1+sum2+f(b))
    end subroutine simpson

end module 