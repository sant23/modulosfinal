module derivadas

real(8) :: h

contains

function d1_centrada(f,x)

    real(8)     :: x
    interface
        function f(x)
            real(8) :: x
            real(8) :: f
        end function
    end interface

    real(8) :: d1_centrada

    d1_centrada = (f(x+h)-f(x-h))/(2.d0*h)

end function

function d1_progresiva(f,x)
    real(8) :: x 
    interface 
        function f(x)
            real (8) :: x 
            real (8) :: f
        end function 
        end interface 
    real(8):: d1_progresiva 
    
    d1_progresiva = (f(x+h)-f(x))/(1.d0*h)  
    end function 

function d1_regresiva(f,x)
    real(8) :: x 
    interface 
        function f(x)
            real (8) :: x 
            real (8) :: f
        end function 
        end interface 
    real(8) :: d1_regresiva 
    d1_regresiva = (f(x)-f(x-h))/h  
    end function

function d2_progresiva(f,x)
    real(8) :: x 
    interface 
        function f(x)
            real (8) :: x 
            real (8) :: f
        end function 
        end interface 
        real(8):: d2_progresiva 
    d2_progresiva = (f(x+2.d0*h) + f(x)-2.d0*f(x+h))/(h**2)  

    end function 

    
function d2_regresiva(f,x)
        real(8) :: x 
        interface 
            function f(x)
                real (8) :: x 
                real (8) :: f
            end function 
            end interface
        real(8)::d2_regresiva  
        d2_regresiva = (f(x-2.d0*h) + f(x)-2.d0*f(x-h))/(h**2)  
    
end function 

function d2_centrada(f,x)
    real(8) :: x 
    interface 
        function f(x)
            real (8) :: x 
            real (8) :: f
        end function 
        end interface 
        real(8)::d2_centrada 
    d2_centrada = (f(x+h) + f(x-h)-2.d0*f(x))/(h**2)  

    end function 

    

end module