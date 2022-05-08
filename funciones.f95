module funciones

contains

function f1(x)

    real(8)     :: x
    real(8)     :: f1

    f1 = x**2

end function

function f2(x)

    real(8)     :: x
    real(8)     :: f2

    f2 = cos(x)

end function

function f4(x)

    real(8)     :: x
    real(8)     :: f4 

    f4=sin(x)
end function 

function f3(x)

    real(8)     :: x
    real(8)     :: f3

    f3 = exp(x**2)

end function

function df1(x)

    real(8)     :: x
    real(8)     :: df1

    df1 = 2*x

end function

end module