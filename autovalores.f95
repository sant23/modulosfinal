module autovalores
    use lineal 

    contains

    subroutine potencia_max(A,evalue,evector,n)
    real(8), intent(in)::A(n,n)
    integer, intent(in)::n 
    real(8), intent(out):: evalue 
    real(8), intent(out)::evector(n)

    integer:: iter 
    real(8):: evector_aux(n), evalue_aux
    !inicializacion del autovalor y autovector
    evector = 1.d0/sqrt(n*1.d0)
    evalue= 0.5d0

    do iter=1,1000
        evector_aux = matmul(A,evector)
        evalue_aux = dot_product(evector,evector_aux)/dot_product(evector,evector)
        !write(*,*) iter, evector, evalue_aux
        !normalizacion del autovector con la norma 2
        evector= evector_aux/norm2(evector_aux)
        !convergencia
        if(abs((evalue_aux-evalue)/evalue)<epsilon(1.d0)) exit 

        evalue= evalue_aux 
    enddo 
end subroutine 

    subroutine potencia_inversa(A,evalue,evector,n)
    
    real(8), intent(in) :: A(n,n)
    integer, intent(in) :: n 
    real(8), intent(out) ::evalue
    real(8), intent(out) :: evector(n)

    !variables locales
    integer :: iter 
    real(8):: evector_aux(n), y(n), baux(n), evalue_aux
    real(8) :: L(n,n), U(n,n), P(n,n)
    call LU_factor(A,L,U,P)

    !inicializacion del autovector y el autovalor 
    evector = 1.d0/sqrt(n*1.d0) !evector iniciado
    evalue = 0.5d0 

    do iter=1,1000 
        !calculo del autovector resolviendo el sistema lineal con LU
        baux = matmul(P,evector)
        call lower_solver(L,baux,y)
        call upper_solver(U,y,evector_aux)
        !calculo del autovalor
        !evalue_aux = dot_product(evector, evector_aux)/dot_product(evector,evector)
        evalue_aux = evector(1)/evector_aux(1)
        !write(*,*) iter,evector, evalue_aux
        evector = evector_aux/norm2(evector_aux)
        !convergencia 
        if (abs((evalue_aux-evalue)/evalue)<epsilon(1.d0))exit 
        evalue = evalue_aux 
    enddo 
end subroutine 

end module 
