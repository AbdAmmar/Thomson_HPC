module parabola_test
implicit none
contains

    subroutine parabola(X,alpha_min)
	real*8, dimension(3,1), intent(in) :: X
	real*8,intent(out) :: alpha_min
	integer :: i,j
	real*8 :: A,B
	
	B=(X(2,1)-X(3,1))/2
	A=((X(2,1)+X(3,1)-2*X(1,1))/2)
	alpha_min = (-B)/(2*A)
	
	endsubroutine
    
end module