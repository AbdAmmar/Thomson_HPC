module print_g
	implicit none
	contains
    
    subroutine prin(N,X,D,X_pos)
	implicit none 
	
	real*16, dimension(:,:),intent(in) :: X 
	integer, intent(in) :: N , D 
	integer :: i, j , coor
	real*16, dimension(N,D) :: X_pos 
    
	
	coor = N*D
	if (D == 1 ) then 
	do i = 1,coor
	X_pos(:,1) = X(1:N,1)
	end do
	else if (D == 2) then
	do i = 1,coor 
	X_pos(:,1) = X(1:N,1)
	X_pos(:,2) = X(N+1:2*N,1)
	end do
	else if (D == 3) then
	do i = 1,coor
	X_pos(:,1) = X(1:N,1)
	X_pos(:,2) = X(N+1:2*N,1)
	X_pos(:,3) = X(2*N+1:3*N,1)
	end do
	end if 
	
	
	write(*,'(a)') ""
	if (D == 3) then 
	write (*,'(a)') "                 =========================================================================="
	write (*,'(a)') "                           X                         Y                         Z"
	write (*,'(a)') "                 =========================================================================="
	else if (D == 2) then 
	write (*,'(a)') "                 ==============================================="
	write (*,'(a)') "                           X                         Y"
	write (*,'(a)') "                 ==============================================="
	else 
	write (*,'(a)') "                 ===================="
	write (*,'(a)') "                           X"
	write (*,'(a)') "                 ===================="
	end if 
	write(*,'(a)') ""
	DO i = 1, N
		write (*,'(I4,a,f26.16,f26.16,f26.16)') i,"      " , X_pos(i, :)
	END DO
	
	endsubroutine
    

    subroutine prin_distance(X,N,D,X_distance,Lx,Ly,Lz)
	implicit none 
	
	real*16, dimension(:,:),intent(in) :: X 
	integer, intent(in) :: N , D
	integer :: i, j
	real*16, dimension(N,N) :: X_distance 
    real*16, dimension(N,N) :: r
    real*16, intent(in) :: Lx, Ly , Lz
	
   
	if (D == 1) then
    
	do i = 1,N
        do j = 1 , N
                X_distance(i,j) = sqrt((X(i,1)-X(j,1))**2)
        end do 
	end do
    
    else if (D == 2 ) then 
    
    do i = 1,N
        do j = 1 , N
                X_distance(i,j) = sqrt((X(i,1)-X(j,1))**2+(X(i+N,1)-X(j+N,1))**2)
        end do 
	end do
    
    else if (D == 3) then 
    do i = 1,N
        do j = 1 , N
                X_distance(i,j) = sqrt((X(i,1)-X(j,1))**2 + (X(i+N,1)-X(j+N,1))**2 + (X(i+2*N,1)-X(j+2*N,1))**2)
        end do 
	end do
    
	end if 
    
    write(*,'(a)') "____________________________________________________________________________________________"
    write(*,'(a)') ""
    write(*,'(a)') '                                    The distance matrix (Geodesic)'
    write(*,'(a)')
	write(*,'(a)') ""
    
	do i = 1, N
            write (*,'(a,I3,a,((1x,1000f26.16)))') "(",i,")" , X_distance(i, :)
	end do 
	
    X_distance = 0.d0
    r = 0.d0 
    
    if (D == 1) then 
    
    do i = 1,N
        do j = 1 , N
            r(i,j) = ((Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1))))))**(0.5d0)
        end do 
	end do
    
    
    else if (D == 2) then 
    
    do i = 1,N
        do j = 1 , N
            r(i,j) = ((Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1))))))**(0.5d0)
        end do 
	end do
    
    
    else if (D == 3) then 
    
    do i = 1,N
        do j = 1 , N
            r(i,j) = ((Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1))))+Lz**(-2)*(2.d0-2.d0*cos(Lz*(X(i+2*N,1)-X(j+2*N,1))))))**(0.5d0)
        end do 
	end do
    
	end if
    
    
    write(*,'(a)') "____________________________________________________________________________________________"
    write(*,'(a)') ""
    write(*,'(a)') '                                    The distance matrix (Euclidean)' 
    write(*,'(a)')
	write(*,'(a)') ""
    
    
    do i = 1, N
            write (*,'(a,I3,a,((1x,1000f26.16)))') "(",i,")" , r(i, :)
	end do
    
    
    
    
    
    
	endsubroutine

    


end module