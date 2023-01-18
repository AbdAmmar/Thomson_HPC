module Thomson_problem
	implicit none
	contains
    
    subroutine energy(N,X,D,Lx,Ly,Lz,E)
	
	implicit none 

	real*16, dimension(:,:), intent(in) :: X
	integer, intent(in) :: N , D
	real*16, intent(in) :: Lx, Ly , Lz
	real*16, dimension(N),intent(out) :: E
	integer :: i,j
	real*16 :: ax , ay , az 
	real*16, parameter :: pi = ACOS(-1.0d0)
	
	
	
	
	E=0.d0

	if (D == 3) then 
   
	do i = 1,N
		do j=1,N
		if (j > i) then
		
		ax = (X(i,1)-X(j,1))
		if (ax > 0.5d0 * Lx*2.d0*pi) then 
			ax = ax - Lx*2.d0*pi
		else if (ax <= -0.5d0 * Lx*2.d0*pi) then
			ax = ax + Lx*2.d0*pi
		end if
		
		ay = (X(i+N,1)-X(j+N,1))
		if (ay > 0.5d0 * Ly*2.d0*pi) then 
			ay = ay - Ly*2.d0*pi
		else if (ay <= -0.5d0 * Ly*2.d0*pi) then
			ay = ay + Ly*2.d0*pi
		end if
		
		az = (X(i+2*N,1)-X(j+2*N,1))
		if (az > 0.5d0 * Lz*2.d0*pi) then 
			az = az - Lz*2.d0*pi
		else if (az <= -0.5d0 * Lz*2.d0*pi) then
			az = az + Lz*2.d0*pi
		end if
		
			E(i)=E(i)+ (Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay))+Lz**(-2)*(2.d0-2.d0*cos(Lz*az)))**(-0.5d0)
			
		end if 
		end do
	end do 
	 
	end if 

	if (D == 2) then 
    
	do i = 1,N
		do j=1,N
		if (j > i) then
		
		ax = (X(i,1)-X(j,1))
		if (ax > 0.5d0 * Lx*2.d0*pi) then 
			ax = ax - Lx*2.d0*pi
		else if (ax <= -0.5d0 * Lx*2.d0*pi) then
			ax = ax + Lx*2.d0*pi
		end if
		
		ay = (X(i+N,1)-X(j+N,1))
		if (ay > 0.5d0 * Ly*2.d0*pi) then 
			ay = ay - Ly*2.d0*pi
		else if (ay <= -0.5d0 * Ly*2.d0*pi) then
			ay = ay + Ly*2.d0*pi
		end if
		
			E(i)=E(i)+ (Lx**(-2.d0)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2.d0)*(2.d0-2.d0*cos(Ly*ay)))**(-0.5d0)
		end if 
		end do
	end do 
     
	end if 

	if (D == 1) then 
    
	do i = 1,N
		do j=1,N
		if (j > i) then
		
		ax = (X(i,1)-X(j,1))
		if (ax > 0.5d0 * Lx*2.d0*pi) then 
			ax = ax - Lx*2.d0*pi
		else if (ax < 0.5d0 * Lx*2.d0*pi) then
			ax = ax + Lx*2.d0*pi
		end if
		
			E(i)=E(i)+(Lx**(-2.d0)*(2.d0-2.d0*cos(Lx*ax)))**(-0.5d0)
		end if 
		end do
	end do 
    
	end if 
	
	endsubroutine
    
    subroutine diff(N,X,D,Lx,Ly,Lz,tdf)
	
	implicit none
	
	real*16, dimension(:,:), intent(in) :: X
	integer, intent(in) :: N ,D 
	real*16 , intent(in) :: Lx,Ly,Lz
	real*16, dimension(3*N,3*N) :: df
	real*16 , dimension(3*N) :: tdf
	integer :: i,j
	real*16 :: ax , ay , az 
	real*16, parameter :: pi = ACOS(-1.0d0)
	
	
	df=0.d0
	tdf=0.d0 
		
	if ( D == 3 ) then 
	
	
	do i=1,N
		do j=1,N
			if (j /= i) then
			
			ax = (X(i,1)-X(j,1))
			if (ax > 0.5d0 * Lx*2.d0*pi) then 
				ax = ax - Lx*2.d0*pi
			else if (ax <= -0.5d0 * Lx*2.d0*pi) then
				ax = ax + Lx*2.d0*pi
			end if
			
			ay = (X(i+N,1)-X(j+N,1))
			if (ay > 0.5d0 * Ly*2.d0*pi) then 
				ay = ay - Ly*2.d0*pi
			else if (ay <= -0.5d0 * Ly*2.d0*pi) then
				ay = ay + Ly*2.d0*pi
			end if
			
			az = (X(i+2*N,1)-X(j+2*N,1))
			if (az > 0.5d0 * Lz*2.d0*pi) then 
				az = az - Lz*2.d0*pi
			else if (az <= -0.5d0 * Lz*2.d0*pi) then
				az = az + Lz*2.d0*pi
			end if
				
			df(i,j) = - (Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay))+Lz**(-2)*(2.d0-2.d0*cos(Lz*az)))**(-1.5d0)*Lx**(-1.d0)*sin(Lx*ax)
			end if 
		end do
	end do

	do i=N+1,2*N
		do j=N+1,2*N
			if (j /= i) then
			
			ax = (X(i-N,1)-X(j-N,1))
			if (ax > 0.5d0 * Lx*2.d0*pi) then 
				ax = ax - Lx*2.d0*pi
			else if (ax <= -0.5d0 * Lx*2.d0*pi) then
				ax = ax + Lx*2.d0*pi
			end if
			
			ay = (X(i,1)-X(j,1))
			if (ay > 0.5d0 * Ly*2.d0*pi) then 
				ay = ay - Ly*2.d0*pi
			else if (ay <= -0.5d0 * Ly*2.d0*pi) then
				ay = ay + Ly*2.d0*pi
			end if
			
			az = (X(i+N,1)-X(j+N,1))
			if (az > 0.5d0 * Lz*2.d0*pi) then 
				az = az - Lz*2.d0*pi
			else if (az <= -0.5d0 * Lz*2.d0*pi) then
				az = az + Lz*2.d0*pi
			end if
			
			df(i,j) =  - (Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay))+Lz**(-2)*(2.d0-2.d0*cos(Lz*az)))**(-1.5d0)*Ly**(-1.d0)*sin(Ly*ay)
			end if 
		end do
	end do

	do i=2*N+1,3*N
		do j=2*N+1,3*N
			if (j /= i) then
			
			ax = (X(i-2*N,1)-X(j-2*N,1))
			if (ax > 0.5d0 * Lx*2.d0*pi) then 
				ax = ax - Lx*2.d0*pi
			else if (ax <= -0.5d0 * Lx*2.d0*pi) then
				ax = ax + Lx*2.d0*pi
			end if
			
			ay = (X(i-N,1)-X(j-N,1))
			if (ay > 0.5d0 * Ly*2.d0*pi) then 
				ay = ay - Ly*2.d0*pi
			else if (ay <= -0.5d0 * Ly*2.d0*pi) then
				ay = ay + Ly*2.d0*pi
			end if
			
			az = (X(i,1)-X(j,1))
			if (az > 0.5d0 * Lz*2.d0*pi) then 
				az = az - Lz*2.d0*pi
			else if (az <= -0.5d0 * Lz*2.d0*pi) then
				az = az + Lz*2.d0*pi
			end if
			
			df(i,j) =  - (Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay))+Lz**(-2)*(2.d0-2.d0*cos(Lz*az)))**(-1.5d0)*Lz**(-1.d0)*sin(Lz*az)
			end if 
		end do
	end do
    
	tdf = sum(df,2)
	end if 

	if (D == 2) then 

	do i=1,N
		do j=1,N
			if (j /= i) then
			
			ax = (X(i,1)-X(j,1))
			if (ax > 0.5d0 * Lx*2.d0*pi) then 
				ax = ax - Lx*2.d0*pi
			else if (ax <= -0.5d0 * Lx*2.d0*pi) then
				ax = ax + Lx*2.d0*pi
			end if
			
			ay = (X(i+N,1)-X(j+N,1))
			if (ay > 0.5d0 * Ly*2.d0*pi) then 
				ay = ay - Ly*2.d0*pi
			else if (ay <= -0.5d0 * Ly*2.d0*pi) then
				ay = ay + Ly*2.d0*pi
			end if
			
			df(i,j) = -Lx**(-1)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay)))**(-1.5d0)*sin(Lx*ax)
			end if 
		end do
	end do
	
	do i=N+1,2*N
		do j=N+1,2*N
			if (j /= i) then
			
			ax = (X(i-N,1)-X(j-N,1))
			if (ax > 0.5d0 * Lx*2.d0*pi) then 
				ax = ax - Lx*2.d0*pi
			else if (ax <= -0.5d0 * Lx*2.d0*pi) then
				ax = ax + Lx*2.d0*pi
			end if
			
			ay = (X(i,1)-X(j,1))
			if (ay > 0.5d0 * Ly*2.d0*pi) then 
				ay = ay - Ly*2.d0*pi
			else if (ay <= -0.5d0 * Ly*2.d0*pi) then
				ay = ay + Ly*2.d0*pi
			end if
			
			df(i,j) = -Ly**(-1)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay)))**(-1.5d0)*sin(Ly*ay)
			end if 
		end do
	end do
	tdf = sum(df,2)
	end if 

	if (D == 1) then 

	do i=1,N
		do j=1,N
			if (j /= i) then
			
			ax = (X(i,1)-X(j,1))
			if (ax > 0.5d0 * Lx*2.d0*pi) then 
				ax = ax - Lx*2.d0*pi
			else if (ax <= -0.5d0 * Lx*2.d0*pi) then
				ax = ax + Lx*2.d0*pi
			end if
			
			df(i,j) = -Lx**(-1.d0)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax)))**(-1.5d0)*sin(Lx*ax)
			end if 
		end do
	end do
	tdf = sum(df,2)
	end if
	
	endsubroutine
    
    subroutine hessien(N,X,D,Lx,Ly,Lz,H)
	
	implicit none
	
	real*16, dimension(:,:), intent(in) :: X
	integer, intent(in) ::N,D
	real*16 , intent(in) :: Lx,Ly,Lz
	real*8 ,dimension (N,N) :: Rx , Ry ,Rz , Rxy, Rxz , Ryz
	real*8, dimension(:,:),allocatable :: H
	integer :: i,j
	real*16 :: ax , ay , az 
	real*16, parameter :: pi = ACOS(-1.0d0)

			Rx=0.d0
			Ry=0.d0
			Rz=0.d0
			Rxy=0.d0
			Rxz=0.d0
			Ryz=0.d0
			H=0.d0
			
	if (D == 3) then

	!! calculate the diagonal
    

	do i=1,3*N
		do j=1,3*N
			if (j /= i .and. j<= N .and. i<=N) then
			
			ax = (X(i,1)-X(j,1))
			if (ax > 0.5d0 * Lx*2.d0*pi) then 
				ax = ax - Lx*2.d0*pi
			else if (ax <= -0.5d0 * Lx*2.d0*pi) then
				ax = ax + Lx*2.d0*pi
			end if
			
			ay = (X(i+N,1)-X(j+N,1))
			if (ay > 0.5d0 * Ly*2.d0*pi) then 
				ay = ay - Ly*2.d0*pi
			else if (ay <= -0.5d0 * Ly*2.d0*pi) then
				ay = ay + Ly*2.d0*pi
			end if
			
			az = (X(i+2*N,1)-X(j+2*N,1))
			if (az > 0.5d0 * Lz*2.d0*pi) then 
				az = az - Lz*2.d0*pi
			else if (az <= -0.5d0 * Lz*2.d0*pi) then
				az = az + Lz*2.d0*pi
			end if
			
			H(i,i) = H(i,i) +3*Lx**(-2.d0)*sin(Lx*ax)**2*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay))+Lz**(-2)*(2.d0-2.d0*cos(Lz*az)))**(-2.5d0)-cos(Lx*ax)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay))+Lz**(-2)*(2.d0-2.d0*cos(Lz*az)))**(-1.5d0)

			else if (j /= i .and. j> N .and. j<=2*N .and. i > N .and. i<=2*N) then
			
			ax = (X(i-N,1)-X(j-N,1))
			if (ax > 0.5d0 * Lx*2.d0*pi) then 
				ax = ax - Lx*2.d0*pi
			else if (ax <= -0.5d0 * Lx*2.d0*pi) then
				ax = ax + Lx*2.d0*pi
			end if
			
			ay = (X(i,1)-X(j,1))
			if (ay > 0.5d0 * Ly*2.d0*pi) then 
				ay = ay - Ly*2.d0*pi
			else if (ay <= -0.5d0 * Ly*2.d0*pi) then
				ay = ay + Ly*2.d0*pi
			end if
			
			az = (X(i+N,1)-X(j+N,1))
			if (az > 0.5d0 * Lz*2.d0*pi) then 
				az = az - Lz*2.d0*pi
			else if (az <= -0.5d0 * Lz*2.d0*pi) then
				az = az + Lz*2.d0*pi
			end if
			
			H(i,i) = H(i,i) +3*Ly**(-2.d0)*sin(Ly*ay)**2*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay))+Lz**(-2)*(2.d0-2.d0*cos(Lz*az)))**(-2.5d0)-cos(Ly*ay)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay))+Lz**(-2)*(2.d0-2.d0*cos(Lz*az)))**(-1.5d0)

			else if (j /= i .and. j> 2*N .and. i > 2*N ) then
			
			ax = (X(i-2*N,1)-X(j-2*N,1))
			if (ax > 0.5d0 * Lx*2.d0*pi) then 
				ax = ax - Lx*2.d0*pi
			else if (ax <= -0.5d0 * Lx*2.d0*pi) then
				ax = ax + Lx*2.d0*pi
			end if
			
			ay = (X(i-N,1)-X(j-N,1))
			if (ay > 0.5d0 * Ly*2.d0*pi) then 
				ay = ay - Ly*2.d0*pi
			else if (ay <= -0.5d0 * Ly*2.d0*pi) then
				ay = ay + Ly*2.d0*pi
			end if
			
			az = (X(i,1)-X(j,1))
			if (az > 0.5d0 * Lz*2.d0*pi) then 
				az = az - Lz*2.d0*pi
			else if (az <= -0.5d0 * Lz*2.d0*pi) then
				az = az + Lz*2.d0*pi
			end if
			
			H(i,i) = H(i,i) +3*Lz**(-2.d0)*sin(Lz*az)**2*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay))+Lz**(-2)*(2.d0-2.d0*cos(Lz*az)))**(-2.5d0)-cos(Lz*az)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay))+Lz**(-2)*(2.d0-2.d0*cos(Lz*az)))**(-1.5d0)

			end if 
		end do
	end do
    
    
	!! calculate R(x) !!


	do i = 1,N
		do j = 1,N
			if (j /= i) then
			
			ax = (X(i,1)-X(j,1))
			if (ax > 0.5d0 * Lx*2.d0*pi) then 
				ax = ax - Lx*2.d0*pi
			else if (ax <= -0.5d0 * Lx*2.d0*pi) then
				ax = ax + Lx*2.d0*pi
			end if
			
			ay = (X(i+N,1)-X(j+N,1))
			if (ay > 0.5d0 * Ly*2.d0*pi) then 
				ay = ay - Ly*2.d0*pi
			else if (ay <= -0.5d0 * Ly*2.d0*pi) then
				ay = ay + Ly*2.d0*pi
			end if
			
			az = (X(i+2*N,1)-X(j+2*N,1))
			if (az > 0.5d0 * Lz*2.d0*pi) then 
				az = az - Lz*2.d0*pi
			else if (az <= -0.5d0 * Lz*2.d0*pi) then
				az = az + Lz*2.d0*pi
			end if
			
			Rx(i,j) = -3.d0*Lx**(-2.d0)*sin(Lx*ax)**2*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay))+Lz**(-2)*(2.d0-2.d0*cos(Lz*az)))**(-2.5d0)+cos(Lx*ax)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay))+Lz**(-2)*(2.d0-2.d0*cos(Lz*az)))**(-1.5d0)

			end if 
		end do 
	end do
	
	!! calculate R(y) !!
    
    
	do i = 1,N
		do j = 1,N
			if (j /= i) then
			
			ax = (X(i,1)-X(j,1))
			if (ax > 0.5d0 * Lx*2.d0*pi) then 
				ax = ax - Lx*2.d0*pi
			else if (ax <= -0.5d0 * Lx*2.d0*pi) then
				ax = ax + Lx*2.d0*pi
			end if
			
			ay = (X(i+N,1)-X(j+N,1))
			if (ay > 0.5d0 * Ly*2.d0*pi) then 
				ay = ay - Ly*2.d0*pi
			else if (ay <= -0.5d0 * Ly*2.d0*pi) then
				ay = ay + Ly*2.d0*pi
			end if
			
			az = (X(i+2*N,1)-X(j+2*N,1))
			if (az > 0.5d0 * Lz*2.d0*pi) then 
				az = az - Lz*2.d0*pi
			else if (az <= -0.5d0 * Lz*2.d0*pi) then
				az = az + Lz*2.d0*pi
			end if
			
			Ry(i,j) = -3.d0*Ly**(-2.d0)*sin(Ly*ay)**2*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay))+Lz**(-2)*(2.d0-2.d0*cos(Lz*az)))**(-2.5d0)+cos(Ly*ay)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay))+Lz**(-2)*(2.d0-2.d0*cos(Lz*az)))**(-1.5d0)

			end if 
		end do 
	end do
    


	!! calculate R(z) !!
    

    
	do i = 1,N
		do j = 1,N
			if (j /= i) then
			
			ax = (X(i,1)-X(j,1))
			if (ax > 0.5d0 * Lx*2.d0*pi) then 
				ax = ax - Lx*2.d0*pi
			else if (ax <= -0.5d0 * Lx*2.d0*pi) then
				ax = ax + Lx*2.d0*pi
			end if
			
			ay = (X(i+N,1)-X(j+N,1))
			if (ay > 0.5d0 * Ly*2.d0*pi) then 
				ay = ay - Ly*2.d0*pi
			else if (ay <= -0.5d0 * Ly*2.d0*pi) then
				ay = ay + Ly*2.d0*pi
			end if
			
			az = (X(i+2*N,1)-X(j+2*N,1))
			if (az > 0.5d0 * Lz*2.d0*pi) then 
				az = az - Lz*2.d0*pi
			else if (az <= -0.5d0 * Lz*2.d0*pi) then
				az = az + Lz*2.d0*pi
			end if
			
			Rz(i,j) = -3.d0*Lz**(-2.d0)*sin(Lz*az)**2*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay))+Lz**(-2)*(2.d0-2.d0*cos(Lz*az)))**(-2.5d0)+cos(Lz*az)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay))+Lz**(-2)*(2.d0-2.d0*cos(Lz*az)))**(-1.5d0)
			
			end if 
		end do 
	end do
	

    
	!! calculate Rxy !!
    

    
	do i = 1,N
		do j = 1,N
			if (j /= i) then
			
			ax = (X(i,1)-X(j,1))
			if (ax > 0.5d0 * Lx*2.d0*pi) then 
				ax = ax - Lx*2.d0*pi
			else if (ax <= -0.5d0 * Lx*2.d0*pi) then
				ax = ax + Lx*2.d0*pi
			end if
			
			ay = (X(i+N,1)-X(j+N,1))
			if (ay > 0.5d0 * Ly*2.d0*pi) then 
				ay = ay - Ly*2.d0*pi
			else if (ay <= -0.5d0 * Ly*2.d0*pi) then
				ay = ay + Ly*2.d0*pi
			end if
			
			az = (X(i+2*N,1)-X(j+2*N,1))
			if (az > 0.5d0 * Lz*2.d0*pi) then 
				az = az - Lz*2.d0*pi
			else if (az <= -0.5d0 * Lz*2.d0*pi) then
				az = az + Lz*2.d0*pi
			end if
			
			    Rxy(i,j) = - 3.d0 *Lx**(-1)*Ly**(-1)* sin(Lx*ax)*sin(Ly*ay)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay))+Lz**(-2)*(2.d0-2.d0*cos(Lz*az)))**(-2.5d0)
			end if 
		end do 
	end do	
	
    
    !!$omp end parallel do
    
	!! calculate Rxz !!
    
    !!$omp parallel do collapse(2)
    
	do i = 1,N
		do j = 1,N
			if (j /= i) then
			
			ax = (X(i,1)-X(j,1))
			if (ax > 0.5d0 * Lx*2.d0*pi) then 
				ax = ax - Lx*2.d0*pi
			else if (ax <= -0.5d0 * Lx*2.d0*pi) then
				ax = ax + Lx*2.d0*pi
			end if
			
			ay = (X(i+N,1)-X(j+N,1))
			if (ay > 0.5d0 * Ly*2.d0*pi) then 
				ay = ay - Ly*2.d0*pi
			else if (ay <= -0.5d0 * Ly*2.d0*pi) then
				ay = ay + Ly*2.d0*pi
			end if
			
			az = (X(i+2*N,1)-X(j+2*N,1))
			if (az > 0.5d0 * Lz*2.d0*pi) then 
				az = az - Lz*2.d0*pi
			else if (az <= -0.5d0 * Lz*2.d0*pi) then
				az = az + Lz*2.d0*pi
			end if
			
			    Rxz(i,j) = - 3.d0 *Lx**(-1)*Lz**(-1)* sin(Lx*ax)*sin(Lz*az)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay))+Lz**(-2)*(2.d0-2.d0*cos(Lz*az)))**(-2.5d0)
			
			end if 
		end do 
	end do
	
    !!$omp end parallel do
    
	!! calculate Ryz !!
    
    !!$omp parallel do collapse(2)
    
	do i = 1,N
		do j = 1,N
			if (j /= i) then
			
			ax = (X(i,1)-X(j,1))
			if (ax > 0.5d0 * Lx*2.d0*pi) then 
				ax = ax - Lx*2.d0*pi
			else if (ax <= -0.5d0 * Lx*2.d0*pi) then
				ax = ax + Lx*2.d0*pi
			end if
			
			ay = (X(i+N,1)-X(j+N,1))
			if (ay > 0.5d0 * Ly*2.d0*pi) then 
				ay = ay - Ly*2.d0*pi
			else if (ay <= -0.5d0 * Ly*2.d0*pi) then
				ay = ay + Ly*2.d0*pi
			end if
			
			az = (X(i+2*N,1)-X(j+2*N,1))
			if (az > 0.5d0 * Lz*2.d0*pi) then 
				az = az - Lz*2.d0*pi
			else if (az <= -0.5d0 * Lz*2.d0*pi) then
				az = az + Lz*2.d0*pi
			end if
			
			    Ryz(i,j) = - 3.d0 *Ly**(-1)*Lz**(-1)* sin(Ly*ay)*sin(Lz*az)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay)+Lz**(-2)*(2.d0-2.d0*cos(Lz*az))))**(-2.5d0)
			end if 
		end do 
	end do
    
    !!$omp end parallel do
    
	!! Rxy(i) !!
    
    !!$omp parallel do collapse(2)
    
	do  i = 1,N
    do j = 1,N
        if (j/=i) then
		
		ax = (X(i,1)-X(j,1))
		if (ax > 0.5d0 * Lx*2.d0*pi) then 
			ax = ax - Lx*2.d0*pi
		else if (ax <= -0.5d0 * Lx*2.d0*pi) then
			ax = ax + Lx*2.d0*pi
		end if
			
		ay = (X(i+N,1)-X(j+N,1))
		if (ay > 0.5d0 * Ly*2.d0*pi) then 
			ay = ay - Ly*2.d0*pi
		else if (ay <= -0.5d0 * Ly*2.d0*pi) then
			ay = ay + Ly*2.d0*pi
		end if
			
		az = (X(i+2*N,1)-X(j+2*N,1))
		if (az > 0.5d0 * Lz*2.d0*pi) then 
			az = az - Lz*2.d0*pi
		else if (az <= -0.5d0 * Lz*2.d0*pi) then
			az = az + Lz*2.d0*pi
		end if
		
        H(i,i+N)=H(i,i+N)+3.d0*Lx**(-1.d0)*Ly**(-1.d0)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay))+Lz**(-2)*(2.d0-2.d0*cos(Lz*az)))**(-2.5d0)*sin(Lx*ax)*sin(Ly*ay)
        H(i+N,i)=H(i+N,i)+3.d0*Lx**(-1.d0)*Ly**(-1.d0)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay))+Lz**(-2)*(2.d0-2.d0*cos(Lz*az)))**(-2.5d0)*sin(Lx*ax)*sin(Ly*ay)
        end if 
    end do
	end do 

    !!$omp end parallel do
    
	!! Rxz(i) !!
    
    !!$omp parallel do collapse(2)
    
	do i = 1,N
		do j = 1,N
        if (j/=i) then
		
		ax = (X(i,1)-X(j,1))
		if (ax > 0.5d0 * Lx*2.d0*pi) then 
			ax = ax - Lx*2.d0*pi
		else if (ax <= -0.5d0 * Lx*2.d0*pi) then
			ax = ax + Lx*2.d0*pi
		end if
			
		ay = (X(i+N,1)-X(j+N,1))
		if (ay > 0.5d0 * Ly*2.d0*pi) then 
			ay = ay - Ly*2.d0*pi
		else if (ay <= -0.5d0 * Ly*2.d0*pi) then
			ay = ay + Ly*2.d0*pi
		end if
			
		az = (X(i+2*N,1)-X(j+2*N,1))
		if (az > 0.5d0 * Lz*2.d0*pi) then 
			az = az - Lz*2.d0*pi
		else if (az <= -0.5d0 * Lz*2.d0*pi) then
			az = az + Lz*2.d0*pi
		end if
		
        H(i,i+2*N)=H(i,i+2*N)+3.d0*Lx**(-1.d0)*Lz**(-1.d0)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax)+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay))+Lz**(-2)*(2.d0-2.d0*cos(Lz*az))))**(-2.5d0)*sin(Lx*ax)*sin(Lz*az)
        H(i+2*N,i)=H(i+2*N,i)+3.d0*Lx**(-1.d0)*Lz**(-1.d0)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax)+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay))+Lz**(-2)*(2.d0-2.d0*cos(Lz*az))))**(-2.5d0)*sin(Lx*ax)*sin(Lz*az)
        end if
		end do
	end do

    !!$omp end parallel do
    
	!! Ryz(i) !!
    
    !!$omp parallel do collapse(2)
    
	do  i = 1,N
		do  j = 1,N
        if (j/=i) then
		
		ax = (X(i,1)-X(j,1))
		if (ax > 0.5d0 * Lx*2.d0*pi) then 
			ax = ax - Lx*2.d0*pi
		else if (ax <= -0.5d0 * Lx*2.d0*pi) then
			ax = ax + Lx*2.d0*pi
		end if
			
		ay = (X(i+N,1)-X(j+N,1))
		if (ay > 0.5d0 * Ly*2.d0*pi) then 
			ay = ay - Ly*2.d0*pi
		else if (ay <= -0.5d0 * Ly*2.d0*pi) then
			ay = ay + Ly*2.d0*pi
		end if
			
		az = (X(i+2*N,1)-X(j+2*N,1))
		if (az > 0.5d0 * Lz*2.d0*pi) then 
			az = az - Lz*2.d0*pi
		else if (az <= -0.5d0 * Lz*2.d0*pi) then
			az = az + Lz*2.d0*pi
		end if
		
        H(i+N,i+2*N)=H(i+N,i+2*N)+3*Ly**(-1)*Lz**(-1)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+(Ly**(-2)*(2.d0-2.d0*cos(Ly*ay))+Lz**(-2)*(2.d0-2.d0*cos(Lz*az))))**(-2.5d0)*sin(Ly*ay)*sin(Lz*az)
        H(i+2*N,i+N)=H(i+2*N,i+N)+3*Ly**(-1)*Lz**(-1)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+(Ly**(-2)*(2.d0-2.d0*cos(Ly*ay))+Lz**(-2)*(2.d0-2.d0*cos(Lz*az))))**(-2.5d0)*sin(Ly*ay)*sin(Lz*az)
        end if 
		end do 
	end do
	
    !!$omp end parallel do
    
	!! add matrix together !!
	!! diagonal       !!
	do i = 1 , 3*N 
		do j = 1,3*N
			if (i <= N .and. j <= N) then
				H(i,j)=H(i,j)+Rx(i,j);
			else if (i > N .and. i <= 2*N .and. j > N .and. j <=2*N) then
				H(i,j)=H(i,j)+Ry(i-N,j-N);
			else if (i > 2*N .and.  j >2*N) then
				H(i,j)=H(i,j)+Rz(i-2*N,j-2*N);
			end if
		end do
	end do 
	!! non diagonal !!
	do  i = 1,3*N 
		do j = 1,3*N
			if (i <= N .and. j > N .and. j <= 2*N) then
			H(i,j)=H(i,j)+Rxy(i,j-N);
			H(j,i)=H(j,i)+Rxy(j-N,i);
			else if (i <=N .and. j > 2*N) then
			H(i,j)=H(i,j)+Rxz(i,j-2*N);
			H(j,i)=H(j,i)+Rxz(j-2*N,i);
			else if (i > N .and. i <= 2*N .and. j > 2*N) then
			H(i,j)=H(i,j)+Ryz(i-N,j-2*N);
			H(j,i)=H(j,i)+Ryz(j-2*N,i-N);
			end if 
		end do
	end do
	
	end if 
	
	if (D == 2) then 

	!! calculate the diagonal
	do i=1,2*N
		do j=1,2*N
			if (j /= i .and. j<= N .and. i<= N) then
			
			ax = (X(i,1)-X(j,1))
			if (ax > 0.5d0 * Lx*2.d0*pi) then 
				ax = ax - Lx*2.d0*pi
			else if (ax <= -0.5d0 * Lx*2.d0*pi) then
				ax = ax + Lx*2.d0*pi
			end if
			
			ay = (X(i+N,1)-X(j+N,1))
			if (ay > 0.5d0 * Ly*2.d0*pi) then 
				ay = ay - Ly*2.d0*pi
			else if (ay <= -0.5d0 * Ly*2.d0*pi) then
				ay = ay + Ly*2.d0*pi
			end if
			
			
			H(i,i) = H(i,i) +3*Lx**(-2.d0)*sin(Lx*ax)**2*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay)))**(-2.5d0)-cos(Lx*ax)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay)))**(-1.5d0)
			
			else if (j /= i .and. j> N .and. j<=2*N .and. i > N .and. i<=2*N) then
			
			ax = (X(i-N,1)-X(j-N,1))
			if (ax > 0.5d0 * Lx*2.d0*pi) then 
				ax = ax - Lx*2.d0*pi
			else if (ax <= -0.5d0 * Lx*2.d0*pi) then
				ax = ax + Lx*2.d0*pi
			end if
			
			ay = (X(i,1)-X(j,1))
			if (ay > 0.5d0 * Ly*2.d0*pi) then 
				ay = ay - Ly*2.d0*pi
			else if (ay <= -0.5d0 * Ly*2.d0*pi) then
				ay = ay + Ly*2.d0*pi
			end if
			
			H(i,i) = H(i,i) +3*Ly**(-2.d0)*sin(Ly*ay)**2*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay)))**(-2.5d0)-cos(Ly*ay)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay)))**(-1.5d0)
			
			end if 
		end do
	end do

	!! calculate R(x) !!
	do i = 1,N
		do j = 1,N
			if (j /= i) then
			
			ax = (X(i,1)-X(j,1))
			if (ax > 0.5d0 * Lx*2.d0*pi) then 
				ax = ax - Lx*2.d0*pi
			else if (ax <= -0.5d0 * Lx*2.d0*pi) then
				ax = ax + Lx*2.d0*pi
			end if
			
			ay = (X(i+N,1)-X(j+N,1))
			if (ay > 0.5d0 * Ly*2.d0*pi) then 
				ay = ay - Ly*2.d0*pi
			else if (ay <= -0.5d0 * Ly*2.d0*pi) then
				ay = ay + Ly*2.d0*pi
			end if
			
			
			Rx(i,j) = -3*Lx**(-2)*sin(Lx*ax)**2*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay)))**(-2.5d0)+cos(Lx*ax)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay)))**(-1.5d0)
			end if 
		end do 
	end do
	
	!! calculate R(y) !!
	do i = 1,N
		do j = 1,N
			if (j /= i) then
			
			ax = (X(i,1)-X(j,1))
			if (ax > 0.5d0 * Lx*2.d0*pi) then 
				ax = ax - Lx*2.d0*pi
			else if (ax <= -0.5d0 * Lx*2.d0*pi) then
				ax = ax + Lx*2.d0*pi
			end if
			
			ay = (X(i+N,1)-X(j+N,1))
			if (ay > 0.5d0 * Ly*2.d0*pi) then 
				ay = ay - Ly*2.d0*pi
			else if (ay <= -0.5d0 * Ly*2.d0*pi) then
				ay = ay + Ly*2.d0*pi
			end if
			
			Ry(i,j) = -3*Ly**(-2)*sin(Ly*ay)**2*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay)))**(-2.5d0)+cos(Ly*ay)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay)))**(-1.5d0)
			end if 
		end do 
	end do
	
	!! calculate Rxy !!
	do i = 1,N
		do j = 1,N
			if (j /= i) then
			
			ax = (X(i,1)-X(j,1))
			if (ax > 0.5d0 * Lx*2.d0*pi) then 
				ax = ax - Lx*2.d0*pi
			else if (ax <= -0.5d0 * Lx*2.d0*pi) then
				ax = ax + Lx*2.d0*pi
			end if
			
			ay = (X(i+N,1)-X(j+N,1))
			if (ay > 0.5d0 * Ly*2.d0*pi) then 
				ay = ay - Ly*2.d0*pi
			else if (ay <= -0.5d0 * Ly*2.d0*pi) then
				ay = ay + Ly*2.d0*pi
			end if
			
			    Rxy(i,j) = - 3 *Lx**(-1.d0)*Ly**(-1.d0)* sin(Lx*ax)*sin(Ly*ay)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay)))**(-2.5d0)
			end if 
		end do 
	end do	
	
	!! Rxy(i) !!
	do  i = 1,N
		do j = 1,N
        if (j/=i) then
		
		ax = (X(i,1)-X(j,1))
		if (ax > 0.5d0 * Lx*2.d0*pi) then 
			ax = ax - Lx*2.d0*pi
		else if (ax <= -0.5d0 * Lx*2.d0*pi) then
			ax = ax + Lx*2.d0*pi
		end if
			
		ay = (X(i+N,1)-X(j+N,1))
		if (ay > 0.5d0 * Ly*2.d0*pi) then 
			ay = ay - Ly*2.d0*pi
		else if (ay <= -0.5d0 * Ly*2.d0*pi) then
			ay = ay + Ly*2.d0*pi
		end if
		
        H(i,i+N)=H(i,i+N)+3*Lx**(-1.d0)*Ly**(-1.d0)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay)))**(-2.5d0)*sin(Lx*ax)*sin(Ly*ay)
        H(i+N,i)=H(i+N,i)+3*Lx**(-1.d0)*Ly**(-1.d0)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*ax))+Ly**(-2)*(2.d0-2.d0*cos(Ly*ay)))**(-2.5d0)*sin(Lx*ax)*sin(Ly*ay)
        end if 
		end do
	end do 
	
	!! add matrix together !!
	!! diagonal       !!
	do i = 1 , 2*N 
		do j = 1,2*N
			if (i <= N .and. j <= N) then
				H(i,j)=H(i,j)+Rx(i,j)
			else if (i > N  .and. j > N ) then
				H(i,j)=H(i,j)+Ry(i-N,j-N)
			else if (i > N .and. j <= N) then
				H(i,j)=H(i,j)+Rxy(i-N,j)
				H(j,i)=H(j,i)+Rxy(j,i-N)
			end if 
		end do
	end do

	end if 

	if (D == 1 ) then 
	!! calculate the diagonal
	do i=1,N
		do j=1,N
			if (j /= i ) then
			
			ax = (X(i,1)-X(j,1))
			if (ax > 0.5d0 * Lx*2.d0*pi) then 
				ax = ax - Lx*2.d0*pi
			else if (ax < 0.5d0 * Lx*2.d0*pi) then
				ax = ax + Lx*2.d0*pi
			end if
			
			H(i,i) = H(i,i) +3.d0*Lx**(-2.d0)*sin(Lx*ax)**2*(Lx**(-2.d0)*(2.d0-2.d0*cos(Lx*ax)))**(-2.5d0)-cos(Lx*ax)*(Lx**(-2.d0)*(2.d0-2.d0*cos(Lx*ax)))**(-1.5d0)
			end if 
		end do
	end do
	
	do i = 1,N
		do j = 1,N
			if (j /= i) then
			
			ax = (X(i,1)-X(j,1))
			if (ax > 0.5d0 * Lx*2.d0*pi) then 
				ax = ax - Lx*2.d0*pi
			else if (ax < 0.5d0 * Lx*2.d0*pi) then
				ax = ax + Lx*2.d0*pi
			end if
			
			H(i,j) = H(i,j) -3.d0*Lx**(-2.d0)*sin(Lx*ax)**2*(Lx**(-2.d0)*(2.d0-2.d0*cos(Lx*ax)))**(-2.5d0)+cos(Lx*ax)*(Lx**(-2.d0)*(2.d0-2.d0*cos(Lx*ax)))**(-1.5d0)
			end if 
		end do 
	end do
	end if
	endsubroutine
    
end module