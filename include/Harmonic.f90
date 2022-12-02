module Harmonic
	implicit none
	contains
    
    subroutine dista(N,X,R0)
	
	real*16, dimension(:,:), intent(in) :: X
	integer, intent(in) :: N 
	real*16, dimension(N,N) :: R0
	integer :: i,j
	
	do i=1,N
		do j=1,N
			R0(i,j) = sqrt((X(i,1)-X(j,1))**2+(X(i+N,1)-X(j+N,1))**2+(X(i+2*N,1)-X(j+2*N,1))**2)
		end do 
	end do 
	
	end subroutine
	
	subroutine EnergyHP(Kon,N,X,D,Lx,Ly,Lz,V,R0)

	real*16, dimension(:,:), intent(in) :: X
	real*8, dimension(8,8), intent(in) :: Kon
	integer, intent(in) :: N , D
	real*16, intent(in) :: Lx, Ly , Lz
	real*16, dimension(N),intent(inout) :: V
	real*16, dimension(N,N) :: r
	real*16, dimension(N,N),intent(in) :: R0
	integer :: i,j

	V=0.d0
 	r=0.d0
	if (D == 3) then 
			do i=1,N
				do j=1,N
					r(i,j) = ((Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1))))+Lz**(-2)*(2.d0-2.d0*cos(Lz*(X(i+2*N,1)-X(j+2*N,1))))))**(0.5d0)
						if (j > i) then
							V(i)=V(i)+Kon(i,j)/2*(r(i,j)-R0(i,j))**2
						end if 
						!write(*,*) i , j , R0(i,j) , r(i,j)
				end do
			end do
  			
	end if 

	endsubroutine
		
	subroutine diffHP(Kon,N,X,D,Lx,Ly,Lz,tdfHP,R0)
	
	implicit none
	
	real*16, dimension(:,:), intent(in) :: X
	real*8, dimension(8,8), intent(in) :: Kon
	integer, intent(in) :: N ,D 
	real*16 , intent(in) :: Lx,Ly,Lz
	real*8, dimension(3*N,3*N) :: df
	real*8, dimension(N,N) :: r
	real*16, dimension(N,N) :: R0
	real*8 :: C
	real*16 , dimension(3*N),intent(out) :: tdfHP
	integer :: i,j

	tdfHP=0.d0
	df=0.d0
	r=0.d0
	if (D == 3) then 
	
	do i=1, N
			do j=1,N
			r(i,j) = ((Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1))))+Lz**(-2)*(2.d0-2.d0*cos(Lz*(X(i+2*N,1)-X(j+2*N,1))))))**(0.5d0)
			end do 
	end do 
	
	
	do i = 1,N 
		do j = 1, N 
			if (j/=i) then 
			df(i,j) = 2.d0*Kon(i,j)*Lx**(-1.d0)*sin(Lx*(X(i,1)-X(j,1)))*(r(i,j)-R0(i,j))
			end if 
		end do 		
	end do 
	
	do i = N+1,2*N 
		do j = N+1, 2*N 
			if (j/=i) then 
			df(i,j) = 2.d0*Kon(i-N,j-N)*Ly**(-1.d0)*sin(Ly*(X(i,1)-X(j,1)))*(r(i-N,j-N)-R0(i-N,j-N))
			end if 
		end do 		
	end do 
	
		do i = 2*N+1,3*N 
		do j = 2*N+1, 3*N 
			if (j/=i) then 
			df(i,j) = 2.d0*Kon(i-2*N,j-2*N)*Lz**(-1.d0)*sin(Lz*(X(i,1)-X(j,1)))*(r(i-2*N,j-2*N)-R0(i-2*N,j-2*N))
			end if 
		end do 		
		end do
	tdfHP = sum(df,2)
	end if 

	endsubroutine
		
	subroutine hessienHP(Kon,N,X,D,Lx,Ly,Lz,H,R0)
	implicit none 
	real*16, dimension(:,:), intent(in) :: X
	real*8, dimension(8,8), intent(in) :: Kon
	integer, intent(in) ::N,D
	real*16 , intent(in) :: Lx,Ly,Lz
	real*8, dimension(N,N) :: r
	real*16, dimension(N,N) :: R0
	real*8 :: C
	real*8 ,dimension (N,N) :: Rx , Ry ,Rz , Rxy, Rxz , Ryz
	real*8, dimension(3*N,3*N) :: H
	integer :: i,j
		
		r=0.d0
		Rx=0.d0
		Ry=0.d0
		Rz=0.d0
		Rxy=0.d0
		Rxz=0.d0
		Ryz=0.d0
		H=0.d0
		
	
	if (D == 3) then 
	
	do i=1, N
			do j=1,N
			r(i,j) = ((Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1))))+Lz**(-2)*(2.d0-2.d0*cos(Lz*(X(i+2*N,1)-X(j+2*N,1))))))**(0.5d0)
			end do 
	end do
	
	
	do i=1,N
		do j=1,N
			if (j /= i) then
			H(i,i)=H(i,i)+Kon(i,j)*((-2*cos(Lx*(X(i,1)-X(j,1))))*(r(i,j)-R0(i,j))+4*Lx**(-2.d0)*(sin(Lx*(X(i,1)-X(j,1))))**2)
			end if 
		end do
	end do
	
	do i=N+1,2*N
		do j=N+1,2*N
			if (j /= i) then 
			H(i,i)=H(i,i)+Kon(i-N,j-N)*((-2*cos(Ly*(X(i,1)-X(j,1))))*(r(i-N,j-N)-R0(i-N,j-N))+4*Ly**(-2.d0)*(sin(Ly*(X(i,1)-X(j,1))))**2)
			end if 
		end do
	end do
	
	
	do i=2*N+1,3*N
		do j=2*N+1,3*N
			if (j /= i) then 
			H(i,i)=H(i,i)+Kon(i-2*N,j-2*N)*((-2*cos(Lz*(X(i,1)-X(j,1))))*(r(i-2*N,j-2*N)-R0(i-2*N,j-2*N))+4*Lz**(-2.d0)*(sin(Ly*(X(i,1)-X(j,1))))**2)
			end if 
		end do
	end do

	do i = 1,N
		do j = 1,N
			if (j /= i) then 
			Rx(i,j) = -Kon(i,j)*((-2*cos(Lx*(X(i,1)-X(j,1))))*(r(i,j)-R0(i,j))+4*Lx**(-2.d0)*(sin(Lx*(X(i,1)-X(j,1))))**2)
			end if 
		end do 
	end do


	do i = 1,N
		do j = 1,N
			if (j /= i) then
			Ry(i,j) = -Kon(i,j)*((-2*cos(Ly*(X(i+N,1)-X(j+N,1))))*(r(i,j)-R0(i,j))+4*Ly**(-2.d0)*(sin(Ly*(X(i+N,1)-X(j+N,1))))**2)
			end if 
		end do 
	end do

	do i = 1,N
		do j = 1,N
			if (j /= i) then
			Rz(i,j) = -Kon(i,j)*((-2*cos(Lz*(X(i+2*N,1)-X(j+2*N,1))))*(r(i,j)-R0(i,j))+4*Lz**(-2.d0)*(sin(Ly*(X(i+2*N,1)-X(j+2*N,1))))**2)
			end if 
		end do 
	end do


	do i = 1,N
		do j = 1,N
			if (j /= i) then
				Rxy(i,j) = -4*Kon(i,j)*Lx**(-1)*Ly**(-1)*sin(Lx*(X(i,1)-X(j,1)))*sin(Ly*(X(i+N,1)-X(j+N,1)))
			end if 
		end do 
	end do
	
	do i = 1,N
		do j = 1,N
			if (j /= i) then
			Rxz(i,j) = -4*Kon(i,j)*Lx**(-1)*Lz**(-1)*sin(Lx*(X(i,1)-X(j,1)))*sin(Lz*(X(i+2*N,1)-X(j+2*N,1)))
			end if 
		end do 
	end do
	
	do i = 1,N
		do j = 1,N
			if (j /= i) then
				Ryz(i,j) = -4*Kon(i,j)*Ly**(-1)*Lz**(-1)*sin(Ly*(X(i+N,1)-X(j+N,1)))*sin(Lz*(X(i+2*N,1)-X(j+2*N,1)))
			end if 
		end do 
	end do
	
	!! Rxy(i) !!
	do i = 1,N
    do j = 1,N
        if (j/=i) then
		H(i,i+N)=H(i,i+N)+4*Kon(i,j)*Lx**(-1)*Ly**(-1)*sin(Lx*(X(i,1)-X(j,1)))*sin(Ly*(X(i+N,1)-X(j+N,1)))
        H(i+N,i)=H(i+N,i)+4*Kon(i,j)*Lx**(-1)*Ly**(-1)*sin(Lx*(X(i,1)-X(j,1)))*sin(Ly*(X(i+N,1)-X(j+N,1)))
        end if 
    end do
	end do
	
	!! Rxz(i) !!
	do i = 1,N
		do j = 1,N
			if (j/=i) then
		H(i,i+2*N)=H(i,i+2*N)+4*Kon(i,j)*Lx**(-1)*Lz**(-1)*sin(Lx*(X(i,1)-X(j,1)))*sin(Lz*(X(i+2*N,1)-X(j+2*N,1)))
        H(i+2*N,i)=H(i+2*N,i)+4*Kon(i,j)*Lx**(-1)*Lz**(-1)*sin(Lx*(X(i,1)-X(j,1)))*sin(Lz*(X(i+2*N,1)-X(j+2*N,1)))
        end if
		end do
	end do

	!! Ryz(i) !!
	do  i = 1,N
		do  j = 1,N
			if (j/=i) then
		H(i+N,i+2*N)=H(i+N,i+2*N)+4*Kon(i,j)*Ly**(-1)*Lz**(-1)*sin(Ly*(X(i+N,1)-X(j+N,1)))*sin(Lz*(X(i+2*N,1)-X(j+2*N,1)))
        H(i+2*N,i+N)=H(i+2*N,i+N)+4*Kon(i,j)*Ly**(-1)*Lz**(-1)*sin(Ly*(X(i+N,1)-X(j+N,1)))*sin(Lz*(X(i+2*N,1)-X(j+2*N,1)))
        end if 
		end do 
	end do
	
 ! add matrix together!!
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
  ! add  non diagonal !!	

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
	
	endsubroutine
	
	subroutine convertFC(FC,Kon)
	real*8 :: FC(28)
	real*8 :: Kon(8,8)
	integer :: i,j
	
	do i = 1, 8
		do j = 1, 8
			Kon(i,j) = 0 
		end do 
	end do
	
	do i = 2, 8
			Kon(1,i) = FC(i-1)
			Kon(i,1) = FC(i-1)
	end do
	
	do i = 3, 8
		Kon(2,i) = FC(i+5)
		Kon(i,2) = FC(i+5)
	end do 
	
	do i = 4, 8
		Kon(3,i) = FC(i+10)
		Kon(i,3) = FC(i+10)
	end do
	
	do i = 5, 8
		Kon(4,i) = FC(i+14)
		Kon(i,4) = FC(i+14)
	end do
	
	do i = 6, 8
		Kon(5,i) = FC(i+17)
		Kon(i,5) = FC(i+17)
	end do
	
	do i = 7, 8
		Kon(6,i) = FC(i+19)
		Kon(i,6) = FC(i+19)
	end do
	
	kon(7,8) = FC(28)
	kon(8,7) = FC(28)
	
	end subroutine
     
end module