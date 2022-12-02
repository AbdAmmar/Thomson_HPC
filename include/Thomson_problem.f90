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
	
	E=0.d0

	if (D == 3) then 
    !!$omp parallel do collapse(2)
	do i = 1,N
		do j=1,N
		if (j > i) then
			E(i)=E(i)+ (Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1))))+Lz**(-2)*(2.d0-2.d0*cos(Lz*(X(i+2*N,1)-X(j+2*N,1)))))**(-0.5d0)
		end if 
		end do
	end do 
	!!$omp end parallel do 
	end if 

	if (D == 2) then 
    !!$omp parallel do collapse(2)
	do i = 1,N
		do j=1,N
		if (j > i) then
			E(i)=E(i)+ (Lx**(-2.d0)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1))))+Ly**(-2.d0)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1)))))**(-0.5d0)
		end if 
		end do
	end do 
    !!$omp end parallel do 
	end if 

	if (D == 1) then 
    !!$omp parallel do collapse(2)
	do i = 1,N
		do j=1,N
		if (j > i) then
			E(i)=E(i)+(Lx**(-2.d0)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1)))))**(-0.5d0)
		end if 
		end do
	end do 
    !!$omp end parallel do 
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
	
	
	df=0.d0
	tdf=0.d0 
		
	if ( D == 3 ) then 
    !!$omp parallel do collapse(2)
	do i=1,N
		do j=1,N
			if (j /= i) then
			df(i,j) = - (Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1))))+Lz**(-2)*(2.d0-2.d0*cos(Lz*(X(i+2*N,1)-X(j+2*N,1)))))**(-1.5d0)*Lx**(-1.d0)*sin(Lx*(X(i,1)-X(j,1)))
			end if 
		end do
	end do
	!!$omp end parallel do
    
    !!$omp parallel do collapse(2)
	do i=N+1,2*N
		do j=N+1,2*N
			if (j /= i) then
			df(i,j) =  - (Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i-N,1)-X(j-N,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i,1)-X(j,1))))+Lz**(-2)*(2.d0-2.d0*cos(Lz*(X(i+N,1)-X(j+N,1)))))**(-1.5d0)*Ly**(-1.d0)*sin(Ly*(X(i,1)-X(j,1)))
			end if 
		end do
	end do
	!!$omp end parallel do
    
    !!$omp parallel do collapse(2)
	do i=2*N+1,3*N
		do j=2*N+1,3*N
			if (j /= i) then
			df(i,j) =  - (Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i-2*N,1)-X(j-2*N,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i-N,1)-X(j-N,1))))+Lz**(-2)*(2.d0-2.d0*cos(Lz*(X(i,1)-X(j,1)))))**(-1.5d0)*Lz**(-1.d0)*sin(Lz*(X(i,1)-X(j,1)))
			end if 
		end do
	end do
    !!$omp end parallel do
	tdf = sum(df,2)
	end if 

	if (D == 2) then 

	do i=1,N
		do j=1,N
			if (j /= i) then
			df(i,j) = -Lx**(-1)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1)))))**(-1.5d0)*sin(Lx*(X(i,1)-X(j,1)))
			end if 
		end do
	end do
	
	do i=N+1,2*N
		do j=N+1,2*N
			if (j /= i) then
			df(i,j) = -Ly**(-1)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i-N,1)-X(j-N,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i,1)-X(j,1)))))**(-1.5d0)*sin(Ly*(X(i,1)-X(j,1)))
			end if 
		end do
	end do
	tdf = sum(df,2)
	end if 

	if (D == 1) then 

	do i=1,N
		do j=1,N
			if (j /= i) then
			df(i,j) = -Lx**(-1.d0)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1)))))**(-1.5d0)*sin(Lx*(X(i,1)-X(j,1)))
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

			Rx=0.d0
			Ry=0.d0
			Rz=0.d0
			Rxy=0.d0
			Rxz=0.d0
			Ryz=0.d0
			H=0.d0
			
	if (D == 3) then

	!! calculate the diagonal
    
    !!$omp parallel do collapse(2)

	do i=1,3*N
		do j=1,3*N
			if (j /= i .and. j<= N .and. i<=N) then
			H(i,i) = H(i,i) +3*Lx**(-2.d0)*sin(Lx*(X(i,1)-X(j,1)))**2*(Lx**(-2)*(2.d0-2.d0*cos(X(i,1)-X(j,1)))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1))))+Lz**(-2)*(2.d0-2.d0*cos(Lz*(X(i+2*N,1)-X(j+2*N,1)))))**(-2.5d0)-cos(Lx*(X(i,1)-X(j,1)))*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1))))+Lz**(-2)*(2.d0-2.d0*cos(Lz*(X(i+2*N,1)-X(j+2*N,1)))))**(-1.5d0)
			else if (j /= i .and. j> N .and. j<=2*N .and. i > N .and. i<=2*N) then
			H(i,i) = H(i,i) +3*Ly**(-2.d0)*sin(Ly*(X(i,1)-X(j,1)))**2*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i- N,1)-X(j- N,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i,1)-X(j,1))))+Lz**(-2)*(2.d0-2.d0*cos(Lz*(X(i+ N,1)-X(j+ N,1)))))**(-2.5d0)-cos(Ly*(X(i,1)-X(j,1)))*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i- N,1)-X(j- N,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i,1)-X(j,1))))+Lz**(-2)*(2-2.d0*cos(Lz*(X(i+ N,1)-X(j+ N,1)))))**(-1.5d0)
			else if (j /= i .and. j> 2*N .and. i > 2*N ) then
			H(i,i) = H(i,i) +3*Lz**(-2.d0)*sin(Lz*(X(i,1)-X(j,1)))**2*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i-2*N,1)-X(j-2*N,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i-N,1)-X(j-N,1))))+Lz**(-2)*(2.d0-2.d0*cos(Lz*(X(i,1)-X(j,1)))))**(-2.5d0)-cos(Lz*(X(i,1)-X(j,1)))*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i-2*N,1)-X(j-2*N,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i-N,1)-X(j-N,1))))+Lz**(-2)*(2.d0-2.d0*cos(Lz*(X(i,1)-X(j,1)))))**(-1.5d0)
			end if 
		end do
	end do
    
    !!$omp end parallel do
    
	!! calculate R(x) !!

    !!$omp parallel do collapse(2)

	do i = 1,N
		do j = 1,N
			if (j /= i) then
			Rx(i,j) = -3.d0*Lx**(-2.d0)*sin(Lx*(X(i,1)-X(j,1)))**2*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1))))+Lz**(-2)*(2.d0-2.d0*cos(Lz*(X(i+2*N,1)-X(j+2*N,1)))))**(-2.5d0)+cos(Lx*(X(i,1)-X(j,1)))*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1))))+Lz**(-2)*(2.d0-2.d0*cos(Lz*(X(i+2*N,1)-X(j+2*N,1)))))**(-1.5d0)
			end if 
		end do 
	end do
	
    !!$omp end parallel do
    
    
    
	!! calculate R(y) !!
    
    !!$omp parallel do collapse(2)
    
	do i = 1,N
		do j = 1,N
			if (j /= i) then
			Ry(i,j) = -3.d0*Ly**(-2.d0)*sin(Ly*(X(i+N,1)-X(j+N,1)))**2*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1))))+Lz**(-2)*(2.d0-2.d0*cos(Lz*(X(i+2*N,1)-X(j+2*N,1)))))**(-2.5d0)+cos(Ly*(X(i+N,1)-X(j+N,1)))*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1))))+Lz**(-2)*(2.d0-2.d0*cos(Lz*(X(i+2*N,1)-X(j+2*N,1)))))**(-1.5d0)
			end if 
		end do 
	end do
    
    !!$omp end parallel do

	!! calculate R(z) !!
    
    !!$omp parallel do collapse(2)
    
	do i = 1,N
		do j = 1,N
			if (j /= i) then
			Rz(i,j) = -3.d0*Lz**(-2.d0)*sin(Lz*(X(i+2*N,1)-X(j+2*N,1)))**2*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1))))+Lz**(-2)*(2.d0-2.d0*cos(Lz*(X(i+2*N,1)-X(j+2*N,1)))))**(-2.5d0)+cos(Lz*(X(i+2*N,1)-X(j+2*N,1)))*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1))))+Lz**(-2)*(2.d0-2.d0*cos(Lz*(X(i+2*N,1)-X(j+2*N,1)))))**(-1.5d0)
			end if 
		end do 
	end do
	
    !!$omp end parallel do
    
	!! calculate Rxy !!
    
    !!$omp parallel do collapse(2)
    
	do i = 1,N
		do j = 1,N
			if (j /= i) then
			    Rxy(i,j) = - 3.d0 *Lx**(-1)*Ly**(-1)* sin(Lx*(X(i,1)-X(j,1)))*sin(Ly*(X(i+N,1)-X(j+N,1)))*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1))))+Lz**(-2)*(2.d0-2.d0*cos(Lz*(X(i+2*N,1)-X(j+2*N,1)))))**(-2.5d0)
			end if 
		end do 
	end do	
	
    
    !!$omp end parallel do
    
	!! calculate Rxz !!
    
    !!$omp parallel do collapse(2)
    
	do i = 1,N
		do j = 1,N
			if (j /= i) then
			     Rxz(i,j) = - 3.d0 *Lx**(-1)*Lz**(-1)* sin(Lx*(X(i,1)-X(j,1)))*sin(Lz*(X(i+2*N,1)-X(j+2*N,1)))*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1))))+Lz**(-2)*(2.d0-2.d0*cos(Lz*(X(i+2*N,1)-X(j+2*N,1)))))**(-2.5d0)
			end if 
		end do 
	end do
	
    !!$omp end parallel do
    
	!! calculate Ryz !!
    
    !!$omp parallel do collapse(2)
    
	do i = 1,N
		do j = 1,N
			if (j /= i) then
			        Ryz(i,j) = - 3.d0 *Ly**(-1)*Lz**(-1)* sin(Ly*(X(i+N,1)-X(j+N,1)))*sin(Lz*(X(i+2*N,1)-X(j+2*N,1)))*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1)))+Lz**(-2)*(2.d0-2.d0*cos(Lz*(X(i+2*N,1)-X(j+2*N,1))))))**(-2.5d0)
			end if 
		end do 
	end do
    
    !!$omp end parallel do
    
	!! Rxy(i) !!
    
    !!$omp parallel do collapse(2)
    
	do  i = 1,N
    do j = 1,N
        if (j/=i) then
        H(i,i+N)=H(i,i+N)+3.d0*Lx**(-1.d0)*Ly**(-1.d0)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1))))+Lz**(-2)*(2.d0-2.d0*cos(Lz*(X(i+2*N,1)-X(j+2*N,1)))))**(-2.5d0)*sin(Lx*(X(i,1)-X(j,1)))*sin(Ly*(X(i+N,1)-X(j+N,1)))
        H(i+N,i)=H(i+N,i)+3*Lx**(-1)*Ly**(-1)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1))))+Lz**(-2)*(2.d0-2.d0*cos(Lz*(X(i+2*N,1)-X(j+2*N,1)))))**(-2.5d0)*sin(Lx*(X(i,1)-X(j,1)))*sin(Ly*(X(i+N,1)-X(j+N,1)))
        end if 
    end do
	end do 

    !!$omp end parallel do
    
	!! Rxz(i) !!
    
    !!$omp parallel do collapse(2)
    
	do i = 1,N
		do j = 1,N
        if (j/=i) then
        H(i,i+2*N)=H(i,i+2*N)+3.d0*Lx**(-1.d0)*Lz**(-1.d0)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1)))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1))))+Lz**(-2)*(2.d0-2.d0*cos(Lz*(X(i+2*N,1)-X(j+2*N,1))))))**(-2.5d0)*sin(Lx*(X(i,1)-X(j,1)))*sin(Lz*(X(i+2*N,1)-X(j+2*N,1)))
        H(i+2*N,i)=H(i+2*N,i)+3.d0*Lx**(-1.d0)*Lz**(-1.d0)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1)))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1))))+Lz**(-2)*(2.d0-2.d0*cos(Lz*(X(i+2*N,1)-X(j+2*N,1))))))**(-2.5d0)*sin(Lx*(X(i,1)-X(j,1)))*sin(Lz*(X(i+2*N,1)-X(j+2*N,1)))
        end if
		end do
	end do

    !!$omp end parallel do
    
	!! Ryz(i) !!
    
    !!$omp parallel do collapse(2)
    
	do  i = 1,N
		do  j = 1,N
        if (j/=i) then
        H(i+N,i+2*N)=H(i+N,i+2*N)+3*Ly**(-1)*Lz**(-1)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1))))+(Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1))))+Lz**(-2)*(2.d0-2.d0*cos(Lz*(X(i+2*N,1)-X(j+2*N,1))))))**(-2.5d0)*sin(Ly*(X(i+N,1)-X(j+N,1)))*sin(Lz*(X(i+2*N,1)-X(j+2*N,1)))
        H(i+2*N,i+N)=H(i+2*N,i+N)+3*Ly**(-1)*Lz**(-1)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1))))+(Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1))))+Lz**(-2)*(2.d0-2.d0*cos(Lz*(X(i+2*N,1)-X(j+2*N,1))))))**(-2.5d0)*sin(Ly*(X(i+N,1)-X(j+N,1)))*sin(Lz*(X(i+2*N,1)-X(j+2*N,1)))
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
			
			H(i,i) = H(i,i) +3*Lx**(-2.d0)*sin(Lx*(X(i,1)-X(j,1)))**2*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1)))))**(-2.5d0)-cos(Lx*(X(i,1)-X(j,1)))*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1)))))**(-1.5d0)
			
			else if (j /= i .and. j> N .and. j<=2*N .and. i > N .and. i<=2*N) then
			
			H(i,i) = H(i,i) +3*Ly**(-2.d0)*sin(Ly*(X(i,1)-X(j,1)))**2*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i- N,1)-X(j- N,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i,1)-X(j,1)))))**(-2.5d0)-cos(Ly*(X(i,1)-X(j,1)))*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i-N,1)-X(j-N,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i,1)-X(j,1)))))**(-1.5d0)
			
			end if 
		end do
	end do

	!! calculate R(x) !!
	do i = 1,N
		do j = 1,N
			if (j /= i) then
			Rx(i,j) = -3*Lx**(-2)*sin(Lx*(X(i,1)-X(j,1)))**2*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1)))))**(-2.5d0)+cos(Lx*(X(i,1)-X(j,1)))*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1)))))**(-1.5d0)
			end if 
		end do 
	end do
	
	!! calculate R(y) !!
	do i = 1,N
		do j = 1,N
			if (j /= i) then
			Ry(i,j) = -3*Ly**(-2)*sin(Ly*(X(i+N,1)-X(j+N,1)))**2*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1)))))**(-2.5d0)+cos(Ly*(X(i+N,1)-X(j+N,1)))*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1)))))**(-1.5d0)
			end if 
		end do 
	end do
	
	!! calculate Rxy !!
	do i = 1,N
		do j = 1,N
			if (j /= i) then
			    Rxy(i,j) = - 3 *Lx**(-1.d0)*Ly**(-1.d0)* sin(Lx*(X(i,1)-X(j,1)))*sin(Ly*(X(i+N,1)-X(j+N,1)))*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1)))))**(-2.5d0)
			end if 
		end do 
	end do	
	
	!! Rxy(i) !!
	do  i = 1,N
		do j = 1,N
        if (j/=i) then
        H(i,i+N)=H(i,i+N)+3*Lx**(-1.d0)*Ly**(-1.d0)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1)))))**(-2.5d0)*sin(Lx*(X(i,1)-X(j,1)))*sin(Ly*(X(i+N,1)-X(j+N,1)))
        H(i+N,i)=H(i+N,i)+3*Lx**(-1)*Ly**(-1)*(Lx**(-2)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1))))+Ly**(-2)*(2.d0-2.d0*cos(Ly*(X(i+N,1)-X(j+N,1)))))**(-2.5d0)*sin(Lx*(X(i,1)-X(j,1)))*sin(Ly*(X(i+N,1)-X(j+N,1)))
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
			H(i,i) = H(i,i) +3.d0*Lx**(-2)*sin(Lx*(X(i,1)-X(j,1)))**2*(Lx**(-2.d0)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1)))))**(-2.5d0)-cos(Lx*(X(i,1)-X(j,1)))*(Lx**(-2.d0)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1)))))**(-1.5d0)
			end if 
		end do
	end do
	
	do i = 1,N
		do j = 1,N
			if (j /= i) then
			    H(i,j) = H(i,j) -3.d0*Lx**(-2.d0)*sin(Lx*(X(i,1)-X(j,1)))**2*(Lx**(-2.d0)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1)))))**(-2.5d0)+cos(Lx*(X(i,1)-X(j,1)))*(Lx**(-2.d0)*(2.d0-2.d0*cos(Lx*(X(i,1)-X(j,1)))))**(-1.5d0)
			end if 
		end do 
	end do
	end if
	endsubroutine
    
end module