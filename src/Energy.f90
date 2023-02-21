subroutine energy(N,D,X,Lx,Ly,Lz,E)
	
  implicit none 
  
  ! ---- ! input  ! ---- ! 
  
  integer               , intent(in)         :: N , D  
  double precision      , intent(in)         :: Lx , Ly , Lz 
  double precision      , intent(in)         :: X(N,D)
  
  ! ---- ! local  ! ---- !  
  
  integer                                    :: i , j 
  double precision , parameter               :: pi = acos(-1.0d0)
  double precision                           :: LLx , LLy , LLz
  double precision                           :: ax  ,  ay , az  
  
  ! ---- ! output ! ---- !  
   
  double precision     , intent(out)         :: E
  
  
  ! ---- ! code  ! ---- !
	
  E = 0.d0

  LLx = 2.d0*pi/Lx
  LLy = 2.d0*pi/Ly
  LLz = 2.d0*pi/Lz


  if (D == 3) then 
   
  do i = 1,N-1
    do j = i+1,N
        
        ax = abs(X(i,1)-X(j,1))         
        ax = ax - nint(ax/Lx)*Lx
        
        ay = abs(X(i,2)-X(j,2))     
        ay = ay - nint(ay/Ly)*Ly
        
        az = abs(X(i,3)-X(j,3)) 
        az = az - nint(az/Lz)*Lz

		
        E = E + (LLx**(-2)*(2.d0-2.d0*cos(LLx*ax))+LLy**(-2)*(2.d0-2.d0*cos(LLy*ay))+LLz**(-2)*(2.d0-2.d0*cos(LLz*az)))**(-0.5d0)
      
    end do
  end do 
	 
  end if 

  if (D == 2) then 
    
  do i = 1,N-1
    do j= i+1,N
		
		
      ax = abs(X(i,1)-X(j,1))         
      ax = ax - nint(ax/Lx)*Lx
        
      ay = abs(X(i,2)-X(j,2))     
      ay = ay - nint(ay/Ly)*Ly
		
      E = E + (LLx**(-2.d0)*(2.d0-2.d0*cos(LLx*ax))+LLy**(-2.d0)*(2.d0-2.d0*cos(LLy*ay)))**(-0.5d0)
      
    
    end do
  end do 
     
  end if 

  if (D == 1) then 
    
  do i = 1,N-1
    do j=i+1,N
		
		
      ax = abs(X(i,1)-X(j,1))         
      ax = ax - nint(ax/Lx)*Lx
        
      E = E +(LLx**(-2.d0)*(2.d0-2.d0*cos(LLx*ax)))**(-0.5d0)
      
      
    end do
  end do 
    
  end if 
	
  end subroutine
  
  subroutine diff(N,D,X,Lx,Ly,Lz,der)
	
  implicit none
  
  ! ---- ! input  ! ---- ! 
  
  integer               , intent(in)         :: N , D  
  double precision      , intent(in)         :: Lx , Ly , Lz 
  double precision      , intent(in)         :: X(N,D)
  
  ! ---- ! local  ! ---- !  
  
  integer                                    :: i , j 
  double precision , parameter               :: pi = acos(-1.0d0)
  double precision                           :: LLx , LLy , LLz
  double precision                           :: ax  ,  ay , az  
  double precision                           :: dist 
  double precision                           :: df(D*N,D*N)
  
  ! ---- ! output ! ---- !  
   
  double precision     , intent(out)         :: der(D*N,1)
  
  
  ! ---- ! code  ! ---- !
  
  der(:,:) = 0.d0
  df (:,:) = 0.d0
  
  LLx = 2.d0*pi/Lx
  LLy = 2.d0*pi/Ly
  LLz = 2.d0*pi/Lz

	
  if ( D == 3 ) then 

  do i=1,N-1
    do j=i+1,N
			
        ax = X(i,1)-X(j,1)         
        
        ay = X(i,2)-X(j,2)     
        
        az = X(i,3)-X(j,3) 
                
        dist = (LLx**(-2)*(2.d0-2.d0*cos(LLx*ax))+LLy**(-2)*(2.d0-2.d0*cos(LLy*ay))+LLz**(-2)*(2.d0-2.d0*cos(LLz*az)))
        
        df(i,j)         =  - dist**(-1.5d0)*LLx**(-1.d0)*sin(LLx*ax) 
        df(i+N,j+N)     =  - dist**(-1.5d0)*LLy**(-1.d0)*sin(LLy*ay)
        df(i+2*N,j+2*N) =  - dist**(-1.5d0)*LLz**(-1.d0)*sin(LLz*az)
        
    end do
  end do
  
  do i = 1 , 3*N 
    do j = 1 , 3*N 
      df(j,i) = - df(i,j)
    end do 
  end do 
     
  der(:,1) = - sum(df,2)
	
  end if 

  if (D == 2) then 
 
  do i=1,N-1
    do j=i+1,N
			  
      ax = X(i,1)-X(j,1)         

      ay = X(i,2)-X(j,2)     
			
      dist = (LLx**(-2)*(2.d0-2.d0*cos(LLx*ax))+LLy**(-2)*(2.d0-2.d0*cos(LLy*ay)))
      
      df(i,j)     = - dist**(-1.5d0)*LLx**(-1.d0)*sin(LLx*ax)
      df(i+N,j+N) = - dist**(-1.5d0)*LLy**(-1.d0)*sin(LLy*ay)
			
    end do
  end do
  
  do i = 1 ,2*N 
    do j = 1 ,2*N 
      df(j,i) = - df(i,j)
    end do 
  end do 

  der(:,1) = - sum(df,2)
  
  end if 
  
  if (D == 1) then 
 
  do i=1,N-1
    do j=i+1,N
			  
      ax = X(i,1)-X(j,1)
        
      dist = LLx**(-2)*(2.d0-2.d0*cos(LLx*ax))
      
      df(i,j)     = - dist**(-1.5d0)*LLx**(-1.d0)*sin(LLx*ax)
			
    end do
  end do
	
  do i = 1 ,N 
    do j = 1 ,N 
      df(j,i) = - df(i,j)
    end do 
  end do 

  der(:,1) = - sum(df,2)
  
  end if 
  

	
  end subroutine
  
  
  subroutine hessian(N,D,X,Lx,Ly,Lz,H)
	
  implicit none
	
  ! ---- ! input  ! ---- ! 
  
  integer               , intent(in)         :: N , D  
  double precision      , intent(in)         :: Lx , Ly , Lz 
  double precision      , intent(in)         :: X(N,D)
  
  ! ---- ! local  ! ---- !  
  
  integer                                    :: i , j 
  double precision , parameter               :: pi = acos(-1.0d0)
  double precision                           :: LLx , LLy , LLz
  double precision                           :: ax  ,  ay , az  
  double precision                           :: dist 
  ! ---- ! output ! ---- !  
   
  double precision     , intent(out)         :: H(D*N,D*N)
  
  
  ! ---- ! code  ! ---- !
  
    LLx = 2.d0*pi/Lx
    LLy = 2.d0*pi/Ly
    LLz = 2.d0*pi/Lz
  
    H(:,:) = 0.d0
			
  if (D == 3) then

  do i = 1,N
    do j = 1,N
        
        if (j .ne. i) then 
        
        ax = X(i,1)-X(j,1)         
        
        ay = X(i,2)-X(j,2)
        
        az = X(i,3)-X(j,3)
        
        
        dist = (LLx**(-2)*(2.d0-2.d0*cos(LLx*ax))+LLy**(-2)*(2.d0-2.d0*cos(LLy*ay))+LLz**(-2)*(2.d0-2.d0*cos(LLz*az)))
        
        ! ---- ! diagonal element    ! ---- !
        
        H(i,i)         = H(i,i)         +3*LLx**(-2.d0)*sin(LLx*ax)**2*dist**(-2.5d0)-cos(LLx*ax)*dist**(-1.5d0)
        H(i+N,i+N)     = H(i+N,i+N)     +3*LLy**(-2.d0)*sin(LLy*ay)**2*dist**(-2.5d0)-cos(LLy*ay)*dist**(-1.5d0)
        H(i+2*N,i+2*N) = H(i+2*N,i+2*N) +3*LLz**(-2.d0)*sin(LLz*az)**2*dist**(-2.5d0)-cos(LLz*az)*dist**(-1.5d0)

        ! ---- ! Rx and Ry and Rz    ! ---- ! 
        
        H(i,j)         =                -3*LLx**(-2.d0)*sin(LLx*ax)**2*dist**(-2.5d0)+cos(LLx*ax)*dist**(-1.5d0)
        H(i+N,j+N)     =                -3*LLy**(-2.d0)*sin(LLy*ay)**2*dist**(-2.5d0)+cos(LLy*ay)*dist**(-1.5d0)
        H(i+2*N,j+2*N) =                -3*LLz**(-2.d0)*sin(LLz*az)**2*dist**(-2.5d0)+cos(LLz*az)*dist**(-1.5d0)
        
        ! ---- ! Rxy and Rxz and Ryz ! ---- ! 
        
        H(i,j+N)       =                -3*LLx**(-1)*LLy**(-1)*sin(LLx*ax)*sin(LLy*ay)*dist**(-2.5d0)
        H(i,j+2*N)     =                -3*LLx**(-1)*LLz**(-1)*sin(LLx*ax)*sin(LLz*az)*dist**(-2.5d0)
        H(i+N,j+2*N)   =                -3*LLy**(-1)*LLz**(-1)*sin(LLy*ay)*sin(LLz*az)*dist**(-2.5d0)
        
        ! ---- ! Rxy(i) and Rxz(i) and Ryz(i) ! ---- !
          
        H(i,i+N)       = H(i,i+N)       +3*LLx**(-1)*LLy**(-1)*sin(LLx*ax)*sin(LLy*ay)*dist**(-2.5d0)  
        H(i,i+2*N)     = H(i,i+2*N)     +3*LLx**(-1)*LLz**(-1)*sin(LLx*ax)*sin(LLz*az)*dist**(-2.5d0)
        H(i+N,i+2*N)   = H(i+N,i+2*N)   +3*LLy**(-1)*LLz**(-1)*sin(LLy*ay)*sin(LLz*az)*dist**(-2.5d0)  
        
        end if 
        
    end do
  end do


  do i = 1 , 3*N
    do j = 1 , 3*N 
      H(j,i) = H(i,j)
    end do 
  end do 
  
  end if 
  
  if (D == 2) then


  do i = 1,N
    do j = 1,N
        
        if (j .ne. i) then 
        
        ax = X(i,1)-X(j,1)         
        
        ay = X(i,2)-X(j,2)     
        
        
        dist = (LLx**(-2)*(2.d0-2.d0*cos(LLx*ax))+LLy**(-2)*(2.d0-2.d0*cos(LLy*ay)))
        
        ! ---- ! diagonal element    ! ---- !
        
        H(i,i)         = H(i,i)         +3*LLx**(-2.d0)*sin(LLx*ax)**2*dist**(-2.5d0)-cos(LLx*ax)*dist**(-1.5d0)
        H(i+N,i+N)     = H(i+N,i+N)     +3*LLy**(-2.d0)*sin(LLy*ay)**2*dist**(-2.5d0)-cos(LLy*ay)*dist**(-1.5d0)

        ! ---- ! Rx and Ry     ! ---- ! 
        
        H(i,j)         =                -3*LLx**(-2.d0)*sin(LLx*ax)**2*dist**(-2.5d0)+cos(LLx*ax)*dist**(-1.5d0)
        H(i+N,j+N)     =                -3*LLy**(-2.d0)*sin(LLy*ay)**2*dist**(-2.5d0)+cos(LLy*ay)*dist**(-1.5d0)
        
        ! ---- ! Rxy  ! ---- ! 
        
        H(i,j+N)       =                -3*LLx**(-1)*LLy**(-1)*sin(LLx*ax)*sin(LLy*ay)*dist**(-2.5d0)
        
        ! ---- ! Rxy(i) ! ---- !
          
        H(i,i+N)       = H(i,i+N)       +3*LLx**(-1)*LLy**(-1)*sin(LLx*ax)*sin(LLy*ay)*dist**(-2.5d0)    
        
        end if 
        
    end do
  end do


  do i = 1 , 2*N
    do j = 1 , 2*N 
      H(j,i) = H(i,j)
    end do 
  end do 
  
  end if 
  
  if (D == 1) then


  do i = 1,N 
    do j = 1,N
        
        if (j .ne. i) then 
        
        ax = X(i,1)-X(j,1)        
        
        dist = LLx**(-2)*(2.d0-2.d0*cos(LLx*ax))
        
        ! ---- ! diagonal element    ! ---- !
        
        H(i,i)         = H(i,i)         +3*LLx**(-2.d0)*sin(LLx*ax)**2*dist**(-2.5d0)-cos(LLx*ax)*dist**(-1.5d0)

        ! ---- ! Rx   ! ---- ! 
        
        H(i,j)         =                -3*LLx**(-2.d0)*sin(LLx*ax)**2*dist**(-2.5d0)+cos(LLx*ax)*dist**(-1.5d0)
         
        end if 
        
    end do
  end do


  do i = 1 , N
    do j = 1 , N 
      H(j,i) = H(i,j)
    end do 
  end do 
  
  end if 
  
  
  endsubroutine
  