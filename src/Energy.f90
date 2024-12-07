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
  double precision                           :: tmp_x, tmp_y, tmp_z
  
  ! ---- ! output ! ---- !  
   
  double precision     , intent(out)         :: E
  
  
  ! ---- ! code  ! ---- !

  E = 0.d0

  LLx = 2.d0*pi/Lx
  LLy = 2.d0*pi/Ly
  LLz = 2.d0*pi/Lz

  tmp_x = 2.d0 / (LLx * LLx)
  tmp_y = 2.d0 / (LLy * LLy)
  tmp_z = 2.d0 / (LLz * LLz)

  if (D == 3) then 
   
  do i = 1,N-1
    do j = i+1,N
        
        ax = abs(X(i,1)-X(j,1))         
        ax = ax - nint(ax/Lx)*Lx
        
        ay = abs(X(i,2)-X(j,2))     
        ay = ay - nint(ay/Ly)*Ly
        
        az = abs(X(i,3)-X(j,3)) 
        az = az - nint(az/Lz)*Lz

        !E = E + (LLx**(-2)*(2.d0-2.d0*cos(LLx*ax))+LLy**(-2)*(2.d0-2.d0*cos(LLy*ay))+LLz**(-2)*(2.d0-2.d0*cos(LLz*az)))**(-0.5d0)
        E = E + 1.d0 / dsqrt(tmp_x * (1.d0 - dcos(LLx*ax)) &
                           + tmp_y * (1.d0 - dcos(LLy*ay)) &
                           + tmp_z * (1.d0 - dcos(LLz*az)) )
      
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

      !E = E + (LLx**(-2.d0)*(2.d0-2.d0*cos(LLx*ax))+LLy**(-2.d0)*(2.d0-2.d0*cos(LLy*ay)))**(-0.5d0)
      
      E = E + 1.d0 / dsqrt(tmp_x * (1.d0 - dcos(LLx*ax)) &
                         + tmp_y * (1.d0 - dcos(LLy*ay)) )
    
    end do
  end do 
     
  end if 

  if (D == 1) then 
    
  do i = 1,N-1
    do j=i+1,N

      ax = abs(X(i,1)-X(j,1))         
      ax = ax - nint(ax/Lx)*Lx
        
      E = E + 1.d0 / dsqrt(tmp_x * (1.d0 - dcos(LLx*ax)))
      !E = E +(LLx**(-2.d0)*(2.d0-2.d0*cos(LLx*ax)))**(-0.5d0)
      
      
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
  double precision                           :: tmp_x, tmp_y, tmp_z, tmp
  
  ! ---- ! output ! ---- !  
   
  double precision     , intent(out)         :: der(D*N,1)
  
  
  ! ---- ! code  ! ---- !
  
  der(:,:) = 0.d0
  df (:,:) = 0.d0
  
  LLx = 2.d0*pi/Lx
  LLy = 2.d0*pi/Ly
  LLz = 2.d0*pi/Lz

  tmp_x = 2.d0 / (LLx * LLx)
  tmp_y = 2.d0 / (LLy * LLy)
  tmp_z = 2.d0 / (LLz * LLz)
        

  if ( D == 3 ) then 

  do i=1,N-1
    do j=i+1,N

        ax = X(i,1)-X(j,1)
        ay = X(i,2)-X(j,2)
        az = X(i,3)-X(j,3)

        dist = tmp_x * (1.d0 - dcos(LLx*ax)) &
             + tmp_y * (1.d0 - dcos(LLy*ay)) &
             + tmp_z * (1.d0 - dcos(LLz*az))

        tmp = 1.d0 / (dist * dsqrt(dist))
        df(i,j)         = - tmp * dsin(LLx*ax) / LLx 
        df(i+N,j+N)     = - tmp * dsin(LLy*ay) / LLy
        df(i+2*N,j+2*N) = - tmp * dsin(LLz*az) / LLz
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

      !dist = (LLx**(-2)*(2.d0-2.d0*cos(LLx*ax))+LLy**(-2)*(2.d0-2.d0*cos(LLy*ay)))
      !df(i,j)     = - dist**(-1.5d0)*LLx**(-1.d0)*sin(LLx*ax)
      !df(i+N,j+N) = - dist**(-1.5d0)*LLy**(-1.d0)*sin(LLy*ay)

      dist = tmp_x * (1.d0 - dcos(LLx*ax)) &
           + tmp_y * (1.d0 - dcos(LLy*ay))
      tmp = 1.d0 / (dist * dsqrt(dist))
      df(i,j)     = - tmp * dsin(LLx*ax) / LLx 
      df(i+N,j+N) = - tmp * dsin(LLy*ay) / LLy
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
        
      !dist = LLx**(-2)*(2.d0-2.d0*cos(LLx*ax))
      !df(i,j) = - dist**(-1.5d0)*LLx**(-1.d0)*sin(LLx*ax)

      dist = tmp_x * (1.d0 - dcos(LLx*ax))
      tmp = 1.d0 / (dist * dsqrt(dist))
      df(i,j) = - tmp * dsin(LLx*ax) / LLx 
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
  
! ---

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
  double precision                           :: tmp_x, tmp_y, tmp_z, tmp
  double precision                           :: cosx, cosy, cosz
  double precision                           :: sinx, siny, sinz

  ! ---- ! output ! ---- !  
   
  double precision     , intent(out)         :: H(D*N,D*N)
  
  
  ! ---- ! code  ! ---- !
  
  LLx = 2.d0*pi/Lx
  LLy = 2.d0*pi/Ly
  LLz = 2.d0*pi/Lz

  tmp_x = 2.d0 / (LLx * LLx)
  tmp_y = 2.d0 / (LLy * LLy)
  tmp_z = 2.d0 / (LLz * LLz)
  
  H(:,:) = 0.d0

  if (D == 3) then

    do j = 1, N
      do i = 1, N
  
        if (j .eq. i) cycle
  
        ax = X(i,1) - X(j,1)         
        ay = X(i,2) - X(j,2)
        az = X(i,3) - X(j,3)

        cosx = dcos(LLx * ax)
        cosy = dcos(LLy * ay)
        cosz = dcos(LLz * az)
  
        dist = tmp_x * (1.d0 - cosx) + tmp_y * (1.d0 - cosy) + tmp_z * (1.d0 - cosz)
  
        tmp = 1.d0 / (dist * dsqrt(dist))
  
        sinx = dsin(LLx * ax)
        siny = dsin(LLy * ay)
        sinz = dsin(LLz * az)
  
        hx = (3.d0 * sinx * sinx / (dist * LLx * LLx) - cosx) * tmp
        hy = (3.d0 * siny * siny / (dist * LLy * LLy) - cosy) * tmp
        hz = (3.d0 * sinz * sinz / (dist * LLz * LLz) - cosz) * tmp
          
        H(i,i)         = H(i,i)         + hx
        H(i+N,i+N)     = H(i+N,i+N)     + hy
        H(i+2*N,i+2*N) = H(i+2*N,i+2*N) + hz
  
        H(i,j)         = -hx
        H(i+N,j+N)     = -hy
        H(i+2*N,j+2*N) = -hz
  
        hx = 3.d0 * sinx * siny * tmp / (dist * LLx * LLy)
        hy = 3.d0 * sinx * sinz * tmp / (dist * LLx * LLz)
        hz = 3.d0 * siny * sinz * tmp / (dist * LLy * LLz)
  
        H(i,j+N)     = -hx
        H(i,j+2*N)   = -hy
        H(i+N,j+2*N) = -hz
  
        H(i,i+N)     = H(i,i+N)     + hx
        H(i,i+2*N)   = H(i,i+2*N)   + hy
        H(i+N,i+2*N) = H(i+N,i+2*N) + hz
      end do
    end do
  
  
    ! TODO
    do i = 1 , 3*N
      do j = 1 , 3*N 
        H(j,i) = H(i,j)
      end do 
    end do 
  
  end if 

  
  if (D == 2) then
    do j = 1, N
      do i = 1, N

        if (j .ne. i) cycle

        ax = X(i,1) - X(j,1)         
        ay = X(i,2) - X(j,2)

        cosx = dcos(LLx * ax)
        cosy = dcos(LLy * ay)

        dist = tmp_x * (1.d0 - cosx) + tmp_y * (1.d0 - cosy)

        tmp = 1.d0 / (dist * dsqrt(dist))

        sinx = dsin(LLx * ax)
        siny = dsin(LLy * ay)

        hx = (3.d0 * sinx * sinx / (dist * LLx * LLx) - cosx) * tmp
        hy = (3.d0 * siny * siny / (dist * LLy * LLy) - cosy) * tmp
          
        H(i,i)         = H(i,i)         + hx
        H(i+N,i+N)     = H(i+N,i+N)     + hy

        H(i,j)         = -hx
        H(i+N,j+N)     = -hy

        hx = 3.d0 * sinx * siny * tmp / (dist * LLx * LLy)
        H(i,j+N) = -hx
        H(i,i+N) = H(i,i+N)     + hx
      end do
    end do

    do i = 1, 2*N
      do j = 1, 2*N
        H(j,i) = H(i,j)
      end do
    end do

  endif 

  
  if (D == 1) then

    do j = 1, N
      do i = 1, N 

        if (j .ne. i) cycle
  
        ax = X(i,1) - X(j,1)         

        cosx = dcos(LLx * ax)
  
        dist = tmp_x * (1.d0 - cosx)
  
        tmp = 1.d0 / (dist * dsqrt(dist))
  
        sinx = dsin(LLx * ax)
  
        hx = (3.d0 * sinx * sinx / (dist * LLx * LLx) - cosx) * tmp
          
        H(i,i) = H(i,i) + hx
        H(i,j) = -hx
      end do
    end do
  
    ! TODO
    do i = 1, N
      do j = 1, N
        H(j,i) = H(i,j)
      end do
    end do
  
  endif 
  
  
end subroutine

! ---

call get_grad_hessian(N, D, X, Lx, Ly, Lz, der, H)

  implicit none

  integer         , intent(in)  :: N , D  
  double precision, intent(in)  :: X(N,D)
  double precision, intent(in)  :: Lx, Ly, Lz 
  double precision, intent(out) :: der(D*N,1)
  double precision, intent(out) :: H(D*N,D*N)
  
  integer                       :: i , j 
  double precision , parameter  :: pi = dacos(-1.0d0)
  double precision              :: LLx , LLy , LLz
  double precision              :: ax  ,  ay , az  
  double precision              :: dist 
  double precision              :: tmp_x, tmp_y, tmp_z, tmp
  double precision              :: cosx, cosy, cosz
  double precision              :: sinx, siny, sinz
  double precision              :: df(D*N,D*N)
  
  
  LLx = 2.d0*pi/Lx
  LLy = 2.d0*pi/Ly
  LLz = 2.d0*pi/Lz

  tmp_x = 2.d0 / (LLx * LLx)
  tmp_y = 2.d0 / (LLy * LLy)
  tmp_z = 2.d0 / (LLz * LLz)
  
  der(:,:) = 0.d0
  df(:,:) = 0.d0
  H(:,:) = 0.d0

  if (D == 3) then

    do i = 1, N-1
      do j = i+1, N
  
        ax = X(i,1) - X(j,1)
        ay = X(i,2) - X(j,2)
        az = X(i,3) - X(j,3)

        cosx = dcos(LLx * ax); sinx = dsin(LLx * ax)
        cosy = dcos(LLy * ay); siny = dsin(LLy * ay)
        cosz = dcos(LLz * az); sinz = dsin(LLz * az)

        dist = tmp_x * (1.d0 - cosx) + tmp_y * (1.d0 - cosy) + tmp_z * (1.d0 - cosz)
        tmp = 1.d0 / (dist * dsqrt(dist))

        df(i,j)         = -tmp * sinx / LLx
        df(i+N,j+N)     = -tmp * siny / LLy
        df(i+2*N,j+2*N) = -tmp * sinz / LLz

        hx = (3.d0 * sinx * sinx / (dist * LLx * LLx) - cosx) * tmp
        hy = (3.d0 * siny * siny / (dist * LLy * LLy) - cosy) * tmp
        hz = (3.d0 * sinz * sinz / (dist * LLz * LLz) - cosz) * tmp
          
        H(i,i)         = H(i,i)         + 2.d0 * hx
        H(i+N,i+N)     = H(i+N,i+N)     + 2.d0 * hy
        H(i+2*N,i+2*N) = H(i+2*N,i+2*N) + 2.d0 * hz
  
        H(i,j)         = -hx
        H(i+N,j+N)     = -hy
        H(i+2*N,j+2*N) = -hz
  
        hx = 3.d0 * sinx * siny * tmp / (dist * LLx * LLy)
        hy = 3.d0 * sinx * sinz * tmp / (dist * LLx * LLz)
        hz = 3.d0 * siny * sinz * tmp / (dist * LLy * LLz)
  
        H(i,j+N)     = -hx
        H(i,j+2*N)   = -hy
        H(i+N,j+2*N) = -hz
  
        H(i,i+N)     = H(i,i+N)     + 2.d0 * hx
        H(i,i+2*N)   = H(i,i+2*N)   + 2.d0 * hy
        H(i+N,i+2*N) = H(i+N,i+2*N) + 2.d0 * hz
      end do
    end do

    do i = 1, N-1
      do j = i+1, N
        df(j,i) = -df(i,j)
        df(j+N,i+N) = -df(i+N,j+N)
        df(j+2*N,i+2*N) = -df(i+2*N,j+2*N)

        H(j,i) = H(i,j)
        H(j+N,i+N) = H(i+N,j+N)
        H(j+2*N,i+2*N) = H(i+2*N,j+2*N)
      end do 
    end do 

    der(:,1) = -sum(df, 2)
  endif ! D = 1

  
  if (D == 2) then

    do i = 1, N-1
      do j = i+1, N

        ax = X(i,1) - X(j,1)
        ay = X(i,2) - X(j,2)

        cosx = dcos(LLx * ax)
        cosy = dcos(LLy * ay)

        dist = tmp_x * (1.d0 - cosx) + tmp_y * (1.d0 - cosy)

        tmp = 1.d0 / (dist * dsqrt(dist))

        sinx = dsin(LLx * ax)
        siny = dsin(LLy * ay)

        df(i,j)     = -tmp * sinx / LLx 
        df(i+N,j+N) = -tmp * siny / LLy

        hx = (3.d0 * sinx * sinx / (dist * LLx * LLx) - cosx) * tmp
        hy = (3.d0 * siny * siny / (dist * LLy * LLy) - cosy) * tmp
          
        H(i,i)     = H(i,i)     + 2.d0 * hx
        H(i+N,i+N) = H(i+N,i+N) + 2.d0 * hy

        H(i,j)     = -hx
        H(i+N,j+N) = -hy

        hx = 3.d0 * sinx * siny * tmp / (dist * LLx * LLy)
        H(i,j+N) = -hx
        H(i,i+N) = H(i,i+N) + 2.d0 * hx
      end do
    end do

    do i = 1, N-1
      do j = i+1, N
        df(j,i) = -df(i,j)
        df(j+N,i+N) = -df(i+N,j+N)

        H(j,i) = H(i,j)
        H(j+N,i+N) = H(i+N,j+N)
      end do
    end do

    der(:,1) = -sum(df, 2)
  endif ! D = 2

  
  if (D == 1) then

    do i = 1, N-1
      do j = i+1, N
  
        ax = X(i,1) - X(j,1)

        cosx = dcos(LLx * ax)
  
        dist = tmp_x * (1.d0 - cosx)
  
        tmp = 1.d0 / (dist * dsqrt(dist))
  
        sinx = dsin(LLx * ax)

        df(i,j) = -tmp * sinx / LLx
  
        hx = (3.d0 * sinx * sinx / (dist * LLx * LLx) - cosx) * tmp
          
        H(i,i) = H(i,i) + 2.d0 * hx
        H(i,j) = -hx
      end do
    end do
  
    do i = 1, N-1
      do j = i+1, N
        df(j,i) = -df(i,j)

        H(j,i) = H(i,j)
      end do
    end do

    der(:,1) = -sum(df, 2)
  endif ! D = 1
  
end

! ---

  

