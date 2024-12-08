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
  double precision                           :: hx, hy, hz
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

subroutine get_grad_hessian(N, D, X, Lx, Ly, Lz, der, H)

  implicit none

  integer         , intent(in)  :: N , D  
  double precision, intent(in)  :: X(N,D)
  double precision, intent(in)  :: Lx, Ly, Lz 
  double precision, intent(out) :: der(D,N,1)
  double precision, intent(out) :: H(D,D,N,N)
  
  integer                       :: i, j 
  double precision              :: LLx, LLy, LLz
  double precision              :: LLx_inv, LLy_inv, LLz_inv
  double precision              :: LLxx_inv, LLyy_inv, LLzz_inv
  double precision              :: LLxy_inv, LLxz_inv, LLyz_inv
  double precision              :: LLxx_inv2, LLyy_inv2, LLzz_inv2
  double precision              :: LLxx_inv3, LLyy_inv3, LLzz_inv3
  double precision              :: LLxy_inv3, LLxz_inv3, LLyz_inv3
  double precision              :: ax, ay, az  
  double precision              :: LLxax, LLyay, LLzaz
  double precision              :: dist, dist_inv
  double precision              :: hxx, hyy, hzz, hxy, hxz, hyz
  double precision              :: tmp
  double precision              :: cosx, cosy, cosz
  double precision              :: sinx, siny, sinz
  double precision              :: pi, pi2
  double precision              :: df(D,N,N)
  

  pi = dacos(-1.0d0)
  pi2 = 2.d0 * pi
  
  LLx = pi2 / Lx
  LLy = pi2 / Ly
  LLz = pi2 / Lz

  LLx_inv = 1.d0 / LLx
  LLy_inv = 1.d0 / LLy
  LLz_inv = 1.d0 / LLz

  LLxx_inv = LLx_inv * LLx_inv
  LLyy_inv = LLy_inv * LLy_inv
  LLzz_inv = LLz_inv * LLz_inv


  LLxy_inv = LLx_inv * LLy_inv
  LLxz_inv = LLx_inv * LLz_inv
  LLyz_inv = LLy_inv * LLz_inv

  LLxx_inv2 = 2.d0 * LLxx_inv
  LLyy_inv2 = 2.d0 * LLyy_inv
  LLzz_inv2 = 2.d0 * LLzz_inv

  LLxx_inv3 = 3.d0 * LLxx_inv
  LLyy_inv3 = 3.d0 * LLyy_inv
  LLzz_inv3 = 3.d0 * LLzz_inv
  LLxy_inv3 = 3.d0 * LLxy_inv
  LLxz_inv3 = 3.d0 * LLxz_inv
  LLyz_inv3 = 3.d0 * LLyz_inv
  
  if (D == 3) then

    !$OMP PARALLEL DEFAULT(NONE)                               &
    !$OMP PRIVATE (i, j, cosx, cosy, cosz, sinx, siny, sinz,   &
    !$OMP          LLxax, LLyay, LLzaz, dist, dist_inv, tmp,   &
    !$OMP          hxx, hyy, hzz, hxy, hxz, hyz)               &
    !$OMP SHARED (N, LLx, LLy, LLz, LLx_inv, LLy_inv, LLz_inv, &
    !$OMP         LLxx_inv, LLyy_inv, LLzz_inv,                &
    !$OMP         LLxy_inv, LLxz_inv, LLyz_inv,                &
    !$OMP         LLxx_inv2, LLyy_inv2, LLzz_inv2,             &
    !$OMP         LLxx_inv3, LLyy_inv3, LLzz_inv3,             &
    !$OMP         LLxy_inv3, LLxz_inv3, LLyz_inv3,             &
    !$OMP         X, df, H)
    !$OMP DO
    do j = 1, N

      df(:,:,j) = 0.d0
      H(:,:,:,j) = 0.d0
      do i = j+1, N
  
        LLxax = LLx * (X(i,1) - X(j,1))
        LLyay = LLy * (X(i,2) - X(j,2))
        LLzaz = LLz * (X(i,3) - X(j,3))

        cosx = dcos(LLxax); sinx = dsin(LLxax);
        cosy = dcos(LLyay); siny = dsin(LLyay);
        cosz = dcos(LLzaz); sinz = dsin(LLzaz);

        dist = LLxx_inv2 * (1.d0 - cosx) + LLyy_inv2 * (1.d0 - cosy) + LLzz_inv2 * (1.d0 - cosz)
        dist_inv = 1.d0 / dist
        tmp = dsqrt(dist) * dist_inv * dist_inv

        df(1,i,j) = -tmp * sinx * LLx_inv
        df(2,i,j) = -tmp * siny * LLy_inv
        df(3,i,j) = -tmp * sinz * LLz_inv

        hxx = (3.d0 * sinx * sinx * dist_inv * LLxx_inv - cosx) * tmp
        hyy = (3.d0 * siny * siny * dist_inv * LLyy_inv - cosy) * tmp
        hzz = (3.d0 * sinz * sinz * dist_inv * LLzz_inv - cosz) * tmp

        hxy = 3.d0 * sinx * siny * dist_inv * LLxy_inv * tmp
        hxz = 3.d0 * sinx * sinz * dist_inv * LLxz_inv * tmp
        hyz = 3.d0 * siny * sinz * dist_inv * LLyz_inv * tmp
          
        H(1,1,i,i) = H(1,1,i,i) + 2.d0 * hxx
        H(2,1,i,i) = H(2,1,i,i) + 2.d0 * hxy
        H(3,1,i,i) = H(3,1,i,i) + 2.d0 * hxz
        H(2,2,i,i) = H(2,2,i,i) + 2.d0 * hyy
        H(3,2,i,i) = H(3,2,i,i) + 2.d0 * hyz
        H(3,3,i,i) = H(3,3,i,i) + 2.d0 * hzz

        H(1,1,i,j) = -hxx
        H(2,1,i,j) = -hxy
        H(3,1,i,j) = -hxz
        H(2,2,i,j) = -hyy
        H(3,2,i,j) = -hyz
        H(3,3,i,j) = -hzz
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP PRIVATE(j, i) &
    !$OMP SHARED(N, df, H)
    !$OMP DO
    do j = 1, N-1
      do i = j+1, N
        df(1,j,i) = -df(1,i,j)
        df(2,j,i) = -df(2,i,j)
        df(3,j,i) = -df(3,i,j)

        H(1,1,j,i) = H(1,1,i,j)
        H(2,1,j,i) = H(2,1,i,j)
        H(3,1,j,i) = H(3,1,i,j)
        H(2,2,j,i) = H(2,2,i,j)
        H(3,2,j,i) = H(3,2,i,j)
        H(3,3,j,i) = H(3,3,i,j)
      end do 
    end do 
    !$OMP END DO
    !$OMP END PARALLEL

    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP PRIVATE(i, j) SHARED(N, H)
    !$OMP DO
    do i = 1, N
      do j = 1, N
        H(1,2,j,i) = H(2,1,j,i)
        H(1,3,j,i) = H(3,1,j,i)
        H(2,3,j,i) = H(3,2,j,i)
      end do 
    end do 
    !$OMP END DO
    !$OMP END PARALLEL

    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP PRIVATE(i, j) SHARED(N, der, df)
    !$OMP DO
    do i = 1, N
      der(1:3,i,1) = 0.d0
      do j = 1, N
        der(1:3,i,1) = der(1:3,i,1) - df(1:3,j,i)
      enddo
    end do 
    !$OMP END DO
    !$OMP END PARALLEL
  endif ! D = 3

  
!  if (D == 2) then
!
!    der = 0.d0
!    df  = 0.d0
!    H   = 0.d0
!
!    do i = 1, N-1
!      do j = i+1, N
!
!        ax = X(i,1) - X(j,1)
!        ay = X(i,2) - X(j,2)
!
!        cosx = dcos(LLx * ax)
!        cosy = dcos(LLy * ay)
!
!        dist = tmp_x * (1.d0 - cosx) + tmp_y * (1.d0 - cosy)
!
!        tmp = 1.d0 / (dist * dsqrt(dist))
!
!        sinx = dsin(LLx * ax)
!        siny = dsin(LLy * ay)
!
!        df(i,j)     = -tmp * sinx / LLx 
!        df(i+N,j+N) = -tmp * siny / LLy
!
!        hx = (3.d0 * sinx * sinx / (dist * LLx * LLx) - cosx) * tmp
!        hy = (3.d0 * siny * siny / (dist * LLy * LLy) - cosy) * tmp
!          
!        H(i,i)     = H(i,i)     + 2.d0 * hx
!        H(i+N,i+N) = H(i+N,i+N) + 2.d0 * hy
!
!        H(i,j)     = -hx
!        H(i+N,j+N) = -hy
!
!        hx = 3.d0 * sinx * siny * tmp / (dist * LLx * LLy)
!        H(i,j+N) = -hx
!        H(i,i+N) = H(i,i+N) + 2.d0 * hx
!      end do
!    end do
!
!    do i = 1, N-1
!      do j = i+1, N
!        df(j,i) = -df(i,j)
!        df(j+N,i+N) = -df(i+N,j+N)
!
!        H(j,i) = H(i,j)
!        H(j+N,i+N) = H(i+N,j+N)
!      end do
!    end do
!
!    der(:,1) = -sum(df, 2)
!  endif ! D = 2
!
!  
!  if (D == 1) then
!
!    der = 0.d0
!    df  = 0.d0
!    H   = 0.d0
!
!    do i = 1, N-1
!      do j = i+1, N
!  
!        ax = X(i,1) - X(j,1)
!
!        cosx = dcos(LLx * ax)
!  
!        dist = tmp_x * (1.d0 - cosx)
!  
!        tmp = 1.d0 / (dist * dsqrt(dist))
!  
!        sinx = dsin(LLx * ax)
!
!        df(i,j) = -tmp * sinx / LLx
!  
!        hx = (3.d0 * sinx * sinx / (dist * LLx * LLx) - cosx) * tmp
!          
!        H(i,i) = H(i,i) + 2.d0 * hx
!        H(i,j) = -hx
!      end do
!    end do
!  
!    do i = 1, N-1
!      do j = i+1, N
!        df(j,i) = -df(i,j)
!
!        H(j,i) = H(i,j)
!      end do
!    end do
!
!    der(:,1) = -sum(df, 2)
!  endif ! D = 1
  
end

! ---

  

