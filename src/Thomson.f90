program Thomson

  implicit none

 ! ----------- !  variable declaration ! ----------- ! 
 
 ! ----   i  ---- ! 
 
 integer                          :: i 
 integer                          :: n_ele
 integer                          :: space 
 integer                          :: nlines
 integer                          :: iter , itermax
  
 ! ----  dp  ---- !
 
 double precision                 :: Lx , Ly , Lz 
 double precision                 :: tol
 double precision                 :: E  
 double precision, allocatable    :: Eign         (:)
 double precision, allocatable    :: geo        (:,:)  
 double precision, allocatable    :: dervative_b(:,:)
 double precision, allocatable    :: dervative_a(:,:)
 double precision, allocatable    :: conj_s     (:,:)
 double precision, allocatable    :: beta       (:,:)
 double precision, allocatable    :: lambda     (:,:)
 double precision, allocatable    :: H          (:,:)
 
 ! ----  ch  ---- !
 
 character (len = 20 )           :: arg , typ , multi 
 character (len = 20 )           :: show , hess , distance , animation
 character (len = 20 )           :: origin 
!-----------------------------------------------------------------------------------!
 
!------------------------------------------------------------------------------------! 
  
  call system('clear')
  if(command_argument_count().eq.0) then
    print *, ""
    write(*,'(a)') "Sorry, you didn't add the input file"
    write(*,'(a)') "----------------------------------------------------"
    write(*,'(a)') "you can run the program using ./Thomson << the name of the input file >>"
    print *, ""
    stop
  else 
    call get_command_argument(1,arg)
  end if 
 
    ! ----- Read from input file  ----- !
  
    call read_f(arg,space,n_ele,tol,itermax,Lx,Ly,Lz,typ,multi,show,hess,distance,animation,origin,nlines)
    
    ! ----- The Title ----- ! 
  
    call title(n_ele,space,itermax,typ,Lx,Ly,Lz)

     
    ! ---- ! memory allocation ! --- !
    
    allocate                   (geo(n_ele,space))
    allocate         (H(n_ele*space,n_ele*space))
    allocate         (dervative_b(n_ele*space,1))
    allocate         (dervative_a(n_ele*space,1))
    allocate              (conj_s(n_ele*space,1))
    allocate                     (beta     (1,1))
    allocate                     (lambda   (1,1))
    allocate                (Eign  (n_ele*space))

    
    ! ---- ! dummy variable ! ---- ! 
    
    iter = 0 
    
    ! ---- ! start the actual code ! ---- ! 
    
    ! ----- Generate first geometry ----- !
    
    call first_geo(arg,n_ele,space,Lx,Ly,Lz,typ,multi,geo,nlines)
    
    write(*,*) ""
    write(*,'(a)') "                                       ((  The starting Geometry  ))  "
    write(*,*) ""
    
    call print_geo(n_ele,space,geo)
    
    call energy(n_ele,space,geo,Lx,Ly,Lz,E)
    
    call diff(n_ele,space,geo,Lx,Ly,Lz,dervative_b)
      
    conj_s(:,1) =  dervative_b(:,1)
    
    ! --------- ! condition if NAN or fall ! -------- ! 
    
    if (isnan(Norm2(dervative_b)) .or. Norm2(dervative_b) > 1e20) then
      write(*,'(a)') "____________________________________________________________________________________________"
      write(*,'(a)') ""
      write (*,'(a)') "you have 2 or more electron at the same position or the energy is infinity"
      Stop
    end if
    
    ! --------- !                        ! ----------- !
    
    write(*,'(a)') "____________________________________________________________________________________________"
    write(*,'(a)') ""
    write(*,'(a,f24.16)') "The energy before optimization = " , E 
    write(*,'(a)') ""
    write(*,'(a,f24.16)') "and the gradiant normalization = " , norm2(dervative_b)
    write(*,'(a)') "____________________________________________________________________________________________"
    write(*,'(a)') ""
    
    ! ----- ! condition for  only one calculation ! ---- !  
    
    if (itermax == 0 ) then 

      call hessian(n_ele,space,geo,Lx,Ly,Lz,H)
     
    end if 
    
    
    if (show /= "show" .and. norm2(dervative_b) > tol) then 
    
      write(*,'(a)') ' Loop   |          Energy           |      Gradient NORM '
      write(*,'(a)') "_______________________________________________________________"
      write(*,'(a)') ""
    end if 
    
    if (animation == "animation") then 
      open(4,file='data_frame.dat',status = 'replace')
      open(8,file='energy.dat'    ,status = 'replace')
    end if 
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%% !
    ! ---- ! the main loop ! ---- ! 
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%% !
    
    do while (iter < itermax ) 
     
    iter = iter + 1 
    
    
    if (norm2(dervative_b) < tol) then 
    
      call energy(n_ele,space,geo,Lx,Ly,Lz,E)
      
      call hessian(n_ele,space,geo,Lx,Ly,Lz,H)
      
      write(*,'(a)') ""
      write(*,'(a,f24.16)') "you are at the minimum, The energy = " , E 
      write(*,'(a)') ""
      write(*,'(a,E20.10,a)') "and the gradiant normalization = " , norm2(dervative_b), "  smaller than the tolerance"
      write(*,'(a)') "____________________________________________________________________________________________"
      write(*,'(a)') ""
      if (animation == "animation") then 
        call anim(n_ele,space,geo,E,iter)
      end if
      exit 
    end if 
    
    ! ---- ! calculate gradiant matrix ! ---- ! 
    
    call diff(n_ele,space,geo,Lx,Ly,Lz,dervative_a)
    
    ! --------------------------------------- !
    
    if (norm2(dervative_a) < tol) then
    
      write(*,'(a)') ""
      write(*,'(a)') "____________________________________________________________________________________________"
      write(*,'(a)') ""
      write(*,'(a,f24.16,a,i5,a)') "Geometry converged at    ",  E ,"   after " ,iter-1  , " Loops"
      write(*,'(a)') ""
      write(*,'(a,E16.10)') "The gradiant norm  =     ", norm2(dervative_a)
      write(*,'(a)') "" 
      exit
    end if 

    ! -------------------------------------- !
    
    ! ---- ! calculate hessian matrix ! ---- ! 
    
    call hessian(n_ele,space,geo,Lx,Ly,Lz,H)
    
    ! -------------------------------------- !
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !
    ! ---- !   conjugated gradiant    ! ---- ! 
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !
    
    beta = matmul(transpose(dervative_a),matmul(H,conj_s))/ matmul(transpose(conj_s),matmul(H,conj_s))

    do i=1,n_ele*space
      conj_s(i,1) =  dervative_a(i,1) + beta(1,1) * conj_s(i,1)
    end do
    
    lambda = matmul(transpose(dervative_a),dervative_a)/matmul(transpose(conj_s),matmul(H,conj_s))
    
    call N_geo(n_ele,space,geo,lambda,conj_s)
    
    call PBC(n_ele,space,geo,Lx,Ly,Lz)
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !
    ! ---- !   conjugated gradiant    ! ---- ! 
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !
    
    call energy(n_ele,space,geo,Lx,Ly,Lz,E)
          
    ! --------------------------------------- !
    
    if (show == "show") then 
    
      write(*,'(a,i7)') "Loop = " , iter
      write(*,'(a)') ""
      write(*,'(a)') "The geometry after the optimization:"
    
      call print_geo(n_ele,space,geo)
    
      write(*,'(a)') "____________________________________________________________________________________________"
      write(*,'(a)') ""
      write(*,'(a,f24.16,a,E8.2)') "The energy after optimization =   ",  E , "     |   The derivative =  ", norm2(dervative_a)
      write(*,'(a)') "____________________________________________________________________________________________"
      write(*,'(a)') ""
    
    else
    
      write(*,'(i7,(1x,f24.16),(3x,f24.16))') iter , E , norm2(dervative_a)
      
    end if
    
    if (animation == "animation") then 
      call anim(n_ele,space,geo,E,iter,origin,Lx,Ly,Lz)
    end if 
    
    end do 
    
    if (show /= "show") then 
    
      write(*,'(a)') ""
      write(*,'(a)') "The final geometry after the optimization:"
      call print_geo(n_ele,space,geo)
            
    end if 
    
    
    if (hess == "hessian") then
    
    write(*,'(a)') "____________________________________________________________________________________________"
    write(*,'(a)') ""
    write(*,'(a)') "                                  The Hessian matrix after the converge "
    write(*,'(a)') ""
    
    do i=1,n_ele*space
      write (*,'(100(3x,f16.8))') H(i,:)
    end do
    
    end if 
    
    if (distance == "distance") then 
      call print_distance(n_ele,space,geo,Lx,Ly,Lz)
    end if
    
    call diagonalize_matrix(n_ele*space,H,Eign)
    
    if (animation == "animation") then 
      close(4)
      close(8)
    end if 
    
    if (animation == "animation") then 
      call plot_anim(n_ele,space,iter,Lx,Ly,Lz)
      call system('gnuplot plot')
      call system('rm -rf plot')
      call system('rm -rf data_frame.dat')
      call system('rm -rf energy.dat')
    end if 
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      call system('rm -rf fort.1')
      call system('rm -rf fort.3')
      call system('rm -rf fort.4')
      call system('rm -rf fort.8')
      call system('rm -rf fort.5')
      call system('rm -rf functions.mod')

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
end program 
