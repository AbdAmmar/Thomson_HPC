program Thomson

  implicit none

#ifdef USE_MPI
    include 'mpif.h'
#endif
 ! ----------- !  variable declaration ! ----------- !
 
 ! ----   i  ---- ! 
 
 integer                          :: i 
 integer                          :: n_ele
 integer                          :: space 
 integer                          :: nlines
 integer                          :: iter , itermax
  
 ! ----  dp  ---- !
 
 double precision                 :: Lx, Ly, Lz 
 double precision                 :: tol
 double precision                 :: E  
 double precision                 :: beta(1,1), lambda(1,1)
 double precision, allocatable    :: Eign(:)
 double precision, allocatable    :: Energy_out(:)
 double precision, allocatable    :: geo(:,:)  
 double precision, allocatable    :: dervative_b(:,:)
 double precision, allocatable    :: dervative_a(:,:)
 double precision, allocatable    :: conj_s(:,:)
 double precision, allocatable    :: H(:,:)
 double precision, allocatable    :: Hconj_s(:,:)

 ! ----  ch  ---- !
 
 character (len = 20 )           :: arg , typ , multi 
 character (len = 20 )           :: show , hess , distance , animation
 character (len = 20 )           :: origin 
!-----------------------------------------------------------------------------------!
! MPI ! 

 integer                         :: ierr
 integer                         :: nprocs         ! total number of processes
 integer                         :: ME             ! CPU index
 character(len=2)                :: Ncore

#ifdef USE_MPI
    call MPI_Init(ierr)

    call MPI_Comm_size(MPI_COMM_WORLD,nprocs,ierr) ! Nombre of procs dispo

    call MPI_Comm_rank(MPI_COMM_WORLD,ME,ierr)     ! Id de chaque proc
#else
    ierr = 0
    nprocs = 1        ! Assume single process execution
    ME = 0            ! Assume single process execution
#endif
  
!-----------------------------------------------------------------------------------!  
  
  if (ME == 0) then
  
  call system('clear')
  if(command_argument_count().eq.0) then
    write(*,'(a)')
    write(*,'(a)') "Sorry, you didn't add the input file"
    write(*,'(a)') "----------------------------------------------------"
    write(*,'(a)') "you can run the program using ./Thomson << the name of the input file >>"
    write(*,'(a)')
    stop
  else 
    call get_command_argument(1,arg)
  end if 
 
    ! ----- Read from input file  ----- !
  
    call read_f(arg,space,n_ele,tol,itermax,Lx,Ly,Lz,typ,multi,show,hess,distance,animation,origin,nlines)
    
  end if 
  
#ifdef USE_MPI
  call MPI_Bcast(space  ,1,mpi_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(n_ele  ,1,mpi_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(itermax,1,mpi_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(nlines ,1,mpi_INTEGER,0,MPI_COMM_WORLD,ierr)
  
  call MPI_Bcast(tol    ,1,mpi_double_precision,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(Lx     ,1,mpi_double_precision,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(Ly     ,1,mpi_double_precision,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(Lz     ,1,mpi_double_precision,0,MPI_COMM_WORLD,ierr)
  
  call MPI_Bcast(arg      ,LEN(arg),mpi_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(typ      ,LEN(arg),mpi_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(multi    ,LEN(arg),mpi_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(show     ,LEN(arg),mpi_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(hess     ,LEN(arg),mpi_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(distance ,LEN(arg),mpi_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(animation,LEN(arg),mpi_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(origin   ,LEN(arg),mpi_CHARACTER,0,MPI_COMM_WORLD,ierr) 
#endif
    ! ---- ! memory allocation ! --- !
    
    allocate(Hconj_s(n_ele*space,1))
    allocate(geo(n_ele,space))
    allocate(H(n_ele*space,n_ele*space))
    allocate(dervative_b(n_ele*space,1))
    allocate(dervative_a(n_ele*space,1))
    allocate(conj_s(n_ele*space,1))
    allocate(Eign(n_ele*space))

    
    ! ---- ! dummy variable ! ---- ! 
    
    iter = 0 
    
    ! ---- ! start the actual code ! ---- ! 
    
    ! ----- Generate first geometry ----- !
    
    call first_geo(arg,n_ele,space,Lx,Ly,Lz,typ,multi,geo,nlines)

#ifdef USE_MPI
    write(Ncore,'(I2)') ME
    open(10,file='rand_'//Ncore//'.out')
#else
    open(10,file='rand.out')
#endif

    ! ----- The Title ----- !
    
    call title(n_ele,space,itermax,typ,Lx,Ly,Lz,10)
    
    write(10,*) ""
    write(10,'(a)') "                                       ((  The starting Geometry  ))  "
    write(10,*) ""
   
    call print_geo(n_ele,space,geo,10)
    
    call energy(n_ele,space,geo,Lx,Ly,Lz,E)
    
    call diff(n_ele,space,geo,Lx,Ly,Lz,dervative_b)
      
    conj_s(:,1) =  dervative_b(:,1)
    
    ! --------- ! condition if NAN or fall ! -------- ! 
    
    if (isnan(Norm2(dervative_b)) .or. Norm2(dervative_b) > 1e20) then
      write(10,'(a)') "____________________________________________________________________________________________"
      write(10,'(a)') ""
      write(10,'(a)') "you have 2 or more electron at the same position or the energy is infinity"
      Stop
    end if
    
    ! --------- !                        ! ----------- !
    
    write(10,'(a)') "____________________________________________________________________________________________"
    write(10,'(a)') ""
    write(10,'(a,f24.16)') "The energy before optimization = " , E 
    write(10,'(a)') ""
    write(10,'(a,f24.16)') "and the gradiant normalization = " , norm2(dervative_b)
    write(10,'(a)') "____________________________________________________________________________________________"
    write(10,'(a)') ""
    
    ! ----- ! condition for  only one calculation ! ---- !  
    
    if (itermax == 0 ) then 

      call hessian(n_ele,space,geo,Lx,Ly,Lz,H)
     
    end if 
    
    
    if (show /= "show" .and. norm2(dervative_b) > tol) then 
    
      write(10,'(a)') ' Loop   |          Energy           |      Gradient NORM '
      write(10,'(a)') "_______________________________________________________________"
      write(10,'(a)') ""
      
    end if 
    
    if (animation == "animation") then 
#ifdef USE_MPI
      open(4,file='data_frame'//Ncore//'.dat',status = 'replace')
      open(8,file='energy'//Ncore//'.dat'    ,status = 'replace')
#else
      open(4,file='data_frame.dat',status = 'replace')
      open(8,file='energy.dat'    ,status = 'replace')
#endif
    end if
    
    
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%% !
  ! ---- ! the main loop ! ---- ! 
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%% !
    
  do while (iter < itermax)

    iter = iter + 1 

    if (norm2(dervative_b) < tol) then 
      call energy(n_ele, space, geo, Lx, Ly, Lz, E)
      call hessian(n_ele, space, geo, Lx, Ly, Lz, H)
      write(10,'(a)') ""
      write(10,'(a,f24.16)') "you are at the minimum, The energy = ", E
      write(10,'(a)') ""
      write(10,'(a,E20.10,a)') "and the gradiant normalization = ", norm2(dervative_b), "  smaller than the tolerance"
      write(10,'(a)') "____________________________________________________________________________________________"
      write(10,'(a)') ""
      if (animation == "animation") then 
        call anim(n_ele, space, geo, E, iter, origin, Lx, Ly, Lz)
      end if
      exit 
    end if 
    
    call get_grad_hessian(n_ele, space, geo(1,1), Lx, Ly, Lz, &
                          dervative_a(1,1), H(1,1))
    
    if (norm2(dervative_a) < tol) then
      write(10,'(a)') ""
      write(10,'(a)') "____________________________________________________________________________________________"
      write(10,'(a)') ""
      write(10,'(a,f24.16,a,i5,a)') "Geometry converged at    ",  E ,"   after " ,iter-1  , " Loops"
      write(10,'(a)') ""
      write(10,'(a,E17.10)') "The gradiant norm  =     ", norm2(dervative_a)
      write(10,'(a)') "" 
      exit
    end if 

    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !
    ! ---- !   conjugated gradiant    ! ---- ! 
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !
    
    call dgemv("N", n_ele*space, n_ele*space, 1.d0, &
               H(1,1), n_ele*space, conj_s(1,1), 1, &
               0.d0, Hconj_s(1,1), 1)

    beta = matmul(transpose(dervative_a), Hconj_s) &
         / matmul(transpose(conj_s), Hconj_s)
    !beta = matmul(transpose(dervative_a), matmul(H, conj_s)) / matmul(transpose(conj_s), matmul(H, conj_s))

    do i = 1, n_ele*space
      conj_s(i,1) = dervative_a(i,1) + beta(1,1) * conj_s(i,1)
    end do
    
    !lambda = matmul(transpose(dervative_a), dervative_a) / matmul(transpose(conj_s), matmul(H, conj_s))
    lambda = matmul(transpose(dervative_a), dervative_a) &
           / matmul(transpose(conj_s), Hconj_s)

    ! TODO
    call N_geo(n_ele, space, geo, lambda, conj_s)
    call PBC(n_ele, space, geo, Lx, Ly, Lz)
    
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !
    ! ---- !   conjugated gradiant    ! ---- ! 
    ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% !

    call energy(n_ele, space, geo, Lx, Ly, Lz, E)

    ! --------------------------------------- !
    
    if (show == "show") then 
    
      write(10,'(a,i7)') "Loop = " , iter
      write(10,'(a)') ""
      write(10,'(a)') "The geometry after the optimization:"
   
      call print_geo(n_ele,space,geo,10)
   
      write(10,'(a)') "____________________________________________________________________________________________"
      write(10,'(a)') ""
      write(10,'(a,f24.16,a,E9.2)') "The energy after optimization =   ",  E , "     |   The derivative =  ", norm2(dervative_a)
      write(10,'(a)') "____________________________________________________________________________________________"
      write(10,'(a)') ""
    
    else
    
      write(10,'(i7,(1x,f24.16),(3x,f24.16))') iter , E , norm2(dervative_a)
      
    end if
    
    if (animation == "animation") then
      call anim(n_ele,space,geo,E,iter,origin,Lx,Ly,Lz)
    end if
    
  end do ! while (iter < itermax)

    
  if (show /= "show") then 
    write(10,'(a)') ""
    write(10,'(a)') "The final geometry after the optimization:"
    call print_geo(n_ele,space,geo,10)
  end if 
    
    
  if (hess == "hessian") then
    write(10,'(a)') "____________________________________________________________________________________________"
    write(10,'(a)') ""
    write(10,'(a)') "                                  The Hessian matrix after the converge "
    write(10,'(a)') ""
    call matout(n_ele*space,n_ele*space,H,10)
  end if 
    
  if (distance == "distance") then 
    call print_distance(n_ele,space,geo,Lx,Ly,Lz,10)
  end if
    
  call diagonalize_matrix(n_ele*space,H,Eign,10)
    
  if (animation == "animation") then 
    close(4)
    close(8)
  end if 
    
  if (animation == "animation") then 
    open (9,file='plot',status = 'replace')
    call plot_anim(n_ele,space,iter,Lx,Ly,Lz,10)
    close(9)
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
  
  close(10)
  
  
  allocate(Energy_out(nprocs))
  
#ifdef USE_MPI
  call MPI_Gather(E, 1, MPI_DOUBLE_PRECISION, Energy_out, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif

#ifdef USE_MPI
  if (ME == 0) then
    open(3,file='energy.dat')
        do i = 1,size(Energy_out)
            write(3,'(I5,f20.12)') i , Energy_out(i)
        end do 
    close(3)
  end if
#endif

#ifdef USE_MPI
  call MPI_Finalize(ierr)
#endif


  deallocate(Hconj_s, geo, H, dervative_b, dervative_a, conj_s, Eign)

end program
