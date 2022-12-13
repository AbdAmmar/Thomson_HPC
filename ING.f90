module functions
implicit none 
contains

subroutine noanswer
character (len=200) :: answer 
integer :: d , M , i , du
real*8  :: r 
real*16, parameter :: pi = ACOS(-1.0d0)
real*16 :: LLx, LLy , LLz
real, allocatable :: x(:),y(:),z(:)

do while (.true.)
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Dimensions number ? (write a number 1,2 or 3)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) d
            if (d <= 3) then
                exit
            else
                write(*,'(a)') "----------------------------------------------"    
                write(*,'(a)') "**** N.P:       Choose number between (1:3)"
            cycle
            end if 
        end do

        write(1,'(a,I1)') "dimension: " , d
        du = d 
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Number of electrons ? (write any number)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) d
        write(1,'(a,I5)') "electron: " , d
        M = d
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Tolerance ? (write any number like this 1e-5)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) r
        write(1,'(a,E10.2)') "tolerance: " , r

        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "The maximum number of iteration ? (write any number)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) d
        write(1,'(a,I10)') "itermax: " , d


        do while (.true.)
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "The size of the box ? (default : 2*pi, 2*pi , 2*pi)"
        write(*,'(a)') "Do you want default size ? (Y/N) "
        write(*,'(a)') "----------------------------------------------"
        read(*,*) answer
            if (answer == "y" .or. answer == "Y") then
                LLx = 2*pi; LLy = 2*pi ;LLz = 2*pi 
                exit
            elseif (answer == "n" .or. answer == "N") then 
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "enter the size of the box (write it as :  *****  ***** *****)"
                write(*,'(a)') "----------------------------------------------"
                read*, LLx , LLy , LLz
                exit
            else
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "your answer is not (Y/N) , kindly try again"
                cycle
            end if 
        end do 

        write(1,'(a,f20.16,f20.16,f20.16)') 'box:' , LLx , LLy , LLz


        write(1,'(a)') "periodic"

        do while (.true.)
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Random geometry for the electron ? (Y/N)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) answer
            if (answer == "y" .or. answer == "Y") then
            write(1,'(a)') "random"
                exit
            elseif (answer == "n" .or. answer == "N") then 
            write(*,'(a)') "----------------------------------------------"
            write(*,'(a)') "*****    you have to write your geometry after the keyword in the input.inp file"
                exit
            else
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "your answer is not (Y/N) , kindly try again"
                cycle
            end if 
        end do

        do while (.true.)
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Multiply the geometry by the size of the box ? (Y/N)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) answer
            if (answer == "y" .or. answer == "Y") then  
            write(1,'(a)') "multiply"
                exit
            elseif (answer == "n" .or. answer == "N") then
                exit
            else
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "your answer is not (Y/N) , kindly try again"
                cycle
            end if 
        end do

        do while (.true.)
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Show all the results (with the geometry in every steps) ? (Y/N)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) answer
            if (answer == "y" .or. answer == "Y") then
            write(1,'(a)') "show"
                exit
            elseif (answer == "n" .or. answer == "N") then 
                exit
            else
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "your answer is not (Y/N) , kindly try again"
                cycle
            end if 
        end do

        do while (.true.)
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Generate video animation ? (Y/N)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) answer
            if (answer == "y" .or. answer == "Y") then
            write(1,'(a)') "animation"
                exit
            elseif (answer == "n" .or. answer == "N") then 
                exit
            else
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "your answer is not (Y/N) , kindly try again"
                cycle
            end if 
        end do


        do while (.true.)
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Show the hessian matrix ? (Y/N)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) answer
            if (answer == "y" .or. answer == "Y") then
            write(1,'(a)') "hessian"
                exit
            elseif (answer == "n" .or. answer == "N") then 
                exit
            else
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "your answer is not (Y/N) , kindly try again"
                cycle
            end if 
        end do
        
        write(1,'(a)') "geometry"
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        allocate(x(M))
        allocate(y(M))
        allocate(z(M))
        
        
        call random_number(x)
        call random_number(y)
        call random_number(z)
        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (du == 3) then 
            do i = 1 , M
                write(1,'(I5,f16.8,f16.8,f16.8)') i , x(i) , y(i) , z(i) 
            end do 
        else if (du == 2) then 
            do i = 1 , M
                write(1,'(I5,f16.8,f16.8)') i , x(i) , y(i)
            end do
        else 
            do i = 1 , M
                write(1,'(I5,f16.8)') i , x(i) 
            end do
        end if 
        write(*,'(a)') "----------------------------------------------------------"
        write(*,'(a)') "***    Now you have your input file with coordinations"
        write(*,'(a)') "***    Please open the input file and modifiy your geometry"
        write(*,'(a)') "----------------------------------------------------------"
end subroutine

subroutine yesanswer
implicit none 

character (len=200) :: answer 
character (len=200) :: FCC , BCC
real*8  :: r , dx, dy
real*16, parameter :: pi = ACOS(-1.0d0)
real*16 :: LLx, LLy , LLz
integer :: i, j, k, n, icont , o
real*16 :: d, d2 , r1 , r2 
real*8, allocatable, dimension(:) :: x, y
integer :: ix, iy, ipoint



do while (.true.)
write(*,'(a)') "----------------------------------------------"
write(*,'(a)') "Ok, Please write one of this option (FCC, BCC or SC in 3D) (HEX , SL in 2D) (ED in 1D)"
write(*,'(a)') "----------------------------------------------"
read(*,*) answer
select case(answer)
    case DEFAULT
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "*****    Kindly review your typing and try again"
    cycle
    
    case("FCC","fcc") 
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Number of unit-cell per dimension ? (write a number) 4*n^3"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) n 
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a,I5)') "The Total number of electron in this case is " , 4*n**3
        
        write(1,'(a)') "dimension: 3"
        write(1,'(a,I5)') "electron: " , 4*n**3
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Tolerance ? (write any number like this 1e-5)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) r
        write(1,'(a,E10.2)') "tolerance: " , r

        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "The maximum number of iteration ? (write any number)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) o
        write(1,'(a,I10)') "itermax: " , o


        do while (.true.)
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "The size of the box ? (default : 2*pi, 2*pi , 2*pi)"
        write(*,'(a)') "Do you want default size ? (Y/N) "
        write(*,'(a)') "----------------------------------------------"
        read(*,*) answer
            if (answer == "y" .or. answer == "Y") then
                LLx = 2*pi; LLy = 2*pi ;LLz = 2*pi 
                exit
            elseif (answer == "n" .or. answer == "N") then 
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "enter the size of the box (write it as :  *****  ***** *****)"
                write(*,'(a)') "----------------------------------------------"
                read*, LLx , LLy , LLz
                exit
            else
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "your answer is not (Y/N) , kindly try again"
                cycle
            end if 
        end do 

        write(1,'(a,f20.16,f20.16,f20.16)') 'box:' , LLx , LLy , LLz


        write(1,'(a)') "periodic"
        write(1,'(a)') "multiply"
        
        
        do while (.true.)
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Show all the results (with the geometry in every steps) ? (Y/N)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) answer
            if (answer == "y" .or. answer == "Y") then
            write(1,'(a)') "show"
                exit
            elseif (answer == "n" .or. answer == "N") then 
                exit
            else
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "your answer is not (Y/N) , kindly try again"
                cycle
            end if 
        end do

        do while (.true.)
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Generate video animation ? (Y/N)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) answer
            if (answer == "y" .or. answer == "Y") then
            write(1,'(a)') "animation"
                exit
            elseif (answer == "n" .or. answer == "N") then 
                exit
            else
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "your answer is not (Y/N) , kindly try again"
                cycle
            end if 
        end do


        do while (.true.)
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Show the hessian matrix ? (Y/N)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) answer
            if (answer == "y" .or. answer == "Y") then
            write(1,'(a)') "hessian"
                exit
            elseif (answer == "n" .or. answer == "N") then 
                exit
            else
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "your answer is not (Y/N) , kindly try again"
                cycle
            end if 
        end do
        
        write(1,'(a)') "geometry"

        d=1.d0/n
        d2=d/2.d0
        icont=0
        
            do i = 0 , n-1 
                do j = 0 , n-1
                    do k = 0 , n-1
                                icont=icont+1
                                                write(1,'(i5,3(2x,f20.18),a)') icont, d*i  ,d*j , d*k 
                    end do		
                end do 
            end do

            do i = 0 , n-1 
                do j = 0 , n-1
                    do k = 0 , n-1
                                icont=icont+1
                                                write(1,'(i5,3(2x,f20.18),a)') icont, d*i+d2  ,d*j+d2 , d*k 
                    end do		
                end do 
            end do


            do i = 0 , n-1 
                do j = 0 , n-1
                    do k = 0 , n-1
                                icont=icont+1
                                                write(1,'(i5,3(2x,f20.18),a)') icont, d*i  ,d*j+d2 , d*k+d2  
                    end do		
                end do 
            end do

            do i = 0 , n-1 
                do j = 0 , n-1
                    do k = 0 , n-1
                                icont=icont+1
                                                write(1,'(i5,3(2x,f20.18),a)') icont, d*i+d2  ,d*j , d*k+d2  
                    end do		
                end do 
            end do
    exit 
    
    case("BCC","bcc")
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Number of unit-cell per dimension ? (write a number) 2*n^3"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) n 
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a,I5)') "The Total number of electron in this case is " , 2*n**3
   
    
        write(1,'(a)') "dimension: 3"
        write(1,'(a,I5)') "electron: " , 2*n**3
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Tolerance ? (write any number like this 1e-5)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) r
        write(1,'(a,E10.2)') "tolerance: " , r

        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "The maximum number of iteration ? (write any number)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) o
        write(1,'(a,I10)') "itermax: " , o


        do while (.true.)
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "The size of the box ? (default : 2*pi, 2*pi , 2*pi)"
        write(*,'(a)') "Do you want default size ? (Y/N) "
        write(*,'(a)') "----------------------------------------------"
        read(*,*) answer
            if (answer == "y" .or. answer == "Y") then
                LLx = 2*pi; LLy = 2*pi ;LLz = 2*pi 
                exit
            elseif (answer == "n" .or. answer == "N") then 
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "enter the size of the box (write it as :  *****  ***** *****)"
                write(*,'(a)') "----------------------------------------------"
                read*, LLx , LLy , LLz
                exit
            else
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "your answer is not (Y/N) , kindly try again"
                cycle
            end if 
        end do 

        write(1,'(a,f20.16,f20.16,f20.16)') 'box:' , LLx , LLy , LLz


        write(1,'(a)') "periodic"
        write(1,'(a)') "multiply"
        
        
        do while (.true.)
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Show all the results (with the geometry in every steps) ? (Y/N)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) answer
            if (answer == "y" .or. answer == "Y") then
            write(1,'(a)') "show"
                exit
            elseif (answer == "n" .or. answer == "N") then 
                exit
            else
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "your answer is not (Y/N) , kindly try again"
                cycle
            end if 
        end do

        do while (.true.)
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Generate video animation ? (Y/N)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) answer
            if (answer == "y" .or. answer == "Y") then
            write(1,'(a)') "animation"
                exit
            elseif (answer == "n" .or. answer == "N") then 
                exit
            else
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "your answer is not (Y/N) , kindly try again"
                cycle
            end if 
        end do


        do while (.true.)
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Show the hessian matrix ? (Y/N)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) answer
            if (answer == "y" .or. answer == "Y") then
            write(1,'(a)') "hessian"
                exit
            elseif (answer == "n" .or. answer == "N") then 
                exit
            else
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "your answer is not (Y/N) , kindly try again"
                cycle
            end if 
        end do
        
        
        
        write(1,'(a)') "geometry"
        
          r1=real(1,8)
          r2=real(2,8)
          d=r1/real(n,8)
          d2=d/r2
          icont=0
            
          do i=0,n-1
             do j=0,n-1
                do k=0,n-1
                   icont=icont+1
                   write (1,'(i5,3(2x,f20.18),a)') icont, d*i, d*j, d*k 
                enddo
             enddo
          enddo
          do i=0,n-1
             do j=0,n-1
                do k=0,n-1
                   icont=icont+1
                   write (1,'(i5,3(2x,f20.18),a)') icont, d*i+d2, d*j+d2, d*k+d2 
                enddo
             enddo
          enddo
          stop
        
         exit
    
        case("SC", "sc")
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Number of unit-cell per dimension ? (write a number) n^3"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) n 
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a,I5)') "The Total number of electron in this case is " , n**3
        write(1,'(a)') "dimension: 3"
        write(1,'(a,I5)') "electron: " , n**3
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Tolerance ? (write any number like this 1e-5)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) r
        write(1,'(a,E10.2)') "tolerance: " , r

        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "The maximum number of iteration ? (write any number)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) o
        write(1,'(a,I10)') "itermax: " , o


        do while (.true.)
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "The size of the box ? (default : 2*pi, 2*pi , 2*pi)"
        write(*,'(a)') "Do you want default size ? (Y/N) "
        write(*,'(a)') "----------------------------------------------"
        read(*,*) answer
            if (answer == "y" .or. answer == "Y") then
                LLx = 2*pi; LLy = 2*pi ;LLz = 2*pi 
                exit
            elseif (answer == "n" .or. answer == "N") then 
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "enter the size of the box (write it as :  *****  ***** *****)"
                write(*,'(a)') "----------------------------------------------"
                read*, LLx , LLy , LLz
                exit
            else
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "your answer is not (Y/N) , kindly try again"
                cycle
            end if 
        end do 

        write(1,'(a,f20.16,f20.16,f20.16)') 'box:' , LLx , LLy , LLz


        write(1,'(a)') "periodic"
        write(1,'(a)') "multiply"
        
        
        do while (.true.)
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Show all the results (with the geometry in every steps) ? (Y/N)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) answer
            if (answer == "y" .or. answer == "Y") then
            write(1,'(a)') "show"
                exit
            elseif (answer == "n" .or. answer == "N") then 
                exit
            else
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "your answer is not (Y/N) , kindly try again"
                cycle
            end if 
        end do

        do while (.true.)
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Generate video animation ? (Y/N)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) answer
            if (answer == "y" .or. answer == "Y") then
            write(1,'(a)') "animation"
                exit
            elseif (answer == "n" .or. answer == "N") then 
                exit
            else
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "your answer is not (Y/N) , kindly try again"
                cycle
            end if 
        end do


        do while (.true.)
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Show the hessian matrix ? (Y/N)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) answer
            if (answer == "y" .or. answer == "Y") then
            write(1,'(a)') "hessian"
                exit
            elseif (answer == "n" .or. answer == "N") then 
                exit
            else
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "your answer is not (Y/N) , kindly try again"
                cycle
            end if 
        end do
        
        write(1,'(a)') "geometry"
        
            d=1.d0/n
            icont=0

            do i = 0 , n-1
                do j = 0 , n-1
                    do k = 0 , n-1
                                icont=icont+1
                                                write(1,'(i3,3(2x,f16.12),a)') icont, d*i  ,d*j , d*k 
                    end do		
                end do 
            end do
    exit
    
    
    case("SL","sl")
    
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Number of electron per dimension ? (write a number) n^2"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) n 
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a,I5)') "The Total number of electron in this case is " , n**2
        write(1,'(a)') "dimension: 2"
        write(1,'(a,I5)') "electron: " , n**2
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Tolerance ? (write any number like this 1e-5)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) r
        write(1,'(a,E10.2)') "tolerance: " , r

        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "The maximum number of iteration ? (write any number)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) o
        write(1,'(a,I10)') "itermax: " , o


        do while (.true.)
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "The size of the box ? (default : 2*pi, 2*pi)"
        write(*,'(a)') "Do you want default size ? (Y/N) "
        write(*,'(a)') "----------------------------------------------"
        read(*,*) answer
            if (answer == "y" .or. answer == "Y") then
                LLx = 2*pi; LLy = 2*pi ;LLz = 2*pi 
                exit
            elseif (answer == "n" .or. answer == "N") then 
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "enter the size of the box (write it as :  *****  ***** )"
                write(*,'(a)') "----------------------------------------------"
                read*, LLx , LLy , LLz
                exit
            else
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "your answer is not (Y/N) , kindly try again"
                cycle
            end if 
        end do 

        write(1,'(a,f20.16,f20.16,f20.16)') 'box:' , LLx , LLy , LLz


        write(1,'(a)') "periodic"
        write(1,'(a)') "multiply"
        
        
        do while (.true.)
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Show all the results (with the geometry in every steps) ? (Y/N)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) answer
            if (answer == "y" .or. answer == "Y") then
            write(1,'(a)') "show"
                exit
            elseif (answer == "n" .or. answer == "N") then 
                exit
            else
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "your answer is not (Y/N) , kindly try again"
                cycle
            end if 
        end do

        do while (.true.)
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Generate video animation ? (Y/N)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) answer
            if (answer == "y" .or. answer == "Y") then
            write(1,'(a)') "animation"
                exit
            elseif (answer == "n" .or. answer == "N") then 
                exit
            else
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "your answer is not (Y/N) , kindly try again"
                cycle
            end if 
        end do


        do while (.true.)
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Show the hessian matrix ? (Y/N)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) answer
            if (answer == "y" .or. answer == "Y") then
            write(1,'(a)') "hessian"
                exit
            elseif (answer == "n" .or. answer == "N") then 
                exit
            else
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "your answer is not (Y/N) , kindly try again"
                cycle
            end if 
        end do
        
        write(1,'(a)') "geometry"
        
        d=1.d0/n
        icont= 0 
        do i=0,n-1
        do j=0,n-1
            icont=icont+1
                write (1,'(i3,2(2x,f16.12),a)') icont, d*i , d*j 
        end do 
        enddo
        exit    
    
    case("HEX", "hex")
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Number of electron on one dimension? (write a number) 4*n^2"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) n 
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a,I5)') "The Total number of electron in this case is " , 4*n**2
        write(1,'(a)') "dimension: 2"
        write(1,'(a,I5)') "electron: " , 4*n**2
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Tolerance ? (write any number like this 1e-5)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) r
        write(1,'(a,E10.2)') "tolerance: " , r

        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "The maximum number of iteration ? (write any number)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) o
        write(1,'(a,I10)') "itermax: " , o


        do while (.true.)
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "The size of the box ? (default : 2*pi, pi*sqrt(3) )"
        write(*,'(a)') "Do you want default size ? (Y/N) "
        write(*,'(a)') "----------------------------------------------"
        read(*,*) answer
            if (answer == "y" .or. answer == "Y") then
                LLx = 2*pi; LLy = 2*pi ;LLz = 2*pi 
                exit
            elseif (answer == "n" .or. answer == "N") then 
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "enter the size of the box (write it as :  *****  *****)"
                write(*,'(a)') "----------------------------------------------"
                read*, LLx , LLy , LLz
                exit
            else
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "your answer is not (Y/N) , kindly try again"
                cycle
            end if 
        end do 

        write(1,'(a,f20.16,f20.16,f20.16)') 'box:' , LLx , LLy*sqrt(3.d0)/2.d0 , LLz


        write(1,'(a)') "periodic"
        
        do while (.true.)
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Show all the results (with the geometry in every steps) ? (Y/N)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) answer
            if (answer == "y" .or. answer == "Y") then
            write(1,'(a)') "show"
                exit
            elseif (answer == "n" .or. answer == "N") then 
                exit
            else
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "your answer is not (Y/N) , kindly try again"
                cycle
            end if 
        end do

        do while (.true.)
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Generate video animation ? (Y/N)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) answer
            if (answer == "y" .or. answer == "Y") then
            write(1,'(a)') "animation"
                exit
            elseif (answer == "n" .or. answer == "N") then 
                exit
            else
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "your answer is not (Y/N) , kindly try again"
                cycle
            end if 
        end do


        do while (.true.)
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Show the hessian matrix ? (Y/N)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) answer
            if (answer == "y" .or. answer == "Y") then
            write(1,'(a)') "hessian"
                exit
            elseif (answer == "n" .or. answer == "N") then 
                exit
            else
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "your answer is not (Y/N) , kindly try again"
                cycle
            end if 
        end do
    
      write(1,'(a)') "geometry"
      allocate (x(0:4*n**2-1),y(0:4*n**2-1))
      d=pi/n
      dx=d
      dy=d*sqrt(3.d0)/2.d0
      ipoint=0
      
      do iy=0,n-1

         do ix=0,2*n-1
            ipoint=ipoint+1
            x(ipoint)=real(ix,8)*dx
            y(ipoint)=real(2*iy,8)*dy
            write (1,'(i3,f20.16,f20.16,a)') ipoint,x(ipoint), y(ipoint) 
         enddo

         do ix=0,2*n-1
            ipoint=ipoint+1
            x(ipoint)=(real(ix,8)+0.5d0)*dx
            y(ipoint)=real(2*iy+1,8)*dy
            write (1,'(i3,f20.16,f20.16,a)') ipoint,x(ipoint), y(ipoint) 
         enddo
      enddo
    
    exit
 
    
    case("ED", "ed")
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Number of electrons? (write a number) n "
        write(*,'(a)') "----------------------------------------------"
        read(*,*) n 
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a,I5)') "The Total number of electron in this case is " , n
        write(1,'(a)') "dimension: 1"
        write(1,'(a,I5)') "electron: " , n
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Tolerance ? (write any number like this 1e-5)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) r
        write(1,'(a,E10.2)') "tolerance: " , r

        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "The maximum number of iteration ? (write any number)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) o
        write(1,'(a,I10)') "itermax: " , o


        do while (.true.)
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "The size of the box ? (default : 2*pi )"
        write(*,'(a)') "Do you want default size ? (Y/N) "
        write(*,'(a)') "----------------------------------------------"
        read(*,*) answer
            if (answer == "y" .or. answer == "Y") then
                LLx = 2*pi; LLy = 2*pi ;LLz = 2*pi 
                exit
            elseif (answer == "n" .or. answer == "N") then 
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "enter the size of the box (write it as :  ***** )"
                write(*,'(a)') "----------------------------------------------"
                read*, LLx , LLy , LLz
                exit
            else
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "your answer is not (Y/N) , kindly try again"
                cycle
            end if 
        end do 

        write(1,'(a,f20.16,f20.16,f20.16)') 'box:' , LLx , LLy , LLz


        write(1,'(a)') "periodic"
        
        
        do while (.true.)
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Show all the results (with the geometry in every steps) ? (Y/N)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) answer
            if (answer == "y" .or. answer == "Y") then
            write(1,'(a)') "show"
                exit
            elseif (answer == "n" .or. answer == "N") then 
                exit
            else
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "your answer is not (Y/N) , kindly try again"
                cycle
            end if 
        end do

        do while (.true.)
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Generate video animation ? (Y/N)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) answer
            if (answer == "y" .or. answer == "Y") then
            write(1,'(a)') "animation"
                exit
            elseif (answer == "n" .or. answer == "N") then 
                exit
            else
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "your answer is not (Y/N) , kindly try again"
                cycle
            end if 
        end do


        do while (.true.)
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "Show the hessian matrix ? (Y/N)"
        write(*,'(a)') "----------------------------------------------"
        read(*,*) answer
            if (answer == "y" .or. answer == "Y") then
            write(1,'(a)') "hessian"
                exit
            elseif (answer == "n" .or. answer == "N") then 
                exit
            else
            write(*,'(a)') "----------------------------------------------"
                write(*,'(a)') "your answer is not (Y/N) , kindly try again"
                cycle
            end if 
        end do
    
      write(1,'(a)') "geometry"
        d=2*pi/n
        icont= 0 
        do i = 1, n
        icont=icont+1
         write (1,'(i4,(2x,f20.18),a)') icont, d*i 
        end do
    
    exit
    
end select
end do 




end subroutine
end module

program ING
use functions
implicit none 

character (len=200) :: answer , answer2 , arg
real*8  :: r 
real*16, parameter :: pi = ACOS(-1.0d0)
real*16 :: LLx, LLy , LLz

call system('clear')
if(command_argument_count().eq.0)then
write(*,'(a)') "je suis désolée, you didn't add the input file"
write(*,'(a)') "----------------------------------------------------"
write(*,'(a)') "you can run the program using ./ING  << the name of the input file >>"
print *, ""
stop
else 
call get_command_argument(1,arg)
end if 

open(1,file=arg,status = 'replace')
do while (.true.)
write(*,'(a)') "----------------------------------------------"
write(*,'(a)') "Do you want to generate known structure ? like (FCC, HEX, ...) (Y/N)"
write(*,'(a)') "----------------------------------------------"
read(*,*) answer
    if (answer == "y" .or. answer == "Y") then
        call yesanswer
        exit 
    elseif (answer == "n" .or. answer == "N") then 
        call noanswer
        exit
    else
        write(*,'(a)') "----------------------------------------------"
        write(*,'(a)') "your answer is not (Y/N) , kindly try again" 
        cycle        
    end if 
end do 

close(1)

call system('rm -rf functions.mod')
end program 