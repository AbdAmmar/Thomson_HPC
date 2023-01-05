program Thomson

!!!!! from other sources !!!
use omp_lib
!!!!! My own codes in include folder !!!!
use Read_f_f
use Thomson_problem
use Print_g
use Diag
use Harmonic
use parabola_test
!!!!!!!!!!!!!!!!!!!!!!!

implicit none

!! integer and strings !!
	
integer :: i ,j , k , space , electron , iter = 1 , itermax , N , coor , nlines
integer, allocatable :: Nu(:)
character (len=200) ::  ranaa , boxa, showa ,arg  ,multia , hara , curva , pera , equaa , anima , showhessiana , distancea
character (len=5),allocatable  :: fixed(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        real *8        !!

real*8 :: test(3,1) , para , tol , f
real*16, parameter :: pi = ACOS(-1.0d0)
real*8, dimension(:),allocatable :: Eign , FC
real*8, dimension(:,:),allocatable ::  Kon, S, beta, lambda,  H_mod, EignV , H , HHP , HHP_mod 

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        real *16       !!

real*16 :: Lx , Ly , Lz , LLx , LLy , LLz 
real*16, dimension(:),allocatable :: TEb, TEa,TEbHP ,TEaHP,TEbT, testt1 , testt2 , testt3  , Gi , Gp , GpHP , GiHP
real*16, dimension(:,:),allocatable :: Gim , Gip , GipHP ,GimHP, X ,Xp , X_pos , X_old , Xb , Xa , GipT  , R0 , y , yp , X_distance

!!!!!! to clean the terminal and option or argument !!!!!

call system('clear')
if(command_argument_count().eq.0)then
write(*,'(a)') "Je suis désolée, you didn't add the input file"
write(*,'(a)') "----------------------------------------------------"
write(*,'(a)') "you can run the program using ./Thomson << the name of the input file >>"
print *, ""
stop
else 
call get_command_argument(1,arg)
end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call omp_set_num_threads(12)
!!!!!!!!!!!  read from the input file !!!!!!!!!!!!!
call read_from_file(arg,nlines,space,N,itermax,tol,LLx,LLy,LLz,showa,ranaa,y,Nu,multia,hara,curva,pera,yp,FC,fixed,anima,showhessiana,distancea)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!      to get Lx,Ly,Lz and the number of coordinates    !!!
Lx=2.d0*pi/LLx ; Ly=2.d0*pi/LLy; Lz=2.d0*pi/LLz 
coor=space*N
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!      The title of the program            !!!!!!

write(*,'(a)') "   _____ _                                                       _     _ "
write(*,'(a)') "  |_   _| |                                                     | |   | |"
write(*,'(a)') "    | | | |__   ___  _ __ ___  ___  ___  _ __    _ __  _ __ ___ | |__ | | ___ _ __ ___"
write(*,'(a)') "    | | | '_ \ / _ \| '_ ` _ \/ __|/ _ \| '_ \  | '_ \| '__/ _ \| '_ \| |/ _ \ '_ ` _ \"
write(*,'(a)') "    | | | | | | (_) | | | | | \__ \ (_) | | | | | |_) | | | (_) | |_) | |  __/ | | | | |"
write(*,'(a)') "    \_/ |_| |_|\___/|_| |_| |_|___/\___/|_| |_| | .__/|_|  \___/|_.__/|_|\___|_| |_| |_|"
write(*,'(a)') "                                                | |"
write(*,'(a)') "                                                |_|"
write(*,'(a)') ""
write(*,'(a)') "                           ,,....................,,.                 "
write(*,'(a)') "                     ,..........,,,,..      ............,*           "
write(*,'(a)') "                 ,.......,,,,,****************,,,,,,.......,,*       "
write(*,'(a)') "              *,.....,,,,***/////((((((((((/////****,,,......,,*     "
write(*,'(a)') "            /,,.....,,**//(((###%%%%%%%%%%%%###(((//**,,,.....,,*/   "
write(*,'(a)') "           /*,,.....,,*/(##%%%%%%%%%%%%%%%%%%%%%%%##(/**,......,,*/  "
write(*,'(a)') "          (/*,,......,*/#%%%,                     %%%#/,,.....,,**/( "
write(*,'(a)') "          #/**,,,......,,                            ,......,,,**/(# "
write(*,'(a)') "          (#(/**,,,,........,                    ,.......,,,***/((#% "
write(*,'(a)') "           %#((//***,,,,,..........................,,,,,,***//((#%%  "
write(*,'(a)') "            *%%#((///*****,,,,,.       ....,,,,,,,*****///(((#%%%%   "
write(*,'(a)') "              #%%%##((((/////******************/////((((###%%%%%     "
write(*,'(a)') "                 %%%%%%%####((((((((((((((((((((####%%%%%%%%%        "
write(*,'(a)') "                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            "
write(*,'(a)') "                            #%%%%%%%%%%%%%%%%%%%%%                   "
write(*,'(a)') ""
write(*,'(a)') " By : A.ALRAKIK                                                Supervision: S.Evangelisti   "
write(*,'(a)') "____________________________________________________________________________________________"
write(*,'(a)') ""
if (ranaa == "random") then
write(*,'(a,i4,a,i1,a,i7,a,a)') "The calculation for " , N ,"  electron in " , space   ,"D dimension , with ", itermax , " loops and random geometry"
else
write(*,'(a,i4,a,i1,a,i7,a,a)') "The calculation for " , N ,"  electron in " , space   ,"D dimension , with ", itermax , " loops and spacific geometry"
end if
write(*,'(a)') "____________________________________________________________________________________________"
write(*,'(a)') ""
write(*,'(a)') "                           ||  The size of the box  ||"
write(*,'(a)') ""
write (*,'(a,f20.16,a,f20.16,a,f20.16,a)') "X = " , LLx ,"   Y = " , LLy , "    Z = " , LLz
write(*,'(a)') ""
write(*,'(a)') "                           || The Pi parameter ||"
write(*,'(a)') ""
write(*,'(a,f20.16,a,f20.16)') "Pi = " , pi , "     2*Pi = " , 2.d0*pi
write(*,'(a)') ""
write(*,'(a)') "                           ||  The ratio between 2*pi and the length of the box  ||"
write(*,'(a)') ""
write (*,'(a,f14.12,a,f14.12,a,f14.12)') "Lx = " , Lx ,"  Ly = " , Ly , "  Lz = " , Lz
write(*,'(a)') ""
write(*,'(a)') "____________________________________________________________________________________________"
write(*,'(a)') "                               _                _             _ "
write(*,'(a)') "                          ___ | |_  __ _  _ __ | |_  ___   __| |"
write(*,'(a)') "                         / __|| __|/ _` || '__|| __|/ _ \ / _` |"
write(*,'(a)') "                         \__ \| |_| (_| || |   | |_|  __/| (_| |"
write(*,'(a)') "                         |___/ \__|\__,_||_|    \__|\___| \__,_|"
write(*,'(a)') "____________________________________________________________________________________________"
write(*,'(a)') ""
!!!   allocate our variable  !!!

allocate (X(coor,1))
allocate (Xp(coor,1))
allocate (X_old(coor,1))
allocate (Xb(coor,1))
allocate (Xa(coor,1))
allocate (S(coor,1))
allocate (Gim(coor,1))
allocate (GimHP(coor,1))
allocate (Gip(coor,1))
allocate (GipHP(coor,1))
allocate (GipT(coor,1))
allocate (X_pos(N,space))
allocate (H_mod(coor,coor))
allocate (HHP_mod(coor,coor))
allocate (Eign(coor))
allocate (EignV(coor,coor))
allocate (lambda(1,1))
allocate (beta(1,1))
allocate (TEbT(N))
allocate (TEbHP(N))
allocate (TEb(N))
allocate (TEa(N))
allocate (TEaHP(N))
allocate (testt1(N))
allocate (testt2(N))
allocate (testt3(N))
allocate (R0(N,N))
allocate (X_distance(N,N))

if (space == 1 ) then 
	allocate (Gp(3*coor))
	allocate (Gi(3*coor))
	allocate (GiHP(3*coor))
	allocate (GpHP(3*coor))
	allocate (H(3*coor,3*coor))
else if (space == 2) then
	allocate (Gp(2*coor))
	allocate (Gi(2*coor))
	allocate (GpHP(2*coor))
	allocate (H(2*coor,2*coor))
elseif (space == 3) then
	allocate (Gp(coor))
	allocate (Gi(coor))
	allocate (GpHP(coor))
	allocate (GiHP(coor))
	allocate (HHP(coor,coor))
	allocate (H(coor,coor))
end if 


!!! condition for random !!!

if (ranaa == "random") then 
	do i=1,coor
		call random_number(X(i,1))
	end do
	else 
	do i=1,coor
	if (i <= N ) then 
	X(i,1)=y(i,1)
	else if (i > N .and. i <= 2*N) then 
    X(i,1)=y(i-N,2)
	else 
	X(i,1)=y(i-2*N,3)
	end if 
	end do
end if


!!! condition for multiplication by L  !!!
if (multia  == "multiply") then 
do i=1,coor
	if (i <= N ) then 
	X(i,1)=LLx*X(i,1)
	else if (i > N .and. i <= 2*N) then 
    X(i,1)=LLy*X(i,1)
	else if (i > 2*N ) then
	X(i,1)=LLz*X(i,1)
end if 
end do 
end if 

!!! get the force constant !!!
if (hara == "harmonic") then
allocate(Kon(8,8))
call convertFC(FC,Kon)
end if 

!!!!!!!!!!!!!!!! print the geometry  !!!!!!!!!!!!!!!!!!!
write(*,'(a)') "The geometry before the Optimization"
call prin(N,X,Space,X_pos)


!!!  get the geodesic distance   !!!
if (hara == "harmonic") then
do i=1,coor
	if (i <= N ) then 
	Xp(i,1)=yp(i,1)
	else if (i > N .and. i <= 2*N) then 
    Xp(i,1)=yp(i-N,2)
	else 
	Xp(i,1)=yp(i-2*N,3)
	end if 
	end do
if (space == 3) then
call dista(N,Xp,R0)
end if 
end if 


!!!!!!!!!!!!!!!! The magic start here !!!!!!!!!!!!!!!!!!

!!!  get the derivative !!!

if (hara == "harmonic") then 
call diffHP(Kon,N,X,space,Lx,Ly,Lz,GpHP,R0)
GpHP=-GpHP
do i = 1 , size(GipHP)
	GipHP(i,1)=GpHP(i)  !! make it as matrix to make transpose 
end do 
do i=1,size(Gip)
	Gip(i,1)=GipHP(i,1)
end do 
else 
call diff(N,X,space,Lx,Ly,Lz,Gp) !! the derivative before the loops
Gp=-Gp
do i = 1 , size(Gip)
Gip(i,1)=Gp(i)  !! make it as matrix to make transpose 
end do
end if 



!!! condition for NAN or Problem !!!

if (isnan(Norm2(GiP)) .or. Norm2(Gip) > 1e20) then
write(*,'(a)') "____________________________________________________________________________________________"
write(*,'(a)') ""
write (*,'(a)') "you have 2 or more electron at the same position or the energy is infinity"
Stop
end if

!!!!! define the first conjucated vector !!!
S=Gip

!!! calculate the energy !!!

if (hara == "harmonic") then 
call energyHP(Kon,N,X,space,Lx,Ly,Lz,TEbHP,R0)	
do i=1, size(TEb)
	TEb(i)=TEbHP(i)
end do
write(*,'(a)') "____________________________________________________________________________________________"
write(*,'(a)') ""
write(*,'(a,f30.16)') "The Total Energy before optimization = " , sum(TEb)
write(*,'(a)') "____________________________________________________________________________________________"
write(*,'(a)') ""
else 
call energy(N,X,space,Lx,Ly,Lz,TEb)
write(*,'(a)') "____________________________________________________________________________________________"
write(*,'(a)') ""
write(*,'(a,f30.16)') "The energy before optimization = " , sum(TEb)
write(*,'(a)') "____________________________________________________________________________________________"
write(*,'(a)') ""
end if 
if (showa /= 'show') then 
write(*,'(a)') ' Loop  |        Energy         |       Energy differance       |         Derivative '
write(*,'(a)') "____________________________________________________________________________________________"
write(*,'(a)') ""
end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!! calculate the second derivative !!!
if (hara == "harmonic") then 
call diffHP(Kon,N,X,space,Lx,Ly,Lz,GiHP,R0)
GiHP=-GiHP
do i = 1 , size(GiHP)
Gim(i,1)=GiHP(i)
end do
else 
call diff(N,X,space,Lx,Ly,Lz,Gi)
Gi=-Gi
do i = 1 , size(Gim)
Gim(i,1)=Gi(i)          
end do
end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! condition for One calculation Only !!!

if (itermax == 0 ) then 
if (hara == "harmonic") then
call hessienHP(Kon,N,X,space,Lx,Ly,Lz,HHP,R0)
	do i=1,coor
			do j = 1,coor
				H_mod(i,j) =   HHP(i,j)	
			end do 
	end do
else 
call hessien(N,X,space,Lx,Ly,Lz,H)
	do i=1,coor
		do j = 1,coor
			H_mod(i,j) =   H(i,j)	
		end do 
	end do
end if
end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (anima == "animation") then 
open(4,file='data_frame.dat',status = 'replace')
open(8,file='energy.dat',status = 'replace')
end if 

do while (iter < itermax + 1)
		if (norm2(Gip) < tol) then
		if (hara == "harmonic") then
			call energyHP(Kon,N,X,space,Lx,Ly,Lz,TEaHP,R0)
			do i=1,size(TEa)
			TEa(i)=TEaHP(i)
			end do 
		else 
			call energy(N,X,space,Lx,Ly,Lz,TEa)
		end if 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(*,'(a,f30.14,a,i9,a)') "you are at the minimum",  sum(TEa) ,"   after " ,iter - 1 , " Loops"
write(*,'(a,E20.10)') "The derivative =     ", norm2(Gim)
write(*,'(a)') "____________________________________________________________________________________________"
write(*,'(a)') ""

if (iter == 1) then 
if (space == 3) then 
write(4,'(10000(1x,f16.8))') X_pos(:, 1),X_pos(:, 2),X_pos(:,3)
elseif (space == 2) then
write(4,'(10000(1x,f16.8))') X_pos(:, 1),X_pos(:, 2)
elseif (space == 1) then
write(4,'(10000(1x,f16.8))') X_pos(:, 1)
end if
end if


if (hara == "harmonic") then
call hessienHP(Kon,N,X,space,Lx,Ly,Lz,HHP,R0)
	do i=1,coor
		do j = 1,coor
			HHP_mod(i,j) =   HHP(i,j)	
		end do 
	end do
	
do i=1,coor
	do j=1,coor
	H_mod(i,j)=HHP_mod(i,j)
	end do 	
end do 
exit
else 
call hessien(N,X,space,Lx,Ly,Lz,H)
	do i=1,coor
		do j = 1,coor
			H_mod(i,j) =   H(i,j)	
		end do 
	end do
	exit
end if 
end if 

if (isnan(norm2(GiP)) .or. norm2(Gip) > 1e20) then
write(*,'(a)') "____________________________________________________________________________________________"
write(*,'(a)') ""
write(*,'(a)')  "you have 2 or more electron at the same position"
exit
end if 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!  condition if reach the minimum !!!!!!!

if (norm2(Gim) < tol) then
if (showa /= 'show') then 
write(*,'(a)') "____________________________________________________________________________________________"
write(*,'(a)') ""
end if 
    write(*,'(a,f20.16,a,i5,a)') "Geometry converged at",  sum(TEa) ,"   after " ,iter-1  , " Loops"
    write(*,'(a)') ""
	write(*,'(a,E20.10)') "The derivative =     ", norm2(Gim)
exit
end if 

if (showa == 'show') then 
write(*,'(a,i7)') "Loop = " , iter
end if 

if (hara == "harmonic") then
call diffHP(Kon,N,X,space,Lx,Ly,Lz,GiHP,R0)
GiHP=-GiHP
do i = 1 , size(GiHP)
GimHP(i,1)=GiHP(i)          !! make it as matrix to make transpose
end do
do i=1,size(Gim)
Gim(i,1)=GimHP(i,1)
end do 
else 
call  diff(N,X,space,Lx,Ly,Lz,Gi)      !! the derivitave in the loops
Gi=-Gi
do i = 1 , size(Gim)
Gim(i,1)=Gi(i)          !! make it as matrix to make transpose
end do
end if 


if (hara == "harmonic") then 
call hessienHP(Kon,N,X,space,Lx,Ly,Lz,HHP,R0) 
do i=1,coor
	do j = 1,coor
		HHP_mod(i,j) =  HHP(i,j)	
	end do 
end do

do i=1,coor
		do j = 1,coor
		H_mod(i,j) = HHP_mod(i,j)	
		end do 
end do
else 
call hessien(N,X,space,Lx,Ly,Lz,H)  !! calculate the hessien matrix

!!! make hessien work for every dimension !!!!!
do i=1,coor
		do j = 1,coor
		H_mod(i,j) = H(i,j)	
		end do 
end do
end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

beta = matmul(transpose(Gim),matmul(H_mod,S))/ matmul(transpose(S),matmul(H_mod,S))
!!!!!!!!    the conjugated vector    !!!!!!!!!!
do i=1,size(S)
S(i,1) =  Gim(i,1) + beta(1,1) * S(i,1)
end do
!!!!!!!! the step for the geometry   !!!!!!!!!!
lambda = matmul(transpose(Gim),Gim)/matmul(transpose(S),matmul(H_mod,S))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,size(X)
	X_old(i,1)=X(i,1)
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! parabola implemintation !!!!!!!!!!!
if (curva == "parabola") then 
	do i=1,size(X)
	Xa(i,1)=X(i,1)+lambda(1,1)*S(i,1)
	end do
	do i=1,size(X)
	Xb(i,1)=X(i,1)-lambda(1,1)*S(i,1)
	end do
	!!---------------------------------------------!!	
	call Energy(N,X_old,space,Lx,Ly,Lz,testt1)
	test(1,1)=sum(testt1)
	call Energy(N,Xa,space,Lx,Ly,Lz,testt2)
	test(2,1)=sum(testt2)
	call Energy(N,Xb,space,Lx,Ly,Lz,testt3)
	test(3,1)=sum(testt3)
	
	call parabola(test,para)
	
	if (para > 0) then 
	do i=1,size(X)
	X(i,1)=X(i,1)+para*S(i,1)
	end do
	end if 
	if (para < 0) then 
	if ( test(1,1) > test(2,1)) then 
	do i=1,size(X)
	X(i,1)=X(i,1)-lambda(1,1)*S(i,1)
	end do
	else if (test(1,1) < test(2,1)) then 
	do i=1,size(X)
	X(i,1)=X(i,1)+lambda(1,1)*S(i,1)
	end do
	end if 
	end if 
end if 
!!!!!!!!!!! parabola implemintation  end !!!!!!!!!!!
if (curva /= "parabola") then 
do i=1,size(X)
	X(i,1)=X(i,1)+lambda(1,1)*S(i,1)
end do
end if 

do i = 1,N
	if (fixed(i) == "f") then 
		X(i,1) = X(i,1)-lambda(1,1)*S(i,1)
		X(i+N,1) = X(i+N,1)-lambda(1,1)*S(i+N,1)
		X(i+2*N,1) = X(i+2*N,1)-lambda(1,1)*S(i+2*N,1)
	end if 
end do 


if (pera == "periodic") then 
	
	if (space == 1 ) then 
	do i = 1,size(X)
	do while (X(i,1) < 0 )
		X(i,1) = X(i,1) + LLx
	end do 
	do while (X(i,1) > LLx)
			X(i,1) = X(i,1) - LLx 
	end do
	if (X(i,1) == 0) then
	X(i,1)=X(i,1)
	end if 
	end do

	else if (space == 2) then 
	do i = 1,(size(X)/2)
	do while (X(i,1) < 0 )
	X(i,1) = X(i,1) + LLx
	end do 
	do while (X(i,1) > LLx)
			X(i,1) = X(i,1) - LLx 
	end do
	if (X(i,1) == 0) then
	X(i,1)=X(i,1)
	end if 
	end do

	do i = (size(X)/2)+1,size(X)
	do while (X(i,1) < 0 )
	X(i,1) = X(i,1) + LLy
	end do 
	do while (X(i,1) > LLy)
			X(i,1) = X(i,1) - LLy
	end do
	if (X(i,1) == 0) then
	X(i,1)=X(i,1)
	end if 
	end do

	else if (space == 3) then 

	do i=1,3*N
	if (i <= N ) then 
	do while (X(i,1) < 0 )
	X(i,1) = X(i,1) + LLx
	end do 
	do while (X(i,1) > LLx)
		X(i,1) = X(i,1) - LLx 
	end do
	if (X(i,1) == 0) then
	X(i,1)=X(i,1)
	end if 
	else if (i > N .and. i <= 2*N) then
	do while (X(i,1) < 0 )
	X(i,1) = X(i,1) + LLy
	end do 
	do while (X(i,1) > LLy)
		X(i,1) = X(i,1) - LLy
	end do
	if (X(i,1) == 0) then
	X(i,1)=X(i,1)
	end if 
    else 
	do while (X(i,1) < 0 )
	X(i,1) = X(i,1) + LLz
	end do 
	do while (X(i,1) > LLz)
		X(i,1) = X(i,1) - LLz
	end do
	if (X(i,1) == 0) then
	X(i,1)=X(i,1)
	end if 
	end if 
	end do 
	end if
	
end if 

!!!!!!!!!!    print the final geometry  !!!!!!!!!!!!
if (showa == "show") then
write(*,'(a)') ""
write(*,'(a)') "The geometry after the optimization:"
call prin(N,X,Space,X_pos)
else 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!! have a look at the print function !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (space == 1 ) then 
	do i = 1,N
	X_pos(:,1) = X(1:N,1)
	end do
	else if (space == 2) then
	do i = 1,2*N 
	X_pos(:,1) = X(1:N,1)
	X_pos(:,2) = X(N+1:2*N,1)
	end do
	else if (space == 3) then
	do i = 1,3*N
	X_pos(:,1) = X(1:N,1)
	X_pos(:,2) = X(N+1:2*N,1)
	X_pos(:,3) = X(2*N+1:3*N,1)
	end do
	end if
end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (hara == "harmonic") then 
call energyHP(Kon,N,X,space,Lx,Ly,Lz,TEaHP,R0)
do i=1,size(TEa)
TEa(i)=TEaHP(i)
end do 
else 
call energy(N,X,space,Lx,Ly,Lz,TEa)
end if 

if (isnan(Norm2(Gim)) .or. Norm2(Gim) > 1e20) then
write(*,'(a)') "____________________________________________________________________________________________"
write(*,'(a)') ""
write (*,'(a)') "you have 2 or more electron at the same position or the energy is infinity" 
Stop
end if

if (showa == 'show') then
write(*,'(a)') "____________________________________________________________________________________________"
write(*,'(a)') ""
write(*,'(a,f16.12,a,E8.2)') "The energy after optimization =  ",  sum(TEa) , "        |   The derivative =  ", norm2(Gim)
write(*,'(a)') "____________________________________________________________________________________________"
write(*,'(a)') ""
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (showa /= 'show') then
call energy(N,X_old,space,Lx,Ly,Lz,TEb)
write(*,'(I4,a,f20.16,a,f20.16,a,E12.5)')  iter ,"    ", sum(TEa) ,"         ", abs(sum(TEa)-sum(TEb)),"              ", Norm2(Gim) 
end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
iter = iter + 1	
if (anima == "animation") then 
if (space == 3) then 
write(4,'(10000(1x,f16.8))') X_pos(:, 1),X_pos(:, 2),X_pos(:,3)
write(8,'(f24.18)') sum(TEa)
elseif (space == 2) then
write(4,'(10000(1x,f16.8))') X_pos(:, 1),X_pos(:, 2)
write(8,'(f24.18)') sum(TEa)
elseif (space == 1) then
write(4,'(10000(1x,f16.8))') X_pos(:, 1)
write(8,'(f24.18)') sum(TEa)
end if 
end if  
end do 

!!!***********************************************!!!
!!!!!!!!!!!!!!! tha magic end :-(  !!!!!!!!!!!!!!!!!!

if (showa /= 'show' ) then
write(*,'(a)') "____________________________________________________________________________________________"
write(*,'(a)') ""
end if

if (showhessiana == "hessian") then
write(*,'(a)') "____________________________________________________________________________________________"
write(*,'(a)') ""
write(*,'(a)') '                                  The Hessian matrix after the converge '
write(*,'(a)')
if (space == 1) then
do i=1,N
write (*,'(100(3x,f16.8))') H_mod(i,:)
end do
elseif  (space == 2) then
do i=1,2*N
write (*,'(100(3x,f16.8))') H_mod(i,:)
end do
elseif  (space == 3) then
do i=1,3*N
write (*,'(100(3x,f16.8))') H_mod(i,:)
end do
end if
end if 


if (distancea == 'distance') then 
    call prin_distance(X,N,space,X_distance,lx,ly,lz)
end if 



call Diagonalization(H_mod,Eign,EignV,coor)
if (showa == 'show') then
write(*,'(a)') "____________________________________________________________________________________________"
write(*,'(a)') ""
print *, "            The EignValues of the hessien"
print *, ""
do i = 1,coor 
write (*,'(f20.12)') Eign(i)
end do
write(*,'(a)') "____________________________________________________________________________________________"
write(*,'(a)') ""
end if 
do i=1,coor
if (abs(Eign(i)) > 1e-6 .and. Eign(i) < 0) then 
write (*,'(a)') "you are at stationary state (Negative eignvalue)"
exit
else 
write (*,'(a)') "all the eignvalue is positive"
exit
end if 
end do 

!!!!!!! print the final geometry out !!!!!!!!!!!!!!!! 
open(2,file='data.out',status = 'replace')
	write(2,'(a,I1)') 'dimension: ' , space
	write(2,'(a,I3)') 'electron: ' , N
	write(2,'(a,I15)') 'itermax:' , itermax
	write(2,'(a,E20.4)') 'tolerance:' , tol
	write(2,'(a,f30.16,f30.16,f30.16)') 'box:' , LLx , LLy , LLz
    if (showa == "show") then 
	write(2,'(a)') 'show'
    end if 
	if (pera == "periodic") then
	write(2,'(a)') 'periodic'
	end if 
	if (hara == "harmonic") then
	write(2,'(a)') 'harmonic'
	end if 
	write(2,'(a)') 'geometry'
	DO i = 1, N 
	if (space == 1) then 
		write(2,'(I5,f24.16,a)') Nu(i) , X_pos(i, :) 
	elseif (space == 2) then 
	    write(2,'(I5,f24.16,f24.16,a)') Nu(i) , X_pos(i, :) 
	elseif (space == 3) then 
	    write(2,'(I5,f24.16,f24.16,f24.16,a)') Nu(i) , X_pos(i, :) 
	end if 
	end do
close(2)

! for animation !

if (anima == "animation") then
if (iter .ne. 1) then
write(*,'(a)') "____________________________________________________________________________________________"
write(*,'(a)') ""
write(*,'(a)') "Please wait for the animation"
write(*,'(a)') "____________________________________________________________________________________________"
write(*,'(a)') ""
open(9,file='plot',status = 'replace')
	write(9,'(a)') 'set term png size 1920,1080'
	write(9,'(a)') 'system "mkdir tmp"'
	write(9,'(a)') "getValue(row,col,filename) = system('awk ''{if (NR == '.row.') print $'.col.'}'' '.filename.'')"
	write(9,'(a,f16.10,a)') 'set xrange [0:',LLx,']'
	write(9,'(a,f16.10,a)') 'set yrange [0:',LLy,']'
	write(9,'(a,f16.10,a)') 'set zrange [0:',LLz,']'
	write(9,'(a)') 'set xtics ("0" 0, "π" 3.14)'
	write(9,'(a)') 'set ytics ("2π = 0 " 0,"π" 3.14,"2π = 0 " 6.28)'
	write(9,'(a)') 'set ztics ("π" 3.14,"2π = 0 " 6.28)'
	if (space == 3) then 
		write(9,'(a)') 'set xtics offset 1,-1,0'
		write(9,'(a)') 'set ytics offset 1,-1,0'
		write(9,'(a)') 'set ztics offset 0,0,0'
		write(9,'(a)') 'set xyplane relative 0'
		write(9,'(a)') 'set view equal xyz'
		write(9,'(a)') 'set border 4095'
	else if (space == 2) then 
		write(9,'(a)') 'set xtics ("0" 0, "π" 3.14 , "2π = 0" 6.28)'
		write(9,'(a)') 'set size ratio -1'
	else if (space == 1) then
		write(9,'(a)') 'set xtics ("0" 0, "π" 3.14 , "2π = 0" 6.28 )'
		write(9,'(a)') 'set size square'
		write(9,'(a)') 'set border 1'
		write(9,'(a)') 'unset ytics'
		write(9,'(a)') 'set xtics nomirror'
	end if
	write(9,'(a)') 'set grid'
	write(9,'(a,I4)') 'N=' , N	
	write(9,'(a,I5,a)') 'do for [i=1:',iter-2,']{'
	write(9,'(a)') 'set output "tmp/image.".i.".png"'
	write(9,'(a)') "x = getValue(i,1,'energy.dat')"
	if (space == 3) then 
		write(9,'(a)') "set label sprintf('Energy = %16.12f', x*1.0) at 9.5,6,7 center font ',18'"
		write(9,'(a,I4,a,I4,a,I4,a)') "splot for [j=1:",N,"] 'data_frame.dat' u (column(j)):(column(",N,"+j)):(column(",2*N,"+j)) every ::i::i w p pt 7 ps 3 notitle"
	else if (space == 2) then 
		write(9,'(a)') "set label sprintf('Energy = %16.12f', x*1.0) at 7.5,6,7 center font ',18'"
		write(9,'(a,I4,a,I4,a)') "p for [j=1:",N,"] 'data_frame.dat' u (column(j)):(column(",N,"+j)) every ::i::i w p pt 7 ps 3 notitle"
	else if (space == 1 ) then 
		write(9,'(a)') "set label sprintf('Energy = %16.12f', x*1.0) at 7.5,6,7 center font ',18'"
		write(9,'(a,I4,a,I4,a)') "p for [j=1:",N,"] 'data_frame.dat' u j:(0) every ::i::i w p pt 7 ps 3 notitle"
	end if 
	write(9,'(a)') "clear"
	write(9,'(a)') "unset label"
	write(9,'(a)') "}"
	write(9,'(a)') "system 'ffmpeg -y -nostats -loglevel 0 -f image2 -r 15.0 -i tmp/image.%d.png -pix_fmt yuv420p -crf 1 animation.mp4'"
	write(9,'(a)') "system 'rm -rf tmp/'"
close(9)
else 
write(*,'(a)') "____________________________________________________________________________________________"
write(*,'(a)') ""
write(*,'(a)') 'I can not provide you animation for one frame ^_^ but i can give you a pic as png'
write(*,'(a)') "____________________________________________________________________________________________"
write(*,'(a)') ""
open(10,file='energy.dat',status = 'replace')
write(10,*)  sum(TEa)
close(10)
open(9,file='plot',status = 'replace')
	write(9,'(a)') 'set term png size 1920,1080'
	write(9,'(a)') 'set xrange [0:6.28]'
	write(9,'(a)') 'set yrange [0:6.28]'
	write(9,'(a)') 'set zrange [0:6.28]'
	write(9,'(a)') 'set xtics ("0" 0, "π" 3.14)'
	write(9,'(a)') 'set ytics ("2π = 0 " 0,"π" 3.14,"2π = 0 " 6.28)'
	write(9,'(a)') 'set ztics ("π" 3.14,"2π = 0 " 6.28)'
	if (space == 3) then 
		write(9,'(a)') 'set xtics offset 1,-1,0'
		write(9,'(a)') 'set ytics offset 1,-1,0'
		write(9,'(a)') 'set ztics offset 0,0,0'
		write(9,'(a)') 'set xyplane relative 0'
		write(9,'(a)') 'set view equal xyz'
	else if (space == 2) then 
		write(9,'(a)') 'set xtics ("0" 0, "π" 3.14 , "2π = 0" 6.28)'
		write(9,'(a)') 'set size square'
	else if (space == 1) then
		write(9,'(a)') 'set xtics ("0" 0, "π" 3.14 , "2π = 0" 6.28)'
		write(9,'(a)') 'set size square'
		write(9,'(a)') 'set border 1'
		write(9,'(a)') 'unset ytics'
		write(9,'(a)') 'set xtics nomirror'
	end if
	write(9,'(a)') 'set grid'
	write(9,'(a,I4)') 'N=' , N	
	write(9,'(a)') 'set output "figure.png"'
	write(9,'(a)') "getValue(row,col,filename) = system('awk ''{if (NR == '.row.') print $'.col.'}'' '.filename.'')"
	write(9,'(a)') "x = getValue(1,1,'energy.dat')"
	if (space == 3) then 
		write(9,'(a)') "set label sprintf('Energy = %16.12f', x*1.0) at 9.5,6,7 center font ',18'"
		write(9,'(a,I4,a,I4,a,I4,a)') "splot for [j=1:",N,"] 'data_frame.dat' u (column(j)):(column(",N,"+j)):(column(",2*N,"+j)) w p pt 7 ps 3 notitle"
	else if (space == 2) then 
		write(9,'(a)') "set label sprintf('Energy = %16.12f', x*1.0) at 7.5,6,7 center font ',18'"
		write(9,'(a,I4,a,I4,a,I4,a)') "plot for [j=1:",N,"] 'data_frame.dat' u (column(j)):(column(",N,"+j)) w p pt 7 ps 3 notitle"
	else if (space == 1 ) then 
	    write(9,'(a)') "set label sprintf('Energy = %16.12f', x*1.0) at 7.5,6,7 center font ',18'"
		write(9,'(a,I4,a,I4,a,I4,a)') "plot for [j=1:",N,"] 'data_frame.dat' u j:(0) w p pt 7 ps 3 notitle"
	end if 
close(9)
end if 
call system('gnuplot plot')
call system('rm -rf plot')
call system('rm -rf data_frame.dat')
call system('rm -rf energy.dat')
end if 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call system('rm -rf fort.7')
call system('rm -rf fort.3')
call system('rm -rf fort.4')
call system('rm -rf fort.8')
call system('rm -rf functions.mod')
end program 