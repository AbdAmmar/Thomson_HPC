module read_f_f
	implicit none
	contains

    subroutine read_from_file(arg,nlines,space,N,itermax,tol,LLx,LLy,LLz,showa,ranaa,&
                              y,Nu,multia,hara,curva,pera,yp,FC,fixed,anima,         &
                              showhessiana,distancea)
	implicit none
	character (len=200),intent(in) :: arg
	integer :: i,j ,nlines,nnlines, k
	character (len=200) :: dima, dimb, elca, elcb, itra, itrb, tola, tolb, ranaa, ranab, geoa, geob,        &
                           boxa, boxb, showa, showb, multib, multia, hara, harb, curva, curvb ,pera, perb , &
                           equab , equaa , fconaa , fconab , anima , animb , showhessiana , showhessianb ,  &
                           distancea , distanceb , fixa , fixb
	character (len=5) , allocatable ::  fixed(:)
	integer, intent(out) :: space , N , itermax
	real*8 , intent(out) :: tol
	real*16, intent(out):: LLx , LLy , LLz
	real*8 ,allocatable :: FC(:) 
	real*16 , dimension(:,:),allocatable :: y , yp
	real*8 , dimension(:,:),allocatable :: eq
	integer, allocatable :: Nu(:) , Nuu(:) , Nuuu(:,:)
	
	dimb  = "dimension:"
	elcb  = "electron:"
	itrb  = "itermax:"
	tolb  = "tolerance:"
	ranab = "random"
	geob  = "geometry"
	boxb  = "box:"
	showb = "show"
    fixb  = "fixed"
	multib= "multiply"
	harb  = "harmonic"
	curvb = "parabola"
	perb  = "periodic"
	equab = "equilibrium"
	fconab = "Fconstant"
	animb ="animation"
	showhessianb="hessian"
    distanceb = "distance"


	nlines = 0
	OPEN (1, file = arg)
	DO
    READ (1,*, END=1)
    nlines = nlines + 1
	END DO
	1 CLOSE (1)
	
	OPEN (8, file ="INPUT_EQUILIBRUIM_DISTANCE.test")
	DO
    READ (8,*, END=2)
    nnlines = nnlines + 1
	END DO
	2 CLOSE (8)
	
	open (3,file = arg)
	do i=1,nlines
	read (3,*,end=3) dima
	if (dima .EQ. dimb) then
	backspace (3)
	read  (3,*) dima , space
	close(3)
	end if
	end do
	3 close(3) 
	
	open (3,file = arg)
	do i=1,nlines
	read (3,*,end=4) elca
	if (elca .EQ. elcb) then
	backspace (3)
	read  (3,*) elca , N
	close(3)
	end if 
	end do 
	4 close(3)
	
	open (3,file = arg)
	do i=1,nlines
	read (3,*,end=5) itra
	if (itra .EQ. itrb) then
	backspace (3)
	read  (3,*) itra , itermax
	close(3)
	end if 
	end do
	5 close(3)

	open (3,file = arg)
	do i=1,nlines
	read (3,*,end=6) tola
	if (tola .EQ. tolb) then
	backspace (3)
	read  (3,*) tola , tol
	close(3)
	end if 
	end do 
	6 close(3)

	open (3,file = arg)
	do i=1,nlines
	read (3,*,end=7) boxa
	if (boxa .EQ. boxb) then
	backspace (3)
	read  (3,*) boxa , LLx , LLy ,LLz
	close(3)
	end if 
	end do
	7 close(3)
	
	open (3,file = arg)
	do i=1,nlines
	read (3,*,end=8) showa
	if (showa .EQ. showb) then
	backspace (3)
	read  (3,*) showa
	close(3)
	end if 
	end do 
	8 close(3)
	
	open (3,file = arg)
	do i=1,nlines
	read (3,*,end=9) ranaa
	if (ranaa .EQ. ranab) then
	backspace (3)
	read  (3,*) ranaa 
	close(3)
	end if 
	end do 
	9 close(3)
    
    open (3,file = arg)
	do i=1,nlines
	read (3,*,end=10) fixa
	if (fixa .EQ. fixb) then
	backspace (3)
	read  (3,*) fixa 
	close(3)
	end if 
	end do 
	10 close(3)
    
    if (fixa == "fixed") then 
	open (3,file = arg)
	do i=1,nlines
	read (3,*,end=11) geoa
	if (geoa .EQ. geob) then
	backspace (3)
	read  (3,*) geoa 
	allocate (y(N,3))
	allocate (Nu(N))
	allocate (fixed(N))
	do k = 1,N
	if (space == 3) then 
	read  (3,*) Nu(k), y(k,1) , y(k,2) , y(k,3) , fixed(k) 
	else if (space == 2 ) then 
	read  (3,*) Nu(k), y(k,1) , y(k,2) , y(k,3) , fixed(k)
	else if (space == 1 ) then
	read  (3,*) Nu(k), y(k,1) , y(k,2) , y(k,3) , fixed(k)
	end if 
	end do 
	close(3)
	end if 
	end do  
	11 close (3)
    
    
    else 
    
    open (3,file = arg)
	do i=1,nlines
	read (3,*,end=12) geoa
	if (geoa .EQ. geob) then
	backspace (3)
	read  (3,*) geoa 
	allocate (y(N,3))
	allocate (Nu(N))
	allocate (fixed(N))
	do k = 1,N
	if (space == 3) then 
	read  (3,*) Nu(k), y(k,1) , y(k,2) , y(k,3)  
	else if (space == 2 ) then 
	read  (3,*) Nu(k), y(k,1) , y(k,2) 
	else if (space == 1 ) then
	read  (3,*) Nu(k), y(k,1) 
	end if 
	end do 
	close(3)
	end if 
	end do  
	12 close (3)
    
    
    end if 

	open (3,file = arg)
	do i=1,nlines
	read (3,*,end=13) multia
	if (multia .EQ. multib) then
	backspace (3)
	read  (3,*) multia
	close(3)
	end if 
	end do 
	13 close(3)

	open (3,file = arg)
	do i=1,nlines
	read (3,*,end=14) hara
	if (hara .EQ. harb) then
	backspace (3)
	read  (3,*) hara
	close(3)
	end if 
	end do 
	14 close(3)


	open (3,file = arg)
	do i=1,nlines
	read (3,*,end=15) curva
	if (curva .EQ. curvb) then
	backspace (3)
	read  (3,*) curva
	close(3)
	end if 
	end do 
	15 close(3)

	open (3,file = arg)
	do i=1,nlines
	read (3,*,end=16) pera
	if (pera .EQ. perb) then
	backspace (3)
	read  (3,*) pera
	close(3)
	end if 
	end do 
	16 close(3)
	
	open (3,file = arg)
	do i=1,nlines
	read (3,*,end=17) anima
	if (anima .EQ. animb) then
	backspace (3)
	read  (3,*) anima
	close(3)
	end if 
	end do 
	17 close(3)
	
	
	if (hara == "harmonic") then
	open (7,file = "INPUT_EQUILIBRUIM_DISTANCE.test")
	do i=1,nnlines
	read (7,*,end=18) equaa
	if (equaa .EQ. equab) then
	backspace (7)
	read  (7,*) equaa
	allocate (yp(N,3))
	allocate (Nuu(N))
	do k = 1,N
	if (space == 3) then 
	read  (7,*) Nuu(k), yp(k,1) , yp(k,2) , yp(k,3) 
	end if 
	end do 
	close(7)
	end if 
	end do  
	18 close (7)
	end if 
	
	if (hara == "harmonic") then
	open (7,file = "INPUT_EQUILIBRUIM_DISTANCE.test")
	do i=1,nnlines
	read (7,*,end=19) fconaa
	if (fconaa .EQ. fconab) then
	backspace (7)
	read  (7,*) fconaa
	allocate(Nuuu(28,2))
	allocate(FC(28))
	do k =1,28
	read  (7,*) Nuuu(k,1) , Nuuu(k,2) , FC(k)
	end do
	close(7)
	end if 
	end do  
	19 close (7)
	end if 

	open (3,file = arg)
	do i=1,nlines
	read (3,*,end=20) showhessiana
	if (showhessiana .EQ. showhessianb) then
	backspace (3)
	read  (3,*) showhessiana
	close(3)
	end if
	end do
	20 close(3)
    
    
    open (3,file = arg)
	do i=1,nlines
	read (3,*,end=21) distancea
	if (distancea .EQ. distanceb) then
	backspace (3)
	read  (3,*) distancea 
	close(3)
	end if 
	end do 
	21 close(3)
    
    
	end subroutine
    
    

end module