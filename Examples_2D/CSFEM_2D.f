
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4     PERIOD)
C     
      INCLUDE 'ABA_PARAM.INC'

       parameter(zero=0.d0, half=0.5, one=1.d0, two=2.d0, three=3.d0, 
     1 four=4.d0, six=6.d0, eight=8.d0, twelve=12.d0,N_element=1,
     2 nsvint=6,ninpt=4)
     
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),
     1     SVARS(NSVARS),ENERGY(8),PROPS(*),COORDS(MCRD,NNODE),
     2     U(NDOFEL),DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),
     3     PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     4     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
     5     JPROPS(*),svint(6,4)

          REAL*8,DIMENSION(:) :: SIG(3),EPS(3)
	  REAL*8,DIMENSION (:,:) ::B(3,NDOFEL),DB(3,NDOFEL),D(3,3),
     1         TB(NDOFEL,3)
	  REAL*8,DIMENSION(:) :: PHI(NNODE),TRI(NNODE,3)  ,CNODES(2,4)
      	  REAL*8,DIMENSION(:) ::VCOORD(2,NNODE+1),X(NNODE),Y(NNODE),
     1       xmid(2)
      	  REAL*8,DIMENSION(:) :: X1(3),X2(3)
      	  REAL*8     AREAP,XC,YC,TEMP1,TEMP2,NX,NY,LSIDE,E1,NU1
	  INTEGER*8, DIMENSION (:,:) :: CTRI(1,3)
          REAL*8,DIMENSION(:) :: FX(NNODE),FY(NNODE)

	  nb = ndofel/2

	  N1 = 3
	  N2 = ndofel

	  ! constitutive matrix...
	  E1 = props(2)
	  nu1 = props(3)	
	
	  ! initialize constitutive matrix...
	  do k1 = 1,n1
	  do k2 = 1,n1
	  D(k2,k2) = 0.
	  enddo
	  enddo
	
	  D1 = E1/(1-nu1*nu1)
	  D(1,1) = D1
	  D(1,2) = D1*nu1
	  D(2,1) = D1*nu1
	  D(2,2) = D1
	  D(3,3) = D1*(1-nu1)*0.5


	  ! initialize the elemental stiffness matrix and right 
	  ! hand side vector
	  do k1=1,ndofel
	  	rhs(k1,NRHS) = 0.
	  	do k2 = 1,ndofel
	  		amatrx(k2,k1) = 0.
	  	enddo
	  enddo


	  ! given an element. first get its coordinates and find its 
	  ! geometric mean. this will be used to triangulate the element. 
	  ! the triangles then form the subcells over which smoothing 
	  ! technique will be employed
	  ! vertex coordinates
	  x = coords(1,:); y = coords(2,:)	
 
	  ! compute area and the geometric center of the current element
	  ! output: A (area), xc,yc (centroid of the subcell)
	  call polygeom(x,y,nb,AREAP,xc,yc)

	  ! initialize the vertex coordinates
	  do i = 1,nb+1
	  	vcoord(1,i) = 0;
	  	vcoord(2,i) = 0;
	  enddo

	  ! update the coordinate list to the existing coordinate list
	  do i = 1,nb
	  	vcoord(1,i) = coords(1,i)
	  	vcoord(2,i) = coords(2,i)
	  enddo
	  vcoord(1,nb+1) = xc
	  vcoord(2,nb+1) = yc

	  ! initialize the triangular subcell connectivity
	  do i = 1,nb
	  	tri(i,1) = 0
	  	tri(i,2) = 0
	  	tri(i,3) = 0
	  enddo

	  ! create triangular subcells...
	  do i = 1,nb
	  	if (i == nb) then
	  		tri(i,1) = i
	  		tri(i,2) = 1
	  		tri(i,3) = nb + 1
	  	else
	  		tri(i,1) = i
	  		tri(i,2) = i+1
	  		tri(i,3) = nb+1
	  	endif
	  enddo

	  !=========================
	  ! loop over the subcells
	  !=========================
          do isub=1,nb
	  	! current sub-cell
	  	ctri(1,:) = tri(isub,:)

	  	! corresponding nodal coordinates
	  	! first initialize the coordinates
	  	do k =1,4
	  		cnodes(1,k) = 0
	  		cnodes(2,k) = 0
	  	enddo

	  	! get the nodal coordinates and add the coordinate
	  	! of the center node to the list
	  	do k=1,3
	  		cnodes(1,k) = vcoord(1,ctri(1,k))
	  		cnodes(2,k) = vcoord(2,ctri(1,k))
	  	enddo
	  	cnodes(1,4) = vcoord(1,ctri(1,1))
	  	cnodes(2,4) = vcoord(2,ctri(1,1))
	  
	  	do k=1,3
	  		x1(k) = 0
	  		x2(k) = 0
	  	enddo
	  
	  	do k = 1,3
	  		x1(k) = cnodes(1,k)
	  		x2(k) = cnodes(2,k)
	  	enddo
		
	  	! compute the area of the subcell
	  	call polygeom(x1,x2,3,Asc,xc,yc)
	  
	  
	  	! initialize the strain-displacement matrix
	  	do ii=1,N1       
	  	do ij=1,ndofel      
	  		B(ii,ij)=0
	  	enddo
	  	enddo

		do ji=1,nb
		fx(ji) = 0
		fy(ji) = 0
		enddo
	
	  	! now loop over the sides of the subcell...
	  	! as we are using triangular subcells, this is always 3
	  	do iside = 1,3
	  		! compute length of the current side
	  		temp1 = (cnodes(1,iside+1)-cnodes(1,iside))**2 
	  		temp2 = (cnodes(2,iside+1)-cnodes(2,iside))**2
	  		lside = sqrt( temp1 + temp2)

	  		! compute the normal of the current side
	  	nx=0.5*(cnodes(2,iside+1)-cnodes(2,iside))/(0.5*lside)
	  	ny=-0.5*(cnodes(1,iside+1)-cnodes(1,iside))/(0.5*lside)
	  	
	  	
	  		! check which side are we integrating on...
	  		if (iside == 1) then
  				! this side is on the boundary. so we use linear 
  				! shape functions. first 
				! initialize the shape function vector...
	  			do k=1,nb
	  				PHI(k) = 0
	  			enddo
	  		
	  			PHI(ctri(1,1)) = 0.5
	  			PHI(ctri(1,2)) = 0.5
	  		else
				! the side is inside the physical element. 
				! we will use wachspress interpolants
				! find the mid-point of the current side....
	  	    	xmid(1)=0.5*(cnodes(1,iside+1)+cnodes(1,iside))
	  	    	xmid(2)=0.5*(cnodes(2,iside+1)+cnodes(2,iside))
	  		
	  			do k=1,nb
	  				PHI(k) = 0
	  			enddo
	  		
	  			! compute wachspress shape function at this point
 	  			call wachspress2d(coords,xmid,nb,PHI)
		
	  		endif

	    do ji = 1,nb
	  	! assemble the equations..
	  	fx(ji) = fx(ji) + nx*PHI(ji)*lside
	  	fy(ji) = fy(ji) + ny*PHI(ji)*lside
	    enddo
	    
	  enddo   ! end loop over sides of the subcell
	  
	  !print*,'printing fx'
	  !print*, 'hehehe'
	  !pause
	  ! now the corrected derivative at the gauss point. note that
	  ! this gauss point is inside the domain, in particular at the
	  ! center of the triangular subcell
	  
	  !***! DEFINE B MATRIX 
	  ! note that B is evaluated at the center of the triangular
	  ! subcell. the coordinates of this are the barycentric
	  ! coordinates of the triangle with weight = Area of the
	  ! triangle
	  do I1=1,nb                         
     	  I2=(I1-1)*2         
	  B(1, 1+I2)= fx(I1)/Asc             
      	  B(2, 2+I2)= fy(I1)/Asc
      	  B(3, 1+I2)= fy(I1)/Asc
      	  B(3, 2+I2)= fx(I1)/Asc
      	  enddo

	  do k2=1,N1
	  	EPS(k2)=0
	    SIG(k2)=0
	  enddo
	   
	  do I=1,N1
	  do J=1,N2	      
		EPS(I)=EPS(I)+B(I, J)*U(J)  	    
	  enddo
	  enddo
	   
	  do I=1,N1
	  	do J=1,N1
	  		SIG(I)=SIG(I)+D(I, J)*EPS(J)        
	  	enddo
	  enddo
	  
	  ! [D]*[B] initialization 
      do I=1,N1
      	do J=1,N2
      		DB(I,J)=0
      	enddo
      enddo

	  ! calculate the value of [D][B] 
      do I=1,N1
      	do J=1,N2
	  	do K2=1,N1
      		DB(I,J)=DB(I,J)+D(I, K2)*B(K2,J)
	  	enddo
	  	enddo
	  enddo
	  
	  ! find the transpose of the strain-displacement matrix
	  do j=1,N1
	  	do i=1,NDOFEL
	  		TB(i,j)=B(j,i)
	  	enddo
	  enddo


       !!!!!!!!!!!print*,'B',B
       !!!!!!!!!!!pause
      ! The residual vector...
  	  do K1=1, N2                 
	    do K4=1,N1
	    	RHS(K1, 1)=RHS(K1, 1)-Asc*B(K4, K1)*SIG(K4)
  	     enddo
  	  enddo

	  !***Compute elemental stiffness matrix
  	  do K1=1,NDOFEL              
	  	do K2=1,NDOFEL
      	         do K3=1,N1       
	    	AMATRX(K1, k2)=AMATRX(K1, k2)+Asc*TB(K1, k3)*DB(k3,k2)
	  	enddo
	  	enddo
	  enddo
          !!!!!!!!!!!print*,'amatrx',amatrx(1,:)
          !!!!!!!!pause
      enddo	! end loop over subcells	                   
      return
      end

!========================================================
!		Wachspress interpolants for
! 		arbitrary polygons
!========================================================

	  subroutine wachspress2d(coords,xpt,n,phi)
      
      implicit none
      
      integer :: im1,modma
      integer :: n,i,j,k
      real*8, dimension(2) :: xpt
      real*8 :: h,wsum,FindDet
      real*8, dimension (2,n) ::  coords
      real*8, dimension(:,:) ::R(n,2),dphi(n,2),un(n,2),p(n,2),mat3(2,2)
      real*8, dimension (:)  :: d(2),w(n),phi(n),mat1(2),mat2(2),phiR(2)


	  !initialize the vectors
      w=0.0; R=0.0d0
      phi=0.0; dphi=0.0; phiR = 0.0d0
      un=0.0; d=0.0
      
      ! compute the unit normal vector for each of the sides
      do i=1,n
        call modinfortran(i,n,modma)
        d=coords(:,modma+1)-coords(:,i)
        un(i,:) = [d(2),-d(1)]/norm2(d)
      enddo
      
      ! compute p = un/h
      do i=1,n
      	h=dot_product(coords(:,i)-xpt,un(i,:))
        p(i,:)=un(i,:)/h
      enddo

      do i=1,n
        call modinfortran(i-2,n,modma)
        im1=modma+1 
      	mat1 =p(im1,:)
      	mat2=p(i,:)
      	do j=1,size(mat1)
        	mat3(1,j)=mat1(j)
        	mat3(2,j)=mat2(j)
        enddo
        
        
      	w(i)=mat3(1,1)*mat3(2,2)-mat3(2,1)*mat3(1,2)
      	R(i,:)=p(im1,:) + p(i,:)
      enddo
      
      ! the shape functions
      wsum=sum(w)
      phi=w/wsum
      
      
      do i=1,size(phi)
      	do j=1,size(R,2)
        phiR(j)=phiR(j)+phi(i)*R(i,j)
        enddo
      enddo

      do k=1,2
      do i=1,size(phi)
      	dphi(i,k)=dphi(i,k)+phi(i)*(R(i,k)-phiR(k))
      enddo
      enddo

      end subroutine wachspress2d



      subroutine modinfortran(a,b,modma)
      implicit none
   
      integer :: a,b
      integer :: modma
      
      if (a*b .lt. 0.0) then
      	modma=mod(a,b)+b
      else
     	modma=mod(a,b)
      endif
      
      end subroutine modinfortran 

!========================================================
!		Compute geometric related information
!	Areap - area of the polygon
!	xc, yc - centroid of the polygon
!========================================================

	  subroutine polygeom(x,y,nb,AREAP,xc,yc)
	  implicit none

	  real*8,dimension(:) :: x(nb),y(nb),dx(nb),dy(nb), xt(nb)
          real*8,dimension(:) :: yt(nb)
	  real*8 :: xm, ym, AREAP, Axc, Ayc
	  real*8 :: xc, yc
	  integer :: nb, i


	  ! initialize the vectors...
	  do i = 1,nb
	  	dx(i) = 0; dy(i) = 0
		xt(i) = 0; yt(i) = 0
	  enddo
	
	  ! initialize the values
	  xm = 0;	ym = 0
	  xc = 0; yc = 0
	  AREAP = 0; Axc= 0; Ayc = 0

	  ! find the mean of the vertices 
	  do i = 1,nb
	  	xm = xm + x(i)
		ym = ym + y(i)
	  enddo

	  xm = xm/nb
	  ym = ym/nb
	

	  ! temporarily shift data to the mean of vertices for improved '
      ! accuracy
	  do i = 1,nb
	  	x(i) = x(i) - xm
		y(i) = y(i) - ym
	  enddo

	  ! re-arrange the vertices...
	  do i = 1,nb-1
		xt(i) = x(i+1)
		yt(i) = y(i+1)
	  enddo
	  xt(nb) = x(1)
	  yt(nb) = y(1)
	

	  ! find the Euclidean distance between the vertices
	  do i = 1,nb-1
		dx(i) = xt(i) - x(i)
		dy(i) = yt(i) - y(i)
	  enddo
	  dx(nb) = xt(nb) - x(nb)
	  dy(nb) = yt(nb) - y(nb)
	
	  ! summation of CW boundary integrals
	  do i = 1,nb
	  AREAP = AREAP + y(i)*dx(i) - x(i)*dy(i)
		
	  Axc=Axc+6*x(i)*y(i)*dx(i)-3*x(i)*x(i)*dy(i)+3*y(i)*dx(i)*dx(i)
     1      + dx(i)*dx(i)*dy(i)
		
	  Ayc=Ayc+3*y(i)*y(i)*dx(i)-6*x(i)*y(i)*dy(i)-3*x(i)*dy(i)*dy(i)
     1      -  dx(i)*dy(i)*dy(i)
	  enddo

	  AREAP = AREAP/2
	  Axc = Axc/12
	  Ayc = Ayc/12
	
	  ! check for CCW vs VCW
	  if (AREAP < 0) then
	  	AREAP = -AREAP
		Axc = -Axc
		Ayc = -Ayc
	  endif
	  
	  ! centroidal moments
	  xc = Axc/AREAP
	  yc = Ayc/AREAP
	
	  ! replace mean of vertices
	  xc = xc + xm
	  yc = yc + ym
      
      return
      end subroutine polygeom
