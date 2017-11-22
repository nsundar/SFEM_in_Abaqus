      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4     PERIOD)
C     
      INCLUDE 'ABA_PARAM.INC'
 
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),
     1     SVARS(NSVARS),ENERGY(8),PROPS(*),COORDS(MCRD,NNODE),
     2     U(NDOFEL),DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),
     3     PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     4     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
     5     JPROPS(*),svint(6,4)
      integer*8 :: element(1,size(coords,2)), numelem, numnode, ndof=3,
     &  totalUnknown
      integer*8 :: iel, nn, i,j,fcnt1,fcnt2,fcnt,iface,pt,nc
      real*8 :: E,nu,D(6,6),cp_vcoord(1,3),xpt(1,3),cp_facecoord(1,3)
      real*8,dimension(:,:),allocatable::kmat, cface_coord,ncoord
      integer*8,dimension(:,:),allocatable::surfplane,g
      integer*8,dimension(:),allocatable::gindex,cface,econ
      integer*8,allocatable:: F_v(:,:), F_h(:,:), G_v(:,:), 
     &   G_h(:,:),itr(:),tri(:,:),aFaces(:,:)
      real*8,dimension(:),allocatable::xvc,yvc,zvc
      real*8,dimension(:,:),allocatable :: aphi,adphi,un
      real*8 ::xsum,ysum,zsum,ctri_coord(3,3),tet_coord(4,3),ws,
     &   node(size(coords,2),3)
      integer*8 ::itri,ctri(3),tetr(1,4),tet_face(4,3),lcface
      real*8,dimension(:,:),allocatable::dNdx(:,:),dNdy(:,:),
     &   dNdz(:,:),Ba(:,:),tri_coord(:,:)	

!     Material Properties
      E=props(2)
      nu=props(3)

      nn=size(node,1)
      if(.not.allocated(econ))allocate(econ(nn))
      if(.not.allocated(g))allocate(g(nn,3))
      if(.not.allocated(ncoord))allocate(ncoord(nn,3))
      if(.not.allocated(aphi))allocate(aphi(nn,1))
      if(.not.allocated(adphi))allocate(adphi(nn,3))
      if(.not.allocated(un))allocate(un((nn/2)+2,3))
      if(.not.allocated(kmat)) allocate(kmat(3*nn,3*nn));kmat=0.0d0
      ncoord=transpose(coords)

!     computing material matrix
      call getmatmtx3d(E,nu,D)

      numelem=size(element,1)
      numnode=size(ncoord,1)
      totalUnknown=ndof*numnode
      ndof=3

      do iel=1,1
      call mean(ncoord,nn,cp_vcoord)
      allocate (gindex(nn*ndof))
      allocate (F_v(2,nn/2))
      allocate (F_h(nn/2,4))
      allocate (aFaces(2+nn/2,10))
      aFaces=0
!     computing local face connectivity
      call LocalFaceConn(nn,F_v,F_h)
      aFaces(1:2,1:size(F_v,2))=F_v(1:2,:)
      aFaces(3:2+nn/2,1:size(F_h,2))=F_h(1:nn/2,:)
      allocate (G_v(2,nn/2))
      allocate (G_h(nn/2,4))
      fcnt1=size(F_v,1)
      fcnt2=size(F_h,1)
      fcnt=fcnt1+fcnt2
!     computing face normals
      call getunitnorm(ncoord,aFaces,nn,un,g)
      if(.not.allocated(surfplane))allocate(surfplane(size(un,1),2))
      call surplane(un,nn,surfplane)

!     Formation of the smoothing cells
!     loop over the element faces
      pt=1
      do iface=1,fcnt
        if (iface.eq.1) then
       allocate (cface(size(F_v,2)))
        cface=F_v(iface,:)
       allocate (cface_coord(size(cface,1),3))
       cface_coord=ncoord(cface,:)
       else if (iface.eq.2) then
       allocate (cface(size(F_v,2)))
       cface=F_v(iface,:)
       allocate (cface_coord(size(cface,1),3))
       cface_coord=ncoord(cface,:)
       else
       allocate (cface(size(F_h,2)))
       cface=F_h(pt,:)
       allocate (cface_coord(size(cface,1),3))
       cface_coord=ncoord(cface,:)
       pt=pt+1
       end if !iface
       lcface=size(cface,1)

       if(.not.allocated(xvc)) allocate(xvc(size(cface_coord,1)))
       if(.not.allocated(yvc)) allocate(yvc(size(cface_coord,1)))
       if(.not.allocated(zvc)) allocate(zvc(size(cface_coord,1)))
       xvc=0.0d0;yvc=0.0d0;zvc=0.0d0
       xvc=cface_coord(:,1)
       yvc=cface_coord(:,2)
       zvc=cface_coord(:,3)

       xsum=sum(xvc)/size(cface_coord,1)
       ysum=sum(yvc)/size(cface_coord,1)
       zsum=sum(zvc)/size(cface_coord,1)
       cp_facecoord(1,1)=xsum
       cp_facecoord(1,2)=ysum
       cp_facecoord(1,3)=zsum

       if(allocated(xvc)) deallocate(xvc)
       if(allocated(yvc)) deallocate(yvc)
       if(allocated(zvc)) deallocate(zvc)
 
      nc=size(cface,1)
      allocate (tri(nc,3))
      call get_tri(nc,tri)
      if ((iface==1) .or. (iface==2))then
      allocate(tri_coord(size(F_v,2)+1,3))
      do i=1,size(F_v,2)
       tri_coord(i,:)=cface_coord(i,:)
      enddo
      tri_coord(size(F_v,2)+1,:)=cp_facecoord(1,:)
      else
      allocate(tri_coord(size(F_h,2)+1,3))
      do i=1,size(F_h,2)
      tri_coord(i,:)=cface_coord(i,:)
      enddo
      tri_coord(size(F_h,2)+1,:)=cp_facecoord(1,:)
      endif

!     Loop over the triangle to form the tetrahedrons smoothing cells
      do itri=1,size(tri,1)
      ctri=tri(itri,:)
      ctri_coord=tri_coord(ctri,:)
      tetr(1,:)=(/1, 2, 3, 4/)
      tet_coord(1,:)= ctri_coord(1,:)
      tet_coord(2,:)= ctri_coord(2,:)
      tet_coord(3,:)= ctri_coord(3,:)
      tet_coord(4,:)= cp_vcoord(1,:)
      tet_face=reshape((/3,1,4,3,1,4,3,4,2,2,2,1/),(/4,3/))
      if (.not.allocated(dNdx)) allocate(dNdx(1,nn))
      if (.not.allocated(dNdy)) allocate(dNdy(1,nn))
      if (.not.allocated(dNdz)) allocate(dNdz(1,nn))

!     compute the smoothed/modified derivative
      call getsmoothderivative(ncoord,nn,tet_coord,tet_face,cface,
     &      lcface,cface_coord,un,g,iface,surfplane,ws,dNdx,dNdy,dNdz)

      if(.not.allocated(Ba)) allocate(Ba(6,3*nn))
      do i=1,1
       Ba=0.0d0
       Ba(1,1:3*nn:3)=dNdx(i,:)
       Ba(2,2:3*nn:3)=dNdy(i,:)
       Ba(3,3:3*nn:3)=dNdz(i,:)
       Ba(4,1:3*nn:3)=dNdy(i,:)
       Ba(4,2:3*nn:3)=dNdx(i,:)
       Ba(5,2:3*nn:3)=dNdz(i,:)
       Ba(5,3:3*nn:3)=dNdy(i,:)
       Ba(6,1:3*nn:3)=dNdz(i,:)
       Ba(6,3:3*nn:3)=dNdx(i,:)
       kmat=kmat+ws*matmul(transpose(Ba),matmul(D,Ba))
       end do !i

      if(allocated(Ba)) deallocate(Ba)
      end do !itri
      deallocate(cface_coord)
      deallocate(cface)
      deallocate(tri)
      deallocate(tri_coord)
      enddo !iface
      enddo !iel

       amatrx(1:ndofel,1:ndofel)=kmat
       RHS(1:ndofel,1)=RHS(1:ndofel,1)-matmul(amatrx,U)
      if(allocated(kmat)) deallocate(kmat)
      end
!        END OF UEL MAIN SUBROUTINE
!=======================================================================
! *********    DEPENDENCIES STARTS FROM HERE    **********
!=======================================================================

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
!      Following subroutine gives material matrix     
       subroutine getmatmtx3d(E,nu,D)
       implicit none
       integer*8:: k1,k2,n1
       real*8 :: E,nu,D(6,6),L, mu
       L=E*nu/((1+nu)*(1-2*nu))
       mu=E/(2*(1+nu))
       n1=6
       ! initialize constitutive matrix...
       do k1=1,n1
       do k2=1,n1
        D(k1,k2)=0.0
       enddo !k1
       enddo !k2
       D(1,1)=L+2*mu
       D(1,2)=L 
       D(1,3)=L
       D(2,1)=L
       D(2,2)=L+2*mu
       D(2,3)=L
       D(3,1)=L
       D(3,2)=L
       D(3,3)=L+2*mu
       D(4,4)=mu
       D(5,5)=mu
       D(6,6)=mu
      end subroutine getmatmtx3d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
!     Following subroutine gives mean of the coordinates 
      subroutine mean(ncoord,nn,cp_vcoord)
       implicit none
       integer*8 ::nn
       real*8 ::ncoord(nn,3)
       real*8::cp_vcoord(1,3), xv(nn),yv(nn),zv(nn),xsum,ysum,zsum
       xv=ncoord(:,1)
       yv=ncoord(:,2)
       zv=ncoord(:,3)
       xsum=sum(xv)/size(ncoord,1)
       ysum=sum(yv)/size(ncoord,1)
       zsum=sum(zv)/size(ncoord,1)
       cp_vcoord(1,1)=xsum
       cp_vcoord(1,2)=ysum
       cp_vcoord(1,3)=zsum
      end subroutine mean

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
!      Following subroutine gives local face connectivity
       subroutine LocalFaceConn(nn,F_v,F_h)
       implicit none
          integer*8 nn, n, F_v(2,nn/2),F_h(nn/2,4)
         n = nn/2
        select case(n)
        case(3)
        F_v=reshape((/3,4,2,5,1,6/),(/2,3/))
        F_h=reshape((/1,3,3,2,2,1,5,5,4,4,6,6/),(/3,4/))

        case(4)
        F_v=reshape((/4,5,3,6,2,7,1,8/),(/2,4/))
        F_h=reshape((/1,2,3,1,2,3,4,5,6,7,8,8,5,6,7,4/),(/4,4/))

        case(5)
        F_v=reshape((/5,6,4,7,3,8,2,9,1,10/),(/2,5/))
        F_h=reshape((/1,2,3,4,1,2,3,4,5,6,7,8,9,10,10,6,7,8,9,5/),
     &       (/5,4/))

        case(6)
        F_v=reshape((/6,7,5,8,4,9,3,10,2,11,1,12/),(/2,6/))
        F_h=reshape((/1,2,3,4,5,1,2,3,4,5,6,7,8,9,10,11,12,12,7,8,9,
     &       10,11,6/),(/6,4/))

        case(7)
        F_v=reshape((/7,8,6,9,5,10,4,11,3,12,2,11,1,14/),(/2,7/))
        F_h=reshape((/1,2,3,4,5,6,7,2,3,4,5,6,7,1,9,10,11,12,13,14,8,8,
     &      9,10,11,12,13,14/),(/7,4/))

        case(8)
        F_v=reshape((/8,9,7,10,6,11,5,12,4,13,3,14,2,15,1,16/),(/2,8/))
        F_h=reshape((/1,2,3,4,5,6,7,8,2,3,4,5,6,7,8,1,10,11,12,13,14,
     &      15,16,9,9,10,11,12,13,14,15,16/),(/8,4/))

        case(9)
        F_v=reshape((/9,10,8,11,7,12,6,13,5,14,4,15,3,16,2,17,1,
     &       18/),(/2,9/))
        F_h=reshape((/1,2,3,4,5,6,7,8,9,2,3,4,5,6,7,8,9,1,11,12,13,14,
     &       15,16,17,18,10,10,11,12,13,14,15,16,17,18/),(/9,4/))
        end select
        end subroutine LocalFaceConn

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
!      Following subroutine gives face/edge normals 
        subroutine getunitnorm(anodef,aFaces,nn,aun,ag)
        implicit none
        integer*8::nn
        real*8,dimension(:,:) :: apc(1,3),anodef(nn,3)
        integer*8,dimension(:,:) :: aFaces((nn/2)+2,10)
        integer*8,dimension(:,:),allocatable :: afacecon
        real*8,dimension(:,:),allocatable :: anodefaces,atn
        real*8,dimension(:) :: crossproduct(3),an(1,3)
        integer :: iface,inn 
        integer,dimension(:) :: dummy(1)
        integer,dimension(:) :: sz2afacecon(1)
        real*8 :: aun((nn/2)+2,3)
        integer*8 ::ag(nn,3)
        integer*8,dimension (:,:),allocatable :: atemp1,atemp,acf
        integer :: atempcount =0
 
        apc(1,:)=(/sum(anodef(:,1))/size(anodef,1),
     & sum(anodef(:,2))/size(anodef,1),sum(anodef(:,3))/size(anodef,1)/)

       do iface=1,size(aFaces,1)
        !if (.not.allocated(afacecon)) allocate(afacecon(1,4))
        !afacecon(1,:)=aFaces(iface,1:4)
        sz2afacecon=minloc(aFaces(iface,:))-1
        !print*,'sz2afacecon',sz2afacecon(1)
        if (.not.allocated(afacecon)) allocate(afacecon(1,sz2afacecon
     &       (1)))
        !print*,'sz1',size(afacecon,2);pause
        afacecon(1,1:size(afacecon,2))=aFaces(iface,1:sz2afacecon(1))
        if (.not.allocated(anodefaces)) allocate(anodefaces(
     &    size(afacecon,2),3))
        if (.not.allocated(atn)) allocate(atn(size(afacecon,2),3)) 
        anodefaces=anodef(afacecon(1,:),:)

        atn(:,1)=anodefaces(:,1)-apc(1,1)
        atn(:,2)=anodefaces(:,2)-apc(1,2)
        atn(:,3)=anodefaces(:,3)-apc(1,3) 

        call cross(atn(2,:)-atn(1,:),atn(3,:)-atn(1,:),crossproduct)
        an(1,:)=crossproduct/(crossproduct(1)**2+crossproduct(2)**2+
     &     crossproduct(3)**2)**0.5

        if(dot_product(an(1,:),atn(1,:)).lt.0.0) then
           afacecon(:,1:size(afacecon,2):1)=afacecon(:,size(afacecon,
     &      2):1:-1)
           anodefaces=anodef(afacecon(1,:),:)
           atn(:,1)=anodefaces(:,1)-apc(1,1)
           atn(:,2)=anodefaces(:,2)-apc(1,2)
           atn(:,3)=anodefaces(:,3)-apc(1,3)
           call cross(atn(2,:)-atn(1,:),atn(3,:)-atn(1,:),crossproduct)
           an(1,:)=crossproduct/(crossproduct(1)**2+crossproduct(2)**2+
     &      crossproduct(3)**2)**0.5
        end if
        aun(iface,:)=an(1,:)
        if (allocated(atn)) deallocate(atn)
        if (allocated(anodefaces)) deallocate(anodefaces)
        if (allocated(afacecon)) deallocate(afacecon)
       end do ! iface
  
        do inn=1,size(anodef,1)
         if (.not.allocated(atemp1)) allocate(atemp1(1,size(aFaces,1)))
            atemp1=0
            do iface=1,size(aFaces,1)
         if (.not.allocated(acf)) allocate(acf(1,size(aFaces(iface,:))))
            acf(1,:)=aFaces(iface,:)
            if(any(acf(1,:)==inn))then
                atempcount=atempcount+1
                atemp1(1,atempcount)=iface
            end if
           if (allocated(acf)) deallocate(acf)
            end do ! iface
            if (.not.allocated(atemp)) allocate(atemp(1,atempcount))
            atemp(1,1:atempcount)=atemp1(1,1:atempcount)
            ag(inn,:)=atemp(1,:)
            if (allocated(atemp1)) deallocate(atemp1)
            if (allocated(atemp)) deallocate(atemp)
            atempcount=0
        end do !inn
      end subroutine getunitnorm

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
!      Following subroutine gives cross product 
      subroutine cross(a, b,crossproduct)
        real*8, DIMENSION(3) :: crossproduct
        real*8, DIMENSION(3) :: a, b

        crossproduct(1) = a(2) * b(3) - a(3) * b(2)
        crossproduct(2) = a(3) * b(1) - a(1) * b(3)
        crossproduct(3) = a(1) * b(2) - a(2) * b(1)
      end subroutine cross

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
!      Following subroutine gives plane information
      subroutine surplane(un,nn,surfplane)
        implicit none
        integer*8 :: i, j, k,n,cnt,q,nn
        real*8 :: un((nn/2)+2,3),uni(3), Ma,Mi
        integer*8,dimension(:,:) :: surfplane(size(un,1),2)
        integer*8, dimension(:,:), allocatable::m

        allocate(m(7,2))
        do i=1,size(un,1)
          uni=un(i,:)
          Ma=maxval(abs(uni))
          Mi=minval(abs(uni))
         cnt=1
         do j=1,size(uni,1)
         if (abs(uni(j)) .eq. Mi) then
            n=j
            m(i,cnt)=n
            cnt=cnt+1
         endif
         if ((abs(uni(j))>Mi) .and. (abs(uni(j))<Ma)) then
            q=j
            m(i,cnt)=q
            cnt=cnt+1
         endif
         surfplane=m
        enddo !j
        enddo !i
        end subroutine surplane

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
!      Following subroutine gives polyhedral shape function 
        subroutine wachspress3d(ncoord,nn,un,g,xpt,aphi,adphi)
         implicit none

          integer*8 :: an,ak,nn,g(nn,3)
          real*8 :: xpt(1,3),ncoord(nn,3)
          real*8, dimension(:,:),allocatable :: aw,ar,ap,apcopy,awloc,
     &           arloc,aphir,temp
          integer*8,dimension(:,:),allocatable :: af
          real*8 :: ah,FindDet,awsum,aphiR1,un((nn/2)+2,3)
          integer :: i,j,d
          real*8,dimension(:,:) :: aphi(nn,1),adphi(nn,3)

          an=size(ncoord,1)
          if(.not.allocated(aw)) allocate(aw(an,1));aw=0.0
          if(.not.allocated(ar)) allocate(ar(an,3));ar=0.0
          do i=1,an
            ak=size(g(i,:),1)
            if(.not.allocated(af)) allocate(af(1,size(g(i,:),1)));
            af(1,:)=g(i,:)
            if(.not.allocated(ap)) allocate(ap(ak,3));
            ap=0.0
            if(.not.allocated(apcopy)) allocate(apcopy(ak,3));
                 do j=1,ak
                  ah=dot_product(ncoord(i,:)-xpt(1,:),un(af(1,j),:))
                  ap(j,:)=un(af(1,j),:)/ah;
                 end do ! for j
                apcopy=ap;
            call determinant(apcopy,size(ap,1),FindDet);apcopy=ap;
            if (FindDet.lt.0) then
             af(1,:)=af(1,size(af,2):1:-1) 
            ak=size(af,2)
            do j=1,ak
                 ah=dot_product(ncoord(i,:)-xpt(1,:),un(af(1,j),:))
                ap(j,:)=un(af(1,j),:)/ah;
            end do ! for j
           end if
          if(allocated(af)) deallocate(af)
          if(.not.allocated(awloc)) allocate(awloc(ak-2,1));awloc=0.0
          if(.not.allocated(arloc)) allocate(arloc(ak-2,3));arloc=0.0
        do j=1,ak-2
          if (.not.allocated(temp)) allocate(temp(3,3));
          temp(1,:)=ap(j,:);temp(2,:)=ap(j+1,:);temp(3,:)=ap(ak,:);
         call determinant(temp,3,FindDet)
         awloc(j,1)=FindDet
         if (allocated(temp)) deallocate(temp);
         arloc(j,:)=ap(j,:)+ap(j+1,:)+ap(ak,:)
       end do ! j
       if(allocated(ap)) deallocate(ap)
       aw(i,1)=sum(awloc)
       if (.not.allocated(temp)) allocate(temp(1,size(arloc,2)))
       temp=matmul(transpose(awloc),arloc)
       if(allocated(awloc)) deallocate(awloc)
       if(allocated(arloc)) deallocate(arloc)
       ar(i,:)=temp(1,:)/aw(i,1);
       if (allocated(temp)) deallocate(temp);
       end do ! i

      awsum=sum(aw)
      aphi=aw/awsum
      if(.not.allocated(aphir)) allocate(aphir(size(aphi,2),size(ar,2)))
      aphir=matmul(transpose(aphi),ar)

      do d=1,3
           aphiR1=aphiR(1,d)
        do i=1,size(adphi,1)
            adphi(i,d) = aphi(i,1)*(ar(i,d)-aphiR1)
        end do ! i
      end do

      if(allocated(aphir)) deallocate(aphir)
      if(allocated(af)) deallocate(af)
      end subroutine wachspress3d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
!      Following subroutine gives determinant 
      Subroutine determinant(matrix, n,FindDet)
        IMPLICIT NONE

       REAL*8, DIMENSION(n,n) :: matrix
       INTEGER*4, INTENT(IN) :: n
       REAL*8 :: m, temp,FindDet
       INTEGER*8 :: i, j, k, l
       LOGICAL :: DetExists = .TRUE.
          l = 1
        !Convert to upper triangular form
        DO k = 1, n-1
             IF (matrix(k,k) == 0) THEN
            DetExists = .FALSE.
            DO i = k+1, n
                IF (matrix(i,k) /= 0) THEN
                    DO j = 1, n
                        temp = matrix(i,j)
                        matrix(i,j)= matrix(k,j)
                        matrix(k,j) = temp
                    END DO

                    DetExists = .TRUE.
                    l=-l
                    EXIT
                ENDIF
            END DO

            IF (DetExists .EQV. .FALSE.) THEN
                FindDet = 0
                return
            END IF
           ENDIF

        DO j = k+1, n
            m = matrix(j,k)/matrix(k,k)
            DO i = k+1, n
                matrix(j,i) = matrix(j,i) - m*matrix(k,i)
            END DO
         END DO
         END DO
    !Calculate determinant by finding product of diagonal elements
       FindDet = l
       DO i = 1, n
        FindDet = FindDet * matrix(i,i)
       END DO
      END SUBROUTINE determinant

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
!      Following subroutine gives norms
      subroutine norm(a, normans)
       real*8, DIMENSION(3) :: a
       real*8 :: normans
       normans=norm2(a) 
      end subroutine norm

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
!      Following subroutine gives tetrahedral shape function 
      subroutine H4shapefunction(point,N,dNdxi)
       implicit none
       real*8 :: N(4),dNdxi(4,3),point(3)
       N=(/1-point(1)-point(2)-point(3), point(1), point(2), point(3)/)
       dNdxi=reshape((/-1.0, 1.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 
     &        0.0, 0.0, 1.0/),(/4,3/))
       end subroutine H4shapefunction 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
!     Following function gives modulus
      integer function mod_ma(a,b)
         implicit none
         integer :: a,b
         if (a*b .lt. 0.0) then
         mod_ma=mod(a,b)+b
        else
        mod_ma=mod(a,b)
        end if
       end function mod_ma

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
!      Following subroutine gives polygonal shape function 
      subroutine wachspress2d(v,x,n,phi)
       integer*4 :: n,im1,mod_ma
       real*8, dimension(2) :: x
       real*8 :: h,wsum,phi(n),FindDet2
       real*8, dimension (n,2) ::  v
       real*8, dimension (:,:),ALLOCATABLE :: R,dphi,un,p,mat3
       real*8, dimension (:),ALLOCATABLE :: d,w,mat1,mat2,phiR

        allocate (R(n,2),w(n),dphi(n,2),un(n,2),mat3(2,2),phiR(2),d(2))
        w=0.0D0
        R=0.0D0
        phi=0.0D0
        dphi=0.0D0
        un=0.0D0

        do i=1,n
           d=v(mod_ma(i,n)+1,:)-v(i,:)
           un(i,:) = [d(2),-d(1)]/norm2(d)
        end do

        allocate (p(n,2))
        do i=1,n
         h=dot_product(v(i,:)-x,un(i,:))
         p(i,:)=un(i,:)/h
        end do

       allocate (mat1(2))
       allocate (mat2(2))
       do i=1,n
         im1=mod_ma(i-2,n)+1
         mat1 =p(im1,:)
         mat2=p(i,:)
        do j=1,size(mat1)
         mat3(1,j)=mat1(j)
         mat3(2,j)=mat2(j)
        end do
        call determinant(mat3,size(mat3,1),FindDet2)
         w(i)=FindDet2
         R(i,:)=p(im1,:) + p(i,:)
        end do
       wsum=sum(w)
       phi=w/wsum
       phiR=0.0D0
       do i=1,size(phi)
         do j=1,size(R,2)
           phiR(j)=phiR(j)+phi(i)*R(i,j)
           end do
       end do

      do k=1,2
      do i=1,size(phi)
          dphi(i,k)=dphi(i,k)+phi(i)*(R(i,k)-phiR(k))
         end do
       end do
       end subroutine

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
!     Following subroutine gives face normal
      subroutine getUnitNormg(anodef,aFaces,anodefrow,anodefcolumn,
     &        aFacesrow,aFacescolumn,aun)
         implicit none

        integer :: anodefrow,anodefcolumn,aFacesrow,aFacescolumn
        integer*8,dimension(aFacesrow,aFacescolumn) :: aFaces
        real*8,dimension(anodefrow,anodefcolumn) :: anodef
        real*8,dimension(aFacesrow,3) :: aun
       real*8,dimension(:,:) :: apc(1,3)
       integer*8,dimension(:,:),allocatable :: afacecon
       real*8,dimension(:,:),allocatable :: anodefaces,atn
       real*8,dimension(:) :: crossproduct(3),an(1,3)
      integer :: iface,inn 
      integer,dimension(:) :: dummy(1)
      integer,dimension(:) :: sz2afacecon(1)

      apc(1,:)=(/sum(anodef(:,1))/size(anodef,1),sum(anodef(:,2))
     &        /size(anodef,1),sum(anodef(:,3))/size(anodef,1)/)

       do iface=1,size(aFaces,1)
        sz2afacecon=3
        if (.not.allocated(afacecon)) allocate(afacecon(1,
     &        sz2afacecon(1)))
        afacecon(1,1:size(afacecon,2))=aFaces(iface,1:sz2afacecon(1))
        if (.not.allocated(anodefaces)) allocate(anodefaces(
     &        size(afacecon,2),3))
        if (.not.allocated(atn)) allocate(atn(size(afacecon,2),3)) 
        anodefaces=anodef(afacecon(1,:),:)
        atn(:,1)=anodefaces(:,1)-apc(1,1)
        atn(:,2)=anodefaces(:,2)-apc(1,2)
        atn(:,3)=anodefaces(:,3)-apc(1,3) 
        call cross(atn(2,:)-atn(1,:),atn(3,:)-atn(1,:),crossproduct)
        an(1,:)=crossproduct/(crossproduct(1)**2+crossproduct(2)**2+
     &         crossproduct(3)**2)**0.5

        if(dot_product(an(1,:),atn(1,:)).lt.0.0) then
           afacecon(:,1:size(afacecon,2):1)=afacecon(:,size(afacecon,2)
     &         :1:-1)
           anodefaces=anodef(afacecon(1,:),:)
           atn(:,1)=anodefaces(:,1)-apc(1,1)
           atn(:,2)=anodefaces(:,2)-apc(1,2)
           atn(:,3)=anodefaces(:,3)-apc(1,3)
           call cross(atn(2,:)-atn(1,:),atn(3,:)-atn(1,:),crossproduct)
           an(1,:)=crossproduct/(crossproduct(1)**2+crossproduct(2)**2+
     &         crossproduct(3)**2)**0.5
        end if
        aun(iface,:)=an(1,:)
       end do ! iface
       end subroutine getUnitNormg

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
!     Following subroutine gives modified derivatives
       subroutine getsmoothderivative(ncoord,nn,tet_coord,tet_face,
     &  cface,lcface,cface_coord,un,g,iface,surfplane,ws,dNdx,dNdy,dNdz)
        implicit none

        integer*8 :: g(nn,3),iface,tet_face(4,3),nn,lcface,cface(lcface)
     &            ,pt,i,j,surfplane((nn/2)+2,2)
        real*8 :: ncoord(nn,3),tet_coord(4,3),ws,cface_coord(lcface,3),
     &             dNdx(1,nn),dNdy(1,nn),dNdz(1,nn)
        real*8 :: w,q(3),vol,coord(4,4),tetcont(3),det,N(4),dNdxi(4,3),
     &       detJ,getJ(3,3),qpt(3),wpt,We,un((nn/2)+2,3)
        real*8 :: v(size(tet_coord,1),size(tet_coord,2)),area,
     &       crossans(3),wt
        integer*8 :: f(size(tet_face,1),size(tet_face,2)),
     &    localsurfconn(size(tet_face,1),size(tet_face,2))
        real*8 :: normals(size(tet_face,1),3),enormal(1,3),
     &   ax1(1,size(ncoord,1))
        real*8 :: fx(1,size(ncoord,1)),fy(1,size(ncoord,1)),
     &      fz(1,size(ncoord,1))
        integer*8 :: plane(1,size(surfplane,2))
        real*8 :: CFnode(size(cface_coord,1),2),xiet(1,2)
        integer*8 :: numver,sztet_coord,ismember
        integer*8, dimension(:,:), allocatable :: eside
        real*8,dimension(:,:),allocatable :: cn
        real*8 :: q1(3,2),xi_pt(1,2),N1d(3,1),xieta_gp(1,3)
        real*8 :: Nwach(nn,1),abc(nn,3),phi(size(CFnode,1))
        integer :: iside,igp
        integer :: gsz1,gsz2
        Nwach=0.0d0
!       weight and gauss point inside a subcell(tetrahedral-CS)
        w=0.166666666666667
        q=(/0.25, 0.25, 0.25/)
        do i=1,4
         tetcont=tet_coord(i,:)
          do j=1,3
           coord(i,j)=tetcont(j)
         enddo !j
        coord(i,4)=1.0
       enddo !i
        call determinant(coord,4,det)
         vol=det/6
        pt=1
    ! map the point into the physical space
        call H4shapefunction(q,N,dNdxi)
        getJ=matmul(transpose(dNdxi),tet_coord)
       call determinant(getJ,3,detJ)
       qpt=matmul(N,tet_coord)
       wpt=w*detJ
       We=wpt
       ws=We
       numver=size(tet_face,1)
       v=tet_coord
       f=tet_face(:,size(tet_face,2):1:-1)
       call getUnitNormg(v,f,size(v,1),size(v,2),size(f,1),size(f,2),
     &           normals);
       fx=0.0d0;fy=0.0d0;fz=0.0d0
       q1(1,1)=0.666666666666667;q1(1,2)=0.166666666666667;
       q1(2,1)=0.166666666666667;q1(2,2)=0.166666666666667;
       q1(3,1)=0.166666666666667;q1(3,2)=0.666666666666667;
       localsurfconn=f;

       do iside=1,numver
        if(.not.allocated(eside)) allocate(eside(1,
     &        size(localsurfconn(iside,:))))
        eside(1,:)=localsurfconn(iside,:)
        if(.not.allocated(cn)) allocate(cn(size(eside,2),3))
        cn=tet_coord(eside(1,:),:)
        enormal(1,:)=normals(iside,:)
        ax1=0.0d0
        call cross(cn(2,:)-cn(1,:),cn(3,:)-cn(1,:),crossans)
        area=0.5d0*norm2(crossans)
       
        do igp=1,size(q1,1)
          xi_pt(1,:)=q1(igp,:)
          N1d(1,1)=1.0d0-xi_pt(1,1)-xi_pt(1,2)
          N1d(2,1)=xi_pt(1,1)
          N1d(3,1)=xi_pt(1,2)
          xieta_gp=matmul(transpose(N1d),cn)
          wt=area/3.0d0

          sztet_coord=size(tet_coord,1);ismember=0
          ismember=minloc(eside(1,:),1,mask=(eside(1,:)).eq.sztet_coord)
           if((any(eside(1,:).eq.sztet_coord))) then
              call wachspress3d(ncoord,nn,un,g,xieta_gp,Nwach,abc)
            else
              plane(1,:)=surfplane(iface,:)
              CFnode(:,1)=cface_coord(:,plane(1,1))
              CFnode(:,2)=cface_coord(:,plane(1,2))
              xiet(:,1)=xieta_gp(:,plane(1,1))
              xiet(:,2)=xieta_gp(:,plane(1,2))
             call wachspress2d(CFnode,xiet,size(CFnode,1),phi)
	     Nwach(cface,1)=phi
            end if
          ax1(1,:)=ax1(1,:)+wt*Nwach(:,1)
         end do ! igp 
        if(allocated(cn)) deallocate(cn)
        if(allocated(eside)) deallocate(eside)
        fx(1,:)=fx(1,:) + enormal(1,1)*ax1(1,:)
        fy(1,:)=fy(1,:) + enormal(1,2)*ax1(1,:)
        fz(1,:)=fz(1,:) + enormal(1,3)*ax1(1,:)
       end do ! iside
         dNdx=fx/We;
         dNdy=fy/We;
         dNdz=fz/We;
       end subroutine getsmoothderivative

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
!     Following subroutine triangle information
       subroutine get_tri(nc,tri)
        implicit none

        integer*8 nc, tri(nc,3)
        select case(nc)
        case(3)
        tri=reshape((/3,2,1/),(/1,3/))

        case(4)
        tri=reshape((/5,5,5,5,2,3,4,1,1,2,3,4/),(/4,3/))

        case(5)
        tri=reshape((/6,6,6,6,6,2,3,4,5,1,1,2,3,4,5/),(/5,3/))

        case(6)
        tri=reshape((/7,7,7,7,7,7,2,3,4,5,6,1,1,2,3,4,5,6/),(/6,3/))

        case(7)
        tri=reshape((/8,8,8,8,8,8,8,2,3,4,5,6,7,1,1,2,3,4,5,6,7/),
     &         (/7,3/))

        case(8)
        tri=reshape((/9,9,9,9,9,9,9,9,2,3,4,5,6,7,8,1,1,2,3,4,5,6,7,8/),
     &       (/8,3/))

        case(9)
        tri=reshape((/10,10,10,10,10,10,10,10,10,2,3,4,5,6,7,8,9,1,1,2,
     &       3,4,5,6,7,8,9/),(/9,3/))
        end select
        end subroutine get_tri
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
