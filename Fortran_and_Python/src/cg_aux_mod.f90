MODULE cg_aux_mod

  use precision, wp => dp

CONTAINS

  !=======================================================================================

  function func(x_in,e,s,r,vp,vs,an,ray) ! Tiempo de viaje
 
    implicit none
    
    real(kind=wp), dimension(:), intent(IN) :: x_in,e,vp,vs,s,r
    real(kind=wp), dimension(:,:), intent(in)::an
    real(kind=wp), dimension(size(x_in)+2)::x
    real(kind=wp) :: func,phi,dx
    character(len=2), intent(in)::ray
    integer::i

    
    x(1)=0.0_wp
    x(2:size(x_in)+1)=x_in
    x(size(x_in)+2)=sqrt((r(1)-s(1))**2+(r(2)-s(2))**2)
  
    
    func=0.0_wp

    do i=2,size(e)+1
       
       dx=abs(x(i-1)-x(i))
       phi=atan(dx/e(i-1))
       
       func=func+dx/(an_vel(phi,vp(i-1),vs(i-1),an(i-1,:),ray)*sin(phi))
   
    end do
  
    
   return
    
  end function func

  !=======================================================================================

  function dfunc(x_in,e,s,r,vp,vs,an,ray)

    !Para calcular la derivada de un punto x(i) se necesita también el punto x(i+1) y
    !el punto x(i-1), por lo tanto se van a calcular las derivadas desde i=2 a i=N-1,
    !siendo x(N) el receptor; x(1) es la fuente. 
    
    implicit none
    
    real(kind=wp), dimension(:), intent(IN) :: x_in,e,vp,vs,s,r
    real(kind=wp), dimension(:,:), intent(in)::an
    real(kind=wp), dimension(size(x_in)+2)::x
    real(kind=wp), dimension(size(x_in))::dfunc
    real(kind=wp)::phi1,phi2,dphi1,dphi2,dx1,dx2,v1,v2,sp1,sp2
    character(len=2), intent(in)::ray
    integer::i

  
    x(1)=0.0_wp
    x(2:size(x_in)+1)=x_in
    x(size(x_in)+2)=sqrt((r(1)-s(1))**2+(r(2)-s(2))**2)
          
    do i=2,size(x)-1
       
       dx1=x(i)-x(i-1)
       dx2=x(i+1)-x(i)
       phi1=atan(dx1/e(i-1))
       phi2=atan(dx2/e(i))
       dphi1=e(i-1)/(e(i-1)**2+dx1**2)
       dphi2=-e(i)/(e(i)**2+dx2**2)
       v1=an_vel(phi1,vp(i-1),vs(i-1),an(i-1,:),ray)
       v2=an_vel(phi2,vp(i),vs(i),an(i,:),ray)
       sp1=sin(phi1)
       sp2=sin(phi2)
       
       dfunc(i-1)=&
            ((sp1*v1-dx1*dphi1*(cos(phi1)*v1+&
            sp1*an_dvel(phi1,vp(i-1),vs(i-1),an(i-1,:),ray)))/&
            (sp1*v1)**2)-& !!!=================================
            ((sp2*v2+dx2*dphi2*(cos(phi2)*v2+&
            sp2*an_dvel(phi2,vp(i),vs(i),an(i,:),ray)))/&
            (sp2*v2)**2)

       
    end do
    
  end function dfunc

  !=======================================================================================
  
  
  function an_vel(phi,vp,vs,an,ray)
    implicit none
    real(kind=wp), intent(in)::phi,vp,vs
    real(kind=wp), dimension(:), intent(in)::an
    real(kind=wp)::an_vel
    character(len=2), intent(in)::ray

    an_vel=0.0_wp
    
    if (ray=='vp' ) then
       an_vel=vp*(1.0_wp+an(2)*(sin(phi)*cos(phi))**2+an(1)*sin(phi)**4)
    elseif(ray=='sv') then
       an_vel=vs*(1.0_wp+(vp/Vs)**2*(an(1)-an(2))*(sin(phi)*cos(phi))**2)
    elseif(ray=='sh') then
       an_vel=vs*(1.0_wp+an(3)*sin(phi)**2)
    end if
    
    return
    
  end function an_vel

   !=======================================================================================

  function an_dvel(phi,vp,vs,an,ray)
    
    implicit none
    real(kind=wp), intent(in)::phi,vp,vs
    real(kind=wp), dimension(:), intent(in)::an
    real(kind=wp)::an_dvel
    character(len=2), intent(in)::ray

    an_dvel=0.0_wp
    
    if (ray=='vp' ) then
       an_dvel=2.0_wp*vp*sin(phi)*cos(phi)*(an(2)*cos(phi)**2+(2.0_wp*an(1)-an(2))*sin(phi)**2)
    elseif(ray=='sv') then
       an_dvel=2.0_wp*vs*(an(1)-an(2))*((vp/vs)**2)*sin(phi)*cos(phi)*(cos(phi)**2-sin(phi)**2)
    elseif(ray=='sh') then
       an_dvel=2.0_wp*vs*an(3)*sin(phi)*cos(phi)
    end if
    
    return
    
  end function an_dvel

 
  !=======================================================================================

END MODULE cg_aux_mod
