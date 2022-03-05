module straight_ray_module

  use precision, wp => dp

contains

  subroutine straight_ray(vp,vs,an,s,r,tp,tsv,tsh)

    !=======================================================================
    !
    ! Esta rutina calcula los tiempos de arrivo de las ondas p,ssv y ssh
    ! para el modelo de weak anisotropy de thompsen.
    !
    !======================================================================
   
    
    implicit none
    real(kind=wp), intent(in)::vp,vs
    real(kind=wp), dimension(:), intent(in)::s,r,an
    
    real(kind=wp)::theta,dist,tp,tsv,tsh

     
    dist=sqrt((r(1)-s(1))**2+(r(2)-s(2))**2+(r(3)-s(3))**2)
    
    theta=atan(sqrt((r(1)-s(1))**2+(r(2)-s(2))**2)/abs(r(3)-s(3)))

    !tiempo de viaje de los rayos
    tp=dist/(vp*(1.0_wp+an(2)*(sin(theta)*cos(theta))**2+an(1)*sin(theta)**4))
    tsv=dist/(vs*(1.0_wp+(vp/vs)**2*(an(1)-an(2))*(sin(theta)*cos(theta))**2))
    tsh=dist/(vs*(1.0_wp+an(3)*sin(theta)**2))



   
    
  end subroutine straight_ray

    
  
end module straight_ray_module
