module isotropic_ray_module

  use precision, wp => dp

contains

  subroutine isotropic_ray(e,v,n_d,s,r,ftol,d_out,t_out)

    !====================================================================
    !
    ! Esta subrutina calcula el rayo en un medio isotropo, utilizando
    ! la técnica propuesta en
    !
    ! "A rapid and accurate two-point ray tracing method in
    ! horizontally layered velocity model" de Tian y Chen  !
    ! IN: e     --> espesores de las capas entre la fuente y el receptor.
    !     v     --> velocidad de la onda, puede ser s o p
    !     n_int --> cantida de disc. entre s y r
    !     s     --> coordenadas de la fuente
    !     r     --> coordenadas del receptor
    !
    ! OUT: d_out --> coordenadas x de las intersecciones entre
    !                el rayo y las disc. del modelo
    !
    !====================================================================


    ! implicit none
    
     integer, intent(in)::n_d
     real(kind=wp), dimension(:), intent(in)::e,v,s,r
     real(kind=wp), dimension(:), intent(inout)::d_out
         

     real(kind=wp)::v_max,e_max,q,delta,delta_f,a,b,ftol,t_out
     integer::i

    
    !averiguo cuantas y cuales son las capas con v=v_max, y calculo el
    !espesor e_max equivalente
    v_max=maxval(v) !velocidad máxima
    e_max=sum(e(pack([(i,i=1,size(v))],v==v_max))) !espesor de v_max equivalente

    !calculo la distancia horizontal entre la fuente y el receptor
    delta=sqrt((s(1)-r(1))**2+(s(2)-r(2))**2)

    !calculo el q inicial
    a=0.0_wp
    b=0.0_wp
    q=0.0_wp
    do i=1,n_d
       a=a+(v(i)/v_max)*e(i)
       if(v(i).ne.v_max) then
          b=b+((v(i)/v_max)*e(i))/sqrt(1-(v(i)/v_max)**2)
       end if
    end do
    a=a/e_max
    if(delta.lt.(a*b)/(a-1))then
       q=delta/a
    elseif (delta.gt.(a*b)/(a-1))then
       q=delta-b
    end if
    
    !hago newton-rhapson para averiguar el parámetro q
    delta_f=F(e,v,v_max,e_max,q) !delta_f inicial
    do while(abs(1.0_wp-delta_f/delta).gt.ftol)
       
       q=q-(F(e,v,v_max,e_max,q)-delta)/F_p(e,v,v_max,e_max,q)
       delta_f=F(e,v,v_max,e_max,q)

    end do
   
    
    !calculo el tiempo de viaje
    t_out=0.0_wp
    do i=1,n_d
       t_out=t_out+e(i)/(v(i)*sqrt(1-((q/(v_max*sqrt(e_max**2+q**2)))*v(i))**2))
    end do

   
    !calculo las distancias parciales hasta el receptor
    !que tambien sirven para calcular las coordenadas de la intersección
    !del rayo con las discontinuidades
    do i=1,n_d-1
       d_out(i)=F(e(1:i),v(1:i),v_max,e_max,q)
    end do
   

  end subroutine isotropic_ray

  !========================================================================
  ! Funciones
  !========================================================================



  real(kind=wp) function F(e,v,v_max,e_max,q)

    !=====================================================================
    ! Esta funcion calcula la distancia delta en función del modelo de
    ! velocidades y del parámetro q
    ! la función devuelve un arreglo con los valores de la distancia
    ! horizontal recorrida por cada segmento del rayo
    !=====================================================================

    implicit none

    real(kind=wp), dimension(:), intent(in)::e,v
    real(kind=wp), intent(in)::v_max,e_max,q

    integer::i

    
    F=0.0_wp
    do i=1,size(e)
       F=F+(v(i)/v_max)*e(i)*q/&
            sqrt(e_max**2+(1-(v(i)/v_max)**2)*q**2)
    end do

    return

  end function F

  ! !========================================================================

  real(kind=wp) function F_p(e,v,v_max,e_max,q)

    !======================================================================
    !Esta funcion es la derivada de la anterior, mas o menos...
    !======================================================================

    implicit none

    real(kind=wp), dimension(:), intent(in)::e,v
    real(kind=wp), intent(in)::v_max,e_max,q
    integer::i

    F_p=0.0_wp
    do i=1,size(e)
       F_p=F_p+(v(i)/v_max)*e(i)/&
            (e_max**2+(1-(v(i)/v_max)**2)*q**2)**(3.0_wp/2.0_wp)
    end do
    F_p=F_p*e_max**2
    return

  end function F_p


end module isotropic_ray_module
