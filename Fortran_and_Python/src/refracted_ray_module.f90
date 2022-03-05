module refracted_ray_module

  use precision, wp => dp

contains

  subroutine refracted_ray(e_in,&
       vp_in,&
       vs_in,&
       an_in,&
       s,&
       r,&
       dp_out,&
       dsv_out,&
       dsh_out,&
       tp,&
       tsv,&
       tsh)

    
    !====================================================================
    !
    ! Ray tracing para medio anisotropo, de placas planas y paralelas.
    ! El rayo anisotropo se consigue por medio de conjugate-gradient
    ! haciendo bending.  Como modelo inicial para el conjugate gradient
    ! se utiliza la solución para el rayo isotropo, obtenida segun:
    !
    ! "A rapid and accurate two-point ray tracing method in
    ! horizontally layered velocity model" de Tian y Chen
    !
    ! La idea es buscar el parámetro de rayo sísmico p tal que el rayo
    ! sea el de menor tiempo y una la fuente y el receptor fijos
    ! Realmente no se busca p, sino q que es algo parecido.
    !
    !
    ! IN: e_in    --> espesores de las capas entre la fuente y el receptor.
    !     vp_in   --> velocidad vertical p 
    !     vs_in   --> velocidad vertical s
    !     ani_in  --> anisotropia: epsilon delta, gamma
    !     s     --> coordenadas de la fuente
    !     r     --> coordenadas del receptor
    !
    ! OUT: dp_out --> intersecciones con las disc del rayo p
    !      dsh_out --> intersecciones con las disc del rayo sh
    !      dsv_out --> intersecciones con las disc del rayo sv
    !      tp  --> tiempo de viaje del rayo p
    !      tsh --> tiempo de viaje del rayo sh
    !      tsv --> tiempo de viaje del rayo sv
    !
    !====================================================================
  

    use isotropic_ray_module
    use cg_aux_mod, only:func,dfunc
    use conjgrad_mod, only:frprmn
    implicit none

    !in
    real(kind=wp), dimension(:), intent(in)::e_in,vp_in,vs_in
    real(kind=wp), dimension(:,:), intent(in)::an_in
    real(kind=wp), dimension(:), intent(in)::s,r

    !out
    real(kind=wp), dimension(size(e_in)-1), intent(inout)::dp_out,&
         dsv_out,dsh_out
    real(kind=wp), intent(out)::tp,tsv,tsh
    !out
    
    !local    
    real(kind=wp), dimension(size(e_in))::e,vp,vs
    real(kind=wp), dimension(size(an_in,1),size(an_in,2))::an
    real(kind=wp), dimension(size(e_in)-1)::ds_out
    real(kind=wp)::tp_out,ts_out,ang_rs,ang_tol
    real(kind=wp)::ftol,fret
    integer::iter,n_e


   
    n_e=size(e)

 
    
    !si la fuente está por debajo del receptor invierto el modelo
    !de espesores y velocidades para que tengan coherencia los indices
    if(s(3).gt.r(3))then
       e(1:n_e)=e_in(n_e:1:-1)
       vp(1:n_e)=vp_in(n_e:1:-1)
       vs(1:n_e)=vs_in(n_e:1:-1)
       an(1:n_e,:)=an_in(n_e:1:-1,:)
    else
       e=e_in
       vp=vp_in
       vs=vs_in
       an=an_in
    end if
 
    !llamo a las subrutinas que van a calcular el rayo para el medio
    !isotropo, para el rayo p y para el rayo s
    ftol=1e-4_wp
    call isotropic_ray(e,vp,n_e,s,r,ftol,dp_out,tp_out)
    call isotropic_ray(e,vs,n_e,s,r,ftol,ds_out,ts_out)

    dsv_out=ds_out
    dsh_out=ds_out
    
    tp=tp_out
    tsv=ts_out
    tsh=ts_out

    !angulo entre fuente y receptor
    ang_rs=&
         (45/atan(1.0))*atan(sqrt((r(1)-s(1))**2+(r(2)-s(2))**2)/abs(r(3)-s(3)))

    ang_tol=1.0 !angulo de tolerancia
    
    !si no es un rayo vertical calculo el anisotropo. Si es casi vertical
    !el anisotropo es igual al isotropo y me quedo con el ya calculado
    if(ang_rs.gt.ang_tol) then
       
       !A partir del rayo obtenido para medio isotropo, ajusto el rayo
       !para medio anisotropo usando bending, ajustando la posición
       !x de las intersecciones del rayo con las disc. (dp_out...)
       !usando conjugate gradient.
       call frprmn(dp_out,ftol,iter,fret,s,r,e,vp,vs,an,'vp')
       call frprmn(dsh_out,ftol,iter,fret,s,r,e,vp,vs,an,'sh')
       call frprmn(dsv_out,ftol,iter,fret,s,r,e,vp,vs,an,'sv')
       
       ! !tiempos de cada rayo en el medio anisotropo
       tp=func(dp_out,e,s,r,vp,vs,an,'vp')
       tsv=func(dsv_out,e,s,r,vp,vs,an,'sv')
       tsh=func(dsh_out,e,s,r,vp,vs,an,'sh')
  
    endif

       
  end subroutine refracted_ray

  
end module refracted_ray_module




