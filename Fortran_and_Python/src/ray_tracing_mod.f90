module ray_tracing_mod

  !==================================================================
  !
  ! Ray tracing para medio anisotropo, de placas planas y paralelas.
  !
  ! 
  ! Ver "A rapid and accurate two-point ray tracing method in
  ! horizontally layered velocity model" de Tian y Chen
  !
  !==================================================================

  use precision, wp => dp
  use refracted_ray_module
  use straight_ray_module

  implicit none
  
contains

  subroutine ray_tracing(z,vp,vs,an,s,r,&
       tp_out,tsh_out,tsv_out)

    !in
    real(kind=wp), intent(in), dimension(:)::z,vp,vs
    real(kind=wp), intent(in), dimension(:,:)::an,s,r

    !out
    real(kind=wp), intent(out), dimension(:,:)::tp_out,tsh_out,tsv_out
    !real(kind=wp), intent(out), dimension(:,:,:,:)::cp_out,csh_out,csv_out
    !integer, intent(out),dimension(:,:)::n_int
  
   
    !local
    real(kind=wp), dimension(size(z)+1)::e !modelos de espesores
    real(kind=wp), dimension(size(s,1),size(r,1),size(z))::dp_out,dsh_out,dsv_out
    !real(kind=wp)::cos_a,sin_a
    integer, dimension(size(s,1),size(r,1))::n_int
    integer, dimension(size(z))::z_ind !indice de las discontinuidades
    integer::i,j,k

    
   
    z_ind=0
    n_int=0
    e=0.0_wp
    tp_out=0.0_wp
    tsh_out=0.0_wp
    tsv_out=0.0_wp
    dp_out=0.0_wp
    dsh_out=0.0_wp
    dsv_out=0.0_wp
    ! cp_out=0.0_wp
    ! csh_out=0.0_wp
    ! csv_out=0.0_wp

    !para cada par fuente-receptor busco los tiempos de viaje
 
    do j=1,size(r,1)
       do i=1,size(s,1)


          !averiguo n_int(i,j): nro de discontinuidades entre cada fuente y cada receptor
          n_int(i,j)=count(z.lt.max(s(i,3),r(j,3)).and.z.gt.min(s(i,3),r(j,3)))

          !si n_int(i,j)=0 la fuente y el receptor están en la misma capa
          !y saco el tiempo de viaje directamente..
          !si n_int(i,j) no es cero hay que estimar los rayos iterativamente
          if(n_int(i,j).eq.0)then


             !me fijo cual es la capa donde están la fuente y el receptor
             z_ind(1)=minval(pack([(k,k=1,size(z))],z.gt.max(s(i,3),r(j,3))))
             call straight_ray(vp(z_ind(1)),&
                  vs(z_ind(1)),&
                  an(z_ind(1),:),&
                  s(i,:),&
                  r(j,:),&
                  tp_out(i,j),&
                  tsv_out(i,j),&
                  tsh_out(i,j))

             
          else

             !==============================================================
             ! el rayo refractado anisotropo se calcula en dos partes:
             !
             ! 1) se calcula el caso de modelo isotropo
             ! 2) usando el resultado de 1) como modelo incial se calcula
             !    el rayo anisotropo
             !
             !==============================================================


             !indexo las discontinuidades entre la fuente y el receptor.
             z_ind(1:n_int(i,j))=pack([(k,k=1,size(z))],&
                  (z.lt.max(s(i,3),r(j,3)).and.z.gt.min(s(i,3),r(j,3))))

             
             !defino el modelo de  espesores 
             e(1)=z(z_ind(1))-min(s(i,3),r(j,3))
             e(2:n_int(i,j))=z(z_ind(2:n_int(i,j)))-z(z_ind(1:n_int(i,j)-1))
             e(n_int(i,j)+1)=max(s(i,3),r(j,3))-z(z_ind(n_int(i,j)))

             !calculo los tiempos de arribo para el rayo refrectado,
             !tambien se sacan las distancias desde la fuente hacia
             !cada una de las intersecciones, en la linea que une
             !la fuente y el receptor
             call refracted_ray(e(1:n_int(i,j)+1),&
                  vp(z_ind(1):z_ind(n_int(i,j))+1),&
                  vs(z_ind(1):z_ind(n_int(i,j))+1),&
                  an(z_ind(1):z_ind(n_int(i,j))+1,:),&
                  s(i,:),&
                  r(j,:),&
                  dp_out(i,j,1:n_int(i,j)),&
                  dsv_out(i,j,1:n_int(i,j)),&
                  dsh_out(i,j,1:n_int(i,j)),&
                  tp_out(i,j),&
                  tsv_out(i,j),&
                  tsh_out(i,j))

       
             ! !guardo la coord xyz de cada segmento del rayo
             ! !-----------------------------------------------------------------
             ! !si la fuente está por debajo de los receptores
             ! !invierto el indice de z_ind para tener la coordenada z
             ! !en el orden correcto
             ! if(s(i,3).gt.r(j,3))then
             !   z_ind(1:n_int(i,j))=z_ind(n_int(i,j):1:-1)
             ! end if
             
             ! !guardo la coordenada xyz de cada interseccion de cada rayo
             ! !con cada discontinuidad
             ! cos_a=(r(j,1)-s(i,1))/sqrt((s(i,1)-r(i,1))**2+(s(i,2)-r(i,2))**2)
             ! sin_a=(r(j,2)-s(i,2))/sqrt((s(i,1)-r(i,1))**2+(s(i,2)-r(i,2))**2)
             ! !rayo p
             ! cp_out(i,j,1:n_int(i,j),1)=dp_out(i,j,1:n_int(i,j))*cos_a   !coord x
             ! cp_out(i,j,1:n_int(i,j),2)=dp_out(i,j,1:n_int(i,j))*sin_a   !coord y
             ! cp_out(i,j,1:n_int(i,j),3)=z(z_ind(1:n_int(i,j)))           !coord z
             ! !rayo sh
             ! csh_out(i,j,1:n_int(i,j),1)=dsh_out(i,j,1:n_int(i,j))*cos_a !coord x
             ! csh_out(i,j,1:n_int(i,j),2)=dsh_out(i,j,1:n_int(i,j))*sin_a !coord y
             ! csh_out(i,j,1:n_int(i,j),3)=z(z_ind(1:n_int(i,j)))          !coord z
             ! !rayo sv
             ! csv_out(i,j,1:n_int(i,j),1)=dsv_out(i,j,1:n_int(i,j))*cos_a !coord x
             ! csv_out(i,j,1:n_int(i,j),2)=dsv_out(i,j,1:n_int(i,j))*sin_a !coord y
             ! csv_out(i,j,1:n_int(i,j),3)=z(z_ind(1:n_int(i,j)))          !coord z

         
             
          end if

       end do  ! i=1,n_s
    end do  ! j=1,n_r

    
    
  end subroutine ray_tracing

end module ray_tracing_mod


!===========================================================

