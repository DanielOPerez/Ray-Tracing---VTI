program ray_tracing_main

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
  use ray_tracing_mod, only:ray_tracing
  use count_file_lines_mod, only:count_file_lines

  implicit none
  
  real(kind=wp), allocatable::z(:),vp(:),vs(:),an(:,:),s(:,:),r(:,:),&
       tp(:,:),tsv(:,:),tsh(:,:),cp(:,:,:,:),csv(:,:,:,:),csh(:,:,:,:)
  integer, allocatable::n_int(:,:)
  integer::i,j,k,n_d,n_s,n_r,n_arg
  character(len=32), allocatable :: arg(:)
  
  !se leen los argumentos por linea de comando
  !modelo, fuente, receptores
  n_arg=command_argument_count()
  if(n_arg.lt.3)then
     write(*,*)'Se usa ./ray modelo fuente receptor'
     stop
  end if

  allocate(arg(n_arg))
  do i = 1, n_arg
     call get_command_argument(i, arg(i))
  end do

  !Se lee el modelo de capas y velocidades
  open(unit=100,file=arg(1),action='read')
  call count_file_lines(arg(1),n_d)
  allocate(z(n_d),&
       vp(n_d),&
       vs(n_d),&
       an(n_d,3))
  do i=1,n_d
     read(100,*)z(i),vp(i),vs(i),an(i,:)
  end do
  close(100)

  !se leen las coordenadas x,y,z de las fuentes
  open(unit=100,file=arg(2),action='read')
  call count_file_lines(arg(2),n_s)
  allocate(s(n_s,3))
  do i=1,n_s
     read(100,*)s(i,:)
  end do
  close(100)
  
  !se leen las coordenadas x,y,z de los receptores
  open(unit=100,file=arg(3),action='read')
  call count_file_lines(arg(3),n_r)
  allocate(r(n_r,3))
  do i=1,n_r
     read(100,*)r(i,:)
  end do
  close(100)


  !alojo variables
  allocate(tp(n_s,n_r),&    !tiempo de viaje de la onda p
       tsv(n_s,n_r),&       !tiempo de viaje de la onda sv
       tsh(n_s,n_r),&       !tiempo de viaje de la onda sh
       cp(n_s,n_r,n_d,3),&      !xyz interseccion de la onda p con las disc
       csh(n_s,n_r,n_d,3),&   !xyz interseccion de la onda sh con las disc
       csv(n_s,n_r,n_d,3),&   !xyz interseccion de la onda sv con las disc
       n_int(n_s,n_r))


  !Llamo al a subrutina que hace el ray tracing
  call ray_tracing(z,vp,vs,an,s,r,&
       tp,tsh,tsv)
  
 

  open(unit=100,file='time',action='write')
  open(unit=110,file='rayos',action='write')

  do i=1,n_s
     do j=1,n_r
        
        write(100,*)i,j,r(j,3),tp(i,j),tsh(i,j),tsv(i,j),&
             tsh(i,j)-tp(i,j),tsv(i,j)-tp(i,j),tsv(i,j)-tsh(i,j)
        
                
        write(110,*)s(i,1),s(i,2),s(i,1),s(i,2),s(i,1),s(i,2),s(i,3)
        do k=1,n_int(i,j)
           write(110,*)s(i,1)+cp(i,j,k,1),&
                s(i,2)+cp(i,j,k,2),&
                s(i,1)+csh(i,j,k,1),&
                s(i,2)+csh(i,j,k,2),&
                s(i,1)+csv(i,j,k,1),&
                s(i,2)+csv(i,j,k,2),&
                cp(i,j,k,3)
        end do
        write(110,*)r(j,1),r(j,2),r(j,1),r(j,2),r(j,1),r(j,2),r(j,3)
        write(110,*)''
        
     end do

      write(110,*)''
      write(110,*)''
      
  end do
  
  close(100)
  close(110)
   

  deallocate(tp,&
       tsh,&
       tsv,&
       cp,&
       csh,&
       csv,&
       n_int)
       

end program ray_tracing_main

!===========================================================

