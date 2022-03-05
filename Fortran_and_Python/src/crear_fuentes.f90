program crear_fuentes

  implicit none

  integer::i,j

  open(unit=100,file='fuentes',action='write')
  
  do i=200,1000,200
     do j=2375,2550,50
        write(100,*)i,0,j
     end do
  end do
  close(100)
  
end program crear_fuentes
