program gen_rec

  implicit none

  integer i,j,prof

  prof=2525
  
  open(unit=100,file='receptores',action='write')
  write(100,*)50
  do i=50,1000,100
     do j=100,1000,100
        write(100,*)i,j,prof
     end do
  end do
  close(100)

end program gen_rec
