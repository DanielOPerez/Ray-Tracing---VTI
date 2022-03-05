module count_file_lines_mod

  use precision, wp=>dp

contains
  
  subroutine count_file_lines(file_in,lines_out)

    !esta subrutina cuenta las lineas de un archivo
    !si el archivo no fue abierto, se abre y se cierra, sino se deja abierto
    
    implicit none
    character(len=*), intent(in)::file_in
    integer, intent(out)::lines_out
    integer::myunit,io
    real(kind=wp)aux

    !veo si el archivo está abierto
    inquire(file=file_in, number=myunit)

    if(myunit.le.0)then !si no está abierto
       lines_out=0
       open(newunit=myunit,file=file_in,action='read')
       do
          read(myunit,*,iostat=io)aux
          if(io/=0)exit
          lines_out=lines_out+1
       end do
       close(myunit)
    else !si está abierto
       lines_out=0
       do
          read(myunit,*,iostat=io)aux
          if(io/=0)exit
          lines_out=lines_out+1
       end do
       rewind(myunit)
    end if
       
    
  end subroutine count_file_lines

  
end module count_file_lines_mod
