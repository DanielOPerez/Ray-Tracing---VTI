module global_variables_mod

  use precision, wp => dp

  implicit none

  real(kind=wp), allocatable::vp(:),vs(:),e(:),an(:)
  integer, allocatable::z_ind(:)
  integer::n_int
  
end module global_variables_mod
