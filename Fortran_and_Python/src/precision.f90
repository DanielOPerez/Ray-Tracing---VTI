  module precision

    ! --------------------------------------------------------------------------
    !  sp: simple precision de la norma IEEE 754
    !  dp: doble precision de la norma IEEE 754
    !
    !  Uso: use precision, wp => sp o use precision, wp => dp
    ! --------------------------------------------------------------------------

    integer,parameter::sp=selected_real_kind(6,37)
    integer,parameter::dp=selected_real_kind(15,307)
    
  end module precision 
    
