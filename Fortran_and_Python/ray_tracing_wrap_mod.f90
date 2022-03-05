module ray_tracing_wrap_mod
  
  use iso_c_binding, only: c_double, c_int
  use ray_tracing_mod, only: ray_tracing

  implicit none

contains

  subroutine ray_tracing_wrap(z_in,&
       vp_in,&
       vs_in,&
       an_in,&
       s_in,&
       r_in,&
       n_d,&
       n_s,&
       n_r,&
       tp_out,&
       tsh_out,&
       tsv_out) bind(c)

    !in
    integer(c_int), intent(in), value::n_d,n_s,n_r
    real(c_double), intent(in)::z_in(n_d),vp_in(n_d),vs_in(n_d),&
         an_in(3,n_d),s_in(3,n_s),r_in(3,n_r)

    !out 
    real(c_double), intent(inout), dimension(n_s,n_r)::tp_out,tsh_out,tsv_out
    integer::i,j

    
    call ray_tracing(&
         z_in,&
         vp_in,&
         vs_in,&
         transpose(an_in),&
         transpose(s_in),&
         transpose(r_in),&
         tp_out,&
         tsh_out,&
         tsv_out)
    



    
    
  end subroutine ray_tracing_wrap
  
end module ray_tracing_wrap_mod
