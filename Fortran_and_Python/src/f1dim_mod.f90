MODULE f1dim_mod
  !Used for communication from linmin to f1dim.
  USE nrtype
  USE cg_aux_mod, only: func
  INTEGER(I4B) :: ncom
  REAL(DP), DIMENSION(:), POINTER :: pcom,xicom
CONTAINS
  FUNCTION f1dim(x,e,s,r,vp,vs,an,ray)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: x
    REAL(DP), DIMENSION(:), INTENT(IN) ::e,vp,vs,s,r
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: an
    CHARACTER(LEN=2), INTENT(IN)::ray
    REAL(DP) :: f1dim
    !Used by linmin as the one-dimensional function passed to mnbrak and brent.
    !INTERFACE
    !   FUNCTION func(x)
    !     USE nrtype
    !     REAL(SP), DIMENSION(:), INTENT(IN) :: x
    !     REAL(SP) :: func
    !   END FUNCTION func
    !END INTERFACE
    REAL(DP), DIMENSION(:), ALLOCATABLE :: xt
    allocate(xt(ncom))
    xt(:)=pcom(:)+x*xicom(:)
    f1dim=func(xt,e,s,r,vp,vs,an,ray)
    deallocate(xt)
  END FUNCTION f1dim
END MODULE f1dim_mod
