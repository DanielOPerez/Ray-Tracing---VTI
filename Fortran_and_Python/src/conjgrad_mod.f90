MODULE conjgrad_mod

  use precision, wp => dp
  use nrtype
 
  
CONTAINS

  SUBROUTINE frprmn(p,ftol,iter,fret,ss,rr,ee,vp,vs,an,ray)
    USE nrutil, ONLY : nrerror
    USE cg_aux_mod, only: func, dfunc

    IMPLICIT NONE
    INTEGER(I4B), INTENT(OUT) :: iter
    REAL(DP), INTENT(IN) :: ftol
    REAL(DP), INTENT(OUT) :: fret
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: p
    real(kind=wp), dimension(:), intent(IN) :: ee,vp,vs,ss,rr
    real(kind=wp), dimension(:,:), intent(in)::an
    character(len=2), intent(in)::ray

    INTEGER(I4B), PARAMETER :: ITMAX=200
    REAL(DP), PARAMETER :: EPS=1.0e-10_DP
    ! Given  a  starting  point  p  that is  a  vector  of  length  N,
    ! Fletcher-Reeves-Polak-Ribiere  minimization  is performed  on  a
    ! function func,  using its  gradient as  calculated by  a routine
    ! dfunc. The convergence tolerance on  the function value is input
    ! as  ftol.  Returned  quantities  are  p  (the  location  of  the
    ! minimum), iter  (the number of iterations  that were performed),
    ! and fret (the minimum value of the function). The routine linmin
    ! is called  to perform line minimizations.   Parameters: ITMAX is
    ! the maximum allowed number of  iterations; EPS is a small number
    ! to  rectify  the special  case  of  converging to  exactly  zero
    ! function value.
    INTEGER(I4B) :: its
    REAL(DP) :: dgg,fp,gam,gg
    REAL(DP), DIMENSION(size(p)) :: g,h,xi

    
    fp=func(p,ee,ss,rr,vp,vs,an,ray)            ! Initializations.
    xi=dfunc(p,ee,ss,rr,vp,vs,an,ray)
    g=-xi
  
    gg=dot_product(g,g)
    if (gg == 0.0) RETURN    ! Unlikely. If gradient is exactly zero
   
    h=g
    xi=h
    
    do its=1,ITMAX        ! Loop over iterations.
       iter=its
     
       call linmin(p,xi,fret,ee,ss,rr,vp,vs,an,ray)   !Next statement is the normal return:
       if (2.0_DP*abs(fret-fp) <= ftol*(abs(fret)+abs(fp)+EPS)) RETURN
       fp=fret
       xi=dfunc(p,ee,ss,rr,vp,vs,an,ray)
       gg=dot_product(g,g)
       dgg=dot_product(xi,xi)   ! This statement for Fletcher-Reeves.
       dgg=dot_product(xi+g,xi) ! This statement for Polak-Ribiere.
       if (gg == 0.0) RETURN    ! Unlikely. If gradient is exactly zero 
       gam=dgg/gg               ! then we are already done.
       g=-xi
       h=g+gam*h
       xi=h
    end do
    call nrerror('frprmn: maximum iterations exceeded')

  END SUBROUTINE frprmn

  SUBROUTINE linmin(p,xi,fret,ee,ss,rr,vp,vs,an,ray)
    USE nrutil, ONLY : assert_eq

    IMPLICIT NONE
    REAL(DP), INTENT(OUT) :: fret
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: p,xi
    REAL(DP), PARAMETER :: TOL=1.0e-4_DP
    real(kind=wp), dimension(:), intent(IN) :: ee,vp,vs,ss,rr
    real(kind=wp), dimension(:,:), intent(in)::an
    character(len=2), intent(in)::ray

    ! Given an N-dimensional  point p and an  N -dimensional direction
    ! xi, both  vectors of length N,  moves and resets p  to where the
    ! fixed-name function func takes on  a minimum along the direction
    ! xi from  p, and  replaces xi by  the actual  vector displacement
    ! that p was moved. Also returns as  fret the value of func at the
    ! returned  location  p.  This  is actually  all  accomplished  by
    ! calling  the routines  mnbrak and  brent.  Parameter:  Tolerance
    ! passed to brent.

    REAL(DP) :: ax,bx,fa,fb,fx,xmin,xx

    ax=0.0                   ! Initial guess for brackets.
    xx=1.0
    call mnbrak(ax,xx,bx,fa,fx,fb,p,xi,ee,ss,rr,vp,vs,an,ray)
    fret=brent(ax,xx,bx,p,xi,TOL,xmin,ee,ss,rr,vp,vs,an,ray)
    
    xi=xmin*xi               ! Construct the vector results to return.
    p=p+xi

  END SUBROUTINE linmin

  !-------------------------------------------------------------------

  FUNCTION brent(ax,bx,cx,pp,xxi,tol,xmin,ee,ss,rr,vp,vs,an,ray)
    USE nrutil, ONLY : nrerror
    USE cg_aux_mod, only: func
    
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: pp,xxi
    REAL(DP), INTENT(IN) :: ax,bx,cx,tol
    REAL(DP), INTENT(OUT) :: xmin
    REAL(DP) :: brent
    real(kind=wp), dimension(:), intent(IN) :: ee,vp,vs,ss,rr
    real(kind=wp), dimension(:,:), intent(in)::an
    character(len=2), intent(in)::ray
    INTEGER(I4B), PARAMETER :: ITMAX=100
    REAL(DP), PARAMETER :: CGOLD=0.3819660_DP,ZEPS=1.0e-3_DP*epsilon(ax)
    ! Given  a  function  func,  and given  a  bracketing  triplet  of
    ! abscissas ax,  bx, cx (such  that bx is  between ax and  cx, and
    ! func(bx) is less than both  func(ax) and func(cx)), this routine
    ! isolates  the minimum  to a  fractional precision  of about  tol
    ! using Brent's method.  !The abscissa of the  minimum is returned
    ! as xmin,  and the minimum  function value is returned  as brent,
    ! the returned function value.  Parameters: Maximum allowed number
    ! of iterations;  golden ratio; and  a small number  that protects
    ! against trying to achieve fractional accuracy for a minimum that
    ! happens to be exactly zero.
    INTEGER(I4B) :: iter
    REAL(DP) :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm

    brent=0.0_dp


    a=min(ax,cx) ! a and b must be in ascending order, though
    b=max(ax,cx) ! the input abscissas need not be.
    v=bx         !  Initializations...
    w=v
    x=v
    e=0.0        ! This will be the distance moved on the step
    fx=func(pp+x*xxi,ee,ss,rr,vp,vs,an,ray)  !before last.
    fv=fx
    fw=fx
    do iter=1,ITMAX ! Main program loop.
       xm=0.5_DP*(a+b)
       tol1=tol*abs(x)+ZEPS
       tol2=2.0_DP*tol1
       if (abs(x-xm) <= (tol2-0.5_DP*(b-a))) then !Test for done here.
          xmin=x ! Arrive here ready to exit with best values.
          brent=fx
          RETURN
       end if
       if (abs(e) > tol1) then ! Construct a trial parabolic fit.
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*r
          q=2.0_DP*(q-r)
          if (q > 0.0) p=-p
          q=abs(q)
          etemp=e
          if (abs(p) >= abs(0.5_DP*q*etemp) .or. &
               p <= q*(a-x) .or. p >= q*(b-x)) then
             ! The above conditions determine the acceptability of the 
             ! parabolic fit. Here it is not o.k., so we take the golden 
             ! section step into the larger of the two segments.
             e=merge(a-x,b-x, x >= xm )
             d=CGOLD*e
          else     ! Take the parabolic step.
             d=p/q
             u=x+d
             if (u-a < tol2 .or. b-u < tol2) d=sign(tol1,xm-x)
          end if
       else        ! Take the golden section step into the larger
          e=merge(a-x,b-x, x >= xm ) ! of the two segments.
          d=CGOLD*e
       end if
       u=merge(x+d,x+sign(tol1,d), abs(d) >= tol1 )
       ! Arrive here with d computed either from parabolic fit, 
       ! or else from golden section.
       fu=func(pp+u*xxi,ee,ss,rr,vp,vs,an,ray)
       ! This is the one function evaluation per iteration.
       if (fu <= fx) then      ! Now we have to decide what to do with our
          if (u >= x) then     ! function evaluation. Housekeeping follows:
             a=x
          else
             b=x
          end if
          call shft(v,w,x,u)
          call shft(fv,fw,fx,fu)
       else
          if (u < x) then
             a=u
          else
             b=u
          end if
          if (fu <= fw .or. w == x) then
             v=w
             fv=fw
             w=u
             fw=fu
          else if (fu <= fv .or. v == x .or. v == w) then
             v=u
             fv=fu
          end if
       end if
    end do !Done with housekeeping. Back for another iteration.
    call nrerror('brent: exceed maximum iterations')

  CONTAINS
    SUBROUTINE shft(a,b,c,d)
      REAL(DP), INTENT(OUT) :: a
      REAL(DP), INTENT(INOUT) :: b,c
      REAL(DP), INTENT(IN) :: d
      a=b
      b=c
      c=d
    END SUBROUTINE shft
  END FUNCTION brent



  SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,pp,xxi,ee,ss,rr,vp,vs,an,ray)

    USE cg_aux_mod, only: func
    
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: pp,xxi
    REAL(DP), INTENT(INOUT) :: ax,bx
    REAL(DP), INTENT(OUT) :: cx,fa,fb,fc
    real(kind=wp), dimension(:), intent(IN) :: ee,vp,vs,ss,rr
    real(kind=wp), dimension(:,:), intent(in)::an
    character(len=2), intent(in)::ray
    REAL(DP), PARAMETER :: GOLD=1.618034_DP,GLIMIT=100.0_DP,TINY=1.0e-20_DP
    ! Given a function func, and  given distinct initial points ax and
    ! bx, this routine searches in  the downhill direction (defined by
    ! the function as evaluated at the initial points) and returns new
    ! points ax, bx,  cx that bracket a minimum of  the function. Also
    ! returned are  the function values  at the three points,  fa, fb,
    ! and  fc.   Parameters:  GOLD  is  the  default  ratio  by  which
    ! successive  intervals  are  magnified;  GLIMIT  is  the  maximum
    ! magnification allowed for a parabolic-fit step.
    REAL(DP) :: fu,q,r,u,ulim
    INTEGER :: k

   
    fa=func(pp+ax*xxi,ee,ss,rr,vp,vs,an,ray)
    fb=func(pp+bx*xxi,ee,ss,rr,vp,vs,an,ray)

    
    if (fb > fa) then            ! Switch roles of a and b so that we
       call swap(ax,bx)          ! can go downhill in the direction
       call swap(fa,fb)          ! from a to b.
    end if

    cx=bx+GOLD*(bx-ax)           ! First guess for c.
    fc=func(pp+cx*xxi,ee,ss,rr,vp,vs,an,ray)
    k=0
   
    do                           ! Do-while-loop: Keep returning here    
       k=k+1     
      
       if (fb < fc) RETURN       ! until we bracket.
       ! Compute u by parabolic extrapolation from a, b, c. 
       ! TINY is used to prevent any possible division by zero.
       r=(bx-ax)*(fb-fc)
       q=(bx-cx)*(fb-fa)
       u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0_DP*sign(max(abs(q-r),TINY),q-r))
       ulim=bx+GLIMIT*(cx-bx) 
       ! We won't go farther than this. Test various possibilities:

       if ((bx-u)*(u-cx) > 0.0) then   ! Parabolic u is between b and c: try it
          fu=func(pp+u*xxi,ee,ss,rr,vp,vs,an,ray)
       
          if (fu < fc) then            ! Got a minimum between b and c.
             ax=bx
             fa=fb
             bx=u
             fb=fu
           
             RETURN
          else if (fu > fb) then       ! Got a minimum between a and u.
             cx=u
             fc=fu
       
             RETURN
          end if
          u=cx+GOLD*(cx-bx)            ! Parabolic fit was no use. Use default
          fu=func(pp+u*xxi,ee,ss,rr,vp,vs,an,ray)                   ! magnification.
          
       else if ((cx-u)*(u-ulim) > 0.0_DP) then ! Parabolic fit is between c and its al-
          fu=func(pp+u*xxi,ee,ss,rr,vp,vs,an,ray)
          ! lowed limit.
      
          if (fu < fc) then
             bx=cx
             cx=u
             u=cx+GOLD*(cx-bx)
        
             call shft(fb,fc,fu,func(pp+u*xxi,ee,ss,rr,vp,vs,an,ray))
            
          end if
       else if ((u-ulim)*(ulim-cx) >= 0.0_DP) then ! Limit parabolic u to maximum allowed value.
         
          u=ulim
          fu=func(pp+u*xxi,ee,ss,rr,vp,vs,an,ray)
       else
         
          u=cx+GOLD*(cx-bx)                     ! magnification.
          fu=func(pp+u*xxi,ee,ss,rr,vp,vs,an,ray)
       end if
       call shft(ax,bx,cx,u)
       call shft(fa,fb,fc,fu)                   ! Eliminate oldest point and continue.
    end do
  

  CONTAINS
    SUBROUTINE shft(a,b,c,d)
      REAL(DP), INTENT(OUT) :: a
      REAL(DP), INTENT(INOUT) :: b,c
      REAL(DP), INTENT(IN) :: d
      a=b
      b=c
      c=d
    END SUBROUTINE shft
    SUBROUTINE swap(a,b)
      REAL(DP) :: a,b,dum
      dum=a
      a=b
      b=dum
    END SUBROUTINE swap


  END SUBROUTINE mnbrak



END MODULE conjgrad_mod
