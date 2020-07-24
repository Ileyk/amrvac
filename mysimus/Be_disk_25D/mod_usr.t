module mod_usr
  use mod_hd
  implicit none

contains

  subroutine usr_init()
    call set_coordinate_system('spherical_2.5D')

    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_source        => pt_grav_source
    usr_get_dt        => get_dt_pt_grav
    usr_special_bc    => specialbound_usr

    call hd_activate()

  end subroutine usr_init

  subroutine initglobaldata_usr

    hd_adiab=1.d-2

  end subroutine initglobaldata_usr

  subroutine initonegrid_usr(ixG^L,ix^L,w,x)

    ! initialize one grid

    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    w(ix^S,rho_)   =  (x(ix^S,1)*dsin(x(ix^S,2)))**(-3.5d0) * dexp(-(1.d0/hd_adiab)*(1.d0/dsin(x(ix^S,2))-1.d0)/x(ix^S,1))
    w(ix^S,mom(1)) =  0.d0
    w(ix^S,mom(2)) =  0.d0
    where (w(ix^S,rho_)>small_density)
      w(ix^S,mom(3)) =  w(ix^S,rho_)*dsqrt(1.d0/(x(ix^S,1)*dsin(x(ix^S,2)))-3.5d0*hd_adiab)
    elsewhere
      w(ix^S,rho_)   =  small_density
      w(ix^S,mom(3)) =  0.d0
    end where

  end subroutine initonegrid_usr

  subroutine pt_grav_source(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)

    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    where (wCT(ixO^S,rho_)>10.d0*small_density)
      w(ixO^S,mom(1))=w(ixO^S,mom(1)) - &
           qdt * wCT(ixO^S,rho_)    / x(ixO^S,1)**2.d0
    endwhere
    ! w(ixO^S,e_ )=w(ixO^S,e_ ) - &
    !      qdt * wCT(ixO^S,mom(1))  / x(ixO^S,1)**2.d0

  end subroutine pt_grav_source

  subroutine get_dt_pt_grav(w,ixG^L,ix^L,dtnew,dx^D,x)

    integer, intent(in) :: ixG^L, ix^L
    double precision, intent(in) :: dx^D, x(ixG^S,1:ndim)
    double precision, intent(in) :: w(ixG^S,1:nw)
    double precision, intent(inout) :: dtnew

    dtnew=bigdouble
    dtnew=min(dtnew,minval(dsqrt(block%dx(ix^S,1)/(1.d0/x(ix^S,1)**2.d0))))

  end subroutine get_dt_pt_grav

  subroutine specialbound_usr(qt,ixG^L,ixB^L,iB,w,x)

    ! special boundary types, user defined

    integer, intent(in) :: ixG^L, ixB^L, iB
    double precision, intent(in) :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)
    double precision :: rho(ixG^S)

    select case(iB)
    case(1)

    w(ixB^S,rho_)   =  (x(ixB^S,1)*dsin(x(ixB^S,2)))**(-3.5d0) * dexp(-(1.d0/hd_adiab)*(1.d0/dsin(x(ixB^S,2))-1.d0)/x(ixB^S,1))
    w(ixB^S,mom(1)) =  0.d0
    w(ixB^S,mom(2)) =  0.d0
    where (w(ixB^S,rho_)>small_density)
      w(ixB^S,mom(3)) =  w(ixB^S,rho_)*dsqrt(1.d0/(x(ixB^S,1)*dsin(x(ixB^S,2)))-3.5d0*hd_adiab)
    elsewhere
      w(ixB^S,rho_)   = small_density
      w(ixB^S,mom(3)) = 0.d0
    end where

    case(2)

    w(ixB^S,rho_)   =  (x(ixB^S,1)*dsin(x(ixB^S,2)))**(-3.5d0) * dexp(-(1.d0/hd_adiab)*(1.d0/dsin(x(ixB^S,2))-1.d0)/x(ixB^S,1))
    w(ixB^S,mom(1)) =  0.d0
    w(ixB^S,mom(2)) =  0.d0
    where (w(ixB^S,rho_)>small_density)
      w(ixB^S,mom(3)) =  w(ixB^S,rho_)*dsqrt(1.d0/(x(ixB^S,1)*dsin(x(ixB^S,2)))-3.5d0*hd_adiab)
    elsewhere
      w(ixB^S,rho_)   = small_density
      w(ixB^S,mom(3)) = 0.d0
    end where

    case default
       call mpistop("Special boundary is not defined for this region")
    end select

  end subroutine specialbound_usr

end module mod_usr
