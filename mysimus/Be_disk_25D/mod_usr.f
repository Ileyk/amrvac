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

  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,&
     ixmax1,ixmax2,w,x)

    ! initialize one grid

    integer, intent(in) :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,ixmin2,&
       ixmax1,ixmax2
    double precision, intent(in) :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)

    w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)   =  (x(ixmin1:ixmax1,ixmin2:ixmax2,&
       1)*dsin(x(ixmin1:ixmax1,ixmin2:ixmax2,&
       2)))**(-3.5d0) * dexp(-(1.d0/hd_adiab)*(1.d0/dsin(x(ixmin1:ixmax1,&
       ixmin2:ixmax2,2))-1.d0)/x(ixmin1:ixmax1,ixmin2:ixmax2,1))
    w(ixmin1:ixmax1,ixmin2:ixmax2,mom(1)) =  0.d0
    w(ixmin1:ixmax1,ixmin2:ixmax2,mom(2)) =  0.d0
    where (w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)>small_density)
      w(ixmin1:ixmax1,ixmin2:ixmax2,mom(3)) =  w(ixmin1:ixmax1,ixmin2:ixmax2,&
         rho_)*(1.d0/dsqrt(x(ixmin1:ixmax1,ixmin2:ixmax2,&
         1)*dsin(x(ixmin1:ixmax1,ixmin2:ixmax2,2)))-3.5d0*hd_adiab)
    elsewhere
      w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)   =  small_density
      w(ixmin1:ixmax1,ixmin2:ixmax2,mom(3)) =  0.d0
    end where

  end subroutine initonegrid_usr

  subroutine pt_grav_source(qdt,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,&
     ixOmin2,ixOmax1,ixOmax2,iwmin,iwmax,qtC,wCT,qt,w,x)

    integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2, iwmin,iwmax
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw),&
        x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)

    where (wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)>10.d0*small_density)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(1))=w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(1)) - qdt * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         rho_)    / x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)**2.d0
    endwhere
    ! w(ixO^S,e_ )=w(ixO^S,e_ ) - &
    !      qdt * wCT(ixO^S,mom(1))  / x(ixO^S,1)**2.d0

  end subroutine pt_grav_source

  subroutine get_dt_pt_grav(w,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,&
     ixmax1,ixmax2,dtnew,dx1,dx2,x)

    integer, intent(in) :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,ixmin2,&
       ixmax1,ixmax2
    double precision, intent(in) :: dx1,dx2, x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       1:ndim)
    double precision, intent(in) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)
    double precision, intent(inout) :: dtnew

    dtnew=bigdouble
    dtnew=min(dtnew,minval(dsqrt(block%dx(ixmin1:ixmax1,ixmin2:ixmax2,&
       1)/(1.d0/x(ixmin1:ixmax1,ixmin2:ixmax2,1)**2.d0))))

  end subroutine get_dt_pt_grav

  subroutine specialbound_usr(qt,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixBmin1,&
     ixBmin2,ixBmax1,ixBmax2,iB,w,x)

    ! special boundary types, user defined

    integer, intent(in) :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixBmin1,ixBmin2,&
       ixBmax1,ixBmax2, iB
    double precision, intent(in) :: qt, x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)
    double precision :: rho(ixGmin1:ixGmax1,ixGmin2:ixGmax2)

    select case(iB)
    case(1)

    w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,rho_)   =  (x(ixBmin1:ixBmax1,&
       ixBmin2:ixBmax2,1)*dsin(x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
       2)))**(-3.5d0) * dexp(-(1.d0/hd_adiab)*(1.d0/dsin(x(ixBmin1:ixBmax1,&
       ixBmin2:ixBmax2,2))-1.d0)/x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,1))
    w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,mom(1)) =  0.d0
    w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,mom(2)) =  0.d0
    where (w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,rho_)>small_density)
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,mom(3)) =  w(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2,rho_)*(1.d0/dsqrt(x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
         1)*dsin(x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,2)))-3.5d0*hd_adiab)
    elsewhere
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,rho_)   = small_density
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,mom(3)) = 0.d0
    end where

    case(2)

    w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,rho_)   =  (x(ixBmin1:ixBmax1,&
       ixBmin2:ixBmax2,1)*dsin(x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
       2)))**(-3.5d0) * dexp(-(1.d0/hd_adiab)*(1.d0/dsin(x(ixBmin1:ixBmax1,&
       ixBmin2:ixBmax2,2))-1.d0)/x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,1))
    w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,mom(1)) =  0.d0
    w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,mom(2)) =  0.d0
    where (w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,rho_)>small_density)
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,mom(3)) =  w(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2,rho_)*(1.d0/dsqrt(x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
         1)*dsin(x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,2)))-3.5d0*hd_adiab)
    elsewhere
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,rho_)   = small_density
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,mom(3)) = 0.d0
    end where

    case default
       call mpistop("Special boundary is not defined for this region")
    end select

  end subroutine specialbound_usr

end module mod_usr
