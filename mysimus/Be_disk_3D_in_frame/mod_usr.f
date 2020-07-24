module mod_usr
  use mod_hd
  implicit none

contains

  subroutine usr_init()
    call set_coordinate_system('spherical_3D')

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

  subroutine initonegrid_usr(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,w,x)

    ! initialize one grid

    integer, intent(in) :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
        ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in) :: x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)
    integer :: Ncells, i, j, i0
    double precision, allocatable :: xcells(:,:), wcells(:,:)

    open(1,file='ini.dat')
    ! skip header
    read(1,*)
    read(1,*) Ncells
    read(1,*)
    allocate(xcells(Ncells,2),wcells(Ncells,4))
    do i=1,Ncells
      read(1,*) xcells(i,1),xcells(i,2),wcells(i,1),wcells(i,2),wcells(i,3),&
         wcells(i,4)
    enddo
    close(1)
    
    do i=ixmin1,ixmax1
      do j=ixmin2,ixmax2
        i0=minloc((xcells(:,1)-x(i,j,3,1))**2.d0+(xcells(:,2)-x(i,j,3,&
           2))**2.d0,DIM=1)
        w(i,j,:,rho_)=wcells(i0,1)
        w(i,j,:,mom(1))=wcells(i0,2)
        w(i,j,:,mom(2))=wcells(i0,3)
        w(i,j,:,mom(3))=wcells(i0,4)
      enddo
    enddo

    where (w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_)<small_density)
      w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,rho_)   =  small_density
      w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mom(1)) =  0.d0
      w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mom(2)) =  0.d0
      w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,mom(3)) =  0.d0
    end where

    deallocate(xcells,wcells)

  end subroutine initonegrid_usr

  subroutine pt_grav_source(qdt,ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,&
     ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,iwmin,iwmax,qtC,&
     wCT,qt,w,x)

    integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
        ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3, iwmin,iwmax
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw), x(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndim)
    double precision, intent(inout) :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:nw)

    where (wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_)>10.d0*small_density)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(1))=w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(1)) - qdt * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_)    / x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)**2.d0
    endwhere
    ! w(ixO^S,e_ )=w(ixO^S,e_ ) - &
    !      qdt * wCT(ixO^S,mom(1))  / x(ixO^S,1)**2.d0

  end subroutine pt_grav_source

  subroutine get_dt_pt_grav(w,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
     ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,dtnew,dx1,dx2,dx3,x)

    integer, intent(in) :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
        ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in) :: dx1,dx2,dx3, x(ixGmin1:ixGmax1,&
       ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:ndim)
    double precision, intent(in) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)
    double precision, intent(inout) :: dtnew

    dtnew=bigdouble
    dtnew=min(dtnew,minval(dsqrt(block%dx(ixmin1:ixmax1,ixmin2:ixmax2,&
       ixmin3:ixmax3,1)/(1.d0/x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       1)**2.d0))))

  end subroutine get_dt_pt_grav

  subroutine specialbound_usr(qt,ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,&
     ixGmax3,ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,ixBmax3,iB,w,x)

    ! special boundary types, user defined

    integer, intent(in) :: ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,&
        ixBmin1,ixBmin2,ixBmin3,ixBmax1,ixBmax2,ixBmax3, iB
    double precision, intent(in) :: qt, x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    double precision, intent(inout) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw)
    double precision :: rho(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3)

    select case(iB)
    case(1)

    w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
       rho_)   =  (x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
       1)*dsin(x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
       2)))**(-3.5d0) * dexp(-(1.d0/hd_adiab)*(1.d0/dsin(x(ixBmin1:ixBmax1,&
       ixBmin2:ixBmax2,ixBmin3:ixBmax3,2))-1.d0)/x(ixBmin1:ixBmax1,&
       ixBmin2:ixBmax2,ixBmin3:ixBmax3,1))
    w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,mom(1)) =  0.d0
    w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,mom(2)) =  0.d0
    where (w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
       rho_)>small_density)
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
         mom(3)) =  w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
         rho_)*(1.d0/dsqrt(x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
         1)*dsin(x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
         2)))-3.5d0*hd_adiab)
    elsewhere
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
         rho_)   = small_density
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,mom(3)) = 0.d0
    end where

    case(2)

    w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
       rho_)   =  (x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
       1)*dsin(x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
       2)))**(-3.5d0) * dexp(-(1.d0/hd_adiab)*(1.d0/dsin(x(ixBmin1:ixBmax1,&
       ixBmin2:ixBmax2,ixBmin3:ixBmax3,2))-1.d0)/x(ixBmin1:ixBmax1,&
       ixBmin2:ixBmax2,ixBmin3:ixBmax3,1))
    w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,mom(1)) =  0.d0
    w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,mom(2)) =  0.d0
    where (w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
       rho_)>small_density)
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
         mom(3)) =  w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
         rho_)*(1.d0/dsqrt(x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
         1)*dsin(x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
         2)))-3.5d0*hd_adiab)
    elsewhere
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
         rho_)   = small_density
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,mom(3)) = 0.d0
    end where

    case default
       call mpistop("Special boundary is not defined for this region")
    end select

  end subroutine specialbound_usr

end module mod_usr
