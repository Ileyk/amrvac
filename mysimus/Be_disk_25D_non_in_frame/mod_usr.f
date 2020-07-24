module mod_usr
  use mod_hd
  implicit none
  double precision :: q, f
  double precision :: a, Egg, Om0

contains

  subroutine usr_init()
    call set_coordinate_system('spherical_2.5D')

    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_source        => pt_grav_source
    usr_get_dt        => get_dt_pt_grav
    usr_special_bc    => specialbound_usr
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output

    call params_read(par_files)

    call hd_activate()

  end subroutine usr_init

  !> Read parameters from a file
  subroutine params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /my_list/ q, f

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       read(unitpar, my_list, end=111)
111    close(unitpar)
    end do

  end subroutine params_read

  subroutine initglobaldata_usr

    Egg = ( 0.49*q**(2./3.) ) / ( 0.6*q**(2./3.) + dlog(1.+q**(1./3.)) )
    a  = 1.d0/(f*Egg) ! in units of Rstar

    ! ang. speed of frame in units of sqrt(GM/R3)
    Om0 = dsqrt((f*Egg)**3.d0*(1.d0+q)/q)
    print*, 'Orb. sep. :', a
    print*, 'Om0 :', Om0

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
         rho_)*dsqrt(1.d0/(x(ixmin1:ixmax1,ixmin2:ixmax2,&
         1)*dsin(x(ixmin1:ixmax1,ixmin2:ixmax2,&
         2)))-3.5d0*hd_adiab)-w(ixmin1:ixmax1,ixmin2:ixmax2,&
         rho_)*x(ixmin1:ixmax1,ixmin2:ixmax2,1)*dsin(x(ixmin1:ixmax1,&
         ixmin2:ixmax2,2))*Om0
    elsewhere
      w(ixmin1:ixmax1,ixmin2:ixmax2,rho_)   =  small_density
      w(ixmin1:ixmax1,ixmin2:ixmax2,mom(3)) =  -w(ixmin1:ixmax1,ixmin2:ixmax2,&
         rho_)*x(ixmin1:ixmax1,ixmin2:ixmax2,1)*dsin(x(ixmin1:ixmax1,&
         ixmin2:ixmax2,2))*Om0
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
    double precision :: cent_str(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir),&
        cor(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir)

    call get_cntfgl_str(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
       ixOmax1,ixOmax2,x,cent_str)
    call get_coriolis(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
       ixOmax2,x,wCT,cor)
    where (wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)>10.d0*small_density)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(1)) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(1)) - qdt * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         rho_)    / x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)**2.d0
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(1)) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(1)) + qdt * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         rho_  ) * (cent_str(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         1)+cor(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1))
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(2)) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(2)) + qdt * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         rho_  ) * (cent_str(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         2)+cor(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,mom(3)) = w(ixOmin1:ixOmax1,&
         ixOmin2:ixOmax2,mom(3)) + qdt * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         rho_  ) * (cent_str(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
         3)+cor(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3))
    endwhere

  end subroutine pt_grav_source

  subroutine get_dt_pt_grav(w,ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,&
     ixmax1,ixmax2,dtnew,dx1,dx2,x)

    integer, intent(in) :: ixGmin1,ixGmin2,ixGmax1,ixGmax2, ixmin1,ixmin2,&
       ixmax1,ixmax2
    double precision, intent(in) :: dx1,dx2, x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       1:ndim)
    double precision, intent(in) :: w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:nw)
    double precision, intent(inout) :: dtnew
    double precision :: force(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:ndir),&
        cent_str(ixGmin1:ixGmax1,ixGmin2:ixGmax2,1:ndir), cor(ixGmin1:ixGmax1,&
       ixGmin2:ixGmax2,1:ndir)

    dtnew=bigdouble
    force=0.d0
    force(ixmin1:ixmax1,ixmin2:ixmax2,1)=1.d0/x(ixmin1:ixmax1,ixmin2:ixmax2,&
       1)**2.d0
    call get_cntfgl_str(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,ixmax1,&
       ixmax2,x,cent_str)
    call get_coriolis(ixGmin1,ixGmin2,ixGmax1,ixGmax2,ixmin1,ixmin2,ixmax1,&
       ixmax2,x,w,cor)
    force(ixmin1:ixmax1,ixmin2:ixmax2,1:ndim)=force(ixmin1:ixmax1,&
       ixmin2:ixmax2,1:ndim)+cor(ixmin1:ixmax1,ixmin2:ixmax2,&
       1:ndim)+cent_str(ixmin1:ixmax1,ixmin2:ixmax2,1:ndim)
    ! dtnew=min(dtnew,minval(dsqrt(block%dx(ix^S,1)/(1.d0/x(ix^S,1)**2.d0))))
    dtnew=min(dtnew,minval(dsqrt(block%dx(ixmin1:ixmax1,ixmin2:ixmax2,&
       1)/dabs(force(ixmin1:ixmax1,ixmin2:ixmax2,1)))),&
       minval(dsqrt(block%dx(ixmin1:ixmax1,ixmin2:ixmax2,&
       1)*block%dx(ixmin1:ixmax1,ixmin2:ixmax2,2)/dabs(force(ixmin1:ixmax1,&
       ixmin2:ixmax2,2)))))

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
         ixBmin2:ixBmax2,rho_)*dsqrt(1.d0/(x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
         1)*dsin(x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
         2)))-3.5d0*hd_adiab)-w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
         rho_)*x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,1)*dsin(x(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2,2))*Om0
    elsewhere
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,rho_)   = small_density
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,mom(3)) = -w(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2,rho_)*x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
         1)*dsin(x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,2))*Om0
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
         ixBmin2:ixBmax2,rho_)*dsqrt(1.d0/(x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
         1)*dsin(x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
         2)))-3.5d0*hd_adiab)-w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
         rho_)*x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,1)*dsin(x(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2,2))*Om0
    elsewhere
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,rho_)   = small_density
      w(ixBmin1:ixBmax1,ixBmin2:ixBmax2,mom(3)) = -w(ixBmin1:ixBmax1,&
         ixBmin2:ixBmax2,rho_)*x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,&
         1)*dsin(x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,2))*Om0
    end where

    case default
       call mpistop("Special boundary is not defined for this region")
    end select

  end subroutine specialbound_usr

  ! -------------------------------------------------------------------------------
  ! Center @ center of masses
  ! x axis along line joining 2 bodies, from donor to accretor
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine get_cntfgl_str(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,x,F_Roche)
  use mod_global_parameters
  integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2
  double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim)
  double precision, intent(out) :: F_Roche(ixImin1:ixImax1,ixImin2:ixImax2,&
     1:ndir)
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! centrifugal centered on star rotating @ Omega orbital so as the star does not spin
  ! in the inertial frame (but it does wobble around the center of mass, hence the need
  ! for the centrifugal term in the next subroutine)
  F_Roche(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)= Om0**2.d0 * x(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,1)*(dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)))**2.d0
  F_Roche(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)= Om0**2.d0 * x(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,1)*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     2))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))
  F_Roche(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3)= 0.d0
  end subroutine get_cntfgl_str
  ! -------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------
  ! Rotation axis along +z, w/ y axis pointing upward in slices of the orbital plane
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine get_coriolis(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,x,wCT,F_Coriolis)
  use mod_global_parameters
  integer, intent(in) :: ixImin1,ixImin2,ixImax1,ixImax2, ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2
  double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2,1:ndim),&
      wCT(ixImin1:ixImax1,ixImin2:ixImax2,1:nw)
  double precision, intent(out) :: F_Coriolis(ixImin1:ixImax1,ixImin2:ixImax2,&
     1:ndir)
  integer :: i
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  F_Coriolis(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     1) = - 2. * Om0 * (-wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     mom(3))*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)))
  F_Coriolis(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     2) = - 2. * Om0 * (-wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     mom(3))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)))
  F_Coriolis(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     3) = - 2. * Om0 * ( wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     mom(2))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))+wCT(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,mom(1))*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)) )
  do i=1,ndir
    F_Coriolis(ixOmin1:ixOmax1,ixOmin2:ixOmax2,i) = F_Coriolis(ixOmin1:ixOmax1,&
       ixOmin2:ixOmax2,i) / wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_  )
  enddo
  end subroutine get_coriolis
  ! -------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------
  subroutine specialvar_output(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
     ixOmax1,ixOmax2,w,x,normconv)
  integer, intent(in)                :: ixImin1,ixImin2,ixImax1,ixImax2,&
     ixOmin1,ixOmin2,ixOmax1,ixOmax2
  double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
     1:ndim)
  double precision                   :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
     nw+nwauxio)
  double precision                   :: normconv(0:nw+nwauxio)
  double precision :: cent_str(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir),&
      cor(ixImin1:ixImax1,ixImin2:ixImax2,1:ndir), pthermal(ixImin1:ixImax1,&
     ixImin2:ixImax2), gradP(ixImin1:ixImax1,ixImin2:ixImax2)

   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+1) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      mom(1))/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+2) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      mom(2))/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+3) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      mom(3))/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+4) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      mom(3))/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,rho_)+x(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,1)*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2))*Om0
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+5) = dsqrt(w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,nw+1)**2.d0+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      nw+2)**2.d0+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+3)**2.d0)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+6) = - 1.d0 / x(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,1)**2.d0
   call get_cntfgl_str(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
      ixOmax2,x,cent_str)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+7) = cent_str(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,1)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+8) = cent_str(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,2)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+9) = cent_str(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,3)
   call get_coriolis(ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,ixOmax1,&
      ixOmax2,x,w,cor)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+10) = cor(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,1)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+11) = cor(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,2)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+12) = cor(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,3)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+13) = cent_str(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,1)+cor(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+14) = cent_str(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,2)+cor(ixOmin1:ixOmax1,ixOmin2:ixOmax2,2)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+15) = cent_str(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,3)+cor(ixOmin1:ixOmax1,ixOmin2:ixOmax2,3)
   call hd_get_pthermal(w,x,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
      ixOmax1,ixOmax2,pthermal)
   call gradient(pthermal,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
      ixOmax1,ixOmax2,1,gradP)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+16) = -gradP(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_)
   call gradient(pthermal,ixImin1,ixImin2,ixImax1,ixImax2,ixOmin1,ixOmin2,&
      ixOmax1,ixOmax2,2,gradP)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+17) = -gradP(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2)/w(ixOmin1:ixOmax1,ixOmin2:ixOmax2, rho_)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+18) = ((w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2, mom(2))**2.d0+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       mom(3))**2.d0) / w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
       rho_)**2.d0) / x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+19) = w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,nw+18) - 1.d0 / x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      1)**2.d0 + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+16)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+20) = w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,nw+19) + w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+13)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,nw+21) = w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,nw+18)+cent_str(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      1)+cor(ixOmin1:ixOmax1,ixOmin2:ixOmax2,1)
   end subroutine specialvar_output
  ! -------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------
  subroutine specialvarnames_output(varnames)
    character(len=*) :: varnames
    varnames='vr vt vp vp_inrtl v'
    varnames=trim(varnames)//' g ctrr ctrt ctrp crr'
    varnames=trim(varnames)//' crt crp nninr nnint nninp'
    varnames=trim(varnames)//' _gdP_r_ovr_rh _gdP_t_ovr_rh'
    varnames=trim(varnames)//' cnt_gm_r all_bt_nn_in_r all v2nnin'
  end subroutine specialvarnames_output
  ! -------------------------------------------------------------------------------

end module mod_usr
