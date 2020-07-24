module mod_usr
  use mod_hd
  implicit none
  double precision :: q, f
  double precision :: a, Egg, Om0, gm2

contains

  subroutine usr_init()
    call set_coordinate_system('spherical_3D')

    usr_set_parameters=> initglobaldata_usr
    usr_init_one_grid => initonegrid_usr
    usr_source        => pt_grav_source
    usr_get_dt        => get_dt_pt_grav
    usr_special_bc    => specialbound_usr
    usr_aux_output    => specialvar_output
    usr_add_aux_names => specialvarnames_output
    usr_refine_grid   => specialrefine_grid

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

    gm2 = 1.d0/q

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

    open(1,file='ini_f50.dat')
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
        w(i,j,:,mom(3))=wcells(i0,4) !-w(i,j,:,rho_)*x(i,j,3,1)*dsin(x(i,j,3,2))*Om0
      enddo
    enddo

    ! where (w(ix^S,rho_)<small_density)
    !   w(ix^S,rho_)   =  small_density
    !   w(ix^S,mom(1)) =  0.d0
    !   w(ix^S,mom(2)) =  0.d0
    !   w(ix^S,mom(3)) =  0.d0
    ! end where

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
    double precision :: cent_str(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir), cor(ixImin1:ixImax1,ixImin2:ixImax2,&
       ixImin3:ixImax3,1:ndir)

    ! call get_cntfgl_str(ixI^L,ixO^L,x,cent_str)
    call get_roche_else(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
       ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,cent_str)
    call get_coriolis(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,ixOmin1,&
       ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,wCT,cor)
    where (wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       rho_)>10.d0*small_density)
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(1)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(1)) - qdt * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_) / x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1)**2.d0
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(1)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(1)) + qdt * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_) * (cent_str(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         1)+cor(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,1))
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(2)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(2)) + qdt * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_) * (cent_str(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         2)+cor(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))
      w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(3)) = w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         mom(3)) + qdt * wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         rho_) * (cent_str(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
         3)+cor(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))
    endwhere

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

    double precision :: force(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       1:ndir), cent_str(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,&
       1:ndir), cor(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3,1:ndir)

    dtnew=bigdouble
    force=0.d0
    force(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)=1.d0/x(ixmin1:ixmax1,&
       ixmin2:ixmax2,ixmin3:ixmax3,1)**2.d0
    ! call get_cntfgl_str(ixG^L,ix^L,x,cent_str)
    call get_roche_else(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,ixmin1,&
       ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,x,cent_str)
    call get_coriolis(ixGmin1,ixGmin2,ixGmin3,ixGmax1,ixGmax2,ixGmax3,ixmin1,&
       ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,x,w,cor)
    force(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       1:ndim)=force(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       1:ndim)+cor(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       1:ndim)+cent_str(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1:ndim)
    ! dtnew=min(dtnew,minval(dsqrt(block%dx(ix^S,1)/(1.d0/x(ix^S,1)**2.d0))))
    dtnew=min(dtnew,minval(dsqrt(block%dx(ixmin1:ixmax1,ixmin2:ixmax2,&
       ixmin3:ixmax3,1)/dabs(force(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       1)))),minval(dsqrt(block%dx(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       1)*block%dx(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       2)/dabs(force(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,2)))),&
       minval(dsqrt(block%dx(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       1)*dsin(block%dx(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       2))*block%dx(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       3)/dabs(force(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,3)))))

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
         rho_)*dsqrt(1.d0/(x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
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
         rho_)*dsqrt(1.d0/(x(ixBmin1:ixBmax1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,&
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


  ! -------------------------------------------------------------------------------
  ! Center @ center of masses
  ! x axis along line joining 2 bodies, from donor to accretor
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine get_cntfgl_str(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,F_Roche)
  use mod_global_parameters
  integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
  double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:ndim)
  double precision, intent(out) :: F_Roche(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:ndir)
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! centrifugal centered on star rotating @ Omega orbital so as the star does not spin
  ! in the inertial frame (but it does wobble around the center of mass, hence the need
  ! for the centrifugal term in the next subroutine)
  F_Roche(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     1)= Om0**2.d0 * x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     1)*(dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)))**2.d0
  F_Roche(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     2)= Om0**2.d0 * x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     1)*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     2))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))
  F_Roche(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)= 0.d0
  end subroutine get_cntfgl_str
  ! -------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------
  ! Rotation axis along +z, w/ y axis pointing upward in slices of the orbital plane
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine get_coriolis(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,wCT,F_Coriolis)
  use mod_global_parameters
  integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
  double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:ndim), wCT(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:nw)
  double precision, intent(out) :: F_Coriolis(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:ndir)
  integer :: i
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  F_Coriolis(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     1) = - 2. * Om0 * (-wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     mom(3))*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)))
  F_Coriolis(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     2) = - 2. * Om0 * (-wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     mom(3))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)))
  F_Coriolis(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     3) = - 2. * Om0 * ( wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     mom(2))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     2))+wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     mom(1))*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)) )
  do i=1,ndir
    F_Coriolis(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       i) = F_Coriolis(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
       i) / wCT(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_  )
  enddo
  end subroutine get_coriolis
  ! -------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------
  ! Center @ center of masses
  ! x axis along line joining 2 bodies, from donor to accretor
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine get_roche_else(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,x,F_Roche)
  use mod_global_parameters
  integer, intent(in) :: ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
      ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
  double precision, intent(in)  :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:ndim)
  double precision, intent(out) :: F_Roche(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:ndim)
  double precision :: xG
  integer :: i
  double precision :: RtoB(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3),&
      dr_smooth

  dr_smooth=0.15d0

 ! RtoB(ixO^S) = dsqrt( x(ixO^S,1)**2.d0 + a**2.d0 - 2.d0*x(ixO^S,1)*a*dsin(x(ixO^S,2))*dcos(x(ixO^S,3)) + 3.d0*block%dx(ixO^S,1)**2.d0 ) ! (a*(x(3,3,4,3)-x(3,3,3,3)))**2.d0 )
  RtoB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = dsqrt( &
     x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     1)**2.d0 + a**2.d0 - 2.d0*x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     ixOmin3:ixOmax3,1)*a*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     ixOmin3:ixOmax3,2))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     ixOmin3:ixOmax3,3)) + 3.d0*dr_smooth**2.d0 )

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ! For centrifugal, w.r.t CM in xG=a*(1/(1+q))
  xG=a/(1.d0+q)

  F_Roche=0.d0
  ! grav of secondary + centrifugal
  F_Roche(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     1)= ( - (gm2*(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     1)-a*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     2))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     3)))) / RtoB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     ixOmin3:ixOmax3)**3.d0 + Om0**2.d0 * (x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     ixOmin3:ixOmax3,1)-xG*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     ixOmin3:ixOmax3,2))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     ixOmin3:ixOmax3,3))-x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     1)*(dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2)))**2.d0) )
  F_Roche(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     2)= ( - (gm2*(          -a*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     ixOmin3:ixOmax3,2))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     ixOmin3:ixOmax3,3)))) / RtoB(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
     ixOmin3:ixOmax3)**3.d0 + Om0**2.d0 * (          &
     -xG*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     2))*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     3))+x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     1)*dcos(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     2))*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,2))  ) )
  F_Roche(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
     3)= ( - (gm2*(           a                 *dsin(x(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,ixOmin3:ixOmax3,3)))) / RtoB(ixOmin1:ixOmax1,&
     ixOmin2:ixOmax2,ixOmin3:ixOmax3)**3.d0 + Om0**2.d0 * (           xG       &
               *dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,3))) )


  end subroutine get_roche_else
  ! -------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------
  subroutine specialvar_output(ixImin1,ixImin2,ixImin3,ixImax1,ixImax2,ixImax3,&
     ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3,w,x,normconv)
  integer, intent(in)                :: ixImin1,ixImin2,ixImin3,ixImax1,&
     ixImax2,ixImax3,ixOmin1,ixOmin2,ixOmin3,ixOmax1,ixOmax2,ixOmax3
  double precision, intent(in)       :: x(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,1:ndim)
  double precision                   :: w(ixImin1:ixImax1,ixImin2:ixImax2,&
     ixImin3:ixImax3,nw+nwauxio)
  double precision                   :: normconv(0:nw+nwauxio)

   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+1) = w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(1))/w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+2) = w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(2))/w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+3) = w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(3))/w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,nw+4) = w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(3))/w(ixOmin1:ixOmax1,&
      ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)+x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      ixOmin3:ixOmax3,1)*dsin(x(ixOmin1:ixOmax1,ixOmin2:ixOmax2,&
      ixOmin3:ixOmax3,2))*Om0
   w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
      nw+5) = dsqrt(w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
      nw+1)**2.d0+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
      nw+2)**2.d0+w(ixOmin1:ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3,&
      nw+3)**2.d0)
  end subroutine specialvar_output
  ! -------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------
  subroutine specialvarnames_output(varnames)
    character(len=*) :: varnames
    varnames='vr vt vp vp_inertial vmag'
  end subroutine specialvarnames_output
  ! -------------------------------------------------------------------------------

  ! -------------------------------------------------------------------------------
  subroutine specialrefine_grid(igrid,level,ixGmin1,ixGmin2,ixGmin3,ixGmax1,&
     ixGmax2,ixGmax3,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3,qt,w,x,refine,&
     coarsen)
    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    ! you must set consistent values for integers refine/coarsen:
    ! refine = -1 enforce to not refine
    ! refine =  0 doesn't enforce anything
    ! refine =  1 enforce refinement
    ! coarsen = -1 enforce to not coarsen
    ! coarsen =  0 doesn't enforce anything
    ! coarsen =  1 enforce coarsen
    use mod_global_parameters
    integer, intent(in) :: igrid, level, ixGmin1,ixGmin2,ixGmin3,ixGmax1,&
       ixGmax2,ixGmax3, ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,ixmax3
    double precision, intent(in) :: qt, w(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:nw), x(ixGmin1:ixGmax1,ixGmin2:ixGmax2,&
       ixGmin3:ixGmax3,1:ndim)
    integer, intent(inout) :: refine, coarsen
    double precision :: RtoB(ixGmin1:ixGmax1,ixGmin2:ixGmax2,ixGmin3:ixGmax3),&
        xrefine
    integer :: ixPmin1,ixPmin2,ixPmin3,ixPmax1,ixPmax2,ixPmax3

    ! Set the center of refinement in xrefine,
    ! a bit ahead of the accretor (orb. sep. minus the accretion radius typically)
    xrefine=a
    RtoB(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3) = dsqrt( x(ixmin1:ixmax1,&
       ixmin2:ixmax2,ixmin3:ixmax3,1)**2.d0 + xrefine**2.d0 - &
       2.d0*x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       1)*xrefine*dsin(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,&
       2))*dcos(x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,3)) )

    ixPmin1=ixmin1+kr(1,1);ixPmin2=ixmin2+kr(1,2);ixPmin3=ixmin3+kr(1,3)
    ixPmax1=ixmax1+kr(1,1);ixPmax2=ixmax2+kr(1,2);ixPmax3=ixmax3+kr(1,3);
    ! Force max refinement for grid which contains the point (xrefine,0,0)
    if (any(RtoB(ixmin1:ixmax1,ixmin2:ixmax2,&
       ixmin3:ixmax3)<2.d0*(x(ixPmin1:ixPmax1,ixPmin2:ixPmax2,ixPmin3:ixPmax3,&
       1)-x(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,1)))) then
      refine = 1
    endif

    ! Never refine @ poles
    if (node(pig2_,igrid)<2**level) refine = -1

  end subroutine specialrefine_grid
  ! -------------------------------------------------------------------------------

end module mod_usr
