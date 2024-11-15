module FSFortranLib
use,intrinsic::iso_c_binding
implicit none
double precision, dimension(:), allocatable, target :: h
contains
  type(c_ptr) function fastscapeFortran(nx,ny,xl,yl,dt_max,dt_n,nstep,vz,topo, &
  kfSet,kfsedSet,m,n,kdSet,kdsedSet,g, gsed, p, setMarine, sealevel, poro_silt, poro_sand, & 
  zporo_silt, zporo_sand, ratio, Lsolve, kds_silt, kds_sand, ibc) bind(c,NAME="fastscapeFortran") 
    ! all boundaries are at base level
    ! initial random topography
    implicit none 
    integer :: nx, ny, istep, nstep, i, j, ind, ind2, ibc, setMarine
    double precision :: xl, yl, dt_max, dt_n, kfsed, m, n, kdsed, g, gsed, p, kdsedSet, kdSet, kfsedSet, kfSet
    double precision :: sealevel, poro_silt, poro_sand, zporo_silt, zporo_sand, ratio, Lsolve, kds_silt, kds_sand;
    double precision, dimension(:), allocatable :: u, chi, kf, kd
 !   double precision, dimension(:), allocatable, target :: h
    integer, parameter :: MAX_NODE = 1000000
    double precision, bind(c) :: vz(MAX_NODE), topo(MAX_NODE)
  ! double precision, dimension(:), allocatable :: vz, topo
  ! interface
  !   type(c_ptr) function dynamic_array2(length,n)bind(c,name='d_array')
  !   import
  !   implicit none
  !   integer(c_int),intent(in),value::length
  !   integer(c_int),intent(out)::n
  !   end function
  ! end interface

    fastscapeFortran = c_null_ptr

    ! initialize FastScape
    call FastScape_Init ()

    ! set grid size
  !  WRITE(*,*) "nx: ",nx,",ny: ",ny
    call FastScape_Set_NX_NY (nx,ny)

    ! allocate memory
    call FastScape_Setup ()

    ! set model dimensions
  !  WRITE(*,*) "xl: ",xl,",yl: ",yl
    call FastScape_Set_XL_YL (xl,yl)

    ! set time step
  !  WRITE(*,*) "dt_max: ",dt_max ,", nsteps: ",nstep, ", dt_n :", dt_n
    call FastScape_Set_DT (dt_max)

    ! set random initial topography
    allocate (h(nx*ny))
  !  allocate (topo(nx*ny))
  !  topo = MALLOC(nx*ny)

    call random_number (h)

    outer_loop_topo: do j = 0,ny-1
      inner_loop_topo: do i = 0,nx-1
      ind = j*nx+i
  !    PRINT *, "topo[",j,"]","[",i,"]",":", topo(ind2)
      h(ind)= topo(ind)
 !     PRINT *, "topo[",j,"]","[",i,"]",":", h(ind) 
      end do inner_loop_topo
    end do outer_loop_topo

    call FastScape_Init_H (h)

    ! set erosional parameters
    allocate (kf(nx*ny),kd(nx*ny))
    kf = kfSet
    kfsed = kfsedSet
    kd = kdSet
    kdsed = kdsedSet
    call FastScape_Set_Erosional_Parameters (kf, kfsed, m, n, kd, kdsed, g, gsed, p)

    ! set marine transport parameters
    if (setMarine == 1) then
      call FastScape_Set_Marine_Parameters &
         (sealevel, poro_silt, poro_sand, zporo_silt, zporo_sand, ratio, Lsolve, kds_silt, kds_sand)
    endif
    
    ! set uplift rate (uniform while keeping boundaries at base level)
    allocate (u(nx*ny))
  !  allocate (vz(nx*ny))
  !  u = 1.d-3
 
  ! cpp按行储存，fortran按列储存，这个也要核对一下
 

    outer_loop_vz: do j = 0,ny-1
      inner_loop_vz: do i = 0,nx-1
        ind = j*nx+i
     !   PRINT *, "vz[",j,"]","[",i,"]",":", topo(ind2)
        u(ind)= vz(ind)
   !     PRINT *, "vz[",j,"]","[",i,"]",":", u(ind) 
      end do inner_loop_vz
    end do outer_loop_vz

    u(1:nx)=0.d0
    u(nx:nx*ny:nx)=0.d0
    u(1:nx*ny:nx)=0.d0
    u(nx*(ny-1)+1:nx*ny)=0.d0
    call FastScape_Set_U (u)
    ! U in left and right boundary = 0;  hight in upper and lower boundary = 0.

    ! set boundary conditions
    call FastScape_Set_BC (ibc)

    ! set number of time steps and initialize counter istep
    call FastScape_Get_Step (istep)
  ! 看一下istep等于0还是1
    !allocate memory to extract chi
    allocate (chi(nx*ny))

    ! loop on time stepping
    do while (istep<nstep)
      if(istep==nstep .and. dt_n>0.d0) then
        call FastScape_Set_DT (dt_n)
        ! execute step
        call FastScape_Execute_Step()
        ! get value of time step counter
        call FastScape_Get_Step (istep)
        ! extract solution
        call FastScape_Copy_Chi (chi)
        ! create VTK file
!        call FastScape_VTK (chi, 2.d0)
        ! outputs h values
        call FastScape_Copy_h (h)
  !      print*, "dt", dt_n
   !     print*,'step',istep
    !    print*,'h range:',minval(h),sum(h)/(nx*ny),maxval(h)
      endif
      ! execute step
      call FastScape_Execute_Step()
      ! get value of time step counter
      call FastScape_Get_Step (istep)
      ! extract solution
      call FastScape_Copy_Chi (chi)
      ! create VTK file
 !     call FastScape_VTK (chi, 2.d0)
      ! outputs h values
      call FastScape_Copy_h (h)
!      print*, "dt", dt_max
 !     print*,'step',istep
  !    print*,'h range:',minval(h),sum(h)/(nx*ny),maxval(h)
    enddo

 !   outer_loop_print: do j = 1,ny
  !    inner_loop_print: do i = 1,nx
  !    ind = j*nx+i
 !     ind2 = (j-1)*nx + i    
 !       PRINT *, "topo_after_fortran[",j,"]","[",i,"]",":", h(ind) 
  !    end do inner_loop_print
  !  end do outer_loop_print

    fastscapeFortran = c_loc(h(1))

    ! output timing
    call FastScape_Debug()

    ! end FastScape run
    call FastScape_Destroy ()

    deallocate (u,kf,kd,chi)
    return
  end function

  subroutine clear_array()bind(c,name="clearArray")
    implicit none
    if(allocated(h))then
        deallocate(h)
   !     write(*,*)"Clear successfully"
    endif
  end subroutine

end module FSFortranLib
