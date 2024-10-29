module FSFortranLib
use,intrinsic::iso_c_binding
implicit none
double precision, dimension(:), allocatable, target :: h
contains
  type(c_ptr) function fastscapeFortran(nx,ny,xl,yl,dt_max,dt_n,nstep,vz,topo) bind(c,NAME="fastscapeFortran") 
    ! all boundaries are at base level
    ! initial random topography
    implicit none 
    integer :: nx, ny, istep, nstep, i, j, ind, ind2
    double precision :: xl, yl, dt_max, dt_n, kfsed, m, n, kdsed, g
    double precision, dimension(:), allocatable :: u, chi, kf, kd
 !   double precision, dimension(:), allocatable, target :: h
    double precision, bind(c) :: vz(1089), topo(1089)
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
    call FastScape_Init_H (h)



    outer_loop_topo: do j = 1,ny
      inner_loop_topo: do i = 1,nx
      ind = j*ny+i
      ind2 = (j-1)*ny + i
  !    PRINT *, "topo[",j,"]","[",i,"]",":", topo(ind2)
      h(ind)= topo(ind2)
  !    PRINT *, "topo[",j,"]","[",i,"]",":", h(ind) 
      end do inner_loop_topo
    end do outer_loop_topo

    ! set erosional parameters
    allocate (kf(nx*ny),kd(nx*ny))
    kf = 2.d-6
    kfsed = -1.d0
    m = 0.6d0
    n = 1.5d0
    kd = 1.d-1
    kdsed = -1.d0
    g = 0.d0
    call FastScape_Set_Erosional_Parameters (kf, kfsed, m, n, kd, kdsed, g, g, -2.d0)

    ! set uplift rate (uniform while keeping boundaries at base level)
    allocate (u(nx*ny))
  !  allocate (vz(nx*ny))
    u = 1.d-3
  ! 确认一下读取cpp的数组，是要从0开始还是1
  ! cpp按行储存，fortran按列储存，这个也要核对一下
    outer_loop_fs: do j = 1,ny
      inner_loop_fs: do i = 1,nx
      ind = j*ny+i
      ind2 = (j-1)*ny + i
  !    PRINT *, "vz[",j,"]","[",i,"]",":", vz(ind2)
      u(ind)= vz(ind2)
  !    PRINT *, "u[",j,"]","[",i,"]",":", u(ind) 
      end do inner_loop_fs
    end do outer_loop_fs

    u(1:nx)=0.d0
    u(nx:nx*ny:nx)=0.d0
    u(1:nx*ny:nx)=0.d0
    u(nx*(ny-1)+1:nx*ny)=0.d0
    call FastScape_Set_U (u)
    ! U in left and right boundary = 0;  hight in upper and lower boundary = 0.

    ! set boundary conditions
    call FastScape_Set_BC (1111)

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
        call FastScape_VTK (chi, 2.d0)
        ! outputs h values
        call FastScape_Copy_h (h)
        print*,'step',istep
        print*,'h range:',minval(h),sum(h)/(nx*ny),maxval(h)
      endif
      ! execute step
      call FastScape_Execute_Step()
      ! get value of time step counter
      call FastScape_Get_Step (istep)
      ! extract solution
      call FastScape_Copy_Chi (chi)
      ! create VTK file
      call FastScape_VTK (chi, 2.d0)
      ! outputs h values
      call FastScape_Copy_h (h)
      print*,'step',istep
      print*,'h range:',minval(h),sum(h)/(nx*ny),maxval(h)
    enddo

    outer_loop_print: do j = 1,ny
      inner_loop_print: do i = 1,nx
      ind = j*nx+i
 !     ind2 = (j-1)*ny + i    
 !       PRINT *, "topo[",j,"]","[",i,"]",":", h(ind) 
      end do inner_loop_print
    end do outer_loop_print

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
        write(*,*)"Clear successfully"
    endif
  end subroutine

end module FSFortranLib