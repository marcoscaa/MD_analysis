program z_orient
 implicit none
! Program to determine the reorentational correlation function
! starting from the dipole correlations 

 integer :: nwat               ! number of water molecules
 
 real(8), dimension(128,3) :: oxy, hy1, hy2  ! oxygen and hydrogens coordinates
 real(8), dimension(128,3) :: dh1o, dh2o     ! H1-O and H2-O distances
 real(8), dimension(128,3) :: dipole, u, u0  ! dipole and orientational vectors
 real(8), parameter :: b2ang =0.529177249
 !real(8), allocatable :: L(:)                 ! cell dimension for cubic geometry
 real(8), dimension(128) :: norm             ! array with norms                  

 real(8) :: L                                 ! cell dimension for cubic geometry
 real(8) :: su                                ! sum to build the correlation 
 real(8) :: z_reo                             ! z_reorentational
 integer :: nstep                             ! number of configurations 

 integer :: i, j, k, t
 character(len=3):: at

 open(30, file='hdl.xyz', status='old')
 open(40, file='correlation.dat', status='replace')

 nstep = 96
 !allocate(L(nstep))
 L = 29.6296*b2ang

 oxy=0.0
 hy1=0.0
 hy2=0.0
 do t = 1, nstep

    ! Reading coordinates 
    read(30,*) nwat
    do i = 1, nwat
       read(30,*)at, (oxy(i,j),j=1,3) 
       read(30,*)at, (hy1(i,j),j=1,3)
       read(30,*)at, (hy2(i,j),j=1,3)
    enddo

    ! Calculate the dipoles
    do i = 1, nwat
       do j = 1, 3
          if(oxy(i,j).gt.L/2.0) oxy(i,j)=oxy(i,j)-L
          if(hy1(i,j).gt.L/2.0) hy1(i,j)=hy1(i,j)-L
          if(hy2(i,j).gt.L/2.0) hy2(i,j)=hy2(i,j)-L
       enddo
       ! Calculate H1-O and H2-O distances
       do j = 1, 3
            dh1o(i,j) = hy1(i,j) - oxy(i,j)
            dh2o(i,j) = hy2(i,j) - oxy(i,j)
       enddo
       ! Apply minimum image convention
       dh1o(i,:) = dh1o(i,:) - L * anint(dh1o(i,:)/L)
       dh2o(i,:) = dh2o(i,:) - L * anint(dh2o(i,:)/L)
       ! Calculate dipole vector
       do j = 1, 3
          dipole(i,j) = (0.52/0.2082)*(dh1o(i,j)+dh2o(i,j))
       enddo
       ! Calculate norm
       norm(i) = sqrt(sum(dipole(i,:)**2.0))
       ! Calculate reorentational unit vector
       do j = 1, 3
          u(i,j) = dipole(i,j) / norm(i) 
       enddo
       if(t.eq.1) u0 = u
    enddo
    
    ! Calculate correlations between u(0) and u(t=nstep)
    su = 0.0
    do i = 1, nwat
       do j = 1, 3
          su = su +  u(i,j)*u0(i,j)  
       enddo
    enddo
    z_reo = su/real(nwat)
    write(40,*) t*0.00048378, z_reo

 enddo

 rewind(30); rewind(40)
 close(30); close(40)

7 format(4x,E21.14,4x,E21.14,4x,E21.14)

 stop
end program z_orient
