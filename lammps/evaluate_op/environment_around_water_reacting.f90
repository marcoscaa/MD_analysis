MODULE analysis
  IMPLICIT NONE
  REAL*8 ::  sw_t(1000000)
  REAL*8 ::  chargew(1000000,2)
  REAL*8, ALLOCATABLE :: dist_mat(:,:)
  INTEGER, ALLOCATABLE :: water_ion_index(:,:)
  INTEGER :: index_frame(1000000) 
  INTEGER :: index_water(1000000,2) 
  INTEGER :: timestep, last_frame
END MODULE analysis

PROGRAM ENVIRONMENT
  USE analysis, only : index_frame, last_frame
  IMPLICIT NONE
  integer :: iostatus, frame, framer

  call initialize 
  iostatus=0
  frame=1
  framer=1

  call sweep_and_find_ionization_events
  
  do
    call read_atom_reduced2(iostatus)
    if ((iostatus/=0).or.(frame>last_frame)) exit
    if (index_frame(framer)==frame) then  
      call compute_dist_matrix
      call find_reacting_water_pair(framer)
      call print_water_environment 
      framer=framer+1
    end if
    frame=frame+1
  end do

  call finish

END PROGRAM ENVIRONMENT

SUBROUTINE initialize
  USE analysis
  IMPLICIT NONE
  CHARACTER(100) :: pos_file, ind_file

  CALL getarg(1, pos_file) 
  CALL getarg(2, ind_file) 

  timestep=0 

  OPEN(unit = 1,file = pos_file)
  OPEN(unit = 2,file = ind_file)
  OPEN(unit = 3,file = "pos_water.xyz")

END SUBROUTINE initialize

SUBROUTINE read_atom_reduced2(iostatus) 
  !Read LAMMPS atom file
  USE parameters, ONLY : pos, natoms, nwater, atype, box
  USE analysis, ONLY : timestep
  IMPLICIT NONE
  INTEGER                    :: iat, ind, atp, iostatus
  REAL*8                     :: box_tmp(2), box_min(3)
 
  READ(1,*,iostat=iostatus)
  READ(1,*) 
  READ(1,*)
  READ(1,*) natoms
  READ(1,*)

  timestep=timestep+1
  nwater=natoms/3 !Assuming only water is present
  CALL allocate_variables

  !Read Box
  DO iat=1,3
    READ(1,*) box_tmp(1), box_tmp(2)
    box(iat) = box_tmp(2)-box_tmp(1)
    box_min(iat)=box_tmp(1)
  END DO

  READ(1,*)

  DO iat = 1, natoms
    READ(1,*,iostat=iostatus) ind, atype(ind), pos(1,ind), pos(2,ind), pos(3,ind)
    pos(:,ind) = pos(:,ind) * box + box_min 
  END DO

END SUBROUTINE READ_ATOM_REDUCED2

SUBROUTINE finish
  IMPLICIT NONE
  CLOSE(1)
  CLOSE(2)
  CLOSE(3)
END SUBROUTINE finish

SUBROUTINE allocate_variables
  USE parameters, only : pos, atype, natoms, nwater
  USE analysis, only : dist_mat
  IMPLICIT NONE
  
  if(allocated(pos)) then
    deallocate(pos,atype,dist_mat)
  end if 

  allocate(pos(3,natoms))
  allocate(atype(natoms))
  allocate(dist_mat(natoms-nwater,nwater))
    
END SUBROUTINE allocate_variables

SUBROUTINE sweep_and_find_ionization_events 
  USE analysis
  IMPLICIT NONE
  INTEGER :: iframe, iostatus, indwater(2), initial
  REAL*8 :: sw, trash

  index_water=0
  iframe=1
  initial=1
  DO

    READ(2,*,iostat=iostatus) index_frame(iframe), sw_t(iframe), trash, &
                              chargew(iframe,1), chargew(iframe,2), &
                              indwater

    IF (iostatus/=0) EXIT

    IF (iframe>1) THEN
      IF ((sw_t(iframe-1)==0.d0).and.(sw_t(iframe)>0.d0)) then
        index_water(initial:iframe,1) = indwater(1)
        index_water(initial:iframe,2) = indwater(2)
      ELSE IF ((sw_t(iframe-1)>0.d0).and.(sw_t(iframe)==0.d0)) then
        initial=iframe
      ELSE 
        index_water(iframe,1) = indwater(1)
        index_water(iframe,2) = indwater(2)
      END IF
    END IF
    iframe=iframe+1


  END DO

  last_frame=iframe

  OPEN(100)
  DO iframe=1,last_frame
    WRITE(100,*) index_frame(iframe), sw_t(iframe), index_water(iframe,:)
  END DO
  CLOSE(100)

END SUBROUTINE sweep_and_find_ionization_events

SUBROUTINE FIND_REACTING_WATER_PAIR(frame)
  !Assigning H atoms for O atoms pertaining to
  !reacting water species
  USE parameters, only : natoms,nwater
  USE analysis
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: frame
  INTEGER :: i, ind1(natoms-nwater), ind2(natoms-nwater)
  INTEGER :: nh1, nh2
  REAL*8 :: d1(natoms-nwater), d2(natoms-nwater), rotmat(3,3)

  IF (.not.allocated(water_ion_index)) then
    allocate(water_ion_index(4,2))
  END IF

  water_ion_index=0

  !Indexes of H atoms
  do i=1,natoms-nwater
    ind1(i)=nwater+i
    ind2(i)=nwater+i
  end do

  !Sort O-H distances for each reacting Ow
  nh1=int(chargew(frame,1))+2
  nh2=int(chargew(frame,2))+2
  d1=Dist_Mat(:,index_water(frame,1))
  d2=Dist_Mat(:,index_water(frame,2))
  CALL Bubble_Sort(d1,ind1,nh1,natoms-nwater)
  CALL Bubble_Sort(d2,ind2,nh2,natoms-nwater)

  !Store the indexes of O and H atoms only for 
  !reacting water
  water_ion_index(1,1)=index_water(frame,1)
  water_ion_index(1,2)=index_water(frame,2)
  water_ion_index(2:nh1+1,1)=ind1(1:nh1)
  water_ion_index(2:nh2+1,2)=ind2(1:nh2)
 
END SUBROUTINE FIND_REACTING_WATER_PAIR

SUBROUTINE PRINT_WATER_ENVIRONMENT
  USE parameters, only : pos,natoms
  USE analysis
  IMPLICIT NONE
  INTEGER :: iat, c, indat(natoms)
  REAL*8, parameter :: rcut=5.0d0
  REAL*8 :: d1(4),d2(4), water_ion(3,natoms)

  CALL translate_water_pair_to_molecular_frame(water_ion,indat)

  c=7
  do iat=1,natoms

    if (.not.any(water_ion_index==iat)) then

      CALL Distance_Vector(pos(:,water_ion_index(1,1)),pos(:,iat),d1)
      CALL Distance_Vector(pos(:,water_ion_index(1,2)),pos(:,iat),d2)
  
      if ((d1(4)<rcut).or.(d2(4)<rcut)) then
        water_ion(:,c) = d1(1:3) !Frame w.r.t. reactive Ow index 1
        indat(c)=iat
        c=c+1
      end if

    end if

  end do

  call project_to_local_frame(water_ion,c-1)

  call print_atoms(water_ion,indat,c-1)

END SUBROUTINE PRINT_WATER_ENVIRONMENT

SUBROUTINE translate_water_pair_to_molecular_frame(water_ion,indat)
  USE parameters, only: pos, natoms
  USE analysis, only: water_ion_index
  IMPLICIT NONE
  REAL*8 :: water_ion(3,natoms)
  INTEGER :: indat(natoms)
  INTEGER :: iat, iw, c
  REAL*8 :: d(4)

  water_ion(:,1) = (/ 0.d0,0.d0,0.d0 /)
  indat(1)=water_ion_index(1,1) !O atom at the center
  indat(6)=water_ion_index(1,2) !O atom
  
  CALL Distance_Vector(pos(:,water_ion_index(1,1)),&
                       pos(:,water_ion_index(1,2)),d)

  water_ion(:,6)=d(1:3)
  
  c=1
  do iw=1,2
    do iat=1,3
      if (water_ion_index(iat+1,iw).ne.0) then
        c=c+1
        CALL Distance_Vector(pos(:,water_ion_index(1,1)),& 
                             pos(:,water_ion_index(iat+1,iw)),d)
        water_ion(:,c) = d(1:3)
        indat(c) = water_ion_index(iat+1,iw)
      end if
    end do
  end do

  if (c>5) then
    print *, 'More Hs than expected for water pair'
    print *, 'STOP'
    STOP
  end if

END SUBROUTINE translate_water_pair_to_molecular_frame

SUBROUTINE print_atoms(ion_pair,indat,n)
  USE parameters, only: nwater
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n, indat(n)
  REAL*8, INTENT(IN) :: ion_pair(3,n)
  INTEGER :: iw

  write(3,fmt="(I5)") n
  write(3,*) 
  do iw=1,n
    IF (indat(iw).le.nwater) then
      write(3,fmt="(A3, 3(F14.8,3X))") "O", ion_pair(:,iw) 
    ELSE
      write(3,fmt="(A3, 3(F14.8,3X))") "H", ion_pair(:,iw) 
    END IF
  end do

END SUBROUTINE print_atoms

SUBROUTINE PROJECT_TO_LOCAL_FRAME (ion_pair, n)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  REAL*8 :: ion_pair(3,n), d,dmax
  REAL*8 :: rotmat(3,3), x(3), y(3), z(3)
  INTEGER :: i,j, indmaxh

  !First, find H farthest from O1
  dmax=0.d0
  do i=2,5
    d=norm2(ion_pair(:,i))
    if (d>dmax) then
      dmax=d
      indmaxh=i
    end if
  end do

  !X axis along O-O
  x = ion_pair(:,6)/norm2(ion_pair(:,6)) 
  !Y is perpendicular to the O(H3O)-O(oh)-H(oh) plane
  CALL CROSSPROD(ion_pair(:,6),ion_pair(:,indmaxh)-ion_pair(:,6),y)
  y = y/norm2(y)
  !Z is perpendicular to X and Y
  CALL CROSSPROD(x,y,z)
  z = z/norm2(z)
  
  DO i = 1,3
    rotmat(1,i)=x(i)
    rotmat(2,i)=y(i)
    rotmat(3,i)=z(i)
  END DO

  do i=1,n
    ion_pair(:,i) = MATMUL(rotmat,ion_pair(:,i))
  end do

END SUBROUTINE PROJECT_TO_LOCAL_FRAME

SUBROUTINE COMPUTE_DIST_MATRIX
  USE parameters, only : natoms, nwater
  USE analysis, only : dist_mat
  IMPLICIT NONE
  REAL*8                     :: Dist
  INTEGER                    :: iw, ih

  DO ih=nwater+1,natoms
    DO iw=1,nwater
      dist_mat(ih-nwater,iw) = Dist(ih,iw)
    END DO
  END DO

END SUBROUTINE COMPUTE_DIST_MATRIX
