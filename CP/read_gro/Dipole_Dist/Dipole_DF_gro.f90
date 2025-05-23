!Z resolved dipole distribution. Length units in nm

MODULE parameters
  IMPLICIT NONE 
  INTEGER, PARAMETER                       :: nbin_h = 30, nlayers = 50
  REAL*8, DIMENSION(:,:,:), ALLOCATABLE    :: pos
  REAL*8, DIMENSION(3)                     :: box
  INTEGER*8, DIMENSION(nlayers,nbin_h)     :: hist
  INTEGER*8, DIMENSION(nlayers)            :: neq
  INTEGER                                  :: natoms, nframes, nequil
  INTEGER                                  :: nwater

END MODULE parameters

PROGRAM DIP_DIST
  USE parameters
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE

  !Remove equilibration part
  DO frame = 1,nequil*(natoms+1)
    READ(1, *)
  END DO
  
  DO frame = nequil+1,nframes
    CALL READ_FRAME
    CALL MAKE_HISTOGRAM
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM DIP_DIST

SUBROUTINE INITIALIZE
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_name, index_file

  !User should provide filename and index_file as arguments 
  CALL getarg(1, file_name)
  CALL getarg(2, index_file)

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)

  !First 2 lines of index file are number of atoms and box
  READ(1, *), natoms, nframes, nequil

  ALLOCATE(pos(3,3,natoms))

  CLOSE(1)

  OPEN(unit  =  1,file  =  file_name)
 
  hist = 0
  neq = 0

END SUBROUTINE INITIALIZE

SUBROUTINE READ_FRAME (frame)
  !Subroutine reads gro files
  USE parameters, ONLY : nwater, natoms, pos, box
  IMPLICIT NONE
  INTEGER                    :: i, junk1, frame, io, index_
  REAL*8                     :: vtmp(3), z_pbc, z
  CHARACTER(5)               :: moltype,atomtype

  io = 0
  !Discarding first two lines
  READ(1,*)
  READ(1,*)

  i = 1
  DO WHILE ( i <= natoms )
    READ(1, fmt="(i5,2a5,i5,3f8.3)"), junk1, moltype, atomtype, junk1, vtmp(:)
    if ( adjustl(trim(moltype)).eq.'SOL' ) then
      io = io + 1
      !Oxygem coordinates
      pos(:,1,io) = vtmp(:)
      !Hydrogen coordinates
      READ(1, fmt="(i5,2a5,i5,3f8.3)"), junk1, moltype, atomtype,junk1,pos(:,2,io)
      READ(1, fmt="(i5,2a5,i5,3f8.3)"), junk1, moltype, atomtype,junk1,pos(:,3,io)
      i = i + 2
    endif
    i = i + 1
  END DO

  nwater = io

  READ(1,*), box(1), box(2), box(3)

END SUBROUTINE

SUBROUTINE MAKE_HISTOGRAM
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i, j, ind_hist, ind_layer
  INTEGER                    :: get_index_layer
  INTEGER, DIMENSION(2)      :: ind_H
  REAL*8                     :: theta, Dipole_z 

  DO i = 1,nwater

    !Index for z layers histogram
    ind_layer = get_index_layer( pos(3,1,i), box(3), nlayers ) 

    !Not consider the particles outised the bins range
    IF ( ind_layer == 0 ) CYCLE

    theta = Dipole_z(pos,i)
    ind_hist = int( (theta + 1.) * DBLE(nbin_h) / 2. ) + 1
    hist(ind_layer,ind_hist) = hist(ind_layer,ind_hist) + 1
    neq(ind_layer) = neq(ind_layer) + 1 ! number of equivalent molecules to be averaged 

  END DO
 
END SUBROUTINE MAKE_HISTOGRAM

SUBROUTINE PRINT_RESULTS
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i, j
  REAL*8                     :: aver(nlayers)

  DO j = 1,nbin_h
    aver = float( hist(:,j) ) / float( neq(:) ) 

    WRITE(*,fmt = "(F10.5, *(3X, F12.7))"), ( float(j) - 0.5 ) *  & 
                                    2. / float(nbin_h) - 1., aver(:)

  END DO

END SUBROUTINE

REAL*8 FUNCTION Dipole_z(pos,ind1)
  !Angle of water dipole with the z axis 
  !ind1: OW index
  USE parameters, ONLY : natoms, box
  IMPLICIT NONE
  REAL*8, INTENT(IN)         :: pos(3, 3, natoms)
  INTEGER, INTENT(IN)        :: ind1
  REAL*8                     :: v1(3), v2(3)
  INTEGER                    :: i

  !Taking the average of the vectors connecting O-H's
  !PBC is now properly considered
  DO i =1,3
    v1(i) = pos(i,2,ind1) - pos(i,1,ind1)
    v1(i) = v1(i) - nint( v1(i) / box(i) ) * box(i)

    v2(i) = pos(i,3,ind1) - pos(i,1,ind1)
    v2(i) = v2(i) - nint( v2(i) / box(i) ) * box(i)
  END DO

  v1 = ( v1 + v2 ) /2.
  v2 = (/ 0., 0., 1. /)

  Dipole_z = dot_product(v1,v2) !/ ( norm2(v1) * norm2(v2) )

END FUNCTION Dipole_z 

INTEGER FUNCTION get_index_layer( z, box_z, nlayers )
  IMPLICIT NONE
  REAL*8, INTENT(IN)         :: z, box_z
  INTEGER, INTENT(IN)        :: nlayers
  INTEGER                    :: i
  REAL*8                     :: z_pbc

  !Applying PBC
  z_pbc = z  - nint(z/box_z)*box_z 
  z_pbc = z_pbc + box_z/2.

  get_index_layer = int( nlayers * z_pbc / box_z ) + 1

END FUNCTION get_index_layer

