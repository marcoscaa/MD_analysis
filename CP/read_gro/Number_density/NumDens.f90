!Z resolved dipole distribution. Length units in nm

MODULE parameters
  IMPLICIT NONE 
  INTEGER, PARAMETER                       :: nlayers = 800
  REAL*8, DIMENSION(:,:,:), ALLOCATABLE    :: pos
  REAL*8, DIMENSION(3)                     :: box
  INTEGER*8, DIMENSION(nlayers)            :: hist
  INTEGER                                  :: natoms, nframes, nequil
  INTEGER                                  :: nwater

END MODULE parameters

PROGRAM ZDF
  USE parameters
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE

  !Remove equilibration part
  DO frame = 1,nequil*(natoms+3)
    READ(1, *)
  END DO
  
  DO frame = nequil+1,nframes
    CALL READ_FRAME
    CALL MAKE_HISTOGRAM
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM ZDF

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
  REAL*8                     :: theta, Angle

  DO i = 1,nwater

    !Index for z layers histogram
    ind_layer = get_index_layer( pos(3,1,i), box(3), nlayers ) 

    !Not consider the particles outised the bins range
    IF ( ind_layer == 0 ) CYCLE

    hist(ind_layer) = hist(ind_layer) + 1

  END DO
 
END SUBROUTINE MAKE_HISTOGRAM

SUBROUTINE PRINT_RESULTS
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i, j
  REAL*8                     :: aver, layer, bin

  OPEN(unit = 2,file = "Number_density.dat")

  bin = box(1)*box(2)*box(3) / float(nlayers) 
  bin = bin * float(nlayers-nequil)

  DO i = 1,nlayers

    layer = box(3) * ( float(i) - 0.5 ) / float(nlayers)
    aver = float(hist(i)) / bin
    
    WRITE(2, fmt = "(F10.5, 3X, F12.7)"), layer, aver 

  END DO

  CLOSE(2)

END SUBROUTINE

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

