!Computes the average coordination number of 
!atoms along the axis normal to the surface (z)

MODULE parameters
  IMPLICIT NONE 
  INTEGER                            :: natoms, ntype, nlayers
  INTEGER                            :: nframes, nequil, stride
  INTEGER                            :: atype1, atype2
  INTEGER, PARAMETER                 :: nhist=12 !Max coordination number
  INTEGER, ALLOCATABLE               :: atype(:), hist_coord(:) 
  INTEGER, ALLOCATABLE, DEVICE       :: atyped(:) 
  REAL*8                             :: box(3), d_cut 
  REAL*8, ALLOCATABLE                :: pos(:,:)
  REAL*8, ALLOCATABLE, DEVICE        :: posd(:,:)
  REAL*8, DEVICE                     :: boxd(3)
  !INTEGER                            :: pbc(3) !Apply (1) pbc or not (0)
END MODULE parameters

PROGRAM Coordination
  USE parameters, ONLY : nframes, stride
  IMPLICIT NONE
  INTEGER :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = 1,nframes
    CALL READ_ATOM
    if (MOD(frame,stride)==0) CALL MAKE_HISTOGRAM
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM Coordination

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, ntype, pos, atype, &
                         nhist, hist_coord, nframes, nequil, &
                         d_cut, box, atype1, atype2, stride, &
                         posd, boxd, atyped
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  !User should provide filename and index_file as arguments 
  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, ntype, nframes, nequil, stride
  READ(1,*) atype1, atype2, d_cut !User defined distance cutoff for coordination number
  CLOSE(1)

  !Allocate module arrays
  ALLOCATE(pos(natoms,3))
  ALLOCATE(atype(natoms))
  ALLOCATE(hist_coord(nhist));hist_coord=0
  ALLOCATE(posd(natoms,3))
  ALLOCATE(atyped(natoms))

  OPEN(unit = 1,file = pos_file)

END SUBROUTINE INITIALIZE

SUBROUTINE REMOVE_EQUIL
  USE parameters, ONLY : natoms, nframes, nequil
  IMPLICIT NONE
  INTEGER                    :: iframe

  DO iframe = 1,nequil*(natoms+9)
    READ(1,*)
  END DO

  nframes = nframes - nequil

END SUBROUTINE REMOVE_EQUIL

SUBROUTINE READ_ATOM
  !Read LAMMPS atom file
  USE parameters, ONLY : pos, natoms, atype, box, &
                         posd, boxd, atyped
  IMPLICIT NONE
  INTEGER                    :: iat, ind
  REAL*8                     :: box_tmp(2)
  CHARACTER(16)              :: jnk
 
  DO iat=1,5
    READ(1,*)
  END DO

  !Read Box
  DO iat=1,3
    READ(1,*) box_tmp(1), box_tmp(2)
    box(iat) = box_tmp(2)-box_tmp(1)
  END DO

  READ(1,*)

  !Read atomic positions and indexes
  DO iat = 1, natoms
    READ(1,*) ind, atype(ind), pos(ind,1), pos(ind,2), pos(ind,3)
    pos(ind,:) = pos(ind,:) * box !Non-reduced coordinates
  END DO
 
  posd=pos
  boxd=box
  atyped=atype

END SUBROUTINE READ_ATOM

SUBROUTINE MAKE_HISTOGRAM
  USE parameters, ONLY : natoms, hist_coord, atype, atype1
  IMPLICIT NONE
  INTEGER                     :: iat, jat
  INTEGER                     :: n_coord 

  DO iat = 1, natoms
    IF (atype(iat)==atype1) THEN
      CALL COORD_NUMBER(iat,n_coord)
      hist_coord ( n_coord+1 ) = hist_coord ( n_coord+1 ) + 1
    END IF
  END DO

END SUBROUTINE MAKE_HISTOGRAM

SUBROUTINE PRINT_RESULTS
  USE parameters, ONLY : nhist, hist_coord
  IMPLICIT NONE
  INTEGER                    :: ic
  
  OPEN(unit = 2,file = "hist_coord.dat")

  write(2,*) '# ', 'Coord_number', 'Number'

  DO ic = 1, nhist

    WRITE(2,*) ic-1, hist_coord(ic) 

  END DO

  CLOSE(2)

END SUBROUTINE

SUBROUTINE COORD_NUMBER( iat, n_coord )
  !Compute the coordination number of atom with index iat
  USE parameters, ONLY : natoms, atyped, d_cut, atype2, posd, boxd
  USE cudafor
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: iat 
  INTEGER                    :: iat2
  INTEGER                    :: n_coord, i
  REAL*8                     :: d
  REAL*8, DEVICE             :: xyz(3)

  n_coord = 0

!$cuf kernel do(1)<<<*,*>>>
  DO iat2 = 1, natoms
    IF ( (iat2.ne.iat).and.(atyped(iat2).eq.atype2) ) THEN
      DO i = 1,3
        xyz(i) = posd(iat,i) - posd(iat2,i)
        xyz(i) = xyz(i) - nint( xyz(i)/boxd(i) ) * boxd(i) 
      END DO
      d = SQRT( SUM(xyz*xyz) )
      IF ( d < d_cut )  n_coord = n_coord + 1
    END IF 
  END DO

END SUBROUTINE COORD_NUMBER

attributes(device) REAL*8 FUNCTION Dist(ind1,ind2)
  ! Distance between two points including pbc
  USE parameters, ONLY : posd, boxd 
  IMPLICIT NONE
  REAL*8                     :: xyz(3)
  INTEGER                    :: i, ind1, ind2

  DO i = 1,3

    xyz(i) = posd(ind1,i) - posd(ind2,i)
    xyz(i) = xyz(i) - nint( xyz(i)/boxd(i) ) * boxd(i) 

  END DO

  Dist = SQRT( SUM(xyz*xyz) )

END FUNCTION Dist
