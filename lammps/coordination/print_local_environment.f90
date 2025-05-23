!Computes the average coordination number of 
!atoms along the axis normal to the surface (z)

MODULE parameters
  IMPLICIT NONE 
  INTEGER                            :: natoms, ntype, nlayers
  INTEGER                            :: nframes, nequil, stride
  INTEGER                            :: atype1, atype2
  INTEGER, PARAMETER                 :: nhist=11 !Max coordination number
  INTEGER, ALLOCATABLE               :: atype(:), hist_coord(:) 
  REAL*8                             :: local_coord(nhist,3)
  REAL*8                             :: box(3), d_cut 
  REAL*8, ALLOCATABLE                :: pos(:,:)
  INTEGER                            :: pbc(3) !Apply (1) pbc or not (0)
END MODULE parameters

PROGRAM Coordination
  USE parameters, ONLY : nframes, stride
  IMPLICIT NONE
  INTEGER :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = 1,nframes
    CALL READ_ATOM
    if (MOD(frame,stride)==0) CALL SELECT_LOCAL_ENVIRONMENT
  END DO

  CALL PRINT_RESULTS

  CLOSE(1);CLOSE(2)

END PROGRAM Coordination

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, ntype, pos, atype, &
                         nhist, hist_coord, nframes, nequil, &
                         d_cut, box, atype1, atype2, stride
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

  OPEN(unit = 1,file = pos_file)
  OPEN(unit = 2,file = 'local_structure.xyz')

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
  USE parameters, ONLY : pos, natoms, atype, box, pbc
  IMPLICIT NONE
  INTEGER                    :: iat, ind
  REAL*8                     :: box_tmp(2)
  CHARACTER(2)               :: apply_pbc(3)
  CHARACTER(16)              :: jnk
 
  DO iat=1,4
    READ(1,*)
  END DO

  !Check if we should apply pbc or not
  !pbc=0
  pbc=1
  READ(1, fmt = "(A16,3(1X,A2))") jnk, apply_pbc
  !DO iat=1,3
  !  IF (apply_pbc(iat).eq.'pp') pbc(iat)=1
  !END DO

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

END SUBROUTINE READ_ATOM

SUBROUTINE SELECT_LOCAL_ENVIRONMENT
  USE parameters, ONLY : natoms, hist_coord, atype, atype1
  IMPLICIT NONE
  INTEGER                     :: iat, jat
  INTEGER                     :: n_coord 

  DO iat = 1, natoms
    IF (atype(iat)==atype1) THEN
      CALL COORD_NUMBER(iat,n_coord)
      hist_coord ( n_coord+1 ) = hist_coord ( n_coord+1 ) + 1
      if (n_coord==4) CALL PRINT_LOCAL_STRUCTURE(n_coord)
    END IF
  END DO

END SUBROUTINE SELECT_LOCAL_ENVIRONMENT

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
  USE parameters, ONLY : natoms, pos, atype,atype2, d_cut, local_coord
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: iat 
  INTEGER                    :: iat2
  INTEGER                    :: n_coord
  REAL*8                     :: d(4)

  n_coord = 0
  local_coord=0

  DO iat2 = 1, natoms
    IF ( (iat2.ne.iat).and.(atype(iat2).eq.atype2) ) THEN
      CALL DISTANCE_VECTOR( iat, iat2, d )
      IF ( d(4) < d_cut ) then
        n_coord = n_coord + 1
        local_coord(n_coord,:) = d(1:3)
      END IF 
    END IF 
  END DO

END SUBROUTINE COORD_NUMBER

SUBROUTINE PRINT_LOCAL_STRUCTURE(n)
  USE parameters, ONLY: pos, local_coord
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: n
  INTEGER                    :: iat

  WRITE(2,fmt='(I2)') n+1
  WRITE(2,*) ''
  WRITE(2,fmt = '(A2, 3(3X,F12.8))') 'Ti', 0.0, 0.0, 0.0

  DO iat = 1,n
    WRITE(2,fmt = '(A2, 3(3X,F12.8))'), 'O', local_coord(iat,1), &
                            local_coord(iat,2), local_coord(iat,3)
  END DO

END SUBROUTINE PRINT_LOCAL_STRUCTURE

REAL*8 FUNCTION Dist(ind1,ind2)
  ! Distance between two points including pbc
  USE parameters, ONLY : pos, box, pbc
  IMPLICIT NONE
  REAL*8                     :: xyz(3)
  INTEGER                    :: i, ind1, ind2

  DO i = 1,3

    xyz(i) = pos(ind1,i) - pos(ind2,i)
    xyz(i) = xyz(i) - nint( xyz(i)/box(i) ) * box(i) * pbc(i)

  END DO

  Dist = SQRT( SUM(xyz*xyz) )

END FUNCTION Dist

SUBROUTINE DISTANCE_VECTOR( iat1, iat2, d ) 
  USE parameters, ONLY : box, pos
  IMPLICIT NONE
  INTEGER, INTENT(IN)          :: iat1, iat2
  DOUBLE PRECISION             :: d(4)
  INTEGER                      :: ipol
  
  DO ipol = 1,3
    d(ipol) = pos(iat2,ipol) - pos(iat1,ipol)
    d(ipol) = d(ipol) - nint(d(ipol)/box(ipol))*box(ipol)
  END DO

  d(4) = SQRT( d(1)*d(1) + d(2)*d(2) + d(3)*d(3) )

END SUBROUTINE DISTANCE_VECTOR
