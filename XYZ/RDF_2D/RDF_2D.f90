!Computes the average coordination number of 
!atoms along the axis normal to the surface (z)

MODULE parameters
  IMPLICIT NONE 
  INTEGER                            :: natoms, stride
  INTEGER                            :: nframes, nequil, nhist
  INTEGER                            :: ind_rdf(2), n
  INTEGER, ALLOCATABLE               :: atype(:), H(:,:) 
  REAL*8                             :: box(3), maxr, number_density
  REAL*8, ALLOCATABLE                :: pos(:,:)
  REAL*8                             :: fixed_coord(3)
END MODULE parameters

PROGRAM RDF
  USE parameters, ONLY : nframes, stride, n
  IMPLICIT NONE
  INTEGER :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = 1,nframes
    CALL READ_XYZ
    if (MOD(frame,stride)==0) THEN
      CALL FIND_CENTER_CNT
      CALL COMPUTE_RDF
      n=n+1
    end if
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM RDF

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms,pos, atype, &
                         nhist, H, nframes, nequil, n, &
                         ind_rdf, maxr, stride
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  !User should provide filename and index_file as arguments 
  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nframes, nequil,stride, nhist, maxr
  READ(*,*) ind_rdf !Atomic species to be computed by RDF
                    !0 means all species

  !Allocate module arrays
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(H(nhist,nhist)); H=0

  n=0 !Number of configurations to average the RDF over

  CLOSE(1)
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

SUBROUTINE READ_XYZ
  !Read XYZ file
  USE parameters, ONLY : pos, natoms, atype, box
  IMPLICIT NONE
  INTEGER                    :: iat 
  CHARACTER(4)               :: typ
 
  READ(1,*)
  READ(1,*) box

  DO iat = 1, natoms
    READ(1,*) typ, pos(iat,1), pos(iat,2), pos(iat,3)
    SELECT CASE (TRIM(typ))
    CASE ("C") 
      atype(iat) = 1
    CASE ("O")
      atype(iat) = 2
    CASE ("H")  
      atype(iat) = 3
    CASE ("X")  
      atype(iat) = 4
    END SELECT 
  END DO

END SUBROUTINE READ_XYZ

SUBROUTINE FIND_CENTER_CNT
  USE parameters, only : pos, box, natoms, fixed_coord, atype
  IMPLICIT NONE
  REAL*8 :: sum_coord(3),d(3)
  INTEGER :: i, nc

  nc=0
  DO i = 1,natoms
    IF (atype(i)==1) THEN
      CALL DISTANCE_VECTOR(pos(:,i),fixed_coord,d)
      sum_coord = sum_coord + d(1:3)
      nc = nc+1
    END IF
  END DO

  fixed_coord = fixed_coord + sum_coord/float(nc)

END SUBROUTINE FIND_CENTER_CNT

SUBROUTINE Vector_Distance_2D(a,b,d_out)
  !Vector distance between vectors a and b. 
  USE parameters, only : box
  IMPLICIT NONE
  REAL*8, INTENT(in) :: a(3), b(3)
  REAL*8 :: d_out(2)
  INTEGER :: pol

  DO pol=1,2
    d_out(pol) = a(pol)-b(pol)
    d_out(pol) = d_out(pol) - nint(d_out(pol)/box(pol))*box(pol)
  END DO

END SUBROUTINE Vector_Distance_2D

SUBROUTINE COMPUTE_RDF
  USE parameters, ONLY : natoms, pos, ind_rdf, box, fixed_coord, &
                         nhist, H, atype, maxr
  IMPLICIT NONE
  INTEGER                     :: iat, iat2
  INTEGER                     :: ind1,ind2 
  REAL*8                      :: d(2), central_atom(3), rotmat(2,2)
  REAL*8                      :: v1(2),v2(2)

    central_atom=fixed_coord !pos(:,iat)
    CALL FIND_LOCAL_FRAME_VECTORS(central_atom,pos,natoms,v1,v2)
    CALL ROTATION_MATRIX_2D(v1,v2,rotmat)
    iat=0

    DO iat2 = 1, natoms

      IF ((( atype(iat2) == ind_rdf(2)) .or. ( ind_rdf(2) == 0 ))) THEN

        if (iat==0) then
                iat=1
                CYCLE
        end if
        CALL Vector_Distance_2D(pos(:,iat2),central_atom,d)
        d = MATMUL(rotmat,d)
        d = d + maxr/2. ! avoid negative d
        ind1 = int(float(nhist)*d(1)/maxr) + 1
        ind2 = int(float(nhist)*d(2)/maxr) + 1
        IF ((ind1<=nhist).and.(ind2<=nhist)) then
          H(ind1,ind2) = H(ind1,ind2) + 1
        END IF

      END IF

    END DO

END SUBROUTINE COMPUTE_RDF

SUBROUTINE FIND_LOCAL_FRAME_VECTORS(central_atom,pos,natoms,v1,v2)
  USE parameters, only : atype,ind_rdf
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: natoms
  REAL*8, INTENT(IN) :: central_atom(3),pos(3,natoms)
  REAL*8 :: v1(2),v2(2)
  INTEGER :: iat

  !Take the first C atom as reference
  DO iat=1,natoms
    if (atype(iat)==ind_rdf(1)) then
      CALL VECTOR_DISTANCE_2D(pos(:,iat),central_atom,v1)
      v1=v1/norm2(v1)
      EXIT
    end if
  END DO
  v2(1) = v1(2)
  v2(2) = -v1(1)

END SUBROUTINE FIND_LOCAL_FRAME_VECTORS

SUBROUTINE ROTATION_MATRIX_2D(v1,v2,rotmat)
  IMPLICIT NONE
  REAL*8, INTENT(IN)         :: v1(2),v2(2)
  REAL*8                     :: rotmat(2,2),invmat(2,2)

  invmat(:,1) = v1
  invmat(:,2) = v2
  
  rotmat=invmat !idenpotent
  !call matinv2(invmat,rotmat)

END SUbROUTINE ROTATION_MATRIX_2D

subroutine matinv2(A,B)
  !! Performs a direct calculation of the inverse of a 2 matrix.
  real*8, intent(in) :: A(2,2)   !! Matrix
  real*8             :: B(2,2)   !! Inverse matrix
  real*8             :: detinv

  ! Calculate the inverse determinant of the matrix
  detinv = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))

  ! Calculate the inverse of the matrix
  B(1,1) = +detinv * A(2,2)
  B(2,1) = -detinv * A(2,1)
  B(1,2) = -detinv * A(1,2)
  B(2,2) = +detinv * A(1,1)
end subroutine


SUBROUTINE PRINT_RESULTS
  USE parameters, ONLY : nhist, H, maxr, n, &
                         natoms, box
  IMPLICIT NONE
  REAL*8                     :: bin, aver, vol
  INTEGER                    :: ihist1,ihist2 
  
  OPEN(unit = 2,file = "rdf2D.dat")

  write(2,*) '#', 'X', 'Y', 'RDF'

  bin =  maxr/float(nhist)
  vol =  bin*bin * box(3) !Assumed square

  DO ihist1 = 1, nhist
    DO ihist2 = 1, nhist

      aver = float(H(ihist1,ihist2)) / vol / float(n) 
      WRITE(2,fmt='(3(F12.8,3X))') &
              bin*(float(ihist1)-0.5)-maxr/2.,&
              bin*(float(ihist2)-0.5)-maxr/2.,&
              aver 
    END DO
  END DO

  CLOSE(2)

END SUBROUTINE

SUBROUTINE DISTANCE_VECTOR( coord, center, d ) 
  USE parameters, ONLY : box
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: coord(3), center(3)
  DOUBLE PRECISION             :: d(3)
  INTEGER                      :: ipol
  
  DO ipol = 1,3
    d(ipol) = coord(ipol) - center(ipol)
    d(ipol) = d(ipol) - nint(d(ipol)/box(ipol))*box(ipol)
  END DO

END SUBROUTINE DISTANCE_VECTOR

REAL*8 FUNCTION Dist(x,y)
  ! Distance between two points along x and y axis including pbc
  USE parameters, ONLY : box
  IMPLICIT NONE
  REAL*8, INTENT(IN)         :: x(3),y(3) 
  REAL*8                     :: xyz(2)
  INTEGER                    :: i

  DO i = 1,2

    xyz(i) = x (i)- y(i)
    xyz(i) = xyz(i) - nint( xyz(i)/box(i) ) * box(i)

  END DO

  Dist = SQRT( SUM(xyz*xyz) )

END FUNCTION Dist
