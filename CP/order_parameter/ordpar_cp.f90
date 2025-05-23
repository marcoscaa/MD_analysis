MODULE parameters
  IMPLICIT NONE 
  INTEGER, PARAMETER                       :: nbins = 250, number_NN=4
  DOUBLE PRECISION, ALLOCATABLE            :: pos(:,:), dist_matrix(:,:)
  DOUBLE PRECISION                         :: box(9), hist(nbins)
  INTEGER                                  :: counter(nbins)
  INTEGER, ALLOCATABLE                     :: ind_atom(:)
  INTEGER, ALLOCATABLE                     :: ind_nn(:,:)
  INTEGER                                  :: natoms, nframes, ntype, nequil

END MODULE parameters

PROGRAM ORDER
  USE parameters
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = nequil+1,nframes
    CALL READ_FRAME
    CALL COMPUTE_DISTANCE_MATRIX
    CALL SORT_NEAREST_NEIGHBORS
    CALL COMPUTE_ORDER_PARAMETER
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM ORDER

SUBROUTINE INITIALIZE
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_name, index_file

  CALL getarg(1, file_name)
  CALL getarg(2, index_file)
  ntype = 4

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)

  !First 2 lines of index file are number of atoms and box
  READ(1, *), natoms, nframes, nequil
  READ(1, *), box
  
  ALLOCATE(pos(natoms,3))
  ALLOCATE(ind_atom(natoms))
  ALLOCATE(dist_matrix(natoms,natoms)); dist_matrix=0.d0
  ALLOCATE(ind_nn(natoms,number_NN)); ind_nn = 0
  
  hist = 0.d0
  counter = 0

  !Read the index file
  DO i = 1,natoms
    READ(1,fmt = "(I2)"), ind_atom(i)  
  END DO
 
  CLOSE(1)

  OPEN(unit  =  1,file  =  file_name)
 
END SUBROUTINE INITIALIZE

SUBROUTINE READ_FRAME
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i

  READ(1,*)

  DO i = 1,natoms
    READ(1, fmt='(3(3X,E22.14))'), pos(i,1), pos(i,2), pos(i,3)
  END DO

  !Bohr to angstrom
  pos = pos * 0.529177

END SUBROUTINE

SUBROUTINE REMOVE_EQUIL
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i

  DO i = 1,nequil*(natoms+1)
    READ(1,*)
  END DO

END SUBROUTINE REMOVE_EQUIL

SUBROUTINE PRINT_RESULTS
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i
  DOUBLE PRECISION           :: bin, average

  OPEN(unit = 2,file = "order.dat")

  bin = box(9) / dble(nbins) 

  DO i = 1,nbins

    IF (counter(i) .ne. 0) THEN
      average = hist(i) / dble(counter(i))
      WRITE(2,fmt = "(E12.5, E15.7)"), bin/2.d0 + DBLE(i-1)*bin, average 
    END IF

  END DO

  CLOSE(2)

END SUBROUTINE

SUBROUTINE COMPUTE_ORDER_PARAMETER
  !Compute the average and standard deviation of an order parameter for a given
  !frame 
  USE parameters, ONLY : pos, ind_atom, hist, natoms, nbins, box, counter
  IMPLICIT NONE
  INTEGER                    :: iat, ind_hist
  INTEGER                    :: get_index_hist
  DOUBLE PRECISION           :: tetrahedral_order, op

  DO iat = 1, natoms
  
    !Selecting only water oxygen atoms
    IF ( ind_atom(iat) .ne. 4 ) CYCLE

    op = tetrahedral_order(iat)

    ind_hist = get_index_hist(pos(iat,3),box(9),nbins) 
    
    counter(ind_hist) = counter(ind_hist) + 1 
    hist(ind_hist) = hist(ind_hist) + op

  END DO

END SUBROUTINE COMPUTE_ORDER_PARAMETER

INTEGER FUNCTION get_index_hist(z,boxz,nbins)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: z, boxz
  DOUBLE PRECISION             :: z_tmp
  INTEGER, INTENT(IN)          :: nbins

  z_tmp = z - nint(z/boxz)*boxz + boxz/2.d0
  
  get_index_hist =  int( dble(nbins) * z_tmp / boxz ) + 1

END FUNCTION get_index_hist

DOUBLE PRECISION FUNCTION tetrahedral_order(index_atom)
  !Compute the tetrahedral order parameter of water oxygen atoms
  USE parameters, ONLY : pos, ind_atom, ind_nn, number_NN
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: index_atom
  INTEGER                    :: in1, in2
  DOUBLE PRECISION           :: sum_angle, Angle

  sum_angle = 0.d0

  DO in1 = 1, number_NN-1
    DO in2 = in1+1, number_NN

      sum_angle = sum_angle &
                + ( Angle(index_atom,ind_nn(index_atom,in1),ind_nn(index_atom,in2)) &
                + 1.d0/3.d0 ) ** 2

    END DO
  END DO

  tetrahedral_order = 1.d0 - 3.d0/8.d0 * sum_angle

END FUNCTION tetrahedral_order

SUBROUTINE COMPUTE_DISTANCE_MATRIX
  !Builds a symmetrical distance matrix. In this case we restrict the matrix 
  !only for water oxygen atoms
  USE parameters, ONLY : dist_matrix, natoms, ind_atom
  IMPLICIT NONE
  INTEGER                              :: iat1, iat2
  DOUBLE PRECISION                     :: Dist
  
  DO iat1 = 1, natoms
  
    !Oxygen atoms
    IF ( ( ind_atom(iat1) .ne. 4 ) .and. ( ind_atom(iat1) .ne. 3 ) ) CYCLE

    DO iat2 = iat1+1, natoms
 
      !Water oxygen and TiO2 oxygen
      IF ( ( ind_atom(iat2) == 4 ) .or. ( ind_atom(iat2) == 3 ) ) THEN

        dist_matrix(iat1,iat2) = Dist(iat1,iat2)
        dist_matrix(iat2,iat1) = dist_matrix(iat1,iat2)
 
      END IF

    END DO

  END DO

END SUBROUTINE COMPUTE_DISTANCE_MATRIX

SUBROUTINE SORT_NEAREST_NEIGHBORS
  !Returns the indexes of the number_NN nearest neighbors of all atoms
  !In this case, we are selecting nearest neighbors of water O only
  USE parameters, ONLY : dist_matrix, natoms, ind_atom, number_NN, ind_nn
  IMPLICIT NONE
  INTEGER                    :: iat1, iat2, ind_tmp(number_NN)
  DOUBLE PRECISION           :: dist_nn(number_NN)
 
  DO iat1 = 1, natoms

    dist_nn = 1000.d0
    ind_tmp = 0 

    !Selecting only water oxygen atoms
    IF ( ind_atom(iat1) .ne. 4 ) CYCLE

    DO iat2 = 1, natoms !Could be more efficient
  
      !Selecting only oxygen atoms
      IF ( ( ind_atom(iat2) .ne. 4 ) .and. ( ind_atom(iat2) .ne. 3 ) ) CYCLE
      IF ( iat1 == iat2 ) CYCLE

      IF ( dist_matrix(iat1, iat2) < dist_nn(4) ) THEN
      
        dist_nn(4) = dist_matrix(iat1, iat2)
        ind_tmp(4) = iat2
        CALL BUBBLE_SORT(dist_nn, ind_tmp, number_NN) 

      END IF

    END DO

    ind_nn(iat1,:) = ind_tmp

  END DO

END SUBROUTINE SORT_NEAREST_NEIGHBORS

SUBROUTINE BUBBLE_SORT(array,indexes,length)
  !Bubble sort an array of length length
  !Ruturns array with smallest value on index 1 
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: length
  DOUBLE PRECISION           :: array(length)
  INTEGER                    :: indexes(length)
  INTEGER                    :: ind, indtmp, i, j
  DOUBLE PRECISION           :: valuetmp

  DO i = 2, length
    DO j = i, length

      IF ( array(j-1) > array(j) ) THEN
        indtmp = indexes(j)
        valuetmp = array(j)
        array(j) = array(j-1)
        indexes(j) = indexes(j-1)
        indexes(j-1) = indtmp
        array(j-1) = valuetmp
      END IF

    END DO
  END DO

END SUBROUTINE BUBBLE_SORT

DOUBLE PRECISION FUNCTION Dist(k, j)
  ! Distance between two points with indexes k and j. Including pbc 
  USE parameters, ONLY : pos, box
  IMPLICIT NONE
  INTEGER, INTENT(IN)                  :: k,j
  INTEGER                              :: i
  DOUBLE PRECISION                     :: xyz(3)

  DO i = 1,3

    xyz(i) = pos(k,i) - pos(j,i)
    xyz(i) = xyz(i) - nint( xyz(i)/box(4*i-3) ) * box(4*i-3)

  END DO

  Dist = SQRT( SUM(xyz*xyz) )

END FUNCTION Dist

DOUBLE PRECISION FUNCTION Angle(ind1, ind2, ind3)
  !Angle between particles with index ind1-ind2 and ind1-ind3
  !Return value in cos(Angle)
  USE parameters, ONLY : natoms, pos, box
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: ind1, ind2, ind3
  INTEGER                    :: ipol
  DOUBLE PRECISION           :: v1(3), v2(3)

  DO ipol = 1,3

    v1(ipol) = pos(ind2,ipol) - pos(ind1,ipol)
    v1(ipol) = v1(ipol) - nint( v1(ipol) / box(4*ipol-3) ) * box(4*ipol-3)
    v2(ipol) = pos(ind3,ipol) - pos(ind1,ipol)
    v2(ipol) = v2(ipol) - nint( v2(ipol) / box(4*ipol-3) ) * box(4*ipol-3)

  END DO

  Angle = dot_product(v1,v2) / ( norm2(v1) * norm2(v2) )

END FUNCTION Angle
