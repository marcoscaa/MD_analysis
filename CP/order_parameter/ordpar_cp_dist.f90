!Distribution of a order parameter in specific (hard-coded) layers 
!close to the surface. 

MODULE histogram 
  IMPLICIT NONE 
  INTEGER, PARAMETER                       :: nbins = 200, number_NN=4
  DOUBLE PRECISION, ALLOCATABLE            :: dist_matrix(:,:)
  DOUBLE PRECISION                         :: hist(8,nbins)
  INTEGER                                  :: counter(8,nbins)
  INTEGER, ALLOCATABLE                     :: ind_nn(:,:)

END MODULE histogram 

PROGRAM ORDER
  USE parameters, ONLY : nframes
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = 1,nframes
    CALL READ_POS
    CALL COMPUTE_DISTANCE_MATRIX
    CALL SORT_NEAREST_NEIGHBORS
    CALL COMPUTE_ORDER_PARAMETER
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM ORDER

SUBROUTINE INITIALIZE
  USE histogram,  ONLY : dist_matrix, ind_nn, hist, counter, &
                         number_NN
  USE parameters, ONLY : natoms, pos
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_name, index_file

  CALL getarg(1, file_name)
  CALL getarg(2, index_file)

  CALL READ_INDEX( index_file )

  ALLOCATE(pos(natoms,3))
  ALLOCATE(dist_matrix(natoms,natoms)); dist_matrix=0.d0
  ALLOCATE(ind_nn(natoms,number_NN)); ind_nn = 0
  
  hist = 0.d0
  counter = 0

  OPEN(unit  =  1,file  =  file_name)
 
END SUBROUTINE INITIALIZE

SUBROUTINE PRINT_RESULTS
  USE histogram,  ONLY : nbins, counter
  USE parameters, ONLY : nlayers
  IMPLICIT NONE
  INTEGER                    :: i,j
  DOUBLE PRECISION           :: bin, average(nlayers), &
                                tot_layer(nlayers)

  OPEN(unit = 2,file = "order_layer.dat")

  bin = 4.d0 / dble(nbins)

  DO j = 1,nlayers
    tot_layer(j) = sum(counter(j,:))
  END DO

  DO i = 1,nbins

    average = dble(counter(:,i)) / tot_layer(:)

    WRITE(2,fmt = "(E12.5, 8(3X, E15.7))"), -3.d0 + bin/2. + dble(i)*bin, average(:) 

  END DO

  CLOSE(2)

END SUBROUTINE

SUBROUTINE COMPUTE_ORDER_PARAMETER
  !Compute the average and standard deviation of an order parameter for a given
  !frame 
  USE histogram,  ONLY : hist, nbins, counter
  USE parameters, ONLY : pos, ind_atom, natoms, box
  IMPLICIT NONE
  INTEGER                    :: iat, ind_hist, ind_layer
  INTEGER                    :: get_index_hist, layer_index
  DOUBLE PRECISION           :: tetrahedral_order, op

  DO iat = 1, natoms
  
    !Selecting only water oxygen atoms
    IF ( ind_atom(iat) .eq. 4 ) THEN 

      op = tetrahedral_order(iat)

      ind_hist = get_index_hist(op,nbins) 
      ind_layer = layer_index(pos(iat,3))
      
      counter(ind_layer, ind_hist) = counter(ind_layer, ind_hist) + 1 
      hist(ind_layer, ind_hist) = hist(ind_layer, ind_hist) + op

    END IF

  END DO

END SUBROUTINE COMPUTE_ORDER_PARAMETER

INTEGER FUNCTION get_index_hist(op,nbins)
  !Histogram of order parameter values. Assuming it is restricted to [-3,1]
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: op 
  INTEGER, INTENT(IN)          :: nbins

  get_index_hist =  int( dble(nbins) * ( op + 3.d0 ) / 4.d0 ) + 1

END FUNCTION get_index_hist

SUBROUTINE COMPUTE_DISTANCE_MATRIX
  !Builds a symmetrical distance matrix. In this case we restrict the matrix 
  !only for water oxygen atoms
  USE histogram,  ONLY : dist_matrix
  USE parameters, ONLY : natoms, ind_atom
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
  USE histogram , ONLY : dist_matrix, number_NN, ind_nn
  USE parameters, ONLY : natoms, ind_atom
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

DOUBLE PRECISION FUNCTION tetrahedral_order(index_atom)
  !Compute the tetrahedral order parameter of water oxygen atoms
  USE histogram,  ONLY : ind_nn, number_NN
  USE parameters, ONLY : pos, ind_atom
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
  !tetrahedral_order = 1.d0 - 3.d0/32.d0 * sum_angle

END FUNCTION tetrahedral_order
