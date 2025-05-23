!Distribution of a order parameter in specific (hard-coded) layers 
!close to the surface. 

MODULE steinhardt 
  IMPLICIT NONE 
  INTEGER, PARAMETER                       :: nbins = 200, number_NN=4, l=6
  INTEGER, PARAMETER                       :: excl_nn=0
  DOUBLE PRECISION, ALLOCATABLE            :: dist_matrix(:,:,:), hist(:)
  INTEGER, ALLOCATABLE                     :: counter(:,:), ind_nn(:,:)

END MODULE steinhardt 

PROGRAM ORDER
  USE parameters, ONLY : nframes
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = 1,nframes
    CALL READ_POS_XYZ
    CALL COMPUTE_DISTANCE_MATRIX
    CALL SORT_NEAREST_NEIGHBORS
    CALL COMPUTE_ORDER_PARAMETER
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM ORDER

SUBROUTINE INITIALIZE
  USE steinhardt,  ONLY : dist_matrix, ind_nn, hist, &
                         counter, number_NN, nbins
  USE parameters, ONLY : natoms, pos, nlayers, nequil, &
                         nframes
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_name, index_file

  CALL getarg(1, file_name)
  CALL getarg(2, index_file)

  !Read index file 
  OPEN(unit = 1,file = index_file)
  READ(1, *), natoms, nframes, nequil, nlayers
  CLOSE(1)

  ALLOCATE(pos(natoms,3))
  ALLOCATE(dist_matrix(natoms,natoms,4)); dist_matrix=0.d0
  ALLOCATE(ind_nn(natoms,number_NN)); ind_nn = 0
  ALLOCATE(counter(nlayers,nbins)); counter = 0
  ALLOCATE(hist(nlayers)); hist = 0 
  
  OPEN(unit  =  1,file  =  file_name)
 
END SUBROUTINE INITIALIZE

SUBROUTINE PRINT_RESULTS
  USE steinhardt,  ONLY : nbins, counter, hist
  USE parameters, ONLY : nlayers, box
  IMPLICIT NONE
  INTEGER                    :: i,j
  INTEGER                    :: tot_layer(nlayers)
  DOUBLE PRECISION           :: bin1,bin2,average(nlayers)
  DOUBLE PRECISION           :: z 
  
  OPEN(unit = 2,file = "order_layer.dat")
  OPEN(unit = 3,file = "average_order_layer.dat")

  bin1 = box(3,3) / dble(nlayers)
  bin2 = 1.d0 / dble(nbins)

  DO j = 1,nlayers
    tot_layer(j) = sum(counter(j,:))
    if (tot_layer(j).eq.0) tot_layer(j)=1

    average(1) = hist(j) / dble(tot_layer(j)) 
    z = bin1/2. + dble(j)*bin1
    z = z - nint(z/box(3,3))*box(3,3) + box(3,3)/2.
    WRITE(3,fmt = "(E12.5, *(3X, E15.7))"), z, average(1) 
  END DO

  DO i = 1,nbins
    average = dble(counter(:,i)) / dble(tot_layer(:))
    WRITE(2,fmt = "(E12.5, *(3X, E15.7))"), bin2/2. + dble(i)*bin2, average(:) 
  END DO

  CLOSE(2)
  CLOSE(3)

END SUBROUTINE

SUBROUTINE COMPUTE_ORDER_PARAMETER
  !Compute the average and standard deviation of an order parameter for a given
  !frame 
  USE steinhardt,  ONLY : hist, nbins, counter
  USE parameters, ONLY : pos, natoms, box
  IMPLICIT NONE
  INTEGER                    :: iat, ind_hist, ind_layer
  INTEGER                    :: get_index_hist, layer_index_equal
  DOUBLE PRECISION           :: Q_local, op

!$omp parallel do private(ind_hist,ind_layer,op)
  DO iat = 1, natoms
  
    op = Q_local(iat)
    ind_hist = get_index_hist(op,nbins) 
    ind_layer = layer_index_equal(pos(iat,3))
    
    counter(ind_layer, ind_hist) = counter(ind_layer, ind_hist) + 1 
    hist(ind_layer) = hist(ind_layer) + op

  END DO
!$omp end parallel do

END SUBROUTINE COMPUTE_ORDER_PARAMETER

INTEGER FUNCTION get_index_hist(op,nbins)
  !Histogram of order parameter values. Assuming it is restricted to [0,1]
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: op 
  INTEGER, INTENT(IN)          :: nbins

  get_index_hist =  int( dble(nbins) * op ) + 1

END FUNCTION get_index_hist

SUBROUTINE COMPUTE_DISTANCE_MATRIX
  !Builds a symmetrical distance matrix
  USE steinhardt,  ONLY : dist_matrix
  USE parameters, ONLY : natoms, ind_atom
  IMPLICIT NONE
  INTEGER                              :: iat1, iat2, ipol
  DOUBLE PRECISION                     :: Dist, Dist_cart
  
  DO iat1 = 1, natoms
  
    DO iat2 = iat1+1, natoms
 
      DO ipol = 1,3
        dist_matrix(iat1,iat2,ipol) = Dist_cart(iat1,iat2,ipol)
        dist_matrix(iat2,iat1,ipol) = dist_matrix(iat1,iat2,ipol)
      END DO

      dist_matrix(iat1,iat2,4) = Dist(iat1,iat2)
      dist_matrix(iat2,iat1,4) = dist_matrix(iat1,iat2,4)
 
    END DO

  END DO

END SUBROUTINE COMPUTE_DISTANCE_MATRIX

SUBROUTINE SORT_NEAREST_NEIGHBORS
  !Returns the indexes of the `number_NN nearest` neighbors for all atoms
  USE steinhardt , ONLY : dist_matrix, number_NN, ind_nn
  USE parameters, ONLY : natoms, ind_atom
  IMPLICIT NONE
  INTEGER                    :: iat1, iat2, ind_tmp(number_NN)
  DOUBLE PRECISION           :: dist_nn(number_NN)
 
  DO iat1 = 1, natoms

    dist_nn = 12.d0 !should be the second rdf minima
    ind_tmp = 0 

    DO iat2 = 1, natoms 
  
      IF (iat2==iat1) CYCLE
      IF ( dist_matrix(iat1, iat2,4) < dist_nn(number_NN) ) THEN
      
        dist_nn(number_NN) = dist_matrix(iat1, iat2,4)
        ind_tmp(number_NN) = iat2
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

! ==================================================================================
!  Importing from Fausto's code
! ==================================================================================

DOUBLE PRECISION FUNCTION Q_local(i)

  USE parameters, ONLY : natoms, pi
  USE steinhardt,  ONLY : l
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: i
  INTEGER :: m

  DOUBLE PRECISION :: rey(-l:l), imy(-l:l)
  DOUBLE PRECISION :: sum_loc
  DOUBLE PRECISION :: Q_tmp

  CALL bond_order(i,rey,imy)
  sum_loc = 0.0d0
  DO m = -l, l
     sum_loc = sum_loc + (rey(m)**2 + imy(m)**2)
  END DO
  Q_tmp =  DSQRT(4*pi*sum_loc / (2*l+1))

  Q_local = Q_tmp

END FUNCTION Q_local

! ----------------------------------------------------------------------------------

SUBROUTINE bond_order(i,rey,imy)
  
  USE parameters, ONLY: natoms, pi
  USE steinhardt,  ONLY: number_NN, ind_NN, dist_matrix, excl_nn, l
  IMPLICIT NONE
  
  INTEGER :: i, j, k, m, ipol
  
  DOUBLE PRECISION :: phi, cosmphi, sinmphi
  DOUBLE PRECISION :: rey(-l:l), imy(-l:l)
  DOUBLE PRECISION :: re, im
  DOUBLE PRECISION :: rij(3) 
  DOUBLE PRECISION :: fact, legendre

  EXTERNAL fact, legendre

  DO m = -l, l
     rey(m) = 0
     imy(m) = 0
  END DO
  
  DO k = excl_nn+1, number_NN
     j = ind_NN(i, k)

     DO ipol = 1,3
       rij(ipol) = dist_matrix(j,i,ipol) / dist_matrix(j,i,4)
     END DO
    
     phi  = ATAN2(rij(2), rij(1))

     DO m = 0, l
        cosmphi = DCOS(m*phi)
        sinmphi = DSIN(m*phi)
        
        re = DSQRT((2*l+1)*fact(l-m) / (4*pi*fact(l+m)))*legendre(m, rij(3))*cosmphi
        im = DSQRT((2*l+1)*fact(l-m) / (4*pi*fact(l+m)))*legendre(m, rij(3))*sinmphi
        rey(m) = rey(m) + re
        imy(m) = imy(m) + im
        
        IF (m .NE. 0) THEN
           rey(-m) = rey(-m) + (-1)**m * re
           imy(-m) = imy(-m) - (-1)**m * im
        END IF
        
     END DO 
  END DO                   
  
  DO m = -l, l
     rey(m) = rey(m) / DBLE(number_NN-excl_nn)
     imy(m) = imy(m) / DBLE(number_NN-excl_nn)
  END DO
  
  RETURN
END SUBROUTINE bond_order

! ----------------------------------------------------------------------
! Function to calculate legendre

FUNCTION legendre(m, x)

  USE steinhardt, ONLY : l
  IMPLICIT NONE
  
  INTEGER :: k, m, mtmp, m1
  
  DOUBLE PRECISION :: x, legendre
  DOUBLE PRECISION :: suml
  DOUBLE PRECISION :: fact
  
  EXTERNAL fact
  
  mtmp = ABS(m)
  m1   = INT((l-mtmp)/2)
  
  suml  = 0.0d0
  
  DO k = 0, m1
     suml = suml + (-1)**k*fact(2*(l-k))*x**(l-2*k-mtmp) / (2**l*fact(k)*fact(l-k)*fact(l-2*k-mtmp))
  END DO
  
  legendre = ((-1)**mtmp)*((1-x**2)**(1.0d0*mtmp/2))*suml
  
  RETURN
END FUNCTION legendre

! --------------------------------------------------
! Function to calculate factorial

FUNCTION fact(k)
  IMPLICIT NONE
  
  INTEGER :: i, k
  DOUBLE PRECISION :: fact
  
  IF (k .LE. 1) THEN
     fact = 1.0d0
     RETURN
  END IF
  
  fact = 1
  
  DO i = 1, k
     fact = fact * i
  END DO
  
  RETURN

END FUNCTION fact

