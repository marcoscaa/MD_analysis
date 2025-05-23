!Distribution of an order parameter

MODULE steinhardt 
  IMPLICIT NONE 
  INTEGER, PARAMETER                       :: nbins = 200, number_NN=4, l=6
  INTEGER, PARAMETER                       :: excl_nn=0
  REAL*8, ALLOCATABLE                      :: dist_matrix(:,:,:), hist(:)
  INTEGER, ALLOCATABLE                     :: ind_nn(:,:)

END MODULE steinhardt 

PROGRAM ORDER
  USE parameters, ONLY : nframes
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = 1,nframes
    CALL READ_ATOM_REDUCED
    CALL COMPUTE_DISTANCE_MATRIX
    CALL SORT_NEAREST_NEIGHBORS
    CALL COMPUTE_ORDER_PARAMETER
    !CALL WRITE_XYZ
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM ORDER

SUBROUTINE INITIALIZE
  USE steinhardt,  ONLY : dist_matrix, ind_nn, hist, &
                         number_NN, nbins
  USE parameters, ONLY : natoms, pos, nlayers, nequil, &
                         nframes, atype
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_name, index_file

  CALL getarg(1, file_name)
  CALL getarg(2, index_file)

  !Read index file 
  OPEN(unit = 1,file = index_file)
  READ(1, *) natoms, nframes, nequil
  CLOSE(1)

  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(dist_matrix(4,natoms,natoms)); dist_matrix=0.d0
  ALLOCATE(ind_nn(natoms,number_NN)); ind_nn = 0
  ALLOCATE(hist(nbins)); hist = 0 
  
  OPEN(unit  =  1,file  =  file_name)
  OPEN(unit  =  2,file  =  'time_dependent_histogram.dat') 
 
END SUBROUTINE INITIALIZE

SUBROUTINE PRINT_RESULTS
  USE steinhardt,  ONLY : nbins, hist
  IMPLICIT NONE
  INTEGER                    :: j
  REAL*8                     :: bin, z
  
  OPEN(unit = 3,file = "order_histogram.dat")

  bin = 1.d0 / dble(nbins)
  hist = hist / sum(hist)

  DO j = 1,nbins

    z = bin/2. + dble(j)*bin
    WRITE(3,fmt = "(E12.5, 3X, E15.7)") z, hist(j) 
  END DO

  CLOSE(3)

END SUBROUTINE

SUBROUTINE COMPUTE_ORDER_PARAMETER
  !Compute the average and standard deviation of an order parameter for a given
  !frame 
  USE steinhardt,  ONLY : hist, nbins 
  USE parameters, ONLY : pos, natoms, box
  IMPLICIT NONE
  INTEGER                    :: iat, ind_hist, hist_local(nbins)
  INTEGER                    :: get_index_hist, layer_index_equal
  REAL*8                     :: Q_local, op

  hist_local = 0

  DO iat = 1, natoms
  
    op = Q_local(iat)
    !Qatom(iat) = op
    ind_hist = get_index_hist(op,nbins) 
    hist_local(ind_hist) = hist_local(ind_hist) + 1
    hist(ind_hist) = hist(ind_hist) + 1.d0

  END DO
   
  CALL PRINT_HISTOGRAM(hist_local,nbins)

END SUBROUTINE COMPUTE_ORDER_PARAMETER

SUBROUTINE PRINT_HISTOGRAM(hist,nbins)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nbins, hist(nbins)
  INTEGER             :: i, sum_hist
  REAL*8              :: bin, prob

  sum_hist = sum(hist)
  WRITE(2,*) '# OP Probability'

  DO i = 1, nbins
    prob = float(hist(i)) / float(sum_hist)
    bin = (float(i)-0.5)/float(nbins)
    WRITE(2,fmt="(F12.8,3X,F12.8)") bin, prob
  END DO

END SUBROUTINE PRINT_HISTOGRAM

INTEGER FUNCTION get_index_hist(op,nbins)
  !Histogram of order parameter values. Assuming it is restricted to [0,1]
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: op 
  INTEGER, INTENT(IN)          :: nbins

  get_index_hist =  int( dble(nbins) * op ) + 1

END FUNCTION get_index_hist

SUBROUTINE COMPUTE_DISTANCE_MATRIX
  !Builds a symmetrical distance matrix
  USE steinhardt,  ONLY : dist_matrix
  USE parameters, ONLY : natoms, pos 
  IMPLICIT NONE
  INTEGER                              :: iat1, iat2
  
  DO iat1 = 1, natoms
  
    DO iat2 = iat1+1, natoms
 
      CALL DISTANCE_VECTOR( pos(:,iat1), pos(:,iat2), dist_matrix(:,iat1,iat2) )
      dist_matrix(1:3,iat2,iat1) = -dist_matrix(1:3,iat1,iat2)
      dist_matrix(4,iat2,iat1) = dist_matrix(4,iat1,iat2)
 
    END DO

  END DO

END SUBROUTINE COMPUTE_DISTANCE_MATRIX

SUBROUTINE SORT_NEAREST_NEIGHBORS
  !Returns the indexes of the `number_NN nearest` neighbors for all atoms
  USE steinhardt , ONLY : dist_matrix, number_NN, ind_nn
  USE parameters, ONLY : natoms, ind_atom
  IMPLICIT NONE
  INTEGER                    :: iat1, iat2, ind_tmp(number_NN)
  REAL*8                     :: dist_nn(number_NN)
 
  DO iat1 = 1, natoms

    dist_nn = 6.d0 !should be the second rdf minima
    ind_tmp = 0 

    DO iat2 = 1, natoms 
  
      IF (iat2==iat1) CYCLE
      IF ( dist_matrix(4, iat1, iat2) < dist_nn(number_NN) ) THEN
      
        dist_nn(number_NN) = dist_matrix(4, iat1, iat2)
        ind_tmp(number_NN) = iat2
        CALL BUBBLE_SORT(dist_nn, ind_tmp, number_NN, number_NN) 

      END IF

    END DO

    ind_nn(iat1,:) = ind_tmp

  END DO

END SUBROUTINE SORT_NEAREST_NEIGHBORS

! ==================================================================================
!  Importing from Fausto's code
! ==================================================================================

REAL*8 FUNCTION Q_local(i)

  USE parameters, ONLY : natoms, pi
  USE steinhardt,  ONLY : l
  IMPLICIT NONE

  INTEGER,INTENT(IN) :: i
  INTEGER :: m

  REAL*8 :: rey(-l:l), imy(-l:l)
  REAL*8 :: sum_loc
  REAL*8 :: Q_tmp

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
  
  REAL*8 :: phi, cosmphi, sinmphi
  REAL*8 :: rey(-l:l), imy(-l:l)
  REAL*8 :: re, im
  REAL*8 :: rij(3) 
  REAL*8 :: fact, legendre

  EXTERNAL fact, legendre

  DO m = -l, l
     rey(m) = 0
     imy(m) = 0
  END DO
  
  DO k = excl_nn+1, number_NN
     j = ind_NN(i, k)

     DO ipol = 1,3
       rij(ipol) = dist_matrix(ipol,j,i) / dist_matrix(4,j,i)
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
  
  REAL*8 :: x, legendre
  REAL*8 :: suml
  REAL*8 :: fact
  
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
  REAL*8 :: fact
  
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

