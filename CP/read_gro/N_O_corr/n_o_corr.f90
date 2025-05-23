!Set of subroutines used to analyse the lipid trajectory

MODULE parameters
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: pos_o(:,:,:), pos_n(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: n_o_corr(:,:,:)
  DOUBLE PRECISION                         :: box(3)
  INTEGER, ALLOCATABLE                     :: nn(:,:,:)
  INTEGER,PARAMETER                        :: nframes = 1000
  REAL*8,PARAMETER                         :: dt = 2. ! in ps
  INTEGER                                  :: natoms, natO, natN

END MODULE parameters

PROGRAM Hbond
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE

  DO frame = 1,nframes
    CALL READ_FRAME (frame)
    CALL NEAREST_NEIGHBORS (frame)
  END DO

  CALL CORRELATION_FUNCTION
  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM Hbond

SUBROUTINE INITIALIZE
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: file_name, index_file
  CHARACTER(5)               :: moltype,atomtype
  INTEGER                    :: junk1

  ! init
  natO = 0
  natN = 0

  !User should provide filename as argument 
  CALL getarg(1, file_name)

  !Assuming constant number of particles during the simulation
  !Getting the number of atoms only once
  OPEN(unit = 1, file = file_name)
  READ(1,*)
  READ(1,*), natoms

  !get the number of O atoms and N atoms 
  DO i=1,natoms
    READ(1, fmt="(i5,2a5,i5,3f8.3)"), junk1, moltype, atomtype

    if (adjustl(trim(atomtype)).eq.'OH2') THEN 
      natO = natO + 1
    elseif (adjustl(trim(atomtype)).eq.'N') THEN
      natN = natN + 1
    endif

  ENDDO
  READ(1,*), box(1), box(2), box(3)
  CLOSE(1)

  !Allocating all the arrays to be used
  ALLOCATE(pos_o(3,natO,nframes))
  ALLOCATE(pos_n(3,natO,nframes))
  ALLOCATE(nn(nframes,natN,4)); nn = 0
  ALLOCATE(n_o_corr(nframes,natN,3)); n_o_corr = 0

  OPEN(unit = 1,file = file_name)
  
END SUBROUTINE INITIALIZE

SUBROUTINE READ_FRAME (frame)
  !Subroutine reads gro files
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i, junk1, frame, io, intg
  REAL*8                     :: vtmp(3), z_pbc, z
  CHARACTER(5)               :: moltype,atomtype

  io = 0
  intg = 0
  !Discarding first two lines
  READ(1,*)
  READ(1,*)

  !In this case, I need only coordinates of Ow and N atoms
  DO i = 1,natoms

    READ(1, fmt="(i5,2a5,i5,3f8.3)"), junk1, moltype, atomtype, junk1, vtmp(:)

    if ( adjustl(trim(atomtype)).eq.'OH2' ) then
      io = io + 1
      pos_o(:,io,frame) = vtmp(:)

    elseif ( adjustl(trim(atomtype)).eq.'N' ) then
      intg = intg + 1
      pos_n(:,intg,frame) = vtmp(:)

    endif

  END DO

  READ(1,*), box(1), box(2), box(3)

END SUBROUTINE

SUBROUTINE NEAREST_NEIGHBORS (frame)
  !Get the nearest Ow neighbors to lipd N atoms
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i, j, k, frame
  INTEGER                    :: ind_nn(4)
  REAL*8                     :: dist_ij, Dist
  REAL*8                     :: dist_nn(4) 


  DO i = 1, natN

    !Setting nn distance to a huge number
    dist_nn(:) = 1000.
    ind_nn = 0

    DO j = 1, natO

       dist_ij = Dist(pos_o(:,j,frame),pos_n(:,i,frame)) 

       !Cutoff: N = 0.6 nm; P = 0.45 nm (based on RDF)
       !getting the index if dist_ij is lower than the previous ones
       IF ( ( dist_ij < 0.45) .and. ( dist_ij < dist_nn(4)) ) THEN

         dist_nn(4) =  dist_ij 
         ind_nn(4) = j
         CALL SORT_ARRAY( dist_nn, ind_nn )

       ENDIF

    ENDDO        
 
    nn(frame, i, :) = ind_nn

  ENDDO
  

END SUBROUTINE NEAREST_NEIGHBORS 

SUBROUTINE CORRELATION_FUNCTION 
  !Correlate displacement of N with its Ow nearest neighbors
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: frame1, frame2, ind_N, ind_O, i_corr
  INTEGER                    :: counter, counter_tmp,i 
  DOUBLE PRECISION           :: corr_tmp(3), corr(3) 
  LOGICAL                    :: Remain_nn

  pos_n(3, :, :) = pos_n(3, :, :) + 100.
  pos_o(3, :, :) = pos_o(3, :, :) + 100.

  DO ind_N = 1,natN

    !frame1 is the duration of the interval 
    DO frame1 = 1,nframes
   
      counter = 0
      corr = 0.D0
   
      !Change initial condition, look for combinations with same interval
      DO frame2 = 1,1!nframes-frame1
   
        counter_tmp = 0
        corr_tmp = 0.D0

!$omp parallel do reduction(+:counter_tmp,corr_tmp)
        DO ind_O = 1,natO
   
          !If the particle remains nearest neighbor during interval frame1 
          IF ( Remain_nn(ind_N, ind_O, frame2, frame2+frame1 ) ) THEN
   
            counter_tmp = counter_tmp + 1

            DO i_corr = 1,3
              corr_tmp(i_corr) =  corr_tmp(i_corr) + &
                                  & pos_n(i_corr, ind_N, frame2) * &
                                  & pos_o(i_corr, ind_O, frame2+frame1)  
            ENDDO
   
          ENDIF
        ENDDO
!$omp end parallel do

        counter = counter + counter_tmp
        corr = corr + corr_tmp
   
      ENDDO
 
      IF ( counter > 0 ) n_o_corr(frame1,ind_N,:) = corr / dble( counter )
   
    ENDDO
   
  ENDDO

  !Average over all N
  DO i =2, natN
    n_o_corr (:,1,:) = n_o_corr (:,1,:) + n_o_corr(:,i,:)
  END DO
  
  n_o_corr(:,1,1) = n_o_corr(:,1,1) / n_o_corr(1,1,1)
  n_o_corr(:,1,2) = n_o_corr(:,1,2) / n_o_corr(1,1,2)
  n_o_corr(:,1,3) = n_o_corr(:,1,3) / n_o_corr(1,1,3)

END SUBROUTINE CORRELATION_FUNCTION 

SUBROUTINE SORT_ARRAY( dist_nn, ind_nn )
  !Sort a simple array with respect to nn distance
  !One needs only a single buble sort iteration, as only
  !the last index of array is not sorted
  IMPLICIT NONE
  REAL*8, INTENT(INOUT)      :: dist_nn(4)
  REAL*8                     ::  temp_d
  INTEGER, INTENT(INOUT)     :: ind_nn(4)
  INTEGER                    :: temp_i, i

  DO i=4,2,-1

    IF ( dist_nn(i) < dist_nn(i-1) ) THEN
  
      temp_d = dist_nn(i-1)
      dist_nn(i-1) = dist_nn(i)
      dist_nn(i) = temp_d

      temp_i = ind_nn(i-1)
      ind_nn(i-1) = ind_nn(i)
      ind_nn(i) = temp_i

    ENDIF

  ENDDO

END SUBROUTINE SORT_ARRAY

SUBROUTINE PRINT_RESULTS (mode)
  !Print to stdout
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: i, j
  CHARACTER(10),INTENT(IN)   :: mode


  DO i = 1,nframes-1

    write(*,fmt = '(f8.3,3X,3(f15.6,3X))'), float(i)*dt, n_o_corr(i,1,:) 

  ENDDO

END SUBROUTINE PRINT_RESULTS

REAL*8 FUNCTION Dist(a, b)
  ! Distance between two points including pbc
  USE parameters
  IMPLICIT NONE
  REAL*8, DIMENSION(3)       :: a,b,xyz
  INTEGER                    :: i

  DO i = 1,3

    xyz(i) = a(i) - b(i)
    xyz(i) = xyz(i) - nint( xyz(i)/box(i) ) * box(i)

  END DO

  Dist = SQRT( SUM(xyz*xyz) )

END FUNCTION Dist

LOGICAL FUNCTION Remain_nn(ind_N, ind_O, ini_frame, end_frame)
  !Check if a particlular Ow remains nearest neighbor to a N atom
  !during the interval ini_frame : end_frame
  USE parameters
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: ind_N, ind_O, ini_frame, end_frame
  INTEGER                    :: i, j, n

  Remain_nn = .True.

  !This function will only return true if ind_o is a nn during 
  ! all the interval  
  DO i = ini_frame, end_frame

    n = 0

    DO j = 1,4

      IF ( nn(i, ind_N, j)  == ind_O )  n = n + 1
    
    ENDDO

    IF ( n == 0 ) THEN

      Remain_nn = .False.
      EXIT

    ENDIF

  ENDDO

END FUNCTION Remain_nn

REAL*8 FUNCTION Dist2(a, b)
  ! Distance between two points including pbc
  USE parameters
  IMPLICIT NONE
  REAL*8, DIMENSION(3)       :: a,b,xyz
  INTEGER                    :: i

  DO i = 1,3

    xyz(i) = a(i) - b(i)
    xyz(i) = xyz(i) - nint( xyz(i)/box(i) ) * box(i)

  END DO

  Dist2 = SUM(xyz*xyz)

END FUNCTION Dist2

REAL*8 FUNCTION Dist2xy(a, b)
  ! Distance between two points, in xy only, including pbc
  USE parameters
  IMPLICIT NONE
  REAL*8, DIMENSION(3)       :: a,b,xyz
  INTEGER                    :: i

  DO i = 1,2

    xyz(i) = a(i) - b(i)
    xyz(i) = xyz(i) - nint( xyz(i)/box(i) ) * box(i)

  END DO

  Dist2xy = xyz(1)*xyz(1) + xyz(2)*xyz(2)

END FUNCTION Dist2xy

REAL*8 FUNCTION Angle(C, a1, a2)
  !Angle between C-a1 and C-a2
  !Return value in cos(Angle)
  USE parameters
  IMPLICIT NONE
  REAL*8, DIMENSION(3)       :: C, a1, a2
  REAL*8                     :: Dist2, Dist

  Angle = dot_product((C - a1),(C - a2)) / ( norm2(C - a1) * norm2(C - a2) )
  !Angle = Dist2(C,a1) + Dist2(C,a2) - Dist2(a1,a2)
  !Angle = Angle / ( 2 * Dist(C,a1) * Dist(C,a2) )

END FUNCTION Angle

REAL*8 FUNCTION z_pbc(z, box_z)
  !Apply PBC for the z coordinate olny
  IMPLICIT NONE
  REAL*8                     :: z, box_z

  IF ( z > 0 ) THEN
    z_pbc = z - int( z/box_z ) * box_z 
  ELSE
    z_pbc = z - int( -1 + z/box_z ) * box_z
  ENDIF

END FUNCTION z_pbc

DOUBLE PRECISION FUNCTION sum_hist( hist, len )
  !Idiotic sum of an array
  IMPLICIT NONE
  INTEGER                    :: i, len
  DOUBLE PRECISION           :: sum_out
  DOUBLE PRECISION           :: hist(1,len) 

  sum_out = 0.

  DO i=1,len

    sum_out = sum_out + hist(1,i)

  ENDDO

  sum_hist = sum_out
  
END FUNCTION sum_hist
