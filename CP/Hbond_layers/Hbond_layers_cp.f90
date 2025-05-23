!Computes the average and instantaneous Hbond
!Returns these values with respect to layers in the z axis

MODULE histogram
  IMPLICIT NONE 
  INTEGER*8, ALLOCATABLE             :: n_inlayer
  INTEGER, PARAMETER                 :: nhist=100
  DOUBLE PRECISION, ALLOCATABLE      :: hist_hblen(:,:) 
  DOUBLE PRECISION, ALLOCATABLE      :: aver_orient(:,:) 
  DOUBLE PRECISION, ALLOCATABLE      :: aver_hb(:,:) 

END MODULE histogram

PROGRAM Hbond
  USE parameters, ONLY : nframes
  IMPLICIT NONE
  INTEGER :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = 1,nframes
    CALL READ_POS
    CALL MAKE_HISTOGRAM
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM Hbond

SUBROUTINE INITIALIZE
  USE histogram , ONLY : nhist, hist_hblen, aver_orient, aver_hb
  USE parameters, ONLY : natoms, pos, nlayers
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  !User should provide filename and index_file as arguments 
  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  CALL READ_INDEX( index_file )

  ALLOCATE(pos(natoms,3))
  ALLOCATE(hist_hblen(nlayers,nhist)); hist_hblen=0.
  ALLOCATE(aver_orient(nlayers,nhist)); aver_orient=0.
  ALLOCATE(aver_hb(nlayers,2)); aver_hb=0.

  OPEN(unit = 1,file = pos_file)

END SUBROUTINE INITIALIZE

SUBROUTINE MAKE_HISTOGRAM
  USE histogram , ONLY : nhist, hist_hblen, aver_orient, aver_hb
  USE parameters, ONLY : natoms, pos, nlayers
  IMPLICIT NONE
  INTEGER                     :: i, j, layer, layer2
  INTEGER                     :: n_Hb(nlayers), n_in_layer(nlayers)
  INTEGER                     :: layer_index, ind_acceptor
  INTEGER                     :: ind_H(2), ih, indHblen 
  INTEGER                     :: hbtmp(nlayers,nhist)
  INTEGER                     :: nhb(nlayers,2)
  LOGICAL                     :: is_oxygen
  LOGICAL                     :: found_hb_acceptor 
  DOUBLE PRECISION            :: len_Hb, oh(3) 
  DOUBLE PRECISION            :: OHtmp(nlayers,nhist) 

  n_Hb = 0
  hbtmp = 0
  OHtmp = 0
  n_in_layer = 0
  nhb = 0

  DO i = 1,natoms
        
    IF ( is_oxygen(i) ) THEN

      layer = layer_index(pos(i,3))
      n_in_layer(layer) = n_in_layer(layer) + 1 
      CALL get_H(i, ind_H) 
 
      DO ih = 1,2
   
        IF ( ind_H(ih) .ne. 0 ) THEN

          CALL FIND_HB_ACCEPTOR( i, ind_H(ih), len_Hb, ind_acceptor )

          IF ( ind_acceptor.ne.0 ) THEN
 
            CALL OH_VECT( i, ind_H(ih), oh )
            n_Hb(layer) = n_Hb(layer) + 1
            indHblen = int(nhist*len_Hb/3.d0) + 1
            OHtmp(layer,indHblen) = OHtmp(layer,indHblen) + oh(3)
            hbtmp(layer,indHblen) = hbtmp(layer,indHblen) + 1
           
            layer2 = layer_index(pos(ind_acceptor,3))
            IF ( layer .eq. layer2 ) THEN
              nhb(layer,1) = nhb(layer,1) + 2
            ELSE
              nhb(layer,2) = nhb(layer,2) + 1
              nhb(layer2,2) = nhb(layer2,2) + 1
            END IF

          END IF

        END IF

      END DO

    END IF

  END DO

  DO i = 1, nlayers

    IF ( n_Hb(i) .ne. 0 ) THEN

      aver_orient(i,:) = aver_orient(i,:) + OHtmp(i,:)/dble(n_Hb(i))
      hist_hblen(i,:)  = hist_hblen(i,:) + dble(hbtmp(i,:))/dble(n_Hb(i))
      aver_hb(i,:) = aver_hb(i,:) + dble(nhb(i,:))/dble(n_in_layer(i))

    END IF

  END DO

END SUBROUTINE MAKE_HISTOGRAM

SUBROUTINE PRINT_RESULTS
  USE histogram,  ONLY : hist_hblen, nhist, aver_orient, aver_hb
  USE parameters, ONLY : nframes, nlayers
  IMPLICIT NONE
  DOUBLE PRECISION           :: bin
  INTEGER                    :: i

  OPEN(unit = 3,file = "Hist_Hb_length.dat")
  OPEN(unit = 4,file = "Hist_OH_orientation.dat")
  OPEN(unit = 5,file = "Average_Hbonds_per_layer.dat")

  hist_hblen = hist_hblen / dble(nframes)
  aver_orient = aver_orient / dble(nframes)
  aver_hb = aver_hb / dble(nframes)
  bin = 3.d0/dble(nhist)

  DO i = 1, nhist

    WRITE(3, fmt = "(E15.7, *(3X, E15.7))"), bin*(.5d0+dble(i)), hist_hblen(:,i)
    WRITE(4, fmt = "(E15.7, *(3X, E15.7))"), bin*(.5d0+dble(i)), aver_orient(:,i)

  END DO

  WRITE(5, *) "#Layer N_intra N_inter N_total"
  DO i = 1, nlayers
    WRITE(5, fmt = "(I2, 3(3X, F12.8))"), i, aver_hb(i,1), aver_hb(i,2), &
                                             aver_hb(i,1) + aver_hb(i,2)
  END DO

  CLOSE(3)
  CLOSE(4)
  CLOSE(5)

END SUBROUTINE

SUBROUTINE FIND_HB_ACCEPTOR( indO, indH, len_Hb, ind_acc ) 
  !Look for oxygem accepor of Hbond donated by the OH vector 
  !formed by H and O with indexes indH and indO, respectively
  USE parameters, ONLY : pos, natoms
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: indO, indH
  INTEGER                    :: iat, ind_acc
  DOUBLE PRECISION           :: len_Hb, angle_Hb
  DOUBLE PRECISION           :: Angle, Dist, dist_OO
  LOGICAL                    :: is_oxygen

  DO iat = 1, natoms

    IF ( ( is_oxygen(iat) ) .and. ( iat .ne. indO ) ) THEN

      dist_OO = Dist( indO, iat )
      angle_Hb = Angle( indO, indH, iat )

      !Chandler's criteria
      IF ( ( dist_OO .le. 3.5d0 ) .and. ( angle_Hb > 0.8660d0 ) ) THEN

        ind_acc = iat
        len_Hb = Dist( indH, iat )
        RETURN

      END IF

    END IF

  END DO

  ind_acc = 0

END SUBROUTINE FIND_HB_ACCEPTOR
