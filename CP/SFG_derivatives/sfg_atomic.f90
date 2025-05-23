!Set of subroutines used to analyse the lipid trajectory

MODULE histogram
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: oh_dist(:,:,:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: firstder(:,:) 
  DOUBLE PRECISION, ALLOCATABLE            :: secder(:,:) 
  DOUBLE PRECISION, ALLOCATABLE            :: aver_1(:,:,:), aver_2(:,:,:) 
END MODULE histogram

PROGRAM SFG
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE

  DO frame = 1,nframes
    CALL READ_POS
    CALL READ_NUMERICAL_DERIVATIVES
    CALL CHECK_SUM_RULE
    CALL COARSE_GRAIN_POS (frame)
    CALL DERIVATIVES_DISTRIBUTION (frame)
  END DO

  CALL PRINT_RESULTS

  CLOSE(1);CLOSE(2);CLOSE(3)
  CLOSE(51);CLOSE(52)

END PROGRAM SFG

SUBROUTINE INITIALIZE
  USE histogram,  ONLY : oh_dist, firstder, secder, &
                         aver_1, aver_2
  USE parameters, ONLY : natoms, nframes, nwater, &
                         pos, coarse, nlayers
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file
  CHARACTER(100)             :: first_der, second_der 

  CALL getarg(1, pos_file)
  CALL getarg(2, first_der)
  CALL getarg(3, second_der)
  CALL getarg(4, index_file)

  CALL READ_INDEX ( index_file )

  ALLOCATE(pos(natoms,3)) ! coordinates
  ALLOCATE(coarse(nframes,nwater)); coarse=1
  ALLOCATE(firstder(natoms,3)); firstder = 0.d0
  ALLOCATE(secder(natoms,3)); secder = 0.d0
  ALLOCATE(aver_1(3,2,nlayers)); aver_1 = 0.d0
  ALLOCATE(aver_2(3,2,nlayers)); aver_2 = 0.d0

  OPEN(unit = 1,file = pos_file)
  OPEN(unit = 2,file = first_der)
  OPEN(unit = 3,file = second_der)
  OPEN(unit = 51,file= "dump_derivatives")
  OPEN(unit = 52,file= "sum_rule")
  
END SUBROUTINE INITIALIZE

SUBROUTINE READ_NUMERICAL_DERIVATIVES
  USE histogram, ONLY : firstder, secder
  USE parameters,ONLY : natoms
  IMPLICIT NONE
  INTEGER                    :: iat

  READ(2,*)
  READ(3,*)

  DO iat = 1, natoms
    READ(2, *), firstder(iat,1),firstder(iat,2),firstder(iat,3)
    READ(3, *), secder(iat,1),secder(iat,2),secder(iat,3)
  END DO

END SUBROUTINE

SUBROUTINE CHECK_SUM_RULE
  USE histogram, ONLY : firstder, secder
  USE parameters,ONLY : natoms
  IMPLICIT NONE
  INTEGER                    :: iat
  DOUBLE PRECISION           :: sr1(3), sr2(3)

  sr1 = 0
  sr2 = 0

  DO iat = 1, natoms
  
    sr1 = sr1 + firstder(iat,:)
    sr2 = sr2 + secder(iat,:)

  END DO

  WRITE(52,fmt="(6(3X,E12.5))") sr1, sr2

END SUBROUTINE CHECK_SUM_RULE

SUBROUTINE DERIVATIVES_DISTRIBUTION (frame)
  !Generate the OH distribution for each water in a specific frame
  USE histogram,  ONLY : firstder, secder, aver_1, aver_2
  USE parameters, ONLY : natoms, coarse, nlayers, pos
  IMPLICIT NONE
  INTEGER                    :: frame, i, ih
  INTEGER                    :: layer, ind_O, ind_H(2)
  INTEGER                    :: counter(nlayers)
  LOGICAL                    :: is_water_oxygen
  DOUBLE PRECISION           :: tmp1(3,2,nlayers), tmp2(3,2,nlayers), oh(3)

  !Water oxygen atom index
  ind_O = 0
  counter = 0 
  tmp1=0.d0;tmp2=0.d0

  DO i = 1,natoms

    !Assigning water molecules to their respective layers
    IF ( is_water_oxygen(i) ) THEN

      ind_O = ind_O + 1
      layer = coarse(frame,ind_O)
      tmp1(:,1,layer) = tmp1(:,1,layer) + firstder(i,:) 
      tmp2(:,1,layer) = tmp2(:,1,layer) + secder(i,:) 
      CALL get_H(i, ind_H)
      counter(layer) = counter(layer) + 1
   
      !Two H per O in water. OW defines the layer for HW 
      DO ih = 1,2

        CALL OH_VECT( i, ind_H(ih), oh )

        tmp1(:,2,layer) = tmp1(:,2,layer) + firstder(ind_H(ih),:)/2.d0 
        tmp2(:,2,layer) = tmp2(:,2,layer) + secder(ind_H(ih),:)/2.d0 

        WRITE(51,fmt="(2(i2,3X),7(3X,F10.5))"), layer, 2, &
             firstder(ind_H(ih),:), secder(ind_H(ih),:), oh(3) 
         
      END DO

      WRITE(51,fmt="(2(i2,3X),7(3X,F10.5))"), layer, 1, &
           firstder(i,:), secder(i,:), pos(i,3)

    END IF

  ENDDO

  DO i = 1,nlayers
    aver_1(:,:,i) = aver_1(:,:,i) + tmp1(:,:,i)/dble(counter(i))
    aver_2(:,:,i) = aver_2(:,:,i) + tmp2(:,:,i)/dble(counter(i))
  END DO

END SUBROUTINE DERIVATIVES_DISTRIBUTION

SUBROUTINE PRINT_RESULTS 
  USE histogram, ONLY : aver_1, aver_2
  USE parameters,ONLY : nlayers, nframes
  IMPLICIT NONE
  INTEGER                    :: il, itype
  DOUBLE PRECISION           :: out1(3), out2(3)

  OPEN(unit = 50,file = "derivatives_atomic")

  DO il = 1, nlayers
  
    DO itype = 1,2

      out1 = aver_1(:,itype,il) / dble(nframes)
      out2 = aver_2(:,itype,il) / dble(nframes)

      WRITE(50,fmt="(2(i2,2X),6(3X,F10.5))") il, itype, out1, out2 

    END DO

  END DO

  CLOSE(50)

END SUBROUTINE PRINT_RESULTS

SUBROUTINE OH_VECT(O, H, OH)
  ! Normalized vector from atom H to atom O
  USE parameters, ONLY : box, pos
  IMPLICIT NONE
  INTEGER, INTENT(IN)                  :: O, H
  DOUBLE PRECISION, INTENT(INOUT)      :: OH(3)
  INTEGER                              :: ipol

  DO ipol=1,3

    OH(ipol) = pos(H,ipol) - pos(O,ipol)
    OH(ipol) = OH(ipol) - nint( OH(ipol)/box(ipol,ipol) ) * box(ipol,ipol)

  END DO

  OH = OH / norm2(OH)

END SUBROUTINE oh_vect
