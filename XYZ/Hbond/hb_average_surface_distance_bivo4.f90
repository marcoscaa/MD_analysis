MODULE histogram2 
  IMPLICIT NONE 
  INTEGER, ALLOCATABLE         :: NC_coordnumb(:), H_id(:)
  INTEGER, ALLOCATABLE         :: hb_counter(:,:), n_in_layer(:)
  REAL*8, ALLOCATABLE          :: hb_length(:,:) 
  CHARACTER(5)                 :: ind_hb(2)
  LOGICAL, ALLOCATABLE         :: is_OW(:)
END MODULE histogram2

PROGRAM ZDF
  USE parameters, ONLY : nequil, nframes, stride
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_EXTXYZ
    IF (frame==1) CALL IDENTIFY_OW
    if (MOD(frame,stride)==0) THEN
      CALL ANALYSIS
    end if
  END DO

  CALL PRINT_RESULTS

END PROGRAM ZDF

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, nframes, nequil, pos, &
                         atype, stride, nhist, zoffset 
  USE histogram2, ONLY : ind_hb, hb_counter, hb_length, n_in_layer
  IMPLICIT NONE
  INTEGER                    :: i,j
  CHARACTER(100)             :: pos_file, index_file

  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)
  READ(1, *) natoms, nframes, nequil, stride, nhist 
  READ(1, *) zoffset
  READ(*,*) ind_hb
  
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(hb_counter(2,nhist)); hb_counter = 0.0
  ALLOCATE(hb_length(2,nhist)); hb_length = 0.d0
  ALLOCATE(n_in_layer(nhist)); n_in_layer = 0

  CLOSE(1)

  OPEN(unit = 1, file =  pos_file)
 
END SUBROUTINE INITIALIZE

SUBROUTINE ANALYSIS
  USE parameters, ONLY : pos, atype, natoms, nhist 
  USE histogram2, ONLY : ind_hb, hb_counter, hb_length, n_in_layer
  IMPLICIT NONE
  INTEGER                    :: i, j, k
  INTEGER                    :: ind_H(4)
  REAL*8                     :: d_OOw, Dist, Angle, ang_OHO

  DO i = 1,natoms

    IF ( trim(atype(i)) == trim(ind_hb(1)) ) THEN

      ind_H = 0
      CALL get_H(i, ind_H)
      CALL UPDATE_COUNTER(i,n_in_layer)

      DO j = 1,natoms

        IF ( ( trim(atype(j)) == trim(ind_hb(2)) ) .and. (i.ne.j) ) THEN

          d_OOw = Dist(i, j)

          DO k = 1,4

            IF (ind_H(k) /= 0) THEN

              ang_OHO = Angle(pos(:,i), pos(:,j), pos(:,ind_H(k)))

              !Definition of H-bond - Luzar and Chandler
              IF ( (d_OOw < 3.5) .and. (ang_OHO > 0.8660) ) THEN

                CALL UPDATE_HBLEN_HISTOGRAM(i,j,ind_H(k),hb_length)
                CALL UPDATE_HBNUMBER_COUNTER(i,j,hb_counter)

              END IF !d_OOw

            END IF ! ind_H

          END DO ! k

        END IF ! atype
 
      END DO !j

    END IF !ind_atom

  END DO !i 

END SUBROUTINE ANALYSIS

SUBROUTINE UPDATE_COUNTER(i,counter)
  USE parameters, only : nhist
  USE histogram2, only : is_OW
  IMPLICIT NONE 
  INTEGER, INTENT(IN) :: i
  INTEGER             :: indz, get_z_index, counter(nhist)

  IF (is_OW(i)) THEN
    indz = get_z_index(i)
    counter(indz) = counter(indz) + 1
  END IF

END SUBROUTINE UPDATE_COUNTER 

INTEGER FUNCTION get_z_index(iatom)
  USE parameters, only : nhist, box, pos, zoffset
  INTEGER, INTENT(IN) :: iatom 
  INTEGER             :: indz 
  REAL*8              :: posz,boxz

  boxz=box(3,3)
  posz = pos(3,iatom) + zoffset
  posz = posz - nint(posz/boxz)*boxz
  posz = posz + boxz/2.
  indz = int(posz/boxz*float(nhist)) + 1 

  get_z_index = indz

END FUNCTION get_z_index

SUBROUTINE UPDATE_HBNUMBER_COUNTER(i,j,counter)
  !i: Hbond donor, j: Hbond acceptor
  USE parameters, ONLY : nhist
  USE histogram2, only : is_OW
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i,j
  INTEGER             :: indz, get_z_index, counter(2,nhist)

  IF (is_OW(i)) THEN
    indz = get_z_index(i)
    counter(1,indz) = counter(1,indz)+1
  END IF
  IF (is_OW(j)) THEN
    indz = get_z_index(j)
    counter(2,indz) = counter(2,indz)+1
  END IF

END SUBROUTINE UPDATE_HBNUMBER_COUNTER

SUBROUTINE get_H(ind_O, ind_H)
  USE parameters, ONLY : atype, pos, cutoff_OH, natoms
  IMPLICIT NONE
  INTEGER, DIMENSION(4)      :: ind_H
  INTEGER                    :: i, k, ind_O
  REAL*8                     :: Dist

  i=1;k=1

  DO WHILE (( i <= natoms) .and. ( k <= 4))

    !Should be the index for H
    IF ( trim(atype(i)) == "H" ) THEN

      IF ( Dist(ind_O,i) < cutoff_OH ) THEN
  
        ind_H(k) = i
        k = k+1
 
      END IF
 
    END IF

    i = i + 1 

  END DO

END SUBROUTINE get_H

SUBROUTINE UPDATE_HBLEN_HISTOGRAM(Donor, Acceptor, H, hb_length)
  USE parameters, only : nhist
  USE histogram2, only : is_OW
  IMPLICIT NONE
  INTEGER, INTENT(IN)         :: Donor, Acceptor, H
  REAL*8                      :: Dist, d1, hb_length(2,nhist)
  INTEGER                     :: ind, get_z_index

  d1 = Dist(Acceptor, H)
  IF (is_OW(Donor)) THEN
    ind = get_z_index(Donor)
    hb_length(1,ind) = hb_length(1,ind) + d1
  END IF
  IF (is_OW(Acceptor)) THEN
    ind = get_z_index(Acceptor)
    hb_length(2,ind) = hb_length(2,ind) + d1
  END IF

END SUBROUTINE UPDATE_HBLEN_HISTOGRAM

SUBROUTINE IDENTIFY_OW()
  USE parameters, ONLY : atype, pos, cutoff_NH, natoms
  USE histogram2, ONLY : is_OW
  IMPLICIT NONE
  INTEGER                    :: iat, iV
  REAL*8                     :: Dist

  IF (.not. allocated(is_OW)) THEN
      allocate (is_OW(natoms))
  END IF

  is_OW=.false.

  DO iat = 1,natoms
    IF (trim(atype(iat))=="O") THEN
        is_OW(iat)=.true.
        DO iV = 1,natoms
           IF (trim(atype(iV))=="V") THEN
               IF (Dist(iat,iV)<2.6) THEN
                   is_OW(iat)=.false.
                   EXIT
               END IF
           END IF
        END DO
    END IF
  END DO

END SUBROUTINE IDENTIFY_OW

SUBROUTINE PRINT_RESULTS
  USE histogram2, ONLY : hb_counter, hb_length, n_in_layer
  USE parameters, ONLY : nhist, box
  IMPLICIT NONE
  INTEGER :: ihist, j
  REAL*8 :: dZ

  OPEN(unit = 2, file = "HB_Number.dat")
  OPEN(unit = 3, file = "HB_Length.dat")

  WRITE(2, *) "# Z(A) Mean_Hbond(Donated) Mean_Hbond(Accepted)"
  WRITE(3, *) "# Z(A) HBond_Length_Donated(A) HBond_Length_Accepted(A)"

  dZ = box(3,3)/float(nhist)

  !Avoid NaNs
  DO ihist=1,nhist
    IF (n_in_layer(ihist)==0) n_in_layer(ihist)=1
    WRITE(2,fmt='(F12.8,2(3X,F12.8))') dz*float(ihist-1)-box(3,3)/2.d0, &
            float(hb_counter(1,ihist))/float(n_in_layer(ihist)), &
            float(hb_counter(2,ihist))/float(n_in_layer(ihist))
    DO j=1,2
      IF (hb_counter(j,ihist)==0) hb_counter(j,ihist)=1
    END DO
  END DO

  DO ihist=1,nhist
    WRITE(3,fmt='(F12.8,2(3X,F12.8))') dz*float(ihist-1)-box(3,3)/2.d0, &
            hb_length(1,ihist)/float(hb_counter(1,ihist)), &
            hb_length(2,ihist)/float(hb_counter(2,ihist))
  END DO

  CLOSE(2);CLOSE(3)

END SUBROUTINE PRINT_RESULTS
