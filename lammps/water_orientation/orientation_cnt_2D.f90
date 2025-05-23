MODULE histogram 
  IMPLICIT NONE 
  INTEGER, ALLOCATABLE         :: ntot(:)
  INTEGER, PARAMETER           :: nhist2=50 !Histogram of angle cos
  REAL*8, ALLOCATABLE          :: OH_xy(:,:), OH_z(:,:), HOH_xy(:,:), HOH_z(:,:) 
  REAL*8                       :: rmax, center(3)
END MODULE histogram

PROGRAM Orientation
  USE parameters, ONLY : nequil, nframes, stride
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_ATOM
    if (MOD(frame,stride)==0) CALL ANALYSIS
  END DO

  CALL PRINT_RESULTS

END PROGRAM Orientation

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, nframes, nequil, pos, &
                         atype, stride, nhist
  USE histogram, ONLY : OH_xy,OH_z,HOH_xy,HOH_z, ntot, &
                        rmax, center, nhist2
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)
  READ(1, *) natoms, nframes, nequil, stride, nhist, rmax
  READ(1, *) center 
  
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(ntot(nhist)); ntot=0
  ALLOCATE(OH_xy(nhist,nhist2)); OH_xy=0.0
  ALLOCATE(OH_z(nhist,nhist2)); OH_z=0.0
  ALLOCATE(HOH_xy(nhist,nhist2)); HOH_xy=0.0
  ALLOCATE(HOH_z(nhist,nhist2)); HOH_z=0.0

  CLOSE(1)

  OPEN(unit = 1, file =  pos_file)
 
END SUBROUTINE INITIALIZE

SUBROUTINE ANALYSIS
  USE parameters, ONLY : pos, atype, natoms, nhist 
  USE histogram, ONLY : OH_xy,OH_z,HOH_xy,HOH_z, ntot, rmax, center
  IMPLICIT NONE
  INTEGER                    :: i, indri
  INTEGER, DIMENSION(2)      :: ind_H 
  REAL*8                     :: Dist_Center, Angle, n(nhist)
  REAL*8                     :: ri, center_(3), Angle_Bissector 

  DO i = 1,natoms

    !Selecting only Ow
    IF ( atype(i) == 2 )  THEN

      ri = Dist_Center(center,pos(:,i))
      indri = int(nhist*ri/rmax)+1
      ntot(indri) = ntot(indri)+1
      CALL get_H(i, ind_H)

      center_=center
      center_(3) = pos(3,i)
      CALL Histogram_OH_Angle(i,ind_H,indri,center_,OH_xy,HOH_xy)

      center_=center
      center_(1:2) = pos(1:2,i)
      CALL Histogram_OH_Angle(i,ind_H,indri,center_,OH_z,HOH_z)

    END IF !ind_atom

  END DO !i 

END SUBROUTINE ANALYSIS

SUBROUTINE Histogram_OH_Angle(iOx,indH,rind,reference,OH_angle,HOH_angle)
  USE parameters, only : pos,box,nhist
  USE histogram, only : nhist2
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iOx, indH(2), rind
  REAL*8, INTENT(IN) :: reference(3)
  REAL*8 :: OH_angle(nhist,nhist2), HOH_angle(nhist,nhist2)
  REAL*8 :: oh(2), hoh, Angle, Angle_Bissector
  INTEGER :: index_angle, ind

  oh(1) = -Angle(pos(:,iOx),reference,pos(:,indH(1)))
  oh(2) = -Angle(pos(:,iOx),reference,pos(:,indH(2)))
  hoh   = Angle_Bissector(iOx,indH,reference) 

  ind=index_angle(oh(1),nhist2)
  OH_angle(rind,ind) = OH_angle(rind,ind) + 1.0
  ind=index_angle(oh(2),nhist2)
  OH_angle(rind,ind) = OH_angle(rind,ind) + 1.0
  ind=index_angle(hoh,nhist2)
  HOH_angle(rind,ind)  = HOH_angle(rind,ind)  + 1.0

END SUBROUTINE

INTEGER FUNCTION index_angle(angle_cos,n)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: angle_cos
  INTEGER, INTENT(IN) :: n

  index_angle = int((angle_cos+1.)/2.*n) + 1

END FUNCTION index_angle

SUBROUTINE get_H(ind_O, ind_H)
  USE parameters, ONLY : atype, pos, cutoff_OH, natoms
  IMPLICIT NONE
  INTEGER, DIMENSION(2)      :: ind_H
  INTEGER                    :: i, k, ind_O
  REAL*8                     :: Dist

  i=1;k=1

  ind_H = 0
  DO WHILE (( i <= natoms) .and. ( k <= 2))

    !Should be the index for H
    IF ( atype(i) == 3 ) THEN

      IF ( Dist(ind_O,i) < cutoff_OH ) THEN
  
        ind_H(k) = i
        k = k+1
 
      END IF
 
    END IF

    i = i + 1 

  END DO

  IF (k.ne.3) THEN
    print *, 'Found water ion with charge ', k-3
  END IF
END SUBROUTINE get_H

SUBROUTINE PRINT_RESULTS
  USE parameters, ONLY : nhist
  USE histogram, ONLY : nhist2, rmax,OH_xy, OH_z, HOH_xy, HOH_z,ntot
  IMPLICIT NONE
  INTEGER :: ihist, ihist2

  OPEN(unit = 2, file = "Angle_distribution_2D.dat")

  WRITE(2, *) "# R Angle OH_xy OH_z HOH_xy HOH_z"

  DO ihist=1,nhist
    if (ntot(ihist).ne.0) then
      DO ihist2=1,nhist2
        WRITE(2,fmt='(6(F12.8,3X))') &
                (float(ihist-1))*rmax/float(nhist), &
                (float(ihist2-1))*2.0/float(nhist2)-1., &
                OH_xy (ihist,ihist2)/float(ntot(ihist)), &
                OH_z  (ihist,ihist2)/float(ntot(ihist)), &
                HOH_xy(ihist,ihist2)  /float(ntot(ihist)), &
                HOH_z (ihist,ihist2)  /float(ntot(ihist))
      END DO
    end if
  END DO

  CLOSE(2)

END SUBROUTINE PRINT_RESULTS
