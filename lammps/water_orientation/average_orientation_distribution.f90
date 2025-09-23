MODULE histogram 
  IMPLICIT NONE 
  INTEGER, ALLOCATABLE         :: ntot(:)
  REAL*8, ALLOCATABLE          :: OH_z(:,:), HOH_xy(:), HOH_z(:) 
  INTEGER                      :: zdir
END MODULE histogram

PROGRAM Orientation
  USE parameters, ONLY : nequil, nframes, stride
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_ATOM_REDUCED
    if (MOD(frame,stride)==0) CALL ANALYSIS
  END DO

  CALL PRINT_RESULTS

END PROGRAM Orientation

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, nframes, nequil, pos, &
                         atype, stride, nhist, box, zoffset, &
                         atypeO, atypeH
  USE histogram, ONLY : OH_z,HOH_z, ntot, zdir
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)
  READ(1, *) natoms, nframes, nequil, stride, nhist, zdir
  READ(1, *) zoffset, atypeO, atypeH
  
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(ntot(nhist)); ntot=0
  ALLOCATE(OH_z(2,nhist)); OH_z=0.0
  ALLOCATE(HOH_z(nhist)); HOH_z=0.0

  CLOSE(1)

  OPEN(unit = 1, file =  pos_file)
 
END SUBROUTINE INITIALIZE

SUBROUTINE ANALYSIS
  USE parameters, ONLY : pos, atype, natoms, nhist, box, zoffset, atypeO 
  USE histogram, ONLY : OH_z,HOH_z, ntot, zdir
  IMPLICIT NONE
  INTEGER                    :: i, ind
  INTEGER, DIMENSION(2)      :: ind_H 
  REAL*8                     :: center_(3), z

  DO i = 1,natoms

    !Selecting only Ow
    IF ( atype(i) == atypeO )  THEN

      CALL get_H(i, ind_H)

      IF (ind_H(2).ne.0) THEN

        z = pos(zdir,i) + zoffset
        z = z  - nint(z/box(zdir))*box(zdir) 
        z = z + box(zdir)/2. ! for data analysis only

        ind = int( z * float(nhist) / box(zdir) ) + 1
        ntot(ind) = ntot(ind)+1

        CALL Normal_Vector(zdir,center_)
        CALL Histogram_OH_Angle(i,ind_H,ind,center_,OH_z,HOH_z)

      END IF

    END IF !ind_atom

  END DO !i 

END SUBROUTINE ANALYSIS

SUBROUTINE Histogram_OH_Angle(iOx,indH,rind,reference,OH_angle,HOH_angle)
  USE parameters, only : pos,box,nhist
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iOx, indH(2), rind
  REAL*8, INTENT(IN) :: reference(3)
  REAL*8 :: OH_angle(2,nhist), HOH_angle(nhist)
  REAL*8 :: oh(2), hoh, Angle_reference, Angle_Bissector

  oh(1) = -Angle_reference(pos(:,iOx),pos(:,indH(1)),reference)
  oh(2) = -Angle_reference(pos(:,iOx),pos(:,indH(2)),reference)
  hoh   = -Angle_Bissector(iOx,indH,reference) 

  OH_angle(1,rind) = OH_angle(1,rind) + maxval(oh)
  OH_angle(2,rind) = OH_angle(2,rind) + minval(oh)
  HOH_angle(rind)  = HOH_angle(rind)  + hoh

END SUBROUTINE

SUBROUTINE get_H(ind_O, ind_H)
  USE parameters, ONLY : atype, pos, cutoff_OH, natoms, atypeH
  IMPLICIT NONE
  INTEGER, DIMENSION(2)      :: ind_H
  INTEGER                    :: i, k, ind_O
  REAL*8                     :: Dist

  i=1;k=1

  ind_H = 0
  DO WHILE (( i <= natoms) .and. ( k <= 2))

    !Should be the index for H
    IF ( atype(i) == atypeH ) THEN

      IF ( Dist(ind_O,i) < cutoff_OH ) THEN
  
        ind_H(k) = i
        k = k+1
 
      END IF
 
    END IF

    i = i + 1 

  END DO

  !IF (k.ne.3) THEN
  !  print *, 'Found water ion with charge ', k-3
  !END IF
END SUBROUTINE get_H

SUBROUTINE Normal_Vector(zdir,vector)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: zdir
  REAL*8, DIMENSION(3) :: vector

  if (zdir==1) then
    vector=(/1.,0.,0./)
  elseif (zdir==2) then
    vector=(/0.,1.,0./)
  elseif (zdir==3) then
    vector=(/0.,0.,1./)
  else
    PRINT *, "Wrong direction of surface normal"
    STOP
  end if

END SUBROUTINE Normal_Vector

SUBROUTINE PRINT_RESULTS
  USE parameters, ONLY : nhist, box
  USE histogram, ONLY : OH_z, HOH_z,ntot, zdir
  IMPLICIT NONE
  INTEGER :: ihist
  REAL*8  :: z

  OPEN(unit = 2, file = "Angle_distribution.dat")

  WRITE(2, *) "# Z OH_vector(min) OH_vector(max) HOH_bissector"

  DO ihist=1,nhist
    z = (float(ihist-1))*box(zdir)/float(nhist)-box(zdir)/2.
    if (ntot(ihist).gt.100) then
      WRITE(2,fmt='(F12.8,3(3X,F12.8))') z, &
              OH_z  (1,ihist)/float(ntot(ihist)), &
              OH_z  (2,ihist)/float(ntot(ihist)), &
              HOH_z (ihist)  /float(ntot(ihist))
    else
      WRITE(2,fmt='(F12.8,3(3X,F12.8))') z,&
              0.0,0.0,0.0
    end if
  END DO

  CLOSE(2)

END SUBROUTINE PRINT_RESULTS
