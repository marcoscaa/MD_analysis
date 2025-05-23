MODULE histogram 
  IMPLICIT NONE 
  INTEGER, ALLOCATABLE         :: ntot(:)
  REAL*8, ALLOCATABLE          :: OH_z(:,:), HOH_z(:,:) 
  INTEGER                      :: zdir
END MODULE histogram

PROGRAM Orientation
  USE parameters, ONLY : nequil, nframes, stride
  IMPLICIT NONE
  INTEGER                              :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_EXTXYZ
    if (MOD(frame,stride)==0) CALL ANALYSIS
  END DO

  CALL PRINT_RESULTS

END PROGRAM Orientation

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, nframes, nequil, pos, box,&
                         atype, stride, nhist, nlayers, layers, zoffset
  USE histogram, ONLY : OH_z,HOH_z, ntot, zdir
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)
  READ(1, *) natoms, nframes, nequil, stride, nhist, zdir, nlayers
  READ(1, *) zoffset
  ALLOCATE(layers(nlayers+1))
  READ(1, *) layers
  
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(ntot(nlayers)); ntot=0
  ALLOCATE(OH_z(nlayers,nhist)); OH_z=0.0
  ALLOCATE(HOH_z(nlayers,nhist)); HOH_z=0.0

  CLOSE(1)

  OPEN(unit = 1, file =  pos_file)
 
END SUBROUTINE INITIALIZE

SUBROUTINE ANALYSIS
  USE parameters, ONLY : pos, atype, natoms, nhist, box, zoffset 
  USE histogram, ONLY : OH_z,HOH_z, ntot, zdir
  IMPLICIT NONE
  INTEGER                    :: i, ind, layer, layer_index
  INTEGER, DIMENSION(2)      :: ind_H 
  REAL*8                     :: center_(3), z

  DO i = 1,natoms

    !Selecting only Ow
    IF ( TRIM(atype(i)) == "O" )  THEN

      CALL get_H(i, ind_H)

      IF (ind_H(2).ne.0) THEN

        z = pos(zdir,i) + zoffset
        layer=layer_index(z,zdir)
        ntot(layer) = ntot(layer)+1

        CALL Normal_Vector(zdir,center_)
        CALL Histogram_OH_Angle(i,ind_H,layer,center_,OH_z,HOH_z)

      END IF

    END IF !ind_atom

  END DO !i 

END SUBROUTINE ANALYSIS

SUBROUTINE Histogram_OH_Angle(iOx,indH,layer,reference,OH_angle,HOH_angle)
  USE parameters, only : pos,box,nhist,nlayers
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iOx, indH(2), layer
  REAL*8, INTENT(IN) :: reference(3)
  INTEGER :: angle_index, iangle
  REAL*8 :: OH_angle(nlayers,nhist), HOH_angle(nlayers,nhist)
  REAL*8 :: oh(2), hoh, Angle_reference, Angle_Bissector

  oh(1) = -Angle_reference(pos(:,iOx),pos(:,indH(1)),reference)
  iangle=angle_index(oh(1),nhist)
  OH_angle(layer,iangle) = OH_angle(layer,iangle) + 1

  oh(2) = -Angle_reference(pos(:,iOx),pos(:,indH(2)),reference)
  iangle=angle_index(oh(2),nhist)
  OH_angle(layer,iangle) = OH_angle(layer,iangle) + 1

  hoh   = -Angle_Bissector(iOx,indH,reference) 
  iangle=angle_index(hoh,nhist)
  HOH_angle(layer,iangle)  = HOH_angle(layer,iangle) + 1

END SUBROUTINE Histogram_OH_Angle

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
    IF ( TRIM(atype(i)) == "H" ) THEN

      IF ( Dist(ind_O,i) < cutoff_OH ) THEN
  
        ind_H(k) = i
        k = k+1
 
      END IF
 
    END IF

    i = i + 1 

  END DO

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
    PRINT *, "ERROR: Non-supported direction of surface normal"
    STOP
  end if

END SUBROUTINE Normal_Vector

INTEGER FUNCTION angle_index(angle,nhist)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nhist
  REAL*8, INTENT(IN) :: angle

  angle_index = int(nhist*(angle+1.)/2.) +1

END FUNCTION angle_index

SUBROUTINE PRINT_RESULTS
  USE parameters, ONLY : nhist, box
  USE histogram, ONLY : OH_z, HOH_z,ntot, zdir
  IMPLICIT NONE
  INTEGER :: ihist

  OPEN(unit = 2, file = "OH_layer_resolved_distribution.dat")
  OPEN(unit = 3, file = "HOH_layer_resolved_distribution.dat")

  WRITE(2, *) "# Angle_cosine P(angle_cos)"
  WRITE(3, *) "# Angle_cosine P(angle_cos)"

  DO ihist=1,nhist
    if (ntot(ihist).gt.100) then
      WRITE(2,fmt='(F12.8,*(3X,F12.8))') 2*(float(ihist-1))/float(nhist)-1.0, &
              OH_z(:,ihist)/float(ntot)
      WRITE(3,fmt='(F12.8,*(3X,F12.8))') 2*(float(ihist-1))/float(nhist)-1.0, &
              HOH_z(:,ihist)/float(ntot)
    end if
  END DO

  CLOSE(2)
  CLOSE(3)

END SUBROUTINE PRINT_RESULTS
