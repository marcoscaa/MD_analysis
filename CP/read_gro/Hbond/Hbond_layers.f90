!Computes the average and instantaneous Hbond
!Returns these values with respect to layers in the z axis

MODULE parameters
  IMPLICIT NONE
  INTEGER                            :: nframes, natoms, nwater
  INTEGER                            :: nlayers, nequil 
  REAL*8, ALLOCATABLE      :: pos(:,:,:)
  REAL*8                   :: box(3)

END MODULE parameters

MODULE histogram
  IMPLICIT NONE 
  REAL*8, ALLOCATABLE             :: aver_nhb(:)
END MODULE histogram

PROGRAM Hbond
  USE parameters, ONLY : nframes
  IMPLICIT NONE
  INTEGER :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = 1,nframes
    CALL READ_GRO
    CALL MAKE_HISTOGRAM
  END DO

  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM Hbond

SUBROUTINE INITIALIZE
  USE histogram , ONLY : aver_nhb 
  USE parameters, ONLY : nframes, natoms, nwater, pos, nlayers, nequil
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  !User should provide filename and index_file as arguments 
  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nframes, nequil, nlayers
  CLOSE(1)

  nwater = natoms/3 ! Assuming only water is present with no dummy atoms

  ALLOCATE(pos(nwater,3,3))
  ALLOCATE(aver_nhb(nlayers)); aver_nhb=0.

  OPEN(unit = 1,file = pos_file)

END SUBROUTINE INITIALIZE

SUBROUTINE REMOVE_EQUIL
  USE parameters, ONLY : natoms, nequil
  IMPLICIT NONE
  INTEGER                    :: i

  DO i = 1,nequil*(natoms+3)
    READ(1,*)
  END DO

END SUBROUTINE REMOVE_EQUIL

SUBROUTINE READ_GRO (frame)
  !Subroutine reads gro files
  USE parameters, ONLY : nwater, natoms, pos, box
  IMPLICIT NONE
  INTEGER                    :: i, j, k, junk1, frame
  CHARACTER(5)               :: moltype,atomtype

  !Discarding first two lines
  READ(1,*)
  READ(1,*)

  DO i = 1,nwater 
    DO j = 1,3
      READ(1, fmt="(i5,2a5,i5,3f8.3)"), junk1,moltype,atomtype,junk1,(pos(i,j,k),k=1,3)
    END DO
  END DO

  READ(1,*), box(1), box(2), box(3)

END SUBROUTINE READ_GRO

SUBROUTINE MAKE_HISTOGRAM
  USE histogram , ONLY : aver_nhb
  USE parameters, ONLY : nwater, pos, nlayers
  IMPLICIT NONE
  INTEGER                     :: i, j, layer
  INTEGER                     :: n_Hb(nlayers)
  INTEGER                     :: n_inlayer(nlayers)
  INTEGER                     :: layer_index
  LOGICAL                     :: is_hbonded

  n_Hb = 0
  n_inlayer = 0

  DO i = 1,nwater

    layer = layer_index(pos(i,1,3))
    n_inlayer(layer) = n_inlayer(layer) + 1
 
    DO j = 1,nwater
   
      IF (is_hbonded(i,j) .and. (i.ne.j)) THEN
        n_Hb(layer) = n_Hb(layer) + 1 
        !layer = layer_index(pos(j,1,3))
        !n_Hb(layer) = n_Hb(layer) + 1
      END IF

    END DO

  END DO

  DO i = 1, nlayers

    IF ( n_inlayer(i) .ne. 0 ) THEN

      aver_nhb(i) = aver_nhb(i) + float(n_Hb(i))/float(n_inlayer(i))

    END IF

  END DO

END SUBROUTINE MAKE_HISTOGRAM

SUBROUTINE PRINT_RESULTS
  USE histogram,  ONLY : aver_nhb 
  USE parameters, ONLY : nlayers, box, nframes, nequil
  IMPLICIT NONE
  REAL*8           :: bin, z
  INTEGER                    :: i

  OPEN(unit = 2,file = "Aver_Hb_length.dat")

  bin = box(3)/dble(nlayers)
  aver_nhb = aver_nhb/dble(nframes-nequil)

  DO i = 1, nlayers

    z = bin*(.5d0+dble(i))
    !z = z - nint(z/box(3))*box(3) + box(3)/2. !Not necessary
    WRITE(2, fmt = "(E15.7, *(3X, E15.7))"), z, aver_nhb(i)

  END DO

  CLOSE(2)

END SUBROUTINE

REAL*8 FUNCTION Dist( X, Y )
  USE parameters, ONLY : box
  IMPLICIT NONE
  REAL*8, INTENT(in) :: X(3), Y(3)
  REAL*8 :: d(3) 
  INTEGER :: ipol

  DO ipol = 1,3
  
    d(ipol) = X(ipol) - Y(ipol)
    d(ipol) = d(ipol) - nint(d(ipol)/box(ipol))*box(ipol)

  END DO

  Dist = norm2(d)

END FUNCTION Dist

REAL*8 FUNCTION Angle( X, Y, Z )
  !X-Y-Z angle
  USE parameters, ONLY : box
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: X(3), Y(3), Z(3)
  REAL*8 :: d1(3), d2(3) 
  INTEGER :: ipol

  DO ipol = 1,3
  
    d1(ipol) = Y(ipol) - X(ipol)
    d1(ipol) = d1(ipol) - nint(d1(ipol)/box(ipol))*box(ipol)

    d2(ipol) = Z(ipol) - X(ipol)
    d2(ipol) = d2(ipol) - nint(d2(ipol)/box(ipol))*box(ipol)

  END DO

  Angle = sum(d1*d2)/(norm2(d1)*norm2(d2))

END FUNCTION Angle

LOGICAL FUNCTION is_hbonded( ind1, ind2 )
  !Determine if oxygen ind1 H-bonds to ind2 based on 
  !Chandler's definition of H-bond
  USE parameters, ONLY : pos
  IMPLICIT NONE 
  INTEGER, INTENT(IN)                  :: ind1, ind2
  INTEGER                              :: ih, iO 
  REAL*8                     :: Angle, OHO_ang, Dist
  REAL*8, PARAMETER          :: cutoff_OO = 0.35
  REAL*8, PARAMETER          :: cutoff_Hb = 0.23 !0.182 
  
  IF ( Dist(pos(ind1,1,:),pos(ind2,1,:)) < cutoff_OO ) THEN

    DO ih = 2,3

      !H-bond donation by ind1
      !OHO_ang = Angle( pos(ind1,1,:), pos(ind1,ih,:), pos(ind2,1,:) )
      !IF ( OHO_ang > 0.8660d0 ) THEN
      IF ( Dist(pos(ind1,ih,:),pos(ind2,1,:)) < cutoff_Hb ) THEN
        is_hbonded = .true.
        RETURN
      END IF

      !H-bond accepted by ind1
      !OHO_ang = Angle( pos(ind2,1,:), pos(ind2,ih,:), pos(ind1,1,:) )
      !IF ( OHO_ang > 0.8660d0 ) THEN
      IF ( Dist(pos(ind2,ih,:),pos(ind1,1,:)) < cutoff_Hb ) THEN
        is_hbonded = .true.
        RETURN
      END IF

    END DO

  END IF

  is_hbonded = .false. 

END FUNCTION is_hbonded

INTEGER FUNCTION layer_index(pos_z)
  USE parameters, ONLY : box, nlayers
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: pos_z
  REAL*8 :: z
  
  z = pos_z - nint(pos_z/box(3))*box(3) + box(3)/2.
  layer_index = int(nlayers*z/box(3)) + 1

END FUNCTION layer_index
