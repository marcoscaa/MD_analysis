MODULE loc 
  IMPLICIT NONE 
  LOGICAL, ALLOCATABLE               :: is_water_ion(:)
  INTEGER, ALLOCATABLE               :: ind_Ow(:)
  REAL*8, PARAMETER                  :: rcut=6.0d0
END MODULE loc

PROGRAM Print_env 
  USE parameters, ONLY : nframes, stride
  IMPLICIT NONE
  INTEGER :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = 1,nframes
    CALL READ_ATOM_REDUCED
    if (MOD(frame,stride)==0) THEN
      CALL FIND_WATER_ION
      CALL PRINT_ENVIRONMENT 
    end if
  END DO

  CLOSE(1);CLOSE(2)

END PROGRAM Print_env

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms,pos, atype, atypeO, &
                         nframes, nequil, atypeH, &
                         atype, stride
  USE loc,      ONLY : is_water_ion, ind_Ow
  IMPLICIT NONE
  CHARACTER(100)             :: pos_file, index_file

  !User should provide filename and index_file as arguments 
  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nframes, nequil,stride
  READ(1,*) atypeO, atypeH 

  !Allocate module arrays
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(ind_Ow(natoms))
  ALLOCATE(is_water_ion(natoms))

  CLOSE(1)
  OPEN(unit = 1,file = pos_file)
  OPEN(unit = 2,file = "water_ion_env.xyz")

END SUBROUTINE INITIALIZE

SUBROUTINE PRINT_ENVIRONMENT
  USE parameters, ONLY : pos, atype, natoms, atypeO, atypeH 
  USE loc,        ONLY : ind_Ow, rcut, is_water_ion
  IMPLICIT NONE
  INTEGER                    :: i, j, ind_wi, n
  REAL*8                     :: pos_env(3,natoms), Dist
  CHARACTER(1)               :: atype_wi(natoms)

  n=0

  DO i = 1,natoms
    !Selecting only Ow from water ion
    IF ( is_water_ion(i) )  THEN

      ind_wi = i
      CALL SAVE_ATOMS(ind_wi,i,pos_env,atype_wi,n)
      
    END IF !is_water_ion

  END DO !i 

  DO i = 1,natoms

    IF ((atype(i)==atypeO).and.(i.ne.ind_wi)) THEN

      IF (Dist(i,ind_wi)<rcut) CALL SAVE_ATOMS(i,ind_wi,pos_env,atype_wi,n)

    END IF

  END DO

  CALL PRINT_ATOMS(pos_env,atype_wi,n)

END SUBROUTINE PRINT_ENVIRONMENT

SUBROUTINE PRINT_ATOMS(pos_env,atype_wi,n)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  REAL*8 :: pos_env(3,n)
  CHARACTER(1) :: atype_wi(n)
  INTEGER :: i

  WRITE(2,*) n
  WRITE(2,*)

  DO i=1,n
    WRITE(2,fmt='(A1,3(3X,F12.8))') atype_wi(i), pos_env(:,i)
  END DO

END SUBROUTINE PRINT_ATOMS

SUBROUTINE SAVE_ATOMS(atom, center, pos_env, atype_wi, n)
  USE parameters, only : pos, natoms, atypeH, atype
  USE loc, only: ind_Ow
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: atom, center
  REAL*8 :: pos_env(3,natoms)
  CHARACTER(1) :: atype_wi(natoms)
  INTEGER  :: ih, n
  REAL*8 :: d(4), dh(4)

  n=n+1
  CALL DISTANCE_VECTOR_IND(atom,center,d)
  atype_wi(n)="O"
  pos_env(:,n)=d(1:3)
  DO ih=1,natoms
    IF (atype(ih)==atypeH) THEN 
      IF (ind_Ow(ih)==atom) THEN
        n=n+1
        CALL DISTANCE_VECTOR_IND(ih,atom,dh)
        pos_env(:,n)=d(1:3)+dh(1:3)
        atype_wi(n)="H"
      ENDIF
    ENDIF
  END DO

END SUBROUTINE SAVE_ATOMS

SUBROUTINE FIND_WATER_ION
  USE parameters, ONLY : atype, natoms, atypeO, atypeH
  USE loc, ONLY : is_water_ion, ind_Ow
  IMPLICIT NONE
  INTEGER :: iat, cn_ow(natoms)
  INTEGER :: closest_atom

  cn_ow=0
  ind_Ow=0
  is_water_ion=.false.
  DO iat=1,natoms
    if (atype(iat)==atypeH) then !2: H index
      ind_Ow(iat) = closest_atom(iat,atypeO)
      cn_ow(ind_Ow(iat)) = cn_ow(ind_Ow(iat)) + 1
    END IF
  END DO
      
  DO iat=1,natoms
    if (atype(iat)==atypeO) then !1: O index
      if(cn_Ow(iat)/=2) then
        is_water_ion(iat)=.True.
        !DO jat=1,natoms
        !if ((atype(jat)==atypeH).and.(ind_Ow(jat)==iat)) then
        !    is_water_ion(jat)=.True.
        !  end if
        !END DO
      end if
    end if
  END DO

END SUBROUTINE FIND_WATER_ION

INTEGER FUNCTION closest_atom(icenter,atyp)
  USE parameters, ONLY : natoms, atype
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: icenter, atyp
  INTEGER :: iat, iat_min
  REAL*8 :: d, min_d, Dist

  min_d = 100.
  iat_min = 0

  DO iat=1,natoms

    IF (atype(iat)==atyp) THEN
      d = Dist(iat,icenter)
      IF (d < min_d) THEN
        min_d=d
        iat_min=iat
      END IF
    END IF

  END DO

  IF (iat_min==0) THEN
    PRINT *, "Could not find closest atom to atom", icenter
    PRINT *, "STOP!!!!"
    STOP
  END IF

  closest_atom = iat_min

END FUNCTION closest_atom

