
PROGRAM ProtonTransfer
  USE parameters, ONLY : nequil, nframes
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL

  DO frame = nequil+1,nframes
    CALL READ_ATOM_REDUCED
    IF (frame==nequil+1) CALL FIND_O2c
    CALL ANALYSIS(frame)
  END DO

  CLOSE(1)

END PROGRAM ProtonTransfer

SUBROUTINE INITIALIZE
  USE parameters,    ONLY : natoms, nframes, nequil, &
                            pos, atype
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file

  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)

  !Read number of atoms and box size
  OPEN(unit = 1,file = index_file)
  READ(1, *), natoms, nframes, nequil
  CLOSE(1)
  
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))

  OPEN(unit  =  1,file  =  pos_file)
  OPEN(unit  =  2,file  = "proton_wire_tio2.dat")
  WRITE(2,fmt='(A21)'), "# step path_O path_H"
 
END SUBROUTINE INITIALIZE

SUBROUTINE ANALYSIS(step)
  USE parameters, ONLY : nTi5c, index_O2c
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: step
  INTEGER                    :: iat, ttype, index_OwTi(nTi5c)
  REAL*8                     :: path_O, path_H

  WRITE(2,fmt='(A4,3X,I9)'), "####", step

  CALL Find_OwTi(index_OwTi)
  
  DO iat=1,nTi5c 

    IF ( index_OwTi(iat).ne.0 ) THEN 
 
      CALL proton_wire_cvs(index_OwTi(iat),path_O,path_H)

      !Only print if found a connected path
      IF (path_O<100.d0) THEN
        WRITE(2,fmt="(i9,2(3X,E17.5))"), & 
          step, path_O, path_H
      END IF

    END IF

  END DO

END SUBROUTINE ANALYSIS

SUBROUTINE proton_wire_cvs(index_OwTi,path_O,path_H)
  USE parameters, ONLY : cutoff_OH, cutoff_TiO, nTi5c, natoms, atype, index_O2c
  IMPLICIT NONE
  INTEGER, INTENT(IN)         :: index_OwTi
  INTEGER                     :: iat, CoordNumb, iO2c, iO
  REAL*8                      :: path_O,path_H
  LOGICAL                     :: is_terminal_hydroxyl, Hbond_Connected

  !The code will look for the minimum proton wire path connection OwTi and O2c
  path_O=100.d0

  !We first need to know if the OwTi is OH2 or OH (terminal hydroxide)
  IF (CoordNumb(index_OwTi,2,cutoff_OH)==1) THEN
    is_terminal_hydroxyl=.true.
  ELSE
    is_terminal_hydroxyl=.false.
  END IF

  !Look for a reactive proton wire from Ow to O2c with 3 members
  !OwTi is the first member, some O2c is the last
  DO iO=1,natoms
    IF ((atype(iO)==3).and.(CoordNumb(iO,1,cutoff_TiO)==0)) THEN
      IF (HBond_Connected(index_OwTi,iO,is_terminal_hydroxyl)) THEN
        DO iO2c=1,nTi5c
          IF (HBond_Connected(iO,index_O2c(iO2c),is_terminal_hydroxyl)) THEN
            CALL Compute_Proton_Write_CV(index_OwTi,iO,index_O2c(iO2c),path_O,path_H)
          END IF
        END DO
      END IF
    END IF
  END DO

END SUBROUTINE proton_wire_cvs

SUBROUTINE Compute_Proton_Write_CV(index_OwTi,iO,iO2c,path_O,path_H)
  USE parameters, only : atype 
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: index_OwTi,iO,iO2c
  REAL*8                     :: path_O, path_H, dOOO
  REAL*8                     :: d1(2),d2(2), Dist

  dOOO=Dist(index_OwTi,iO)+Dist(iO,iO2c)

  !Only the shortest proton wire will be considered
  IF (dOOO<path_O) THEN
    CALL MIN_H_PATH(index_OwTi,iO,d1(1),d2(1))
    CALL MIN_H_PATH(iO,iO2c,d1(2),d2(2))
    path_O=dOOO
    path_H=(d1(1)-d2(1)+d1(2)-d2(2))/2.d0
  END IF

END SUBROUTINE Compute_Proton_Write_CV

SUBROUTINE MIN_H_PATH(O1,O2,d1,d2)
  USE parameters, ONLY : natoms, atype
  IMPLICIT NONE
  INTEGER,INTENT(IN)         :: O1, O2
  REAL*8                     :: d1,d2,d1_,d2_,Dist
  REAL*8                     :: dsum,dsum_min
  INTEGER                    :: iH

  dsum_min=100.

  DO iH=1,natoms
    IF (atype(iH)==2) THEN
      d1_=Dist(iH,O1)
      d2_=Dist(iH,O2)
      if (dsum<dsum_min) then
        d1=d1_
        d2=d2_
        dsum_min=dsum
      end if
    END IF
  END DO

END SUBROUTINE MIN_H_PATH

LOGICAL FUNCTION Hbond_Connected(O1,O2,O2_donates_Hbond_to_O1)
  USE parameters, only : natoms,atype, cutoff_OH
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: O1, O2
  LOGICAL, INTENT(IN)        :: O2_donates_Hbond_to_O1
  LOGICAL                    :: is_connected
  INTEGER                    :: iH, Hb_donor, Hb_acceptor
  REAL*8                     :: Dist,Anglei
  REAL*8,PARAMETER           :: LZ_angle=0.866d0,LZ_distance=3.5d0

  is_connected=.false.

  IF (Dist(O1,O2)<LZ_distance) THEN
    
    IF (O2_donates_Hbond_to_O1) THEN
      Hb_donor=O2
      Hb_acceptor=O1
    ELSE
      Hb_donor=O1
      Hb_acceptor=O2
    END IF

    DO iH=1,natoms
      IF (atype(iH)==2) THEN
        IF ((Dist(Hb_donor,iH)<cutoff_OH).and. &
             (Anglei(Hb_donor,iH,Hb_acceptor)>LZ_angle)) THEN
          is_connected=.true.
          GO TO 182
        END IF
      END IF
    END DO

  END IF

182  Hbond_Connected=is_connected

END FUNCTION HBond_Connected

SUBROUTINE Find_O2c
  USE parameters, ONLY : natoms, index_O2c, &
                         cutoff_TiO, atype, nTi5c
  IMPLICIT NONE
  INTEGER                    :: iat,iO
  INTEGER                    :: CoordNumb
  INTEGER                    :: bufferO2c(natoms)
  !
  iO = 0
  !
  DO iat=1,natoms
    !
    !Find 2 coordinated O atoms
    IF ( atype(iat) == 3 ) THEN
      !
      IF (CoordNumb(iat,1,cutoff_TiO)==2) THEN
        !
        iO=iO+1
        bufferO2c(iO)=iat
        !
      END IF
      !
    END IF
    !
  END DO
  !
  ALLOCATE(index_O2c(iO));   index_O2c=bufferO2c(1:iO)
  nTi5c=iO
  !
END SUBROUTINE Find_O2c

SUBROUTINE Find_OwTi(index_OwTi)
  USE parameters, ONLY : natoms, cutoff_OH, atype, cutoff_TiO, nTi5c
  IMPLICIT NONE
  INTEGER                    :: index_OwTi(nTi5c)
  INTEGER                    :: iat, iO
  INTEGER                    :: CoordNumb
  !
  iO = 0
  index_OwTi=0
  !
  DO iat=1,natoms
    !
    !Find coordinated OHx bound to Ti5c
    IF ( atype(iat) == 3 ) THEN
      !
      IF ((CoordNumb(iat,1,cutoff_TiO)==1).and.(CoordNumb(iat,2,cutoff_OH)>0)) THEN
        !
        iO=iO+1
        index_OwTi(iO)=iat
        !
      END IF
      !
    END IF
    !
  END DO
  !
END SUBROUTINE Find_OwTi
