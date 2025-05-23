MODULE parameters_pt
  IMPLICIT NONE
  REAL*8, PARAMETER                        :: cutoff_OH=1.2,cutoffOwTi=2.6
  REAL*8, PARAMETER                        :: cutoff_TiS=3.5
  REAL*8, PARAMETER                        :: beta=1.d0/(0.00831446*330)
  REAL*8, DIMENSION(:,:), ALLOCATABLE      :: pos
  REAL*8, DIMENSION(3,3)                   :: box
  INTEGER, DIMENSION(:), ALLOCATABLE       :: ind_atom
  INTEGER                                  :: natoms, nframes, tstep,nequil
  INTEGER                                  :: nTi5c
  INTEGER, ALLOCATABLE                     :: index_ti5c(:),index_S2c(:)
  INTEGER, PARAMETER                       :: typeTi=2,typeS=1,typeO=3,typeH=4 
CONTAINS

SUBROUTINE REMOVE_EQUIL
  IMPLICIT NONE
  INTEGER                    :: i

  DO i = 1,nequil*(natoms+9)
    READ(1,*)
  END DO

END SUBROUTINE REMOVE_EQUIL

SUBROUTINE READ_FRAME
  !Read LAMMPS atom file
  IMPLICIT NONE
  INTEGER                    :: iat, ind
  REAL*8                     :: box_tmp(2), bd(3)
  CHARACTER(10)              :: boxtype, trash

  DO iat=1,4
    READ(1,*)
  END DO

  READ(1,*) trash, trash, trash, boxtype

  box=0.0;bd=0.0
  IF( trim(boxtype) .eq. 'xy' ) THEN
  !Non-Orthorhombic case

    DO iat=1,3
      READ(1,*) box_tmp(1), box_tmp(2), bd(iat)
      box(iat,iat) = box_tmp(2)-box_tmp(1)
    END DO

    box(1,1) = box(1,1) - max(0.0,bd(1),bd(2),bd(1)+bd(2)) + min(0.0,bd(1),bd(2),bd(1)+bd(2))
    box(2,2) = box(2,2) - max(0.0,bd(3)) + min(0.0,bd(3))
    box(2,1) = bd(1)
    box(3,1) = bd(2)
    box(3,2) = bd(3)

  ELSE
  !Orthorhombic case

    DO iat=1,3
      READ(1,*) box_tmp(1), box_tmp(2)
      box(iat,iat) = box_tmp(2)-box_tmp(1)
    END DO

  END IF

  READ(1,*)

  DO iat = 1, natoms
    READ(1,*) ind, ind_atom(ind), pos(:,ind)
    !pos(:,ind) = pos(:,ind) * box !Non-reduced coordinates
  END DO

END SUBROUTINE

SUBROUTINE COMPUTE_CVS(ind_S2c,cv1,cv2)
  IMPLICIT NONE
  INTEGER, INTENT(IN)    :: ind_S2c
  REAL*8, INTENT(INOUT)  :: cv1,cv2
  REAL*8                 :: d_SH, d_HO
  INTEGER                :: iS, iat, i_H, i_O

  CALL mindist_index(ind_S2c,typeH,i_H,d_SH)
  CALL mindist_index(i_H,typeO,i_O,d_HO)

  cv1 = d_SH - d_HO
  cv2 = d_SH + d_HO

END SUBROUTINE COMPUTE_CVS

SUBROUTINE mindist_index(ind_center, type_mindist, ind_mindist, dmin)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ind_center, type_mindist
  INTEGER             :: ind_mindist, iat
  REAL*8              :: dmin, d

  dmin=100.

  DO iat=1,natoms
    IF ( (ind_atom(iat)==type_mindist ) .and. (iat.ne.ind_center) ) THEN
      d=Dist(ind_center, iat) 
      IF ( d < dmin ) THEN
        dmin=d 
        ind_mindist=iat
      ENDIF
    ENDIF
  ENDDO

  IF (dmin>10) THEN
    print *, 'Too large mindist. Stop!'
    STOP
  ENDIF

END SUBROUTINE mindist_index

REAL*8 FUNCTION Dist(inda, indb)
  ! Distance between two points including pbc
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: inda, indb
  REAL*8, DIMENSION(3)       :: xyz
  INTEGER                    :: i

  DO i = 1,3

    xyz(i) = pos(i,inda) - pos(i,indb)
    xyz(i) = xyz(i) - nint( xyz(i) ) 

  END DO

  xyz=matmul(box,xyz)
  Dist = SQRT( SUM(xyz*xyz) )

END FUNCTION Dist

SUBROUTINE Find_Ti5c_and_S2c (coord_file)
  IMPLICIT NONE
  CHARACTER(100), INTENT(IN) :: coord_file 
  INTEGER                    :: iat, iTi, iS, iS2c, S2cbound, buffer
  INTEGER                    :: bufferTi(natoms), bufferS2c(natoms)
  REAL(8)                    :: d
  !
  !Read 1 frame to get box and coordinates
  OPEN(unit  =  1,file  =  coord_file)
  CALL READ_FRAME
  CLOSE(1)
  !
  iS = 0
  iTi= 0
  !
  DO iat=1,natoms
    !
    !Find uncoordinated Ti atoms
    IF ( ind_atom(iat) == typeTi ) THEN
      !
      IF (CordNumb(iat,typeS,cutoff_TiS)/=6) THEN
        !
        iTi=iTi+1
        bufferTi(iTi)=iat
        !
      END IF
      !
    !Find 2-fold coordinated S atoms
    ELSE IF ( ind_atom(iat) == typeS ) THEN
      !
      IF (CordNumb(iat,typeTi,cutoff_TiS)/=3) THEN
        !
        iS=iS+1
        bufferS2c(iS)=iat
        !
      END IF
      !
    END IF
    !
  END DO
  !
  nTi5c=iTi !Hardcoded
  ALLOCATE(index_ti5c(iTi));     index_Ti5c=bufferTi(1:iTi)
  ALLOCATE(index_S2c(iTi));       index_S2c=bufferS2c(1:iTi)
  !
  IF ((iTi.ne.iS).and.(2*iTi.ne.iS)) THEN
    PRINT *, "STOP!!!"
    PRINT *, "Inconsistent assignment of undercoordinated Ti and O", iTi, iS
    STOP
  END IF
  !
END SUBROUTINE Find_Ti5c_and_S2c

INTEGER FUNCTION CordNumb(ind,atype,cutoff)
  ! Number of nearest atoms with type atype 
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: ind, atype
  INTEGER                    :: iat, cn
  REAL*8, INTENT(IN)         :: cutoff !angstrom

  cn=0

  DO iat = 1,natoms
    IF (ind_atom(iat)==atype) THEN
      IF (Dist(ind,iat)<=cutoff) cn=cn+1
    ENDIF
  ENDDO
 
  CordNumb=cn 

END FUNCTION CordNumb

INTEGER FUNCTION Count_Hydroxyl(ind_S2c)
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: ind_S2c(nTi5c)
  INTEGER              :: iat, cn

  cn = 0

  DO iat=1,nTi5c

    !IF ( CordNumb(ind_S2c(iat),2,cutoff_OH)==1 ) cn = cn + 1 
    IF ( CordNumb(ind_S2c(iat),typeH,1.15d0)==1 ) cn = cn + 1 

  END DO

  Count_Hydroxyl = cn

END FUNCTION Count_Hydroxyl

!!!!!!
!!!!!!!!!!!Unused subroutines and functions!!!!!!!!!!!!!
!!!!!

SUBROUTINE find_index_OwTi(ind_owti)
  IMPLICIT NONE
  INTEGER              :: ind_owti(nTi5c), iat
  REAL*8               :: d

  ind_owti=0

  DO iat=1,nTi5c
    CALL mindist_index(index_ti5c(iat),typeO,ind_owti(iat),d)
    IF (d > 3.) ind_owti(iat)=0
  END DO

END SUBROUTINE find_index_OwTi

REAL(8) FUNCTION ContinuousMindist(indat, atype)
  INTEGER, INTENT(in)        :: indat, atype
  INTEGER                    :: iat
  REAL(8)                    :: d
  REAL(8), PARAMETER         :: beta = 500.d0

  d=0.d0

  DO iat = 1, natoms 
  
    IF (ind_atom(iat).eq.atype) d = d + exp(beta/Dist(indat,iat))

  END DO

  ContinuousMindist = beta / ( log(d) )
  
END FUNCTION ContinuousMindist

REAL*8 FUNCTION mindist(ind, type_mindist)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ind, type_mindist
  INTEGER             :: iat
  REAL*8              :: dmin, d

  dmin=100.

  DO iat=1,natoms
    IF ( (ind_atom(iat)==type_mindist ) .and. (iat.ne.ind) ) THEN
      d=Dist(ind, iat) 
      IF ( d < dmin ) THEN
        dmin=d 
      ENDIF
    ENDIF
  ENDDO

  mindist=dmin 

END FUNCTION mindist

REAL*8 FUNCTION CordNumbCont(ind)
  ! Number of nearest H - continuous version
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: ind
  INTEGER                    :: iat
  REAL*8,PARAMETER           :: r0=1.2d0 !angstrom
  REAL*8                     :: d, d32, d42, cordnum

  cordnum=0.

  DO iat = 1,natoms
    IF (ind_atom(iat)==typeH) THEN
      d=Dist(ind,iat)/r0
      !d32=d**2;d32=d32*d32;d32=d32*d32;d32=d32*d32;d32=d32*d32
      d32=d**24!28
      d42=d**48!38
      cordnum = cordnum + (1.d0-d32)/(1.d0-d42)
    ENDIF

  ENDDO
 
  CordNumbCont=cordnum 

END FUNCTION CordNumbCont

REAL*8 FUNCTION Dist2(vec1, vec2)
  ! Distance between two points including pbc
  IMPLICIT NONE
  REAL*8, INTENT(IN)         :: vec1(3), vec2(3)
  REAL*8, DIMENSION(3)       :: xyz
  INTEGER                    :: i

  DO i = 1,3

    xyz(i) = vec1(i) - vec2(i)
    xyz(i) = xyz(i) - nint( xyz(i) ) 

  END DO

  xyz = matmul(box,xyz)
  Dist2 = SQRT( SUM(xyz*xyz) )

END FUNCTION Dist2

SUBROUTINE DistVect(inda, indb, vect)
  ! Vector from a to b including pbc
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: inda, indb
  REAL*8, DIMENSION(3)       :: xyz, vect
  INTEGER                    :: i

  DO i = 1,3

    xyz(i) = pos(i,inda) - pos(i,indb)
    xyz(i) = xyz(i) - nint( xyz(i) ) 

  END DO

  vect = matmul(box,xyz) 

END SUBROUTINE DistVect

END MODULE parameters_pt
