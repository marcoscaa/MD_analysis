MODULE parameters_pt
  IMPLICIT NONE
  REAL*8, PARAMETER                        :: cutoff_OH=1.2,cutoffOwTi=2.6
  REAL*8, PARAMETER                        :: cutoff_TiOti=2.6
  REAL*8, PARAMETER                        :: beta=1.d0/(0.00831446*330)
  REAL*8, DIMENSION(:,:), ALLOCATABLE      :: pos
  REAL*8, DIMENSION(3)                     :: box
  INTEGER, DIMENSION(:), ALLOCATABLE       :: ind_atom
  INTEGER                                  :: natoms, nframes, tstep,nequil
  INTEGER                                  :: nTi5c
  INTEGER, ALLOCATABLE                     :: index_ti5c(:),index_o2c(:)
  INTEGER, ALLOCATABLE                     :: index_o2c_ti(:,:)
  LOGICAL, PARAMETER                       :: biased_dynamics=.False.

CONTAINS

SUBROUTINE READ_FRAME
  !Read LAMMPS atom file
  IMPLICIT NONE
  INTEGER                    :: iat, ind, atp
  REAL*8                     :: box_tmp(2)

  DO iat=1,5
    READ(1,*)
  END DO

  !Read Box
  DO iat=1,3
    READ(1,*) box_tmp(1), box_tmp(2)
    box(iat) = box_tmp(2)-box_tmp(1)
  END DO

  READ(1,*)

  DO iat = 1, natoms
    READ(1,*) ind, atp, pos(:,ind)
    pos(:,ind) = pos(:,ind) * box !Non-reduced coordinates
  END DO

END SUBROUTINE

SUBROUTINE find_index_OwTi(ind_owti)
  IMPLICIT NONE
  INTEGER              :: ind_owti(nTi5c), iat
  REAL*8               :: d(3)

  ind_owti=0

  DO iat=1,nTi5c
    CALL mindist_index_vector(index_ti5c(iat),4,ind_owti(iat),d)
    IF (norm2(d) > 3. ) ind_owti(iat)=0
  END DO

END SUBROUTINE find_index_OwTi

SUBROUTINE min_known_path_cvs(ind_ini,cv_O,cv_H, ttype)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ind_ini
  INTEGER             :: i, iO(2), iH(2)
  INTEGER             :: ttype
  REAL*8              :: dooo
  REAL*8              :: cv_O, cv_H

  iH=0
  iO=0

  !Find minimum OOO path, starting from a O2c 
  !CALL MinSumDist(ind_ini,4,ibr,ind_end,dooo)
  !CALL Hs_along_the_path(ind_ini,ind_end,ibr,iH)
  CALL MinSumDistFull(ind_ini,iO,iH)

  IF ((ANY(iH==0)).or.(ANY(iO==0))) THEN
    !PRINT *, 'Could not find path'
    cv_O=0
    cv_H=0
    ttype=-1
    RETURN
    !STOP
  END IF

  dooo=Dist(ind_ini,iO(1))+Dist(iO(1),iO(2))

  cv_O = 0.5 * dooo
  cv_H = 0.5 * ((Dist(ind_ini,iH(1))-Dist(iO(1),iH(1))) &
              + (Dist(iO(1),iH(2))-Dist(iO(2),iH(2))))

  ttype=0
  if ( ANY(index_O2c==iO(2)) ) ttype=iO(2)

END SUBROUTINE min_known_path_cvs

SUBROUTINE Hs_along_the_path(iini,iend,ibr,iH)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iini, iend,ibr
  INTEGER             :: iat, iH(2)
  REAL*8              :: v1, v2

  CALL MinSumDistHyd(iini,ibr,iend,iH)

  IF ( iH(1).eq.iH(2) ) THEN
    PRINT *, "EQUAL H", iH
    !STOP
  ELSEIF (ANY(iH==0)) THEN
    PRINT *, "Could not find H"
    STOP
  ENDIF

END SUBROUTINE Hs_along_the_path

SUBROUTINE MinSumDist(ifix_1,atype,ind1_min,ind2_min,dmin)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ifix_1, atype
  INTEGER             :: ind1_min, ind2_min, ibr, iend
  REAL*8              :: dmin, d, dfix_br

  dmin=1000.

  DO ibr=1,natoms

    IF (ind_atom(ibr).ne.atype) CYCLE
    dfix_br=Dist(ifix_1,ibr)
    IF (dfix_br>4.d0) CYCLE

    DO iend=1,natoms

      IF ((ind_atom(iend).ne.3).and.(ind_atom(iend).ne.4)) CYCLE
      IF ((iend.eq.ifix_1).or.(iend.eq.ibr)) CYCLE

      d=dfix_br+Dist(ibr,iend)

      IF (d<dmin) THEN
        dmin=d
        ind1_min=ibr
        ind2_min=iend
      END IF

    END DO

  END DO

END SUBROUTINE MinSumDist

SUBROUTINE MinSumDistHyd(ifix_1,ifix_2,ifix_3,indH)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ifix_1, ifix_2, ifix_3
  INTEGER             :: indH(2), ih1, ih2
  REAL*8              :: dmin, d
  REAL*8              :: dh1O2c,dh1Obr
  REAL*8              :: dh2Obr,dh2Oend

  dmin=1000.

  DO ih1=1,natoms

    if (ind_atom(ih1).ne.2) CYCLE
    dh1O2c=Dist(ifix_1,ih1)
    dh1Obr=Dist(ih1,ifix_2)
    if ((dh1O2c>1.3).and.(dh1Obr>1.3)) CYCLE

    DO ih2=1,natoms

      IF ((ind_atom(ih2).ne.2).or.(ih1.eq.ih2)) CYCLE
      dh2Obr=Dist(ifix_2,ih2)
      dh2Oend=Dist(ih2,ifix_3)
      if ((dh2Obr>1.3).and.(dh2Oend>1.3)) CYCLE
      !if ((dh1Obr<1.1).and.(dh2Obr<1.1)) CYCLE

      d=dh1O2c+dh1Obr+dh2Obr+dh2Oend

      IF (d<dmin) THEN
        dmin=d
        indH(1)=ih1
        indH(2)=ih2
      END IF

    END DO

  END DO

END SUBROUTINE MinSumDistHyd

SUBROUTINE MinSumDistFull(ifix_1,indO,indH)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ifix_1
  INTEGER             :: indO(2), indH(2)
  INTEGER             :: ih1, ih2, iend, ibr
  REAL*8              :: dmin, d
  REAL*8              :: dh1O1,dh1O2,dO1O2
  REAL*8              :: dh2O2,dh2O3,dO2O3

  dmin=1000.

  DO ih1=1,natoms

    IF (ind_atom(ih1).ne.2) CYCLE
    dh1O1=Dist(ifix_1,ih1)
    IF (dh1O1>3.) CYCLE

    DO ibr=1,natoms

      IF (ind_atom(ibr).ne.4) CYCLE
      IF (ibr.eq.ifix_1) CYCLE
      dO1O2=Dist(ifix_1,ibr)
      dh1O2=Dist(ibr,ih1)
      IF (dO1O2>4.d0) CYCLE
      IF (dh1O2>3.) CYCLE

      DO ih2=1,natoms

        IF ((ind_atom(ih2).ne.2).or.(ih1.eq.ih2)) CYCLE
        
        dh2O2=Dist(ibr,ih2)
        IF (dh2O2>3.) CYCLE

        DO iend=1,natoms
 
          IF ((ind_atom(iend).ne.3).and.(ind_atom(iend).ne.4)) CYCLE
          IF ((iend.eq.ifix_1).or.(iend.eq.ibr)) CYCLE
          dh2O3=Dist(ih2,iend)
          IF (dh2O3>3.) CYCLE
          IF ( ( abs((dh1O1-dh1O2)-(dh2O2-dh2O3))>0.4 ) .and. &
               ( (dh1O1-dh1O2)*(dh2O2-dh2O3)<0.d0 ) ) CYCLE 
 
          d=dh1O1+dh1O2+dh2O2+dh2O3

          IF (d<dmin) THEN
            dmin=d
            indO(1)=ibr
            indO(2)=iend
            indH(1)=ih1
            indH(2)=ih2
          END IF

        END DO

      END DO
  
    END DO
  
  END DO

END SUBROUTINE MinSumDistFull

SUBROUTINE mindist_index_vector(ind_center, type_mindist, ind_mindist, vec_mindist)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ind_center, type_mindist
  INTEGER             :: ind_mindist, iat
  REAL*8              :: vec_mindist(3)
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

  CALL DistVect(ind_center,ind_mindist, vec_mindist)

END SUBROUTINE mindist_index_vector

REAL*8 FUNCTION Dist(inda, indb)
  ! Distance between two points including pbc
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: inda, indb
  REAL*8, DIMENSION(3)       :: xyz
  INTEGER                    :: i

  DO i = 1,3

    xyz(i) = pos(i,inda) - pos(i,indb)
    xyz(i) = xyz(i) - nint( xyz(i)/box(i) ) * box(i)

  END DO

  Dist = SQRT( SUM(xyz*xyz) )

END FUNCTION Dist

REAL*8 FUNCTION Dist2(vec1, vec2)
  ! Distance between two points including pbc
  IMPLICIT NONE
  REAL*8, INTENT(IN)         :: vec1(3), vec2(3)
  REAL*8, DIMENSION(3)       :: xyz
  INTEGER                    :: i

  DO i = 1,3

    xyz(i) = vec1(i) - vec2(i)
    xyz(i) = xyz(i) - nint( xyz(i)/box(i) ) * box(i)

  END DO

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
    xyz(i) = xyz(i) - nint( xyz(i)/box(i) ) * box(i)

  END DO

  vect = xyz 

END SUBROUTINE DistVect

SUBROUTINE CenterOfLine(inda, indb, vect)
  ! center of line between a to b including pbc
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: inda, indb
  REAL*8, DIMENSION(3)       :: xyz, vect
  INTEGER                    :: i

  DO i = 1,3

    xyz(i) = pos(i,inda) + pos(i,indb)
    xyz(i) = xyz(i) - nint( xyz(i)/box(i) ) * box(i)

  END DO

  vect = xyz 

END SUBROUTINE CenterOfLine

SUBROUTINE Find_Ti5c_and_O2c (coord_file)
  IMPLICIT NONE
  CHARACTER(100), INTENT(IN) :: coord_file 
  INTEGER                    :: iat, iTi, iO, iO2c, o2cbound, buffer
  INTEGER                    :: bufferTi(natoms), bufferO2c(natoms)
  REAL(8)                    :: d
  !
  !Read 1 frame to get box and coordinates
  OPEN(unit  =  1,file  =  coord_file)
  CALL READ_FRAME
  CLOSE(1)
  !
  iO = 0
  iTi= 0
  !
  DO iat=1,natoms
    !
    !Find 5 coordinated Ti atoms
    !IF ( ( ind_atom(iat) == 1 ) .and. ( pos(3,iat) .lt. 15. ) ) THEN
    IF ( ind_atom(iat) == 1 ) THEN
      !
      IF (CordNumb(iat,3,cutoff_TiOti)==5) THEN
        !
        iTi=iTi+1
        bufferTi(iTi)=iat
        !
      END IF
      !
    !Find 2 coordinated O atoms
    !ELSE IF ( ( ind_atom(iat) == 3 ) .and. ( pos(3,iat) .lt. 15. ) ) THEN
    ELSE IF ( ind_atom(iat) == 3 ) THEN
      !
      IF (CordNumb(iat,1,cutoff_TiOti)==2) THEN
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
  iTi=6;iO=6
  nTi5c=iTi !Hardcoded
  ALLOCATE(index_ti5c(iTi));     index_Ti5c=bufferTi(1:iTi)
  index_Ti5c=(/ 1,2,3,31,32,33 /)
  ALLOCATE(index_O2c(iTi));       index_O2c=bufferO2c(1:iTi)
  index_O2c=(/ 243,244,245,303,304,305 /)
  ALLOCATE(index_o2c_ti(5,iTi)); index_o2c_ti=0
  !
  IF (iTi.ne.iO) THEN
    PRINT *, "STOP!!!"
    PRINT *, "Inconsistent assignment of undercoordinated Ti and O", iTi, iO
    STOP
  END IF
  !
  DO iTi=1,nTi5c
    !
    iO2c = 0
    o2cbound=0
    !print *, iTi
    !
    DO iO=1,nTi5c
      !
      d=Dist(index_Ti5c(iTi),index_O2c(iO))
      IF ( d < 4.6 ) THEN
        iO2c=iO2c+1
        index_o2c_ti(iO2c,iTi) = index_O2c(iO)
        IF ( d < 2.5 ) o2cbound = index_O2c(iO)
      END IF
      !
    END DO
    !
    IF ( iO2c .ne. 5 ) THEN
      PRINT *, "Inconsistent assignment of undercoordinated Ti and O", iO2c
      STOP
    END IF
    !
    iO2c = 0
    !
    ! Indexes of two adjacent O2c go first
    DO iO=1,5
      !
      IF ( index_o2c_ti(iO,iTi) .ne. o2cbound ) THEN
        !
        d = Dist( o2cbound, index_o2c_ti(iO,iTi) )
        !
        IF ( d > 4.5 ) THEN 
          !
          iO2c = iO2c + 1
          buffer = index_o2c_ti(iO2c,iTi)
          index_o2c_ti(iO2c, iTi) = index_o2c_ti(iO,iTi)
          index_o2c_ti(iO,iTi) = buffer
          !
        END IF
        !
      END IF
      !
    END DO
    !
    IF ( iO2c .ne. 2 ) THEN
      PRINT *, "Inconsistent assignment of undercoordinated Ti and O"
      STOP
    END IF
    !
    !print *, index_o2c_ti(:,iTi)
    !
  END DO
  !
END SUBROUTINE Find_Ti5c_and_O2c

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

REAL*8 FUNCTION CordNumbCont(ind)
  ! Number of nearest H - continuous version
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: ind
  INTEGER                    :: iat
  REAL*8,PARAMETER           :: r0=1.2d0 !angstrom
  REAL*8                     :: d, d32, d42, cordnum

  cordnum=0.

  DO iat = 1,natoms
    IF (ind_atom(iat)==2) THEN
      d=Dist(ind,iat)/r0
      !d32=d**2;d32=d32*d32;d32=d32*d32;d32=d32*d32;d32=d32*d32
      d32=d**24!28
      d42=d**48!38
      cordnum = cordnum + (1.d0-d32)/(1.d0-d42)
    ENDIF

  ENDDO
 
  CordNumbCont=cordnum 

END FUNCTION CordNumbCont

SUBROUTINE REMOVE_EQUIL
  IMPLICIT NONE
  INTEGER                    :: i

  DO i = 1,nequil*(natoms+9)
    READ(1,*)
  END DO

END SUBROUTINE REMOVE_EQUIL

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

INTEGER FUNCTION AtomClosestToLine(line0,linef,atype)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: line0, linef, atype
  INTEGER             :: iat, ind_min
  REAL*8              :: d, dmin
  REAL*8              :: center_of_line(3)

  dmin = 100.

  CALL CenterOfLine(line0,linef,center_of_line)

  DO iat=1,natoms
    
    IF  ( ind_atom(iat) == atype ) THEN

      d=Dist2(pos(:,iat),center_of_line)

      IF ( d .lt. dmin ) THEN 
        dmin=d
        ind_min=iat
      END IF

    END IF

  END DO

  !IF ( (dmin>1.) .or. (ALL(doh>1.4)) ) ind_min=0
  AtomClosestToLine = ind_min
    
END FUNCTION AtomClosestToLine

INTEGER FUNCTION Count_Hydroxyl(ind_O2c)
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: ind_O2c(nTi5c)
  INTEGER              :: iat, cn

  cn = 0

  DO iat=1,nTi5c

    !IF ( CordNumb(ind_O2c(iat),2,cutoff_OH)==1 ) cn = cn + 1 
    IF ( CordNumb(ind_O2c(iat),2,1.15d0)==1 ) cn = cn + 1 

  END DO

  Count_Hydroxyl = cn

END FUNCTION Count_Hydroxyl

REAL(8) FUNCTION BiasPotential(index_O2c,nsites)
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: nsites, index_O2c(nsites)
  INTEGER                    :: iat, iO
  REAL(8)                    :: mindist, TotalBias

  TotalBias=0.d0

  DO iO=1,nsites
    
    mindist = ContinuousMindist(index_O2c(iO),2)
    TotalBias = TotalBias + Bias(mindist)

  END DO

  BiasPotential = TotalBias

END FUNCTION BiasPotential

REAL(8) FUNCTION Bias(x)
  IMPLICIT NONE
  REAL(8), INTENT(IN)        :: x
  
  Bias = Gauss(x,-18.d0,1.17d0,0.04d0) + Gauss(x,-2.4d0,1.36d0,0.04d0)

END FUNCTION Bias

REAL(8) FUNCTION Gauss(x,A,x0,B)
  IMPLICIT NONE
  REAL(8), INTENT(IN)        :: x, A, x0,B

  Gauss = A * exp( -( (x-x0)**2 )/B )

END FUNCTION Gauss

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

END MODULE parameters_pt
