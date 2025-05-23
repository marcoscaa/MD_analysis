!All units in atomic units

MODULE histogram
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: wfc(:,:), wfc_field(:,:,:) 
  DOUBLE PRECISION, ALLOCATABLE            :: bec(:,:,:) 
  DOUBLE PRECISION, ALLOCATABLE            :: dipole(:,:,:) 
  DOUBLE PRECISION                         :: efield(3)
  INTEGER                                  :: nwfc
  INTEGER, ALLOCATABLE                     :: nwann(:)
  LOGICAL, ALLOCATABLE                     :: is_dipole(:)
  LOGICAL, PARAMETER                       :: include_induced_field=.False.
  LOGICAL, PARAMETER                       :: compute_bec=.False.
END MODULE histogram

PROGRAM SFG
  USE parameters, ONLY : nframes
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL_WFC

  DO frame = 1,nframes
    CALL READ_POS_WFC_CEL
    CALL READ_POL_BEC
    !CALL COARSE_GRAIN_POS (frame)
    CALL ASSIGN_DIPOLES
    !CALL WANNIER_DISPLACEMENT 
    CALL MOLECULAR_DIPOLE_POLARIZABILITY
  END DO

  CALL The_End


END PROGRAM SFG

SUBROUTINE REMOVE_EQUIL_WFC
  USE parameters, ONLY : natoms, nequil, nframes
  USE histogram,  ONLY : nwfc
  IMPLICIT NONE
  INTEGER                    :: iwfc, iat

  DO iwfc = 1, nequil*(nwfc+1)
    READ(2,*)
  END DO

  DO iat = 1, nequil*(natoms+1)
    READ(1,*)
  END DO

  nframes = nframes-nequil

END SUBROUTINE REMOVE_EQUIL_WFC

SUBROUTINE INITIALIZE
  USE histogram,  ONLY : wfc, nwfc, bec, & 
                         wfc_field, dipole, is_dipole, &
                         compute_bec
  USE parameters, ONLY : natoms, nframes, nwater, &
                         pos, coarse, nlayers
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file
  CHARACTER(100)             :: wfc_file, bec_file, pol_file 
  CHARACTER(100)             :: cel_file 

  CALL getarg(1, pos_file)
  CALL getarg(2, wfc_file)
  CALL getarg(3, pol_file)
  CALL getarg(4, cel_file)
  CALL getarg(5, index_file)
  IF (compute_bec) CALL getarg(6, bec_file)

  CALL READ_INDEX_WFC ( index_file )

  ALLOCATE(pos(natoms,3)) ! coordinates
  IF ( nlayers > 1 ) THEN 
    ALLOCATE(coarse(nframes,nwater)); coarse=0
  END IF
  ALLOCATE(wfc(nwfc,3)); wfc=0
  ALLOCATE(wfc_field(nwfc,3,3)); wfc_field=0
  IF (compute_bec) ALLOCATE(bec(natoms,3,3)); bec=0
  ALLOCATE(is_dipole(natoms)); is_dipole=.false.
  ALLOCATE(dipole(3,4,natoms)); dipole=0.d0

  OPEN(unit = 1,file = pos_file)
  OPEN(unit = 2,file = wfc_file)
  OPEN(unit = 3,file = pol_file)
  IF (compute_bec) OPEN(unit = 4,file = bec_file)
  OPEN(unit = 9,file = cel_file)
  OPEN(unit = 5,file = '01.dipole.dat')
  OPEN(unit = 7,file = '03.born_effective_charge.dat')
  OPEN(unit = 8,file = '04.wannier_distance.dat')
  OPEN(unit =50,file = '05.polarizability_labframe.dat')
  
END SUBROUTINE INITIALIZE

SUBROUTINE READ_INDEX_WFC ( index_file ) 
  USE histogram,  ONLY : nwfc, efield, nwann
  USE parameters, ONLY : nlayers, nframes, nequil,&
                         natoms, layers, box, ind_atom, ntype
  IMPLICIT NONE
  CHARACTER(100), INTENT(IN) :: index_file
  INTEGER                    :: i,j
  
  OPEN(2, file=index_file)

  READ (2, *) natoms, ntype, nframes, nequil, nwfc, nlayers
  READ (2, *) efield 
  ALLOCATE (nwann(ntype))
  READ (2, *) nwann
  print *, nwann

  IF ( nlayers > 1 ) THEN
    ALLOCATE( layers(nlayers+1) )
    READ(2, *) layers
  END IF

  ALLOCATE( ind_atom(natoms) ) 

  DO i = 1, natoms
    READ(2,*) ind_atom(i)
  END DO

  CLOSE(2)

END SUBROUTINE READ_INDEX_WFC

SUBROUTINE READ_POS_WFC_CEL
  USE histogram, ONLY : wfc, nwfc 
  USE parameters, ONLY : pos, natoms, box
  IMPLICIT NONE
  INTEGER                    :: i, ipol

  READ(1,*)
  READ(2,*)
  READ(9,*)

  DO i = 1, natoms
    READ(1, *), ( pos(i,ipol), ipol=1,3 )
  END DO

  DO i = 1, nwfc
    READ(2, *), ( wfc(i,ipol), ipol=1,3 )
  END DO

  DO i = 1,3 
    READ(9, *), box(1,i), box(2,i), box(3,i) 
  END DO

END SUBROUTINE

SUBROUTINE READ_POL_BEC
  USE histogram, ONLY : nwfc, wfc_field, bec, &
                        compute_bec 
  USE parameters, ONLY : natoms
  IMPLICIT NONE
  INTEGER                    :: i, ipol1, ipol2

  READ(3,*)

  DO i = 1, nwfc
    READ(3, *), ( ( wfc_field(i,ipol1,ipol2), ipol2=1,3 ), ipol1=1,3)
  END DO
 
  IF (compute_bec) THEN
    READ(4,*)
    DO i = 1, natoms
      READ(4, *), ( ( bec(i,ipol1,ipol2), ipol2=1,3 ), ipol1=1,3)
    END DO
  END IF

END SUBROUTINE

SUBROUTINE ASSIGN_DIPOLES
  !Transform a system of point charges into a system of dipoles
  USE parameters, ONLY : natoms, ind_atom
  USE histogram,  ONLY : dipole, is_dipole, nwann
  IMPLICIT NONE
  INTEGER                    :: iat

  dipole = 0.d0
  is_dipole = .False.

  DO iat=1,natoms

    if (nwann(ind_atom(iat))/=0) THEN
      is_dipole(iat) = .true.
      CALL SUM_NEAREST_WFCs( iat, dipole(:,:,iat), nwann(ind_atom(iat)) )
      if (ind_atom(iat)==4) THEN !Only oxygen atoms
        CALL SUM_WATER_HYDROGENS( iat, dipole(:,:,iat) )
      end if
    end if

  END DO

END SUBROUTINE ASSIGN_DIPOLES

SUBROUTINE SUM_WATER_HYDROGENS(iat,dipole)
  USE parameters,ONLY : pos, ind_atom, natoms
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: iat
  INTEGER                    :: iH, nH
  DOUBLE PRECISION           :: dipole(3,4), d(4)
  DOUBLE PRECISION, PARAMETER :: cutoff_OH=2.26d0 !in bohr

  nH=0
  DO iH=1,natoms

    IF (ind_atom(iH)==2) THEN

      CALL DISTANCE_VECTOR( pos(iat,:), pos(iH,:), d )

      IF ( d(4) < cutoff_OH ) THEN
        dipole(:,1) = dipole(:,1) + d(1:3)
        nH=nH+1
      END IF

    END IF

  END DO

  if (nH==0) THEN
    print *, 'Could not find H for Oxygen ', iat
  end if

END SUBROUTINE SUM_WATER_HYDROGENS

SUBROUTINE SUM_NEAREST_WFCs( iat, dipole, nwf )
  !Sum the coordinates of the WFCs closest to atom with index "iat"
  USE histogram, ONLY : nwfc, wfc, wfc_field, is_dipole
  USE parameters,ONLY : pos, ind_atom
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: iat, nwf
  DOUBLE PRECISION           :: dipole(3,4)
  INTEGER                    :: iwf, ifield, count_wf(4)
  DOUBLE PRECISION           :: d(4) 
  DOUBLE PRECISION           :: cutoff_wfc !in bohr

  cutoff_wfc = 1.2d0 !1.6d0
  count_wf=0

  DO iwf = 1,nwfc

    CALL DISTANCE_VECTOR( pos(iat,:), wfc(iwf,:), d )

    IF ( d(4) < cutoff_wfc ) THEN
      dipole(:,1) = dipole(:,1) - 2.d0*d(1:3)
      count_wf(1) = count_wf(1) + 1
      CALL PRINT_WFC_DISTANCE( ind_atom(iat), pos(iat,3), d(4) )
    END IF
 
    DO ifield = 1,3

      CALL DISTANCE_VECTOR( pos(iat,:), wfc_field(iwf,ifield,:), d )

      IF ( d(4) < cutoff_wfc ) THEN
        dipole(:,ifield+1) = dipole(:,ifield+1) - 2.d0*d(1:3)
        count_wf(ifield+1) = count_wf(ifield+1) + 1
      END IF

    END DO

  END DO

  IF ( ANY(count_wf/=nwf) ) THEN
    PRINT *, "!!!WARNING!!!"
    PRINT *, "Found : ", count_wf, " WFCs. Expecting ", nwf
    is_dipole(iat)=.false.
    !STOP
  END IF

END SUBROUTINE SUM_NEAREST_WFCs

SUBROUTINE PRINT_WFC_DISTANCE( indat, d, z )
  !Output distace d of Wannier center from atom iat
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: d, z
  INTEGER, INTENT(IN) :: indat

  WRITE(8,fmt='(I2, 2(F12.8,3X))')  indat, z, d

END SUBROUTINE PRINT_WFC_DISTANCE

SUBROUTINE MOLECULAR_DIPOLE_POLARIZABILITY
  !Convert the derivatives from the lab frame to the bond frame
  USE parameters, ONLY : natoms
  USE histogram , ONLY : include_induced_field, is_dipole
  IMPLICIT NONE
  INTEGER                    :: frame, i
  DOUBLE PRECISION           :: induced_field(3,3)

  DO i = 1,natoms

    IF ( is_dipole(i) ) THEN !We also include bridging and terminal OH groups

      CALL PRINT_DIPOLE(i) 

      IF (include_induced_field) THEN
        CALL INDUCED_FIELD_EWALD( i, induced_field )
      ELSE
        induced_field = 0.
      END IF

      CALL PRINT_POLARIZABILITY_AND_BEC( i, induced_field )
   
    END IF

  ENDDO
 
END SUBROUTINE MOLECULAR_DIPOLE_POLARIZABILITY

SUBROUTINE PRINT_POLARIZABILITY_AND_BEC( iOx, induced_field )
  !Obtain water polarizability and Born effective charge tensors 
  USE histogram, ONLY : bec, efield, dipole, compute_bec
  USE parameters, ONLY : pos, box
  IMPLICIT NONE
  INTEGER, INTENT(IN)          :: iOx
  DOUBLE PRECISION, INTENT(IN) :: induced_field(3,3)
  DOUBLE PRECISION             :: alpha(3,3)
  DOUBLE PRECISION             :: bec_tmp(3,3)
  INTEGER                      :: iat, iwfc, ipol1, ipol2

  alpha = 0

  DO ipol1=1,3
    DO ipol2=1,3
      alpha(ipol1,ipol2) = dipole(ipol2,ipol1+1,iOx) &
                               - dipole(ipol2,1,iOx)
      alpha(ipol1,ipol2) = alpha(ipol1,ipol2) / &
               ( efield(ipol1) + induced_field(ipol1,ipol2) )
    END DO
  END DO

  write(50, fmt='(9(3X,F8.3))') &
    ( ( alpha(ipol1,ipol2) * (0.529177249**3), ipol2=1,3), ipol1=1,3 )
  
END SUBROUTINE PRINT_POLARIZABILITY_AND_BEC

SUBROUTINE PRINT_DIPOLE( iwat )
  USE parameters, ONLY : ind_atom, pos
  USE histogram, ONLY : dipole
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: iwat
  INTEGER                    :: ind_H(2), n_hydrogens
  DOUBLE PRECISION           :: rotmat(3,3), oh(3,2), dip_orient

  IF (ind_atom(iwat)==4) THEN

    CALL get_H2(iwat, ind_H, n_hydrogens)
    CALL OH_VECTOR( iwat, ind_H, oh )
    CALL ROT_MAT( oh, rotmat ) 
    dip_orient = dipole(3,1,iwat)/norm2(dipole(:,1,iwat))
    dipole(:,1,iwat) = matmul(rotmat,dipole(:,1,iwat))

  END IF

  !Dipole in debye
  WRITE( 5, fmt="(i2, 5(F12.7,3X))" ) ind_atom(iwat), pos(iwat,3), dipole(:,1,iwat)*0.529177249/0.2081943, dip_orient

END SUBROUTINE PRINT_DIPOLE
 
SUBROUTINE The_End
  USE histogram, ONLY : compute_bec

  CLOSE(1);CLOSE(2);CLOSE(3)
  CLOSE(5);CLOSE(6);CLOSE(7);CLOSE(8)
  CLOSE(50)
  
  IF (compute_bec) CLOSE(4)

END SUBROUTINE The_End

SUBROUTINE INDUCED_FIELD_EWALD ( ind_at, induced_field )
  !Calculate the re-irradiated field produced by the change in 
  !polarization after we apply the external electric field
  USE histogram, ONLY : dipole, is_dipole
  USE parameters, ONLY : natoms, box, pos
  IMPLICIT NONE
  INTEGER, INTENT(IN)          :: ind_at
  DOUBLE PRECISION             :: induced_field(3,3)
  DOUBLE PRECISION             :: LR_field(3,3)
  DOUBLE PRECISION             :: SR_field(3,3)
  DOUBLE PRECISION             :: sum_delta_dip(3,3)
  DOUBLE PRECISION             :: delta_dip(3,3,natoms)
  DOUBLE PRECISION             :: boxrec(3), Vcell 
  DOUBLE PRECISION             :: Kvec(3), K2 
  DOUBLE PRECISION             :: r(4,natoms), r2, fact
  DOUBLE PRECISION             :: rLvec(3), rL,rL2
  DOUBLE PRECISION             :: f1, f2
  DOUBLE PRECISION, PARAMETER  :: Pi=3.14159265359,alpha=.1
  DOUBLE PRECISION             :: eps(3,3)
  COMPLEX, PARAMETER           :: i = (0,1.d0)
  INTEGER                      :: iat, ipol, ieps
  INTEGER                      :: ik1, ik2, ik3
  INTEGER,PARAMETER            :: nk=8, nL=0

  Vcell = 1.d0
  DO ipol=1,3
    boxrec(ipol) =  2*Pi/box(ipol,ipol)
    Vcell = Vcell * box(ipol,ipol)
  END DO

  induced_field = 0.d0
  LR_field = 0.d0
  SR_field = 0.d0

  sum_delta_dip = 0.d0
  DO iat=1,natoms
    CALL DISTANCE_VECTOR( pos(iat,:), pos(ind_at,:), r(:,iat) )
    DO ipol=1,3
      delta_dip(:,ipol,iat) = dipole(:,ipol+1,iat) - dipole(:,1,iat)
    END DO
    sum_delta_dip = sum_delta_dip + delta_dip(:,:,iat)
  END DO

  !eps = 0.d0

  !!Fourier part of Ewald sum
  !DO ik1=-nk,nk 
  !  DO ik2=-nk,nk
  !    DO ik3=-nk,nk
  !      Kvec = (/ ik1,ik2,ik3 /) * boxrec
  !      K2 = sum(Kvec*Kvec)
  !      IF (K2==0) CYCLE

  !      f1 = 4.d0*(Pi/K2) * exp(-K2/(4*alpha)) 
  !      DO iat=1,natoms
  !        f2 = f1 * REAL(exp(i*sum(Kvec*r(1:3,iat)))) 

  !        DO ieps = 1,3
  !          eps(:,ieps) = eps(:,ieps) - &
  !            f2 * Kvec * sum(Kvec*delta_dip(:,ieps,iat)) 
  !        END DO

  !      END DO

  !      DO ieps = 1,3
  !        eps(:,ieps) = eps(:,ieps) + &
  !          f1 * Kvec * sum(Kvec*delta_dip(:,ieps,ind_at))
  !      END DO

  !    END DO
  !  END DO
  !END DO
  !
  !!K=0 term 
  !!eps = eps - (4.d0*Pi/3.d0) * ( sum_delta_dip - delta_dip(:,:,ind_at) ) 

  !LR_field = eps/Vcell
  !
  !!Short range part of Ewald sum
  !DO ik1=-nL,nL
  !  DO ik2=-nL,nL
  !    DO ik3=-nL,nL

  !      DO iat=1,natoms

  !        IF ( (iat==ind_at) .or. ( .not. is_dipole(iat)) ) CYCLE

  !        rLvec = r(1:3,iat) + MATMUL((/ ik1, ik2, ik3 /),box)
  !        rL2 = sum(rLvec*rLvec)
  !        IF ( rL2 .lt. 1.D-8 ) CYCLE
  !        rL = sqrt(rL2)

  !        !f1 = 3.d0*erfc(sqrt(alpha)*rL)/(rL2*rL) &
  !        !   + 6.d0*sqrt(alpha/Pi) * Exp(-alpha*rL2)/rL2 &
  !        !   + 4.d0*alpha*sqrt(alpha/Pi) * Exp(-alpha*rL2)
  !        !f1 = f1 / rL2
  !        !     
  !        !f2 = erfc(sqrt(alpha)*rL)/rL &
  !        !   + 2.d0*sqrt(alpha/Pi)*exp(-alpha/rL2)
  !        !f2 = f2 / rL2

  !        f1 = 2.d0*sqrt(alpha/Pi) * rL * exp( -alpha*rL2 ) &
  !           + erfc(sqrt(alpha)*rL)
  !        f2 = 4.d0*alpha*sqrt(alpha/Pi) * exp( -alpha*rL2 ) / rL2

  !        DO ieps = 1,3
  !          !SR_field(:,ieps) = SR_field(:,ieps) &
  !          !                 - f1 * rLvec * sum( rLvec*dipole(:,ieps,iat) ) & 
  !          !                 + f2 * dipole(:,ieps,iat) 
  !          SR_field(:,ieps) = SR_field(:,ieps) & 
  !                           + f1*3.d0*rLvec * sum(rLvec*delta_dip(:,ieps,iat)) & 
  !                           / (rL2*rL2*rL) &
  !                           - f1*delta_dip(:,ieps,iat) / (rL2*rL) &
  !                           + f2* rLvec * sum( rLvec*delta_dip(:,ieps,iat) )
  !       
  !        END DO

  !      END DO

  !    END DO
  !  END DO
  !END DO

  !DO ipol=1,3
  !  DO ieps = 1,3
  !    induced_field(ieps,ipol) = LR_field(ipol,ieps) + SR_field(ipol,ieps)
  !  END DO
  !END DO

  induced_field = 0.0
  induced_field(3,3) = - 4*Pi*sum_delta_dip(3,3) / box(1,1) / box(2,2) / box(3,3)
  !print *, LR_field(:,1:3)
  !print *, SR_field(:,1:3)
  !print *, "########################"
  !print *, LR_field(:,1:3) + SR_field(:,:)
  !print *, 4*Pi*sum_delta_dip(3,3) / box(1,1) / box(2,2) / box(3,3)
  !STOP
    
END SUBROUTINE INDUCED_FIELD_EWALD

SUBROUTINE OH_VECTOR(O, H, OH)
  ! OH vectors of a water molecule 
  USE parameters, ONLY : box, pos
  IMPLICIT NONE
  INTEGER, INTENT(IN)                  :: O, H(2)
  DOUBLE PRECISION                     :: OH(3,2)
  INTEGER                              :: ih, ipol

  OH = 0

  DO ih = 1,2

    IF ( H(ih) .ne. 0 ) THEN

      DO ipol=1,3

        OH(ipol,ih) = pos(H(ih),ipol) - pos(O,ipol)
        OH(ipol,ih) = OH(ipol,ih) - nint( OH(ipol,ih)/box(ipol,ipol) ) * box(ipol,ipol)

      END DO

      !OH(ih,:) = OH(ih,:) / norm2(OH(ih,:))
   
    ELSE

      OH(:,ih) = (/1.,0.,0./) !Arbitrary
      PRINT *, 'WARNING: THERE ARE DISSOCIATED WATER HERE'
 
    END IF

  END DO

END SUBROUTINE OH_VECTOR

SUBROUTINE get_H2(ind_O, ind_H, nhydrogens)
  USE parameters, ONLY : natoms, ind_atom
  IMPLICIT NONE
  INTEGER, DIMENSION(2)                :: ind_H
  INTEGER                              :: i, k, ind_O
  INTEGER                              :: nhydrogens
  DOUBLE PRECISION                     :: Dist, dist_OH, mindist(2)
  DOUBLE PRECISION, PARAMETER          :: cutoff_OH = 3.5 !2.9

  mindist = 100.d0
  ind_H = 0

  DO i = 1,natoms

    !Should be the index for H
    IF ( ind_atom(i) == 2 ) THEN

      dist_OH = Dist(ind_O,i)

      !Select only H's within a threshold
      IF ( dist_OH < cutoff_OH ) THEN
 
        !Very basic sorting of all OH distances within the threshold
        IF ( dist_OH < mindist(1) ) THEN 

          ind_H(2) = ind_H(1)
          ind_H(1) = i
          mindist(2) = mindist(1)
          mindist(1) = dist_OH

        ELSEIF ( dist_OH < mindist(2) ) THEN
 
          ind_H(2) = i
          mindist(2) = dist_OH

        END IF
 
      END IF
 
    END IF

  END DO

  !For consitency: higher index comes first
  IF ( ind_H(1) < ind_H(2) ) THEN
    i = ind_H(2)
    ind_H(2) = ind_H(1)
    ind_H(1) = i
  END IF

  nhydrogens=2
  DO i=1,2
    IF (ind_H(i) == 0) nhydrogens = nhydrogens -1
  END DO

END SUBROUTINE get_H2

SUBROUTINE ROT_MAT( OH, rotmat )
  !Rotation matrix from lab to molecular frame
  !Assuming OH is normalized
  !In the molecular frame: z is along OHs bissector; x is in the H2O plane
  !y is vector poiting out-of-plane
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(INOUT)      :: rotmat(3,3) 
  DOUBLE PRECISION, INTENT(IN)         :: OH(3,2) 
  DOUBLE PRECISION                     :: OH1(3), OH2(3) 
  DOUBLE PRECISION                     :: x(3), y(3), z(3)
  INTEGER                              :: i

  OH1 = OH(:,1)
  OH2 = OH(:,2)

  !Get y
  CALL CROSSPROD(OH1,OH2,y) 
  y = y/norm2(y)

  !Get z
  z = ( OH1 + OH2 ) / norm2( OH1 + OH2 ) 
  
  !Get x
  CALL CROSSPROD(z,y,x)
  x = x/norm2(x)

  DO i = 1,3
    rotmat(1,i)=x(i)
    rotmat(2,i)=y(i)
    rotmat(3,i)=z(i)
  END DO

END SUBROUTINE ROT_MAT 
