!Set of subroutines used to analyse the lipid trajectory

MODULE wannier
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: forces(:,:,:) 
  DOUBLE PRECISION, ALLOCATABLE            :: wfc(:,:) 
  DOUBLE PRECISION, ALLOCATABLE            :: wfc_field(:,:,:) 
  DOUBLE PRECISION, ALLOCATABLE            :: dipole(:,:,:) 
  DOUBLE PRECISION                         :: efield(3) 
  LOGICAL, ALLOCATABLE                     :: is_dipole(:)
  INTEGER                                  :: nwfc
END MODULE wannier

MODULE slab
  IMPLICIT NONE
  LOGICAL,ALLOCATABLE                      :: is_titanium_5c(:)
END MODULE slab

PROGRAM SFG
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE

  DO frame = 1,nframes
    CALL READ_POS_FORCES
    CALL READ_WFC
    !CALL CHECK_SUM_RULE
    CALL ASSIGN_DIPOLES
    CALL BOND_FRAME_DERIVATIVE (frame)
  END DO

  CLOSE(1);CLOSE(2);CLOSE(3)

END PROGRAM SFG

SUBROUTINE INITIALIZE
  USE wannier,  ONLY : forces, wfc, wfc_field, dipole, &
                       nwfc, is_dipole
  USE parameters, ONLY : natoms, nframes, nwater, &
                         pos, coarse, nlayers
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, index_file
  CHARACTER(100)             :: forces_file, wfc_file 
  CHARACTER(100)             :: wfc_field_file 

  CALL getarg(1, pos_file)
  CALL getarg(2, forces_file)
  CALL getarg(3, wfc_file)
  CALL getarg(4, index_file)

  CALL READ_INDEX_WFC ( index_file )

  ALLOCATE(pos(natoms,3)); pos=0 ! coordinates
  ALLOCATE(forces(natoms,7,3)); forces=0
  ALLOCATE(wfc(nwfc,3)); wfc=0
  ALLOCATE(wfc_field(nwfc,3,3)); wfc_field=0
  ALLOCATE(dipole(3,4,natoms)); dipole=0.d0
  ALLOCATE(is_dipole(natoms))
  
  OPEN(unit = 1,file = pos_file)
  OPEN(unit = 2,file = forces_file)
  OPEN(unit = 3,file = wfc_file)
  
END SUBROUTINE INITIALIZE

SUBROUTINE READ_POS_FORCES
  USE wannier, ONLY : forces 
  USE parameters,ONLY : natoms, pos
  IMPLICIT NONE
  INTEGER                    :: iat, ipol, ifield

  READ(1,*)
  READ(2,*)

  DO iat = 1, natoms
    READ(1, *), ( pos(iat,ipol), ipol=1,3 )
    !Each line: nofield, field x, field -x, field y, field -y, field z, field -z
    READ(2, *), ( ( forces(iat,ifield,ipol), ipol=1,3 ), ifield=1,7 )
  END DO

END SUBROUTINE READ_POS_FORCES

SUBROUTINE READ_WFC
  USE wannier, ONLY : nwfc, wfc, wfc_field 
  IMPLICIT NONE
  INTEGER                    :: iwfc, ipol, ifield

  READ(3, *)

  DO iwfc=1,nwfc
    !Each line: nofield, field x, field -x, field y, field -y, field z, field -z
    READ(3, *), ( wfc(iwfc,ipol), ipol=1,3 ), &
                ( ( wfc_field(iwfc,ifield,ipol), ipol=1,3 ), ifield=1,3 )
  END DO
    
END SUBROUTINE READ_WFC

SUBROUTINE BOND_FRAME_DERIVATIVE (frame)
  !Convert the derivatives from the lab frame to the bond frame
  USE parameters, ONLY : natoms, coarse
  IMPLICIT NONE
  INTEGER                    :: frame, i, ih
  INTEGER                    :: ind_O, ind_H(2)
  INTEGER                    :: n_hydrogens
  INTEGER                    :: ind_layer_pf, j
  LOGICAL                    :: is_oxygen
  DOUBLE PRECISION           :: oh(2,3)
  DOUBLE PRECISION           :: der1(3), der2(3)
  DOUBLE PRECISION           :: rotmat(3,3)
  DOUBLE PRECISION           :: induced_field(3,3)

  !Water oxygen atom index
  ind_O = 0

  DO i = 1,natoms

    IF ( is_oxygen(i) ) THEN

      ind_O = ind_O + 1
      ind_H = 0
      CALL get_H2(i, ind_H, n_hydrogens)
  
      IF ( n_hydrogens .eq. 0 ) CYCLE

      CALL OH_VECT( i, ind_H, oh )

      CALL INDUCED_FIELD_EWALD( i, induced_field )

      !Two OH vectors, taking PBC into account 
      DO ih = 1,n_hydrogens

        CALL ROT_MAT( oh, ih, rotmat )

        CALL dF_dEps( i, ind_H(ih), rotmat, induced_field, der1 )

        CALL d2F_dEps2( i, ind_H(ih), rotmat, induced_field, der2 )

        WRITE(*,fmt="(i2,7(3X,E12.5))"), ind_layer_pf(i,n_hydrogens), der1, der2, oh(ih,2)
         
      END DO

    END IF

  ENDDO

END SUBROUTINE BOND_FRAME_DERIVATIVE

SUBROUTINE OH_VECT(O, H, OH)
  ! Normalized OH vectors of a water molecule 
  USE parameters, ONLY : box, pos
  IMPLICIT NONE
  INTEGER, INTENT(IN)                  :: O, H(2)
  DOUBLE PRECISION                     :: OH(2,3), unit_vector(3)
  INTEGER                              :: ih, ipol

  IF ( H(1) .eq. 0 ) THEN
    print *, 'ERROR: First index is zero ', O
    STOP
  END IF

  DO ih = 1,2

    IF ( H(ih) .ne. 0 ) THEN

      DO ipol=1,3

        OH(ih,ipol) = pos(H(ih),ipol) - pos(O,ipol)
        OH(ih,ipol) = OH(ih,ipol) - nint( OH(ih,ipol)/box(ipol,ipol) ) * box(ipol,ipol)

      END DO

    ELSE

      CALL NEAREST_TI(O, unit_vector) 

      OH(ih,:) = unit_vector 
 
    END IF

  END DO

END SUBROUTINE oh_vect

SUBROUTINE NEAREST_TI( iO, d )
  !Find Ti nearest to O with index iO
  USE parameters, ONLY : pos, box, natoms
  INTEGER, INTENT(IN)        :: iO
  DOUBLE PRECISION           :: Dist, d(3), md
  DOUBLE PRECISION           :: mindist 
  INTEGER                    :: iTi, iat, ipol
  LOGICAL                    :: is_titanium

  mindist=1000

  DO iat = 1,natoms

    IF ( is_titanium(iat) ) THEN
  
      md = Dist(iO, iat)

      IF ( md < mindist ) THEN
        iTi = iat
        mindist=md
      END IF

    END IF

  END DO

  DO ipol=1,3
    d(ipol) = pos(iTi,ipol) - pos(iO,ipol)
    d(ipol) = d(ipol) - nint(d(ipol)/box(ipol,ipol))*box(ipol,ipol)
  END DO

  d = d/mindist

END SUBROUTINE NEAREST_TI

SUBROUTINE ROT_MAT( OH, ih, rotmat )
  !Rotation matrix from lab to bond frame
  !Assuming OH1 is normalized
  !In the bond frame: z is along OH; x is in the H2O plane,
  !poiting to nearest OH; y is vector poiting out-of-plane
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(INOUT)      :: rotmat(3,3) 
  DOUBLE PRECISION, INTENT(IN)         :: OH(2,3) 
  DOUBLE PRECISION                     :: OH1(3), OH2(3)
  DOUBLE PRECISION                     :: x(3), y(3)
  INTEGER, INTENT(IN)                  :: ih
  INTEGER                              :: i

  OH1 = OH(ih,:) / norm2(OH(ih,:))
  OH2 = OH(3-ih,:) / norm2(OH(3-ih,:))

  !Get x
  CALL CROSSPROD(OH1,OH2,x) 
  x = x/norm2(x)
  
  !Get y
  CALL CROSSPROD(x,OH1,y)
  y = y/norm2(y)

  rotmat=0
  DO i = 1,3
    rotmat(1,i)=x(i)
    rotmat(2,i)=y(i)
    rotmat(3,i)=OH1(i)
  END DO

END SUBROUTINE ROT_MAT 

SUBROUTINE dF_dEps( iO, iH, rotmat, induced_field, der1 )
  !Numerical derivative of force w.r.t. a uniform electric field
  USE wannier, ONLY : forces, efield
  IMPLICIT NONE
  INTEGER, INTENT(IN)                  :: iO, iH
  DOUBLE PRECISION, INTENT(IN)         :: rotmat(3,3)
  DOUBLE PRECISION, INTENT(IN)         :: induced_field(3,3)
  DOUBLE PRECISION, INTENT(INOUT)      :: der1(3)
  DOUBLE PRECISION                     :: dF(3,3)
  INTEGER                              :: ipol, ipol2, ifield
  
  DO ifield = 1,3
    DO ipol = 1,3
      dF(ifield,ipol) = forces(iH,2*ifield,ipol) - forces(iH,1,ipol) &
                      - forces(iO,2*ifield,ipol) + forces(iO,1,ipol)
    END DO
    !Only the component of the force along the OH bond is considered
    dF(ifield,1) = sum( rotmat(3,:) * dF(ifield,:) ) / & 
                   ( efield(ifield) + induced_field(ifield,ifield) )
  END DO

  !Convert derivative from lab frame to the bond frame
  DO ifield = 1,3
    der1(ifield) = sum( rotmat(ifield,:) * dF(:,1) )
  END DO

END SUBROUTINE dF_dEps

SUBROUTINE d2F_dEps2( iO, iH, rotmat, induced_field, der2 )
  !Numerical derivative of force w.r.t. a uniform electric field
  !Considering only the diagonal terms of the polarizability tensor
  USE wannier, ONLY : forces, efield
  IMPLICIT NONE
  INTEGER, INTENT(IN)                  :: iO, iH
  DOUBLE PRECISION, INTENT(IN)         :: rotmat(3,3)
  DOUBLE PRECISION, INTENT(IN)         :: induced_field(3,3)
  DOUBLE PRECISION, INTENT(INOUT)      :: der2(3)
  DOUBLE PRECISION                     :: dF(3,3)
  INTEGER                              :: ipol, ifield
  
  DO ifield = 1,3
    DO ipol = 1,3
      dF(ifield,ipol) = forces(iH,2*ifield,ipol) + forces(iH,2*ifield+1,ipol) &
                      - 2*forces(iH,1,ipol) &
                      - forces(iO,2*ifield,ipol) - forces(iO,2*ifield+1,ipol) &
                      + 2*forces(iO,1,ipol)
    END DO
    !Only the component of the force along the OH bond is considered
    dF(ifield,1) = sum( rotmat(3,:) * dF(ifield,:) ) / &
                  ( efield(ifield) + induced_field(ifield,ifield) )**2
  END DO

  !Convert derivative from lab frame to the bond frame
  DO ifield = 1,3
    der2(ifield) = sum( rotmat(ifield,:) * dF(:,1) * rotmat(ifield,:) )
  END DO

END SUBROUTINE d2F_dEps2

SUBROUTINE CROSSPROD( x, y, prod )
  !Cross product of 3D vectors x and y
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(3)       :: prod
  DOUBLE PRECISION, INTENT(IN)         :: x(3),y(3)
  
  prod(1) = x(2)*y(3) - x(3)*y(2)
  prod(2) = x(3)*y(1) - x(1)*y(3)
  prod(3) = x(1)*y(2) - x(2)*y(1)

END SUBROUTINE CROSSPROD

SUBROUTINE ASSIGN_DIPOLES
  !Transform a system of point charges into a system of dipoles
  USE parameters, ONLY : natoms
  USE wannier,  ONLY : dipole, is_dipole
  IMPLICIT NONE
  INTEGER                    :: iat, ih, ifield
  INTEGER                    :: ind_H(2), n_hydrogens
  LOGICAL                    :: is_titanium, is_oxygen
  DOUBLE PRECISION           :: oh(3,2), cutoff

  dipole = 0.d0
  is_dipole = .False.

  DO iat=1,natoms

    IF ( is_oxygen(iat) ) THEN 

      CALL get_H2(iat, ind_H, n_hydrogens)

      IF ( n_hydrogens.ne.0 ) THEN 

        CALL OH_VECT( iat, ind_H, oh )

        DO ih=1,n_hydrogens
          DO ifield=1,4
            dipole(:,ifield,iat) = dipole(:,ifield,iat) + oh(:,ih)
          END DO
        END DO

      END IF

      CALL SUM_NEAREST_WFCs( iat, dipole(:,:,iat), 4 )

      is_dipole(iat) = .True.

    ELSE IF ( is_titanium(iat) ) THEN
  
      CALL SUM_NEAREST_WFCs( iat, dipole(:,:,iat), 4 )
     
      is_dipole(iat) = .True.

    END IF

  END DO

END SUBROUTINE ASSIGN_DIPOLES

SUBROUTINE SUM_NEAREST_WFCs( iat, dipole, nwf )
  !Sum the coordinates of the WFCs closest to atom with index "iat"
  USE wannier, ONLY : nwfc, wfc, wfc_field
  USE parameters,ONLY : pos, ind_atom, box
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: iat, nwf
  DOUBLE PRECISION           :: dipole(3,4)
  INTEGER                    :: iwf, ifield, count_wf(4)
  DOUBLE PRECISION           :: d(4) 
  DOUBLE PRECISION           :: cutoff_wfc !in bohr

  cutoff_wfc = 1.6d0
  count_wf=0
  IF ( nwf .eq. 4 ) cutoff_wfc = 1.6d0

  DO iwf = 1,nwfc

    CALL DISTANCE_VECTOR( pos(iat,:), wfc(iwf,:), d )

    IF ( d(4) < cutoff_wfc ) THEN
      dipole(:,1) = dipole(:,1) - 2.d0*d(1:3)
      count_wf(1) = count_wf(1) + 1
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
    !STOP
  END IF

END SUBROUTINE SUM_NEAREST_WFCs

SUBROUTINE INDUCED_FIELD_EWALD ( ind_at, induced_field )
  !Calculate the re-irradiated field produced by the change in 
  !polarization after we apply the external electric field
  USE wannier, ONLY : dipole, is_dipole
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

  eps = 0.d0

  !Fourier part of Ewald sum
  DO ik1=-nk,nk 
    DO ik2=-nk,nk
      DO ik3=-nk,nk
        Kvec = (/ ik1,ik2,ik3 /) * boxrec
        K2 = sum(Kvec*Kvec)
        IF (K2==0) CYCLE

        f1 = 4.d0*(Pi/K2) * exp(-K2/(4*alpha)) 
        DO iat=1,natoms
          f2 = f1 * REAL(exp(i*sum(Kvec*r(1:3,iat)))) 

          DO ieps = 1,3
            eps(:,ieps) = eps(:,ieps) - &
              f2 * Kvec * sum(Kvec*delta_dip(:,ieps,iat)) 
          END DO

        END DO

        DO ieps = 1,3
          eps(:,ieps) = eps(:,ieps) + &
            f1 * Kvec * sum(Kvec*delta_dip(:,ieps,ind_at))
        END DO

      END DO
    END DO
  END DO
  
  !K=0 term 
  !eps = eps - (4.d0*Pi/3.d0) * ( sum_delta_dip - delta_dip(:,:,ind_at) ) 

  LR_field = eps/Vcell
  
  !Short range part of Ewald sum
  DO ik1=-nL,nL
    DO ik2=-nL,nL
      DO ik3=-nL,nL

        DO iat=1,natoms

          IF ( (iat==ind_at) .or. ( .not. is_dipole(iat)) ) CYCLE

          rLvec = r(1:3,iat) + MATMUL((/ ik1, ik2, ik3 /),box)
          rL2 = sum(rLvec*rLvec)
          IF ( rL2 .lt. 1.D-8 ) CYCLE
          rL = sqrt(rL2)

          !f1 = 3.d0*erfc(sqrt(alpha)*rL)/(rL2*rL) &
          !   + 6.d0*sqrt(alpha/Pi) * Exp(-alpha*rL2)/rL2 &
          !   + 4.d0*alpha*sqrt(alpha/Pi) * Exp(-alpha*rL2)
          !f1 = f1 / rL2
          !     
          !f2 = erfc(sqrt(alpha)*rL)/rL &
          !   + 2.d0*sqrt(alpha/Pi)*exp(-alpha/rL2)
          !f2 = f2 / rL2

          f1 = 2.d0*sqrt(alpha/Pi) * rL * exp( -alpha*rL2 ) &
             + erfc(sqrt(alpha)*rL)
          f2 = 4.d0*alpha*sqrt(alpha/Pi) * exp( -alpha*rL2 ) / rL2

          DO ieps = 1,3
            !SR_field(:,ieps) = SR_field(:,ieps) &
            !                 - f1 * rLvec * sum( rLvec*dipole(:,ieps,iat) ) & 
            !                 + f2 * dipole(:,ieps,iat) 
            SR_field(:,ieps) = SR_field(:,ieps) & 
                             + f1*3.d0*rLvec * sum(rLvec*delta_dip(:,ieps,iat)) & 
                             / (rL2*rL2*rL) &
                             - f1*delta_dip(:,ieps,iat) / (rL2*rL) &
                             + f2* rLvec * sum( rLvec*delta_dip(:,ieps,iat) )
         
          END DO

        END DO

      END DO
    END DO
  END DO

  DO ipol=1,3
    DO ieps = 1,3
      induced_field(ieps,ipol) = LR_field(ipol,ieps) + SR_field(ipol,ieps)
    END DO
  END DO

  !print *, LR_field(:,1:3)
  !print *, SR_field(:,1:3)
  !print *, "########################"
  !print *, LR_field(:,1:3) + SR_field(:,:)
  !STOP
    
END SUBROUTINE INDUCED_FIELD_EWALD

SUBROUTINE READ_INDEX_WFC ( index_file ) 
  USE wannier,  ONLY : nwfc, efield
  USE parameters, ONLY : nlayers, nframes, nequil,&
                         natoms, layers, box, ind_atom
  IMPLICIT NONE
  CHARACTER(100), INTENT(IN) :: index_file
  INTEGER                    :: i,j
  
  OPEN(2, file=index_file)

  READ (2, *) natoms, nframes, nequil, nwfc, nlayers
  READ (2, *) ( ( box(j,i), i=1,3 ), j=1,3 )
  READ (2, *) efield 

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

SUBROUTINE get_H2(ind_O, ind_H, nhydrogens)
  USE parameters, ONLY : natoms, ind_atom
  IMPLICIT NONE
  INTEGER, DIMENSION(2)                :: ind_H
  INTEGER                              :: i, k, ind_O
  INTEGER                              :: nhydrogens
  DOUBLE PRECISION                     :: Dist, dist_OH, mindist(2)
  DOUBLE PRECISION, PARAMETER          :: cutoff_OH = 2.456

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

SUBROUTINE FIND_TI5C
  USE slab, ONLY : is_titanium_5c
  USE parameters, ONLY: pos, natoms
  IMPLICIT NONE
  INTEGER                    :: iTi, iOx, ncoord
  DOUBLE PRECISION           :: d, Dist
  DOUBLE PRECISION, PARAMETER :: cutoff=4.4d0 ! in bohr
  LOGICAL                    :: is_titanium
  LOGICAL                    :: is_TiO2_oxygen
  
  ALLOCATE(is_titanium_5c(natoms))
  is_titanium_5c=.false.

  DO iTi = 1, natoms

    IF ( is_titanium(iTi) ) THEN
  
      ncoord = 0
      DO iOx = 1, natoms 

        IF ( is_TiO2_oxygen(iOx) ) THEN

          d = Dist(iTi, iOx)
 
          IF ( d < cutoff ) THEN
   
            ncoord = ncoord + 1

          END IF

        END IF

      END DO

      IF ( ncoord .eq. 5 ) THEN

        is_titanium_5c(iTi) = .true.

      ELSE IF ( ncoord .ne. 6 ) THEN
      
        WRITE(*,*) 'Ti with wrong coordination number', iTi, ncoord
        WRITE(*,*) 'STOP!!!'
        STOP
   
      END IF

    END IF

  END DO

END SUBROUTINE FIND_TI5C

INTEGER FUNCTION ind_layer_pf( ind_O, n_hydrogens )
  !Layer index for prefactor. ind_O is the index of water oxygen in pos
  !This function also identifies terminal and bridging OH groups
  USE parameters, ONLY : pos, box, nlayers, layers
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: ind_O, n_hydrogens
  INTEGER                    :: il
  DOUBLE PRECISION           :: z
  LOGICAL                    :: is_water_oxygen

  IF ( is_water_oxygen( ind_O ) ) THEN 

    IF ( n_hydrogens.eq.2 ) THEN !Molecular water

      z = pos(ind_O,3)
      z = z - nint(z/box(3,3))*box(3,3) + box(3,3)/2.

      DO il = 1,nlayers
     
        if ( z <= layers(il+1) ) THEN
          ind_layer_pf = il
          RETURN
        endif
     
      END DO

      print *, "Could not find layer for ", ind_O, z

    ELSE !Terminal OH

      ind_layer_pf = nlayers + 1
      RETURN

    END IF

  ELSE !Bridging OH

    ind_layer_pf = nlayers + 2
    RETURN

  END IF

  ind_layer_pf = 0

END FUNCTION ind_layer_pf

!SUBROUTINE CHECK_SUM_RULE
!  USE wannier, ONLY : firstder, secder
!  USE parameters,ONLY : natoms
!  IMPLICIT NONE
!  INTEGER                    :: iat
!  DOUBLE PRECISION           :: sr1(3), sr2(3)
!
!  sr1 = 0
!  sr2 = 0
!
!  DO iat = 1, natoms
!  
!    sr1 = sr1 + firstder(iat,:)
!    sr2 = sr2 + secder(iat,:)
!
!  END DO
!
!  WRITE(51,fmt="(6(3X,E12.5))") sr1, sr2
!
!END SUBROUTINE CHECK_SUM_RULE

!SUBROUTINE PRINT_RESULTS 
!  USE wannier, ONLY : counter, dist_1, dist_2
!  USE parameters,ONLY : nlayers 
!  IMPLICIT NONE
!  INTEGER                    :: il
!  DOUBLE PRECISION           :: aver1, aver2
!
!  OPEN(unit = 50,file = "out_average")
!
!  DO il = 1, nlayers
!  
!    aver1 = dist_1(il) / dble(counter(il))
!    aver2 = dist_2(il) / dble(counter(il))
!
!    WRITE(50,fmt="(i1,2(3X,F10.5))") il, aver1, aver2 
!
!  END DO
!
!  CLOSE(50)
!
!END SUBROUTINE PRINT_RESULTS

