!Time correlation function of polarizabilty and polarization
!Both funcitons are expanded in terms of the OH vibrational 
!mode of water. Assuming complete separation between librational 
!and OH vibrational modes, the correlation function can be written
!as the velocity-velocity correlaiton function weighted by the 
!derivatives of bond polarization and polarizability with respect 
!to ROH. The user should parse the precomputed derivatives.


MODULE histogram
  IMPLICIT NONE 
  DOUBLE PRECISION, ALLOCATABLE            :: mu_z(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: alpha_xx(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: autocorr(:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: dMu_dR(:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: dAlpha_dR(:,:,:)
  DOUBLE PRECISION, ALLOCATABLE            :: layerspf(:)
  INTEGER, ALLOCATABLE                     :: ind_H_full(:,:,:) 
  INTEGER, ALLOCATABLE                     :: nhyd(:,:) 
  INTEGER                                  :: nlayerspf
END MODULE histogram

PROGRAM SFG
  USE parameters
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL_VEL

  DO frame = 1,nframes
    CALL READ_POS_VEL
    CALL SET_CENTER_OF_MASS_TO_ZERO
    CALL SET_CENTER_OF_MASS_VELOCITY_TO_ZERO
    CALL IDENTIFY_OH_GROUPS
    !CALL COARSE_GRAIN_POS (frame)
    CALL COMPUTE_MU_AND_ALPHA (frame)
  END DO

  CALL SMOOTH_COARSE_GRAIN
  CALL SSVAC 
  CALL PRINT_RESULTS

  CLOSE(1)

END PROGRAM SFG

SUBROUTINE INITIALIZE
  USE histogram,  ONLY : autocorr, mu_z, alpha_xx, &
                         dMu_dR, dAlpha_dR, ind_H_full, & 
                         nhyd
  USE parameters, ONLY : natoms, nframes, nwater, &
                         pos, vel, coarse, nlayers, &
                         nOxygen, OH_index
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: pos_file, vel_file, index_file,&
                                deriv_files

  CALL getarg(1, pos_file)
  CALL getarg(2, vel_file)
  CALL getarg(3, index_file)
  deriv_files="derivatives"

  CALL READ_INDEX ( index_file )

  ALLOCATE(pos(natoms,3)) ! coordinates
  ALLOCATE(vel(natoms,3)) ! velocities 
  ALLOCATE(coarse(nframes,nOxygen)); coarse=0
  ALLOCATE(OH_index(natoms))
  ALLOCATE(ind_H_full(nframes,nOxygen,2)); ind_H_full=0
  ALLOCATE(autocorr(nlayers,nframes)); autocorr = 0
  ALLOCATE(mu_z(2,nOxygen,nframes)); mu_z = 0
  ALLOCATE(alpha_xx(2,nOxygen,nframes)); alpha_xx = 0
  ALLOCATE(nhyd(nOxygen,nframes)); nhyd=0

  CALL READ_DERIVATIVES ( deriv_files )

  OPEN(unit = 1,file = pos_file)
  OPEN(unit = 2,file = vel_file)
  
END SUBROUTINE INITIALIZE

SUBROUTINE COMPUTE_MU_AND_ALPHA (frame)
  !Obtain OH dipole moment (mu) and polarizability (alpha) as a function of time
  !Using the precomputed derivatives with respect to OH distance
  USE histogram,  ONLY : dMu_dR, dAlpha_dR, mu_z, alpha_xx, ind_H_full, nhyd
  USE parameters, ONLY : natoms, vel, coarse, nwater
  IMPLICIT NONE
  INTEGER                    :: frame, i, ih, layer_pf, layer
  INTEGER                    :: ind_H(2)
  INTEGER                    :: ind_layer, ind_layer_pf
  INTEGER                    :: n_hydrogens, i_Ox, nOH
  INTEGER                    :: ihp 
  LOGICAL                    :: is_oxygen, switched_oh_index
  DOUBLE PRECISION           :: oh(2,3), vOH
  DOUBLE PRECISION           :: rotmat(3,3), alphatmp(3,3)

  !OH atom index
  nOH = 0
  i_Ox = 0 

  DO i = 1,natoms

    !Calculate the OH vector:
    IF ( is_oxygen(i) ) THEN

      i_Ox = i_Ox + 1
      CALL GET_IND_H(i, ind_H, n_hydrogens)

      nhyd(i_Ox,frame) = n_hydrogens

      IF ( n_hydrogens/=0 ) THEN 

        !CALL PRESERVE_OH_INDEX( i_Ox, ind_H, frame ) 

        layer_pf = ind_layer_pf( i, n_hydrogens ) !For prefactor 
        layer = ind_layer( i )

        CALL OH_VECT( i, ind_H, oh)

        coarse( frame, i_Ox ) = layer 
        ind_H_full( frame, i_Ox, :) = ind_H

        DO ih = 1,n_hydrogens

          nOH = nOH + 1

          CALL ROT_MAT( oh,ih, rotmat )
          !Velocity in the bond frame
          vOH = sum( ( vel(ind_H(ih),:) - vel(i,:) ) * rotmat(:,3) )

          IF ( switched_oh_index( i_Ox, ind_H, frame ) ) THEN
            ihp=3-ih
          ELSE
            ihp=ih
          END IF
            
          !Dipole and polarizability in the lab frame
          mu_z(ihp,i_Ox,frame) = sum( dMu_dR(:,layer_pf) * vOH * rotmat(3,:) )
          !Note that this now is alpha_xx
          alphatmp=MATMUL(rotmat,dAlpha_dR(:,:,layer_pf))
          alphatmp=MATMUL(alphatmp,TRANSPOSE(rotmat))
          alpha_xx(ihp,i_Ox,frame) = alphatmp(2,2) * vOH 

        END DO

        IF ( switched_oh_index( i_Ox, ind_H, frame ) ) THEN
          ind_H_full( frame, i_Ox, 1) = ind_H(2)
          ind_H_full( frame, i_Ox, 2) = ind_H(1)
        END IF

      END IF !n_hydrogens

    END IF !is_oxygen

  ENDDO

  IF ( nOH /= 2*nwater ) THEN
    print *, 'WARNING: Wrong number of OH groups ', nOH, ' on step ', frame
    !STOP
  END IF

END SUBROUTINE COMPUTE_MU_AND_ALPHA

SUBROUTINE SSVAC 
  !Surface specific velocity-velocity correlation function
  USE histogram,  ONLY : mu_z, alpha_xx, autocorr, nhyd 
  USE parameters, ONLY : nframes, nOxygen, coarse, nlayers
  IMPLICIT NONE
  DOUBLE PRECISION           :: sum_dip(nlayers), surv_prob(nlayers)
  DOUBLE PRECISION           :: avg_in_layer(nlayers)
  INTEGER                    :: frame1, frame2, iOx, ipol
  INTEGER                    :: counter_in(nlayers), counter_all(nlayers), ih 
  INTEGER                    :: layer, delta_t
  INTEGER                    :: n_in_layer(nlayers)
  INTEGER, PARAMETER         :: nsep=20
  LOGICAL                    :: in_layer, bound_to_same_H

  !frame1 is the duration of the interval 
  DO frame1 = 0,int((nframes-1)/5.)

    n_in_layer = 0
    surv_prob = 0
    delta_t = 0

    !Change initial condition, look for combinations with same interval
    DO frame2 = 1,nframes-frame1,nsep

      counter_in = 0
      counter_all = 0
      sum_dip = 0.D0
      delta_t = delta_t + 1

!$omp parallel do private(layer) reduction(+:counter_in, sum_dip, counter_all)
      DO iOx = 1,nOxygen 

        layer = coarse(frame2,iOx)
        IF ( ( layer .ne. 0 ) .and. ( layer .ne. 2 ) ) THEN
        !IF ( ( layer .ne. 0 ) .and. ( layer .ne. 4 ) ) THEN

          counter_all(layer) = counter_all(layer) + nhyd(iOx, frame2)

          !If the particle dont leave the window during the interval frame1
          IF ( in_layer(iOx,frame2,frame2+frame1,50) ) THEN

            counter_in(layer) = counter_in(layer) + nhyd(iOx, frame2)

            !Autocorrelation
            !IF ( bound_to_same_H(iOx,1,frame2,frame2+frame1,50) ) THEN
              sum_dip(layer) = sum_dip(layer) &
                             + mu_z(1,iOx,frame2) * alpha_xx(1,iOX,frame2+frame1) 
            !END IF
            !IF ( bound_to_same_H(iOx,2,frame2,frame2+frame1,50) ) THEN
              sum_dip(layer) = sum_dip(layer) &
                             + mu_z(2,iOx,frame2) * alpha_xx(2,iOX,frame2+frame1)
            !END IF

            !Intramolecular cross correlation
            !sum_dip(layer) = sum_dip(layer) &
            !               + mu_z(ih,iw,frame2) * alpha_xx(3-ih,iw,frame2+frame1)

          ENDIF

        ENDIF

      ENDDO
!$omp end parallel do

      DO layer = 1,nlayers
        IF ( counter_all(layer) /= 0 ) THEN
          surv_prob(layer) = surv_prob(layer) &
                           + dble(counter_in(layer)) / dble(counter_all(layer))
          !Averaged correlation function
          autocorr(layer,frame1+1) = autocorr(layer,frame1+1) + &
                                sum_dip(layer) / dble( counter_all(layer) ) 
        END IF
      END DO

      n_in_layer = n_in_layer + counter_in

    ENDDO

    surv_prob = surv_prob / dble( delta_t )
    avg_in_layer = dble(n_in_layer) / dble( delta_t )

    DO layer = 1, nlayers
      IF ( surv_prob(layer) /= 0 ) THEN
        !Averaged autocorrelation function
        autocorr(layer,frame1+1) = avg_in_layer(layer) * autocorr(layer,frame1+1) / dble( delta_t ) / surv_prob(layer) 
      END IF
    END DO

  ENDDO

END SUBROUTINE SSVAC 

SUBROUTINE READ_DERIVATIVES ( derivatives_file )
  !Read pre-computed derivatives of polarization and polarizability
  !w.r.t. the rOH vector
  USE histogram,  ONLY : dMu_dR, dAlpha_dR, layerspf, nlayerspf
  USE parameters, ONLY : nlayers
  CHARACTER(100), INTENT(IN) :: derivatives_file 
  INTEGER                    :: ipol1, ipol2, il

  OPEN(unit = 50,file = derivatives_file)

  READ(50,*) nlayerspf

  ALLOCATE(layerspf(nlayerspf+1)); layerspf = 0
  ALLOCATE(dMu_dR(3,nlayerspf+2)); dMu_dR = 0
  ALLOCATE(dAlpha_dR(3,3,nlayerspf+2)); dAlpha_dR = 0

  READ(50,*) layerspf

  DO il = 1, nlayerspf+2 !Last two lines are from terminal and bridging OH

    READ(50,*) ( dMu_dR(ipol1,il), ipol1=1,3 ),& 
               ( ( dAlpha_dR(ipol1,ipol2,il), ipol2=1,3 ), ipol1=1,3 )

  END DO
  
  CLOSE(50)

END SUBROUTINE READ_DERIVATIVES

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE histogram,  ONLY : autocorr
  USE parameters, ONLY : nlayers, nframes, iprint, dt
  IMPLICIT NONE
  INTEGER                    :: i, j
  DOUBLE PRECISION           :: n_const(nlayers)
  DOUBLE PRECISION,PARAMETER :: au_to_fs = 2.418884326509D-2
  DOUBLE PRECISION           :: tc 

  tc = au_to_fs*dt

  OPEN( unit = 2, file = "ssvac.dat" )

  DO i = 0,nframes-1

    write(2,fmt = '(E13.6,3X,*(E14.7,3X))') dble(i)*tc, &
         ( autocorr(j,i+1) , j=1,nlayers) 
         !( autocorr(j,i+1) / abs(autocorr(j,1)) , j=1,nlayers) 

  ENDDO

  CLOSE(2)

END SUBROUTINE PRINT_RESULTS

SUBROUTINE ROT_MAT( OH, ih, rotmat )
  !Rotation matrix from bond to lab frame
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

  OH1 = OH(ih,:)
  OH2 = OH(3-ih,:)

  !Get x
  CALL CROSSPROD(OH1,OH2,x) 
  x = x/norm2(x)
  
  !Get y
  CALL CROSSPROD(x,OH1,y)
  y = y/norm2(y)

  DO i = 1,3
    rotmat(i,1)=x(i)
    rotmat(i,2)=y(i)
    rotmat(i,3)=OH1(i)
  END DO

END SUBROUTINE ROT_MAT 

INTEGER FUNCTION ind_layer_pf( ind_O, n_hydrogens )
  !Layer index for prefactor. ind_O is the index of water oxygen in pos
  !This function also identifies terminal and bridging OH groups
  USE histogram,  ONLY : nlayerspf, layerspf
  USE parameters, ONLY : pos, box
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: ind_O, n_hydrogens
  INTEGER                    :: il
  DOUBLE PRECISION           :: z
  LOGICAL                    :: is_water_oxygen

  IF ( is_water_oxygen( ind_O ) ) THEN 

    IF ( n_hydrogens.eq.2 ) THEN !Molecular water

      z = pos(ind_O,3)
      z = z - nint(z/box(3,3))*box(3,3) + box(3,3)/2.

      DO il = 1,nlayerspf
     
        if ( z <= layerspf(il+1) ) THEN
          ind_layer_pf = il
          RETURN
        endif
     
      END DO

      print *, "Could not find layer for ", ind_O, z

    ELSE !Terminal OH

      ind_layer_pf = nlayerspf + 1
      RETURN

    END IF

  ELSE !Bridging OH

    ind_layer_pf = nlayerspf + 2
    RETURN

  END IF

  ind_layer_pf = 0

END FUNCTION ind_layer_pf

LOGICAL FUNCTION switched_oh_index( i_Ox, ind_H, frame ) 
  !IF the H index switched from frame-1 to frame, switch 
  !it back to avoid a discontinuity on the time-correlation 
  !function
  USE histogram, ONLY : ind_H_full 
  INTEGER, INTENT(IN)        :: i_Ox, frame
  INTEGER, INTENT(INOUT)     :: ind_H(2)
  INTEGER                    :: ih, ih_tmp

  IF ( frame > 1 ) THEN

    DO ih = 1,2

      !IF there is a switch of indexes
      IF ( ind_H_full( frame-1, i_Ox, ih ) .eq. ind_H(3-ih) ) THEN
        switched_oh_index = .true.
      END IF

    END DO

  END IF

  IF ( ind_H(1) .eq. ind_H(2) ) THEN
    print *, 'Equal indexes of H for a same Oxygen ', i_Ox
    STOP
  END IF

  switched_oh_index = .false.

END FUNCTION switched_oh_index 

LOGICAL FUNCTION bound_to_same_H( ind_Ox, ind_H , ti, tf, sep)
  !Check if Oxygen ind_Ox is bound to same H ind_H during the 
  !time interval ti:tf
  USE histogram,  ONLY : ind_H_full
  USE parameters, ONLY : natoms
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: ind_Ox, ind_H, ti, tf, sep
  INTEGER                    :: H_ini, H_end
  
  H_ini = ind_H_full(ti,ind_Ox,ind_H)
  H_end = ind_H_full(tf,ind_Ox,ind_H)

  IF ( ( H_ini .ne. 0 ) .and. ( H_ini .eq. H_end ) ) THEN
    bound_to_same_H = ALL( ind_H_full(ti:tf:sep,ind_Ox,ind_H) == H_ini )
  ELSE
    bound_to_same_H=.false.
  END IF

END FUNCTION bound_to_same_H
