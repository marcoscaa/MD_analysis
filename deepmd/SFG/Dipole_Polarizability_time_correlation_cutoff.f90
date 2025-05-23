
MODULE correlation 
  IMPLICIT NONE 
  INTEGER                         :: tcorr !Length of the time corr. function
  INTEGER                         :: dcart !Cartesian component of dipole
  INTEGER                         :: pcart1, pcart2 !Cartesian of polarizability
  DOUBLE PRECISION                :: layer1down,layer1up
  DOUBLE PRECISION                :: layer2down,layer2up
  DOUBLE PRECISION                :: zcenter 
  DOUBLE PRECISION, ALLOCATABLE   :: polart(:,:), dipolet(:,:)
  DOUBLE PRECISION, ALLOCATABLE   :: intra_corr(:), inter_corr(:) 
  DOUBLE PRECISION, ALLOCATABLE   :: signdipole(:,:) 
  DOUBLE PRECISION, PARAMETER     :: rcut=6.d0
END MODULE correlation

PROGRAM PTC 
  USE parameters, ONLY : nframes, coarse
  USE correlation, ONLY: rcut
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  !CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_RAW_DIPOLE_COMPONENTS(1)
    CALL READ_RAW_POLARIZABILITY_COMPONENTS(1,2)
    CALL READ_RAW_POS_BOX
    CALL ASSIGN_DIPOLES (frame)
    CALL COMPUTE_POLARIZABILITIES (frame) 
    CALL MAKE_NEIGHBOR_LIST (frame,rcut)
    CALL COARSE_GRAIN_POS_DIPOLES (frame)
  END DO

  !CALL SUBTRACT_AVERAGE
  CALL DIPOLE_POLAR_CORR 
  CALL PRINT_RESULTS

  CLOSE(1);CLOSE(2);CLOSE(4);CLOSE(5)

END PROGRAM PTC

SUBROUTINE INITIALIZE
  USE correlation, ONLY : polart, dipolet, tcorr, &
                          intra_corr, inter_corr, &
  !                        dcart, pcart1, pcart2, &
                          layer1down,layer1up,&
                          layer2down,layer2up, &
                          signdipole, zcenter
  USE parameters,  ONLY : pos, polar, natoms, nwater, &
                          nframes, nequil, within_cutoff, &
                          dt, atype, mean_pol, at_surface, &
                          polar, dipole 
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: index_file

  CALL getarg(1, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nframes, nequil, nwater, dt, tcorr
  READ(1,*) mean_pol
  !READ(1,*) dcart, pcart1, pcart2
  READ(1,*) layer1down,layer1up,layer2down,layer2up, zcenter
  CLOSE(1)

  ALLOCATE(pos(3,natoms)) ! coordinates
  ALLOCATE(polar(1,2,nwater)) ! molecular polarizability tensors
  ALLOCATE(dipole(1,nwater)) ! Molecular dipole
  ALLOCATE(polart(nwater,nframes)) ! molecular polarizability 
  ALLOCATE(dipolet(nwater,nframes)) ! Molecular dipole
  ALLOCATE(signdipole(nwater,nframes)) ! Orientation of dipole(z)
  ALLOCATE(within_cutoff(nwater,nwater,nframes))
  ALLOCATE(at_surface(nwater,nframes)) !Coarse grain function
  ALLOCATE(intra_corr(tcorr)); intra_corr=0.d0 
  ALLOCATE(inter_corr(tcorr)); inter_corr=0.d0
  ALLOCATE(atype(natoms))

  !Read atype from type.raw
  OPEN(unit = 1,file = 'type.raw')
  READ(1, *) (atype(i), i=1,natoms)
  CLOSE(1)

  OPEN(unit = 1,file = 'coord.raw')
  OPEN(unit = 2,file = 'box.raw')
  OPEN(unit = 4,file = 'dipole.raw')
  OPEN(unit = 5,file = 'polarizability.raw')
  
END SUBROUTINE INITIALIZE

SUBROUTINE DIPOLE_POLAR_CORR
  !Calculate the lifetime of water in evenly separated bins in the z direction
  USE correlation, ONLY : polart, intra_corr, inter_corr, & 
                          tcorr, dipolet, signdipole
  USE parameters,  ONLY : nframes, nwater, within_cutoff, &
                          at_surface
  IMPLICIT NONE
  DOUBLE PRECISION           :: sum_intra, sum_inter
  DOUBLE PRECISION           :: tracematmul
  INTEGER, PARAMETER         :: sep=1
  INTEGER                    :: frame1, frame2, i_at, j_at, nt0
  INTEGER                    :: frame1p2

  !frame1 is the duration of the interval 
  DO frame1 = 0,tcorr-1

    nt0=0

    !Change initial condition, look for combinations with same interval
    DO frame2 = 1,nframes-frame1,sep

      sum_intra=0.d0;sum_inter=0.d0
      frame1p2=frame1+frame2
      nt0=nt0+1

!$omp parallel do reduction(+:sum_inter,sum_intra)
      DO i_at = 1,nwater

        IF (at_surface(i_at,frame2)) THEN

          sum_intra = sum_intra + signdipole(i_at,frame2) * &
              dipolet(i_at,frame2) * polart(i_at,frame1p2) 

          DO j_at = 1, nwater

            IF ( within_cutoff(j_at, i_at, frame2) ) THEN
!            IF ( at_surface(j_at, frame1p2) ) THEN

              sum_inter = sum_inter + signdipole(i_at,frame2) * &
                dipolet(i_at,frame2) * polart(j_at,frame1p2) 
                

            END IF

          ENDDO

        END IF

      ENDDO
!$omp end parallel do

      intra_corr(frame1+1) = intra_corr(frame1+1) + sum_intra
      inter_corr(frame1+1) = inter_corr(frame1+1) + sum_inter

    ENDDO

    !Time-averaged correlation functions
    intra_corr(frame1+1) = intra_corr(frame1+1) / dble(nt0)
    inter_corr(frame1+1) = inter_corr(frame1+1) / dble(nt0)

  ENDDO

END SUBROUTINE DIPOLE_POLAR_CORR 

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE correlation, ONLY : intra_corr, inter_corr, tcorr
  USE parameters,  ONLY : dt, nframes
  IMPLICIT NONE
  INTEGER                    :: i

  OPEN(unit = 10,file = 'dptc.dat')

  WRITE(10, *), "# Time, Intra-corr, Inter-corr"

  DO i = 1,tcorr

    write(10,fmt = '(F12.8,4(3X,E18.11))'), dble(i-1)*dt, &
          intra_corr(i), inter_corr(i)

  ENDDO

  CLOSE(10)

END SUBROUTINE PRINT_RESULTS

SUBROUTINE ASSIGN_DIPOLES(frame)
  !We assume atomic type to be on the order: O H H O H H
  USE parameters, ONLY : dipole, pos, nwater, atype, natoms, box
  USE correlation , ONLY : dcart, dipolet
  IMPLICIT NONE
  INTEGER, INTENT(IN)                  :: frame
  INTEGER                              :: iat, ih, iw, ind_H(2)
  DOUBLE PRECISION                     :: d(4)
  
  iw = 0

  DO iat=1,natoms

    IF (atype(iat)==0) THEN

      iw = iw + 1 

      !Electronic part
      dipolet(iw,frame) = -2.d0*dipole(1,iw)
      !dipole(:,iw,frame) = 0.d0 

      CALL get_H2(iat, ind_H)

      DO ih=1,2 
        
        CALL DISTANCE_VECTOR(pos(:,iat), pos(:,ind_H(ih)), d)
        dipolet(iw,frame) = dipolet(iw,frame) + d(3)!d(dcart) 

      END DO

      dipolet(iw,frame) = dipolet(iw,frame) - &
        nint(dipolet(iw,frame)/box(3,3))*box(3,3)
        !nint(dipolet(iw,frame)/box(dcart,dcart))*box(dcart,dcart)
    
    END IF

  END DO

END SUBROUTINE ASSIGN_DIPOLES

SUBROUTINE SUBTRACT_AVERAGE
  USE correlation, ONLY : polart
  USE parameters,  ONLY : nframes, nwater
  IMPLICIT NONE 
  INTEGER                    :: it, iw
  DOUBLE PRECISION           :: mean
   
  mean=0.d0

  DO it=1,nframes
    DO iw=1,nwater

      mean = mean + polart(iw,it)

    END DO
  END DO

  polart = polart - mean / dble(nwater*nframes)

END SUBROUTINE SUBTRACT_AVERAGE

SUBROUTINE COMPUTE_POLARIZABILITIES (frame)
  !Compute isotropic and anisotropic polarizabilities
  USE correlation, ONLY : polart, pcart1, pcart2
  USE parameters,  ONLY : polar, nwater
  IMPLICIT NONE
  INTEGER, INTENT(in)        :: frame
  INTEGER                    :: iw, i, j

  DO iw=1,nwater
    
    !DIRTY HACK WARINING: for computational efficiency only. 
    !polart(iw,frame) = ( polar(pcart1,pcart1,iw) + polar(pcart2,pcart2,iw) ) / 2.
    polart(iw,frame) = ( polar(1,1,iw) + polar(1,2,iw) ) / 2.

  END DO

END SUBROUTINE COMPUTE_POLARIZABILITIES

SUBROUTINE COARSE_GRAIN_POS_DIPOLES (frame)
  !Define the surface region in our simulation box
  !We assign the sign of the dipole(z) following the convention:
  !dipole(3) > 0 : pointing away from the surface
  !Layer2 is assumed to be the layer with normal vector 
  !pointing to negative z direction
  USE correlation, ONLY : layer1down,layer1up, layer2down,layer2up, & 
                          dipolet, dcart, signdipole, zcenter
  USE parameters,  ONLY : pos, box, natoms, nwater, atype, &
                          at_surface
  IMPLICIT NONE
  INTEGER, INTENT(IN)        :: frame
  INTEGER                    :: iat, iw
  DOUBLE PRECISION           :: posz

  iw=0
  signdipole(:,frame) = 1

  DO iat=1,natoms
  
    !Should be the index of oxygen
    IF (atype(iat)==0) THEN
  
      iw = iw + 1
      posz = pos(3,iat) - zcenter
      posz = posz - nint( posz / box(3,3) ) * box(3,3)

      IF (( posz>layer1down ).and.( posz<layer1up )) THEN 
        at_surface(iw,frame) = .true.
      ELSEIF (( posz>layer2down ).and.( posz<layer2up )) THEN
        at_surface(iw,frame) = .true.
        !IF (dcart==3) signdipole(iw,frame) = -1
        signdipole(iw,frame) = -1
      ELSE
        at_surface(iw,frame) = .false.
      ENDIF

    ENDIF

  END DO

END SUBROUTINE COARSE_GRAIN_POS_DIPOLES
