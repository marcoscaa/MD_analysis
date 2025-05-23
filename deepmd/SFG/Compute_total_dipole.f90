
MODULE correlation 
  IMPLICIT NONE 
  INTEGER                         :: tcorr !Length of the time corr. function
  INTEGER                         :: dcart !Cartesian component of dipole
  INTEGER                         :: pcart1, pcart2 !Cartesian of polarizability
  DOUBLE PRECISION                :: toplayerdown, bottomlayerup, zcenter
  DOUBLE PRECISION, ALLOCATABLE   :: polart(:,:,:), dipolet(:,:)
  DOUBLE PRECISION, ALLOCATABLE   :: smooth_cutoff(:,:) 
END MODULE correlation

PROGRAM PTC 
  USE parameters, ONLY : nframes, coarse
  IMPLICIT NONE
  INTEGER                    :: frame

  CALL INITIALIZE
  !CALL REMOVE_EQUIL

  DO frame = 1,nframes
    CALL READ_RAW_DIPOLE_COMPONENTS(1)
    CALL READ_RAW_POLARIZABILITY_COMPONENTS(1,2)
    CALL READ_RAW_POS_BOX
    CALL COARSE_GRAIN_POS_DIPOLES (frame)
    CALL ASSIGN_DIPOLES (frame)
    CALL COMPUTE_POLARIZABILITIES (frame) 
  END DO

  CALL PRINT_RESULTS

  CLOSE(1);CLOSE(2);CLOSE(4);CLOSE(5)

END PROGRAM PTC

SUBROUTINE INITIALIZE
  USE correlation
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
  READ(1,*) natoms, nframes, nequil, nwater
  READ(1,*) mean_pol
!  READ(1,*) dcart, pcart1, pcart2
  READ(1,*) toplayerdown, bottomlayerup, zcenter
  CLOSE(1)

  ALLOCATE(pos(3,natoms)) ! coordinates
  ALLOCATE(polar(1,2,nwater)) ! molecular polarizability tensors
  ALLOCATE(dipole(1,nwater)) ! Molecular dipole
  !First index: 1: Layer up; 2: layer down
  ALLOCATE(polart(2,2,nframes)) ! Total polarizability of layer
  !First  index: 1: Layer up  2: layer down
  !Second index: 1: pcart1    2: pcart2
  ALLOCATE(dipolet(2,nframes)) ! Total dipole of layer
  ALLOCATE(smooth_cutoff(2,nwater)) !1: down; 2: up 
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

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE correlation, ONLY : polart, dipolet
  USE parameters,  ONLY : dt, nframes
  IMPLICIT NONE
  INTEGER                    :: i

  OPEN(unit = 10,file = 'total_dipole.raw')
  OPEN(unit = 11,file = 'total_polarizability.raw')
  

  DO i = 1,nframes

    write(10,fmt = '(2(E13.6,3X))'), dipolet(1,i), dipolet(2,i)
    write(11,fmt = '(4(E13.6,3X))'), polart(1,1,i), polart(2,1,i), &
                                     polart(1,2,i), polart(2,2,i)

  ENDDO

  CLOSE(10); CLOSE(11)

END SUBROUTINE PRINT_RESULTS

SUBROUTINE ASSIGN_DIPOLES(frame)
  !We assume atomic type to be on the order: O H H O H H
  USE parameters, ONLY : dipole, pos, nwater, atype, natoms, box
  USE correlation , ONLY : dcart, dipolet, smooth_cutoff
  IMPLICIT NONE
  INTEGER, INTENT(IN)                  :: frame
  INTEGER                              :: iat, ih, iw, ind_H(2)
  DOUBLE PRECISION                     :: d(4), diptmp
  
  iw = 0
  dipolet(:,frame) = 0.d0

  DO iat=1,natoms

    IF (atype(iat)==0) THEN

      iw = iw + 1 

      !Electronic part
      diptmp = -2.d0*dipole(1,iw)

      CALL get_H2(iat, ind_H)

      DO ih=1,2 
        
        CALL DISTANCE_VECTOR(pos(:,iat), pos(:,ind_H(ih)), d)
        diptmp = diptmp + d(3) 

      END DO

      diptmp = diptmp - nint(diptmp/box(3,3))*box(3,3)

      dipolet(:,frame) = dipolet(:,frame) + diptmp*smooth_cutoff(:,iw)
    
    END IF

  END DO

END SUBROUTINE ASSIGN_DIPOLES

SUBROUTINE COMPUTE_POLARIZABILITIES (frame)
  !Compute isotropic and anisotropic polarizabilities
  USE correlation, ONLY : polart, pcart1, pcart2, &
                          smooth_cutoff 
  USE parameters,  ONLY : polar, nwater
  IMPLICIT NONE
  INTEGER, INTENT(in)        :: frame
  INTEGER                    :: iw, i, j

  polart(:,:,frame) = 0.d0

  DO iw=1,nwater
    
    DO i=1,2

      polart(i,1,frame) = polart(i,1,frame) + &
                          polar(1,1,iw)*smooth_cutoff(i,iw) 
      polart(i,2,frame) = polart(i,2,frame) + &
                          polar(1,2,iw)*smooth_cutoff(i,iw)
    ENDDO

  END DO

END SUBROUTINE COMPUTE_POLARIZABILITIES

SUBROUTINE COARSE_GRAIN_POS_DIPOLES
  !Define the molecules in layer_up and layer_down
  !Result is assigned to logical arrays in_layer_up and
  !in_layer_down
  USE correlation, ONLY : toplayerdown, bottomlayerup, & 
                          zcenter, smooth_cutoff
  USE parameters,  ONLY : pos, box, natoms, atype
  IMPLICIT NONE
  INTEGER                    :: iat, iw
  DOUBLE PRECISION           :: posz
  DOUBLE PRECISION, PARAMETER :: damp=1.d0

  iw=0

  DO iat=1,natoms
  
    !Should be the index of oxygen
    IF (atype(iat)==0) THEN
  
      iw = iw + 1
      posz = pos(3,iat) - zcenter
      posz = posz - nint( posz / box(3,3) ) * box(3,3)

      smooth_cutoff(1,iw) = (ERFC((posz-toplayerdown)/damp))/2.d0
      smooth_cutoff(2,iw) = (ERF((posz-bottomlayerup)/damp)+1)/2.d0

    ENDIF

  END DO

END SUBROUTINE COARSE_GRAIN_POS_DIPOLES
