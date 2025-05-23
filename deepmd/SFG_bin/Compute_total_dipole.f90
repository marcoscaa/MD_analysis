
MODULE correlation 
  IMPLICIT NONE 
  REAL                :: toplayerdown, bottomlayerup, zcenter
  REAL, ALLOCATABLE   :: sum_dipole(:,:), sum_polar(:,:,:)
  REAL, ALLOCATABLE   :: smooth_cutoff(:,:,:) 
END MODULE correlation

PROGRAM PTC 
  IMPLICIT NONE

  CALL INITIALIZE
  CALL OPEN_SET 
  CALL COARSE_GRAIN_POS_DIPOLES 
  CALL ASSIGN_DIPOLES 
  CALL COMPUTE_POLARIZABILITIES  
  CALL PRINT_RESULTS
  CALL CLOSE_SET 

END PROGRAM PTC

SUBROUTINE INITIALIZE
  USE correlation, ONLY : toplayerdown, bottomlayerup, zcenter
  USE parameters,  ONLY : natoms, nwater, &
                          atype, mean_pol, pos, nequil
  IMPLICIT NONE
  INTEGER                    :: i
  CHARACTER(100)             :: index_file

  CALL getarg(1, index_file)

  !Read the index file
  OPEN(unit=1, file = index_file)
  READ(1,*) natoms, nequil, nwater
  READ(1,*) mean_pol
  READ(1,*) toplayerdown, bottomlayerup, zcenter
  CLOSE(1)

  ALLOCATE(atype(natoms))
  ALLOCATE(pos(3,natoms))

  !Read atype from type.raw
  OPEN(unit = 1,file = 'type.raw')
  READ(1, *) (atype(i), i=1,natoms)
  CLOSE(1)

  OPEN(unit = 10,file = "total_dipole.raw")
  OPEN(unit = 11,file = "total_polarizability.raw")

END SUBROUTINE INITIALIZE

SUBROUTINE OPEN_SET
  USE correlation, ONLY : smooth_cutoff, sum_dipole, sum_polar
  USE parameters, ONLY : post, boxt, polart, natoms, nwater, &
                         nframes, polart, dipolet
  IMPLICIT NONE
  CHARACTER(100)             :: filename

  !Read the index file
  OPEN(unit=1, file = "nframes")
  READ(1,*) nframes
  CLOSE(1)

  ALLOCATE(post(3,natoms,nframes))
  ALLOCATE(boxt(3,3,nframes))
  ALLOCATE(polart(3,3,nwater,nframes))
  ALLOCATE(dipolet(3,nwater,nframes))
  ALLOCATE(smooth_cutoff(2,nwater,nframes)) !1: down; 2: up 
  ALLOCATE(sum_dipole(2,nframes))    !1: down; 2: up 
  ALLOCATE(sum_polar(2,2,nframes))   !1: down; 2: up 

  OPEN(1, file="coord.bin", action="read", form='unformatted', access='stream')
  OPEN(2, file="box.bin", action="read", form='unformatted', access='stream')
  OPEN(3, file="dipole.bin", action="read", form='unformatted', access='stream')
  OPEN(4, file="polarizability.bin", action="read", form='unformatted', access='stream')
  READ(1) post
  READ(2) boxt
  READ(3) dipolet
  READ(4) polart

  CLOSE(1);CLOSE(2);CLOSE(3);CLOSE(4)

END SUBROUTINE OPEN_SET

SUBROUTINE CLOSE_SET
  USE correlation, ONLY : sum_dipole, sum_polar, smooth_cutoff
  USE parameters, ONLY : post, boxt, polart, dipolet
  IMPLICIT NONE

  DEALLOCATE(post, boxt, polart, dipolet)
  DEALLOCATE(sum_dipole, sum_polar, smooth_cutoff)
  CLOSE(10); CLOSE(11)

END SUBROUTINE CLOSE_SET

SUBROUTINE PRINT_RESULTS 
  !Print to stdout
  USE correlation, ONLY : sum_polar, sum_dipole
  USE parameters,  ONLY : nframes
  IMPLICIT NONE
  INTEGER                    :: i

  DO i = 1,nframes

    write(10,fmt = '(2(E13.6,3X))'), sum_dipole(1,i), sum_dipole(2,i)
    write(11,fmt = '(4(E13.6,3X))'), sum_polar(1,1,i), sum_polar(2,1,i), &
                                     sum_polar(1,2,i), sum_polar(2,2,i)

  ENDDO

END SUBROUTINE PRINT_RESULTS

SUBROUTINE ASSIGN_DIPOLES(frame)
  USE parameters, ONLY : dipolet, pos, post, nwater, atype, &
                         natoms, boxt, nframes, box, boxinv
  USE correlation , ONLY : smooth_cutoff, sum_dipole
  IMPLICIT NONE
  INTEGER, INTENT(IN)                  :: frame
  INTEGER                              :: iat, ih, iw, ind_H(2), t
  REAL                     :: d, diptmp 

  sum_dipole=0.0

  DO t=1,nframes

    iw = 0
    pos=post(:,:,t)
    box=boxt(:,:,t)
    call matinv3(box,boxinv)

    DO iat=1,natoms

      IF (atype(iat)==0) THEN

        iw = iw + 1 

        !Electronic part
        diptmp = -2.d0*dipolet(3,iw,t)

        CALL get_H2(iat, ind_H)

        DO ih=1,2 
          
          d = pos(3,ind_H(ih)) - pos(3,iat)
          d = d - nint(d/box(3,3))*box(3,3)
          diptmp = diptmp + d

        END DO

        sum_dipole(:,t) = sum_dipole(:,t) + diptmp*smooth_cutoff(:,iw,t)
      
      END IF

    END DO

  END DO

END SUBROUTINE ASSIGN_DIPOLES

SUBROUTINE COMPUTE_POLARIZABILITIES
  !Compute isotropic and anisotropic polarizabilities
  USE correlation, ONLY : sum_polar, smooth_cutoff 
  USE parameters,  ONLY : polart, nwater,nframes
  IMPLICIT NONE
  INTEGER                    :: iw, i, t

  sum_polar=0.0

  DO t=1,nframes

    DO iw=1,nwater
      
      DO i=1,2

        sum_polar(i,1,t) = sum_polar(i,1,t) + &
                               polart(1,1,iw,t)*smooth_cutoff(i,iw,t) 
        sum_polar(i,2,t) = sum_polar(i,2,t) + &
                               polart(2,2,iw,t)*smooth_cutoff(i,iw,t)
      ENDDO

    END DO

  END DO

END SUBROUTINE COMPUTE_POLARIZABILITIES

SUBROUTINE COARSE_GRAIN_POS_DIPOLES
  !Define the molecules in layer_up and layer_down
  !Result is assigned to logical arrays in_layer_up and
  !in_layer_down
  USE correlation, ONLY : toplayerdown, bottomlayerup, & 
                          zcenter, smooth_cutoff
  USE parameters,  ONLY : post, boxt, natoms, atype, nframes
  IMPLICIT NONE
  INTEGER                    :: iat, iw, t
  REAL           :: posz
  REAL, PARAMETER :: damp=1.d0

  DO t = 1,nframes

    iw=0

    DO iat=1,natoms
    
      !Should be the index of oxygen
      IF (atype(iat)==0) THEN
    
        iw = iw + 1
        posz = post(3,iat,t) - zcenter
        posz = posz - nint( posz / boxt(3,3,t) ) * boxt(3,3,t)

        smooth_cutoff(1,iw,t) = (ERFC((posz-toplayerdown)/damp))/2.d0
        smooth_cutoff(2,iw,t) = (ERF((posz-bottomlayerup)/damp)+1)/2.d0

      ENDIF

    END DO

  END DO

END SUBROUTINE COARSE_GRAIN_POS_DIPOLES
