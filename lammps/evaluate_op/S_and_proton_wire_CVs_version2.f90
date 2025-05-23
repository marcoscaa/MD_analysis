!Computes the PT CV. This is not a general code!!!
!Only works if atype is ordered with O1 ... ON H1 ... H2N

MODULE OP
  INTEGER,ALLOCATABLE :: min_ind(:,:), wi_index(:,:)
  INTEGER,ALLOCATABLE :: Ox_net(:,:) 
  INTEGER             :: hind_wire, chain_index(9), nchain, nH_i
  REAL*8, ALLOCATABLE :: dist_mat(:,:), S_w(:) 
  REAL*8              :: dOO_total, dOHO
END MODULE OP

PROGRAM Coord_number
  USE parameters, ONLY : nframes
  IMPLICIT NONE
  INTEGER :: frame

  CALL INITIALIZE
  CALL REMOVE_EQUIL
  
  DO frame = 1,nframes
    CALL READ_ATOM_REDUCED
    CALL COMPUTE_DIST_MATRIX
    CALL FIND_OXYGEN_NETWORK (frame)
    CALL PROTON_WIRE_LENGTH(frame)
    CALL COMPUTE_CV
    CALL PRINT_CV(frame)
  END DO

  CLOSE(1); CLOSE(2)

END PROGRAM Coord_number

SUBROUTINE INITIALIZE
  USE parameters, ONLY : natoms, nwater, pos, atype, &
                         nframes, nequil
  USE OP
  IMPLICIT NONE
  INTEGER           :: i
  CHARACTER(100)    :: pos_file, index_file, water_ion_index

  !User should provide filename and index_file as arguments 
  CALL getarg(1, pos_file)
  CALL getarg(2, index_file)
  CALL getarg(3, water_ion_index)

  CALL READ_INDEX_FILES(index_file,water_ion_index)

  !Allocate module arrays
  ALLOCATE(pos(3,natoms))
  ALLOCATE(atype(natoms))
  ALLOCATE(min_ind(2,2*nwater))
  ALLOCATE(dist_mat(natoms,natoms))
  ALLOCATE(Ox_net(6,nwater))

  OPEN(unit = 1,file = pos_file)
  OPEN(unit = 2,file = 'proton_wire_cv.dat')

END SUBROUTINE INITIALIZE

SUBROUTINE READ_INDEX_FILES(ind1,ind2)
  USE parameters, only : natoms,nwater,nframes,nequil        
  USE OP
  IMPLICIT NONE  
  CHARACTER(100),INTENT(IN) :: ind1,ind2
  INTEGER :: iframe,frame

  OPEN(unit=1, file = ind1)
  READ(1,*) natoms, nwater, nframes, nequil
  CLOSE(1)

  ALLOCATE(S_w(nframes))
  ALLOCATE(wi_index(2,nframes))

  OPEN(unit=1, file = ind2)
  DO iframe=1,nframes
    READ(1,*) frame, S_w(iframe), wi_index(:,iframe)
    if (iframe.ne.frame) then
      WRITE(*,*) 'Wrong number of frames in the index file'
      STOP
    endif
  END DO
  CLOSE(1)

END SUBROUTINE READ_INDEX_FILES

SUBROUTINE COMPUTE_DIST_MATRIX
  USE parameters, only : natoms, nwater
  USE OP, ONLY :dist_mat
  IMPLICIT NONE
  REAL*8                     :: Dist
  INTEGER                    :: iat, jat

  DO iat=1,natoms
    DO jat=iat+1,natoms
      dist_mat(iat,jat) = Dist(iat,jat)
      dist_mat(jat,iat) = dist_mat(iat,jat)
    END DO
  END DO

END SUBROUTINE COMPUTE_DIST_MATRIX

SUBROUTINE FIND_OXYGEN_NETWORK (frame)
  USE parameters, only : pos,nwater,natoms
  USE OP, only : min_ind, Ox_net, dist_mat, wi_index, nH_i
  IMPLICIT NONE        
  INTEGER, INTENT(IN) :: frame
  REAL*8 :: oh_dist(nwater), Angle
  INTEGER :: oh_ind(nwater), oh_ind_0(nwater), ih, i, iox(nwater), tmp
  REAL*8, PARAMETER :: min_angle=0.866d0
  !REAL*8, PARAMETER :: min_angle=-1.d0

  iox = 0
  Ox_net=0

  do i=1,nwater
    oh_ind_0(i) = i
  end do

  !For each H, find the two closest O neighbors
  do ih=nwater+1,natoms
    oh_ind = oh_ind_0
    oh_dist=dist_mat(1:nwater,ih)
    call bubble_sort(oh_dist,oh_ind,2,nwater)
    min_ind(:,ih-nwater) = oh_ind(1:2)

    !This insures the two O connected by H are H-bonded. If yes,
    !we consider them to be connected
    if (Angle(pos(:,oh_ind(1)),pos(:,oh_ind(2)),pos(:,ih))>min_angle) then
      iox(oh_ind(1)) = iox(oh_ind(1)) + 1
      !iox(oh_ind(2)) = iox(oh_ind(2)) + 1
      Ox_net(iox(oh_ind(1)),oh_ind(1)) = oh_ind(2)
      !Ox_net(iox(oh_ind(2)),oh_ind(2)) = oh_ind(1)
    end if
  end do

  nH_i = 0
  do ih=1,natoms-nwater
    if (min_ind(1,ih).eq.wi_index(1,frame)) nH_i=nH_i+1
  end do

  if (nH_i==3) then
    tmp = wi_index(1,frame)
    wi_index(1,frame) = wi_index(2,frame)
    wi_index(2,frame)=tmp
  else if (nH_i==1) then
    print *, 'Wront number of H for frame ', frame
  end if

END SUBROUTINE FIND_OXYGEN_NETWORK

SUBROUTINE PROTON_WIRE_LENGTH(frame)
  USE OP, only : Ox_net, dist_mat, chain_index, dOO_total, wi_index
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: frame
  INTEGER :: j(9), io, io1, io2, io3, io4, io5, io6, io7
  REAL*8 :: dchain(8), dmax

  dmax=100.d0
  dchain=0.0d0
  chain_index=0
  j=0

  j(1)=wi_index(2,frame) !Initial O is the H2O(donor)
  chain_index(1)=j(1)
  DO io=1,6 ! 6 is the maximum number of H-bonds H2O will ever make
    j(2)=Ox_net(io,j(1)) !Second O receives H-bond from j(1)
    if (j(2)==0) EXIT
    dchain(1)=dist_mat(j(1),j(2))
    if (j(2)==wi_index(1,frame)) then 
       call check_max_path(j,chain_index,dchain,dmax,2)
    end if

    DO io1=1,6
      j(3)=Ox_net(io1,j(2)) !Third O receives H-bond from j(2)
      if (j(3)==0) EXIT
      dchain(2)=dist_mat(j(2),j(3))
      if (j(3)==wi_index(1,frame)) then 
        call check_max_path(j,chain_index,dchain,dmax,3)
      end if

      DO io2=1,6
        j(4)=Ox_net(io2,j(3)) !Fourth O receives H-bond from j(3)
        if (j(4)==0) EXIT
        dchain(3)=dist_mat(j(3),j(4))
        if (j(4)==wi_index(1,frame)) then
          call check_max_path(j,chain_index,dchain,dmax,4)
        end if

        DO io3=1,6
          j(5)=Ox_net(io3,j(4)) !Fifth O receives H-bond from j(4)
          if (j(5)==0) EXIT
          dchain(4)=dist_mat(j(4),j(5))
          if (j(5)==wi_index(1,frame)) then
            call check_max_path(j,chain_index,dchain,dmax,5)
          end if

          DO io4=1,6
            j(6)=Ox_net(io4,j(5)) !Sixth O receives H-bond from j(5)
            if (j(6)==0) EXIT
            dchain(5)=dist_mat(j(5),j(6))
            if (j(6)==wi_index(1,frame)) then
              call check_max_path(j,chain_index,dchain,dmax,6)
            end if

            DO io5=1,6
              j(7) = Ox_net(io5,j(6))
              if (j(7)==0) EXIT
              dchain(6)=dist_mat(j(6),j(7))
              if (j(7)==wi_index(1,frame)) then
                call check_max_path(j,chain_index,dchain,dmax,7)
              end if

              DO io6=1,6
                j(8) = Ox_net(io6,j(7))
                if (j(8)==0) EXIT
                dchain(7)=dist_mat(j(7),j(8))
                if (j(8)==wi_index(1,frame)) then
                  call check_max_path(j,chain_index,dchain,dmax,8)
                end if

                DO io7=1,6
                  j(9) = Ox_net(io7,j(8))
                  if (j(9)==0) EXIT
                  dchain(8)=dist_mat(j(8),j(9))
                  if (j(9)==wi_index(1,frame)) then
                    call check_max_path(j,chain_index,dchain,dmax,9)
                  end if
                END DO
              END DO
            END DO
          END DO
        END DO
      END DO
    END DO

  END DO

  dOO_total=dmax

  if (dmax==100.d0) then
    WRITE(*,*) 'Could not find path ', frame
  end if

END SUBROUTINE PROTON_WIRE_LENGTH

SUBROUTINE CHECK_MAX_PATH(ind_tmp,ind_chain,dchain,dmax,chain)
  IMPLICIT NONE
  INTEGER :: ind_tmp(9), ind_chain(9), chain
  REAL*8 :: dOO, dchain(8), dmax
  INTEGER :: i

  dOO = sum(dchain(1:chain-1))
  if (dOO<dmax) then
    dmax=dOO
    if (i<9) then
      do i=chain+1,9
        ind_tmp(i)=0
        end do
    end if
    ind_chain=ind_tmp
  end if

END SUBROUTINE CHECK_MAX_PATH

SUBROUTINE COMPUTE_CV
  USE parameters, only : natoms, nwater
  USE OP, only : dOHO, nchain, dOO_total, nchain, chain_index, nH_i, Ox_net, dist_mat
  IMPLICIT NONE
  INTEGER :: io1,io2,ih,ichain
  REAL*8 :: COMPUTE_dOHO

  dOHO=0.d0
  nchain=1
  DO ichain=1,8
    io1=chain_index(ichain)
    io2=chain_index(ichain+1)
    if (io2.ne.0) then
     nchain=nchain+1
     dOHO = dOHO + COMPUTE_dOHO(io1,io2)
    end if
  END DO

  if (nchain>1) then
    dOHO=dOHO/dble(nchain-1)
    dOO_total=dOO_total/dble(nchain-1)
  else
    !If the path is not found, compute CVs using neighbor of proton donor species
    io1=chain_index(1)
    dOHO = COMPUTE_dOHO(io1,Ox_net(1,io1))
    dOO_total = dist_mat(io1,Ox_net(1,io1))
  end if
  if (nH_i.eq.3) dOHO=-dOHO

END SUBROUTINE COMPUTE_CV

REAL*8 FUNCTION COMPUTE_dOHO(io1,io2)
  USE OP, only : min_ind, dist_mat
  USE parameters, only : natoms, nwater
  IMPLICIT NONE
  integer, intent(in) :: io1,io2
  integer :: ih
  real*8 :: dOHO

  DO ih=1,natoms-nwater
    if ((min_ind(1,ih)==io1).and.(min_ind(2,ih)==io2)) then
      dOHO = dist_mat(io1,ih+nwater) - dist_mat(io2,ih+nwater)
      EXIT
    end if
  END DO

  COMPUTE_dOHO = dOHO

END FUNCTION COMPUTE_dOHO

SUBROUTINE PRINT_CV(step)
  USE OP
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: step
  REAL*8              :: Dist

  WRITE(2,fmt='(I10,2(3X,F12.8),3X,I2)') step, dOHO,&
          dOO_total, nchain
  
END SUBROUTINE PRINT_CV
