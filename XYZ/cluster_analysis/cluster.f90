module cluster_mod
  IMPLICIT NONE
  ! Atomic coordinates of atoms of point_type will be stored in points
  REAL*8, ALLOCATABLE        :: points(:,:)
  ! n is the total number of atoms of type point_type
  INTEGER :: n, total_monomer_count, num_clusters_found
  ! point_type is the element symbol to be included in the analysis
  CHARACTER(5) :: point_type
  ! Distance threshold in Angstrom units to construct cluster
  REAL*8 :: dist_threshold
end module

PROGRAM CLUSTER
  ! This code will read a trajectory in extendend xyz format and 
  ! print the total number of monomers in clusters, as well as the 
  ! total number of clusters to the cluster.dat file
  ! Usage:
  ! echo atype dist_thrs | cluster.x trajfile.extxyz 
  ! atype(string):   chemical symbol of the element to be computed
  ! dist_thrs(real): distance in Angstroms to consider two atype 
  !                  atoms to be clustered

  USE parameters, ONLY : nframes, stride
  USE, intrinsic :: iso_fortran_env, Only : iostat_end
  IMPLICIT NONE
  INTEGER :: iostat ! For trajectory reading 

  CALL INITIALIZE
  
  DO 
    CALL READ_EXTXYZ_IO(iostat)
    IF (iostat == iostat_end ) THEN
      EXIT
    END IF
    CALL SELECT_POINT_TYPE
    CALL COUNT_CLUSTER_MONOMERS
    CALL PRINT_RESULTS
  END DO

  CLOSE(1)

END PROGRAM CLUSTER

SUBROUTINE INITIALIZE
  USE cluster_mod, only : point_type,dist_threshold
  IMPLICIT NONE
  CHARACTER(100)             :: pos_file

  !User should provide filename as argument 
  CALL getarg(1, pos_file)
  READ *, point_type, dist_threshold 

  ! Open trajectory file
  OPEN(unit = 1,file = pos_file)
  ! Open the output file and write header
  OPEN (unit = 2, file = 'cluster.dat')

END SUBROUTINE INITIALIZE

SUBROUTINE SELECT_POINT_TYPE
  USE cluster_mod, only : points,point_type,n
  USE parameters, only : pos, atype, natoms
  ! Builds the matrix points using only atomic coordinates with 
  ! elements point_type
  IMPLICIT NONE
  INTEGER                  :: iat
  REAL*8                   :: pos_(3,natoms)

  ! Initialize the total number of elements with type point_type
  n = 0 

  DO iat=1,natoms
      if (trim(atype(iat))==trim(point_type)) then
          n = n+1
          pos_(:,n) = pos(:,iat)
      end if
  END DO

  IF (ALLOCATED(points)) deallocate(points)
  ALLOCATE(points(3,n))

  DO iat=1,n
    points(:,iat) = pos_(:,iat)
  END DO

END SUBROUTINE SELECT_POINT_TYPE

SUBROUTINE COUNT_CLUSTER_MONOMERS
  USE cluster_mod, only : points, n, dist_threshold, &
                          total_monomer_count, num_clusters_found
  ! Description:
  !   This subroutine computes the total number of monomers (individual points)
  !   that belong to any cluster of points in 3D. A cluster is defined by a
  !   distance threshold: two points are considered connected if the distance
  !   between them is less than or equal to 'dist_threshold'. The subroutine
  !   uses a Breadth-First Search (BFS) algorithm to identify connected components
  !   (clusters). For each identified cluster, it prints its size (number of monomers).
  !   It also returns the total number of monomers across all clusters and the
  !   total number of clusters found.

  ! Outputs:
  !   total_monomer_count : INTEGER, INTENT(OUT)
  !                       The total count of all points that are part of any
  !                       identified cluster. This is the sum of all individual
  !                       cluster sizes.
  !   num_clusters_found : INTEGER, INTENT(OUT)
  !                       The total number of distinct clusters identified.

  IMPLICIT NONE

  ! Local Variables:
  LOGICAL, DIMENSION(n) :: visited
  !   'visited' array keeps track of points that have already been processed
  !   or added to a cluster.

  INTEGER, DIMENSION(n) :: queue
  !   'queue' is used for the Breadth-First Search (BFS) traversal. It stores
  !   the indices of points to be processed in the current cluster.

  INTEGER :: q_head, q_tail
  !   'q_head' points to the front of the queue (next point to process).
  !   'q_tail' points to the back of the queue (next available slot).

  INTEGER :: i, j, current_idx
  !   Loop counters and index of the current point being processed.

  REAL*8 :: d
  !   Stores the Euclidean distance between two points.

  REAL*8 :: DistXYZ
  !   Function to compute the Euclidean distance including PBC

  INTEGER :: current_cluster_size
  !   Temporarily stores the size of the cluster currently being processed.

  INTEGER :: cluster_size_distribution(n)
  !   Histogram the number of monomers in clusters

  ! --- Initialization ---
  total_monomer_count = 0
  !   Initialize the total count of monomers across all clusters to zero.

  num_clusters_found = 0
  cluster_size_distribution = 0
  !   Initialize the count of clusters found to zero.

  visited = .FALSE.
  !   Mark all points as unvisited initially.

  ! --- Main Loop: Iterate through all points to find unvisited ones and start new clusters ---
  DO i = 1, n
    ! If the current point 'i' has not been visited yet, it means it either
    ! starts a new cluster or is an isolated point (which itself forms a cluster of size 1).
    IF (.NOT. visited(i)) THEN
      ! Increment the number of clusters found.
      num_clusters_found = num_clusters_found + 1
      current_cluster_size = 0
      !   Reset the size for the new cluster.

      ! --- Start a new Breadth-First Search (BFS) for a cluster ---
      ! Initialize the queue for the new cluster.
      q_head = 1
      q_tail = 1

      ! Add the starting point 'i' to the queue and mark it as visited.
      queue(q_tail) = i
      visited(i) = .TRUE.
      q_tail = q_tail + 1 ! Move tail to next available slot

      ! --- BFS Traversal Loop ---
      ! Continue as long as there are points in the queue to process.
      DO WHILE (q_head < q_tail)
        ! Get the next point from the front of the queue.
        current_idx = queue(q_head)
        q_head = q_head + 1 ! Move head to the next point

        ! This point 'current_idx' is part of the current cluster,
        ! so increment its size.
        current_cluster_size = current_cluster_size + 1

        ! --- Check all other points for proximity to the current_idx point ---
        DO j = 1, n
            ! Only consider points that have not been visited yet to avoid
            ! redundant processing and infinite loops in cyclic graphs.
            IF (.NOT. visited(j)) THEN
                ! Calculate the squared Euclidean distance between current_idx and j.
                d = DistXYZ(points(:,i),points(:,j))  

                ! If the distance is within the squared threshold,
                ! point 'j' is connected to 'current_idx' and belongs to the same cluster.
                IF (d <= dist_threshold) THEN
                    ! Add point 'j' to the queue for future processing and mark it as visited.
                    queue(q_tail) = j
                    visited(j) = .TRUE.
                    q_tail = q_tail + 1 ! Move tail to next available slot
                END IF
            END IF
        END DO
      END DO ! End of BFS Traversal Loop

      ! Add the current cluster's size to the total monomer count.
      total_monomer_count = total_monomer_count + current_cluster_size
      cluster_size_distribution(current_cluster_size) = &
          cluster_size_distribution(current_cluster_size) + 1
    END IF
  END DO ! End of Main Loop (iterating through all points)

  print *, cluster_size_distribution

END SUBROUTINE COUNT_CLUSTER_MONOMERS

SUBROUTINE PRINT_RESULTS 
  USE cluster_mod, only : total_monomer_count,num_clusters_found
  IMPLICIT NONE
  REAL*8 :: ratio

  ratio = float(total_monomer_count)/float(num_clusters_found)
  
  WRITE(2,*) total_monomer_count,num_clusters_found, ratio

END SUBROUTINE

