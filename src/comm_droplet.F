!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Jian Sun - 09/28/2022:                                         !
!    This module is generated following John Dennis's suggestion !
!    to separate the MPI communication of droplets from other    !
!    MPI calls in the existing comm.F code                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define _COUNTPREFIX 1

module comm_droplet_module

implicit none

private

public :: comm_droplet_number,comm_droplet_value
public :: setupIndexPointers, setupDepartDroplet
public :: makeContiguous
public :: CollectDropCount, CollectDropInfo

contains

#ifdef MPI

  subroutine comm_droplet_number(numHoles,Depart,Arrive)

    use input, only : myid,nparcelsLocal,npvals,mynw,mysw,myne, &
                      myse,myeast,mywest,mynorth,mysouth,ierr, &
                      nparcelsLocalActive
    use constants, only: undefined_index,num_nn,inorth,isouth, &
                         iwest,ieast,inw,ine,isw,ise
    use mpi

    implicit none

    integer, intent(in) :: numHoles                                      ! Current number of holes available 

    integer, intent(out) :: Depart(num_nn)                               ! Number of droplets that will enter each 
                                                                         ! nearest neighbor at different directions
    integer, intent(out) :: Arrive(num_nn)                               ! Number of droplets that will enter the 
                                                                         ! current MPI region from each nearest 
                                                                         ! neighbor at different directions
    ! Local variables

    integer :: n_recv,n_send
    integer :: reqs(16)
    integer :: tag_n,tag_s,tag_w,tag_e,tag_nw,tag_ne,tag_sw,tag_se,indx
    integer, dimension(mpi_status_size,16) :: status1
    integer :: i,n

    ! initiate some MPI index and tag values

    tag_n  = 1001
    tag_s  = 1002
    tag_w  = 1003
    tag_e  = 1004
    tag_nw = 1005
    tag_sw = 1006
    tag_ne = 1007
    tag_se = 1008

    ! We are just exchanging some integers, so just use the correct
    ! version that is already present on the host
    ! initiate the MPI non-blocking receive interface to 
    ! know how many droplets from the nearest neighbors

    n_recv = 1
    call mpi_irecv(Arrive(inorth),1,MPI_INT,mynorth,tag_n,MPI_COMM_WORLD, &
                   reqs(n_recv),ierr)

    n_recv = n_recv + 1
    call mpi_irecv(Arrive(isouth),1,MPI_INT,mysouth,tag_s,MPI_COMM_WORLD, &
                   reqs(n_recv),ierr)

    n_recv = n_recv + 1
    call mpi_irecv(Arrive(iwest),1,MPI_INT,mywest,tag_w,MPI_COMM_WORLD, &
                   reqs(n_recv),ierr)

    n_recv = n_recv + 1
    call mpi_irecv(Arrive(ieast),1,MPI_INT,myeast,tag_e,MPI_COMM_WORLD, &
                   reqs(n_recv),ierr)

    n_recv = n_recv + 1
    call mpi_irecv(Arrive(inw),1,MPI_INT,mynw,tag_nw,MPI_COMM_WORLD, &
                   reqs(n_recv),ierr)

    n_recv = n_recv + 1
    call mpi_irecv(Arrive(ine),1,MPI_INT,myne,tag_ne,MPI_COMM_WORLD, &
                   reqs(n_recv),ierr)

    n_recv = n_recv + 1
    call mpi_irecv(Arrive(isw),1,MPI_INT,mysw,tag_sw,MPI_COMM_WORLD, &
                   reqs(n_recv),ierr)

    n_recv = n_recv + 1
    call mpi_irecv(Arrive(ise),1,MPI_INT,myse,tag_se,MPI_COMM_WORLD, &
                   reqs(n_recv),ierr)

    ! initiate the MPI non-blocking send interface to 
    ! send how many droplets entering the nearest neighbor 

    n_send = 9
    call mpi_isend(Depart(inorth),1,MPI_INT,mynorth,tag_s,MPI_COMM_WORLD, &
                   reqs(n_send),ierr)

    n_send = n_send + 1
    call mpi_isend(Depart(isouth),1,MPI_INT,mysouth,tag_n,MPI_COMM_WORLD, &
                   reqs(n_send),ierr)

    n_send = n_send + 1
    call mpi_isend(Depart(iwest),1,MPI_INT,mywest,tag_e,MPI_COMM_WORLD, &
                   reqs(n_send),ierr)

    n_send = n_send + 1
    call mpi_isend(Depart(ieast),1,MPI_INT,myeast,tag_w,MPI_COMM_WORLD, &
                   reqs(n_send),ierr)

    n_send = n_send + 1
    call mpi_isend(Depart(inw),1,MPI_INT,mynw,tag_se,MPI_COMM_WORLD, &
                   reqs(n_send),ierr)

    n_send = n_send + 1
    call mpi_isend(Depart(ine),1,MPI_INT,myne,tag_sw,MPI_COMM_WORLD, &
                   reqs(n_send),ierr)

    n_send = n_send + 1
    call mpi_isend(Depart(isw),1,MPI_INT,mysw,tag_ne,MPI_COMM_WORLD, &
                   reqs(n_send),ierr)

    n_send = n_send + 1
    call mpi_isend(Depart(ise),1,MPI_INT,myse,tag_nw,MPI_COMM_WORLD, &
                   reqs(n_send),ierr)

    ! make sure that all the non-blocking MPI operations are complete

    call mpi_waitall(16,reqs,status1,ierr)

    ! sanity check: if too many droplets enter the current MPI region
    !               and exceed the number of "holes", stop the program
    !               with an error message

    if ( sum(Arrive) .gt. numHoles ) then
        write(*,*) "Too many new droplets will enter the MPI rank: ", myid
        write(*,*) "nparcelsLocalActive,nparcelsLocal",nparcelsLocalActive,nparcelsLocal
        write(*,*) "Arriving, holes to fill: ", sum(Arrive), numHoles
        write(*,*) "Stop the program ..."
        call stopcm1
    end if

#if 0 
    ! sanity check: if the total number of leaving droplets matches
    !               the total number of entering droploets across 
    !               all the MPI ranks

    n_send = sum(Depart)
    if ( myid == 0 ) then
       call mpi_reduce(MPI_IN_PLACE,n_send,1,MPI_INT, &
                       MPI_SUM,0,MPI_COMM_WORLD,ierr)
    else
       call mpi_reduce(n_send,n_send,1,MPI_INT, &
                       MPI_SUM,0,MPI_COMM_WORLD,ierr)
    end if

    n_recv = sum(Arrive)
    if ( myid == 0 ) then
       call mpi_reduce(MPI_IN_PLACE,n_recv,1,MPI_INT, &
                       MPI_SUM,0,MPI_COMM_WORLD,ierr)
    else
       call mpi_reduce(n_recv,n_recv,1,MPI_INT, &
                       MPI_SUM,0,MPI_COMM_WORLD,ierr)
    end if

    if ( myid == 0 ) then
       if ( n_send .ne. n_recv ) then
          write(*,*) "Unmatched total leaving droplets vs. total entering droplets..."
          write(*,*) "Total leaving droplets: ", n_send
          write(*,*) "Total entering droplets: ", n_recv
          write(*,*) "Stop the program ..."
          call stopcm1
       end if
    end if

#endif

  end subroutine comm_droplet_number

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  subroutine comm_droplet_value(holes_ind,Depart_ind,Depart,Arrive, &
                ptrDepart, ptrArrive, pdata)

    use input, only : myid,nparcelsLocal,npvals,mynw,mysw,myne, &
                      myse,myeast,mywest,mynorth,mysouth,ierr, &
                      pract,nparcelsLocalActive, &
                      timestats,mytime
    use constants, only: undefined_index,neg_huge,num_nn,inorth, &
                         isouth,iwest,ieast,inw,ine,isw,ise
    use mpi

    implicit none

    integer, intent(in), dimension(:) :: holes_ind                       ! location index in "pdata" that can 
                                                                         ! add a new droplet;
                                                                         ! not a fixed-size array but we do 
                                                                         ! not need to know its size
    integer, intent(in), dimension(:) :: Depart_ind                      ! location index in "pdata" of droplets
                                                                         ! that  are departing
    integer, intent(in) :: Depart(num_nn)                                ! Number of droplets that will enter each 
                                                                         ! nearest neighbor at different directions
    integer, intent(in) :: Arrive(num_nn)                                ! Number of droplets that will enter the 
                                                                         ! current MPI region from each nearest 
                                                                         ! neighbor at different directions
    integer, intent(in) :: ptrArrive(num_nn)                             ! location index into holes_ind 
    integer, intent(in) :: ptrDepart(num_nn)                             ! location index into Depart_ind 
    real, intent(inout), dimension(nparcelsLocal,npvals) :: pdata        ! droplet information

    ! Local variables

    real, dimension(:,:), allocatable :: droplet_n1, droplet_s1, droplet_w1, &     ! information of droplets that will leave the 
                                         droplet_e1, droplet_nw1, droplet_ne1, &   ! current MPI region and enter the nearest neighbor;
                                         droplet_sw1, droplet_se1, &               ! 1 means sending array;
                                         droplet_n2, droplet_s2, droplet_w2, &     ! 2 means receiving array
                                         droplet_e2, droplet_nw2, droplet_ne2, &
                                         droplet_sw2, droplet_se2

    integer :: n_idx,s_idx,w_idx,e_idx,nw_idx,ne_idx,sw_idx,se_idx
    integer :: reqs(16),tag(num_nn,2)
    integer :: index_n,index_s,index_w,index_e,index_nw,index_ne, &
               index_sw,index_se,n_recv,n_send
    integer, dimension(mpi_status_size,16) :: status1
    integer :: i,j,k,n,indx
    integer :: numDepart,numArrive,num_holes,last_active
    integer :: num_tmp

    ! initialize some MPI related variables

    index_n  = undefined_index
    index_s  = undefined_index
    index_w  = undefined_index
    index_e  = undefined_index
    index_nw = undefined_index
    index_sw = undefined_index
    index_ne = undefined_index
    index_se = undefined_index

    do i = 1, num_nn
       tag(i,1) = undefined_index
       tag(i,2) = undefined_index
    end do

    numDepart = sum(Depart)
    numArrive = sum(Arrive)

     ! allocate spaces to store the droplets that will leave the current MPI region

     if ( Depart(inorth) .ne. 0 ) allocate(droplet_n1(Depart(inorth),npvals))
     if ( Depart(isouth) .ne. 0 ) allocate(droplet_s1(Depart(isouth),npvals))
     if ( Depart(iwest) .ne. 0 ) allocate(droplet_w1(Depart(iwest),npvals))
     if ( Depart(ieast) .ne. 0 ) allocate(droplet_e1(Depart(ieast),npvals))
     if ( Depart(inw) .ne. 0 ) allocate(droplet_nw1(Depart(inw),npvals))
     if ( Depart(ine) .ne. 0 ) allocate(droplet_ne1(Depart(ine),npvals))
     if ( Depart(isw) .ne. 0 ) allocate(droplet_sw1(Depart(isw),npvals))
     if ( Depart(ise) .ne. 0 ) allocate(droplet_se1(Depart(ise),npvals))

     ! allocate temporary arrays if there are new droplets from the nearest neighbor

     if ( Arrive(inorth) .ne. 0 ) allocate(droplet_n2(Arrive(inorth),npvals))
     if ( Arrive(isouth) .ne. 0 ) allocate(droplet_s2(Arrive(isouth),npvals))
     if ( Arrive(iwest) .ne. 0 ) allocate(droplet_w2(Arrive(iwest),npvals))
     if ( Arrive(ieast) .ne. 0 ) allocate(droplet_e2(Arrive(ieast),npvals))
     if ( Arrive(inw) .ne. 0 ) allocate(droplet_nw2(Arrive(inw),npvals))
     if ( Arrive(ine) .ne. 0 ) allocate(droplet_ne2(Arrive(ine),npvals))
     if ( Arrive(isw) .ne. 0 ) allocate(droplet_sw2(Arrive(isw),npvals))
     if ( Arrive(ise) .ne. 0 ) allocate(droplet_se2(Arrive(ise),npvals))

     !$acc data create(droplet_n1,droplet_s1,droplet_w1,droplet_e1,     &
     !$acc             droplet_nw1,droplet_ne1,droplet_sw1,droplet_se1, &
     !$acc             droplet_n2,droplet_s2,droplet_w2,droplet_e2,     &
     !$acc             droplet_nw2,droplet_ne2,droplet_sw2,droplet_se2)


     ! extract the information of droplets that will leave the current MPI region 
     ! and reset the corresponding "pdata" record to a large negative value
     ! pack up droplet send buffers and mark departed droplets in pdata
     if(Depart(inorth) .gt. 0) &
        call packDroplets(inorth, Depart, ptrDepart, Depart_ind, pdata, droplet_n1)

     if(Depart(isouth) .gt. 0) &
        call packDroplets(isouth, Depart, ptrDepart, Depart_ind, pdata, droplet_s1)

     if(Depart(iwest) .gt. 0) &
        call packDroplets(iwest, Depart, ptrDepart, Depart_ind, pdata, droplet_w1)

     if(Depart(ieast) .gt. 0) &
        call packDroplets(ieast, Depart, ptrDepart, Depart_ind, pdata, droplet_e1)

     if(Depart(inw) .gt. 0) &
        call packDroplets(inw, Depart, ptrDepart, Depart_ind, pdata, droplet_nw1)

     if(Depart(ine) .gt. 0) &
        call packDroplets(ine, Depart, ptrDepart, Depart_ind, pdata, droplet_ne1)

     if(Depart(isw) .gt. 0) &
        call packDroplets(isw, Depart,  ptrDepart, Depart_ind, pdata, droplet_sw1)

     if(Depart(ise) .gt. 0) &
        call packDroplets(ise, Depart, ptrDepart, Depart_ind, pdata, droplet_se1)



     ! initiate the MPI non-blocking receive interface to 
     ! obtain new droplet information from the nearest neighbors

     n_recv   = 0
     if ( Arrive(inorth) .ne. 0 ) then
       n_recv = n_recv + 1
       tag(inorth,2)  = 1000 + Arrive(inorth)
       !$acc host_data use_device(droplet_n2)
       call mpi_irecv(droplet_n2,Arrive(inorth)*npvals,MPI_REAL, &
                      mynorth,tag(inorth,2),MPI_COMM_WORLD,reqs(n_recv),ierr)
       !$acc end host_data
       index_n = n_recv
     end if
     if ( Arrive(isouth) .ne. 0 ) then
       n_recv = n_recv + 1
       tag(isouth,2)  = 1000 + Arrive(isouth)
       !$acc host_data use_device(droplet_s2)
       call mpi_irecv(droplet_s2,Arrive(isouth)*npvals,MPI_REAL, &
                      mysouth,tag(isouth,2),MPI_COMM_WORLD,reqs(n_recv),ierr)
       !$acc end host_data
       index_s = n_recv
     end if
     if ( Arrive(iwest) .ne. 0 ) then
       n_recv = n_recv + 1
       tag(iwest,2)  = 1000 + Arrive(iwest)
       !$acc host_data use_device(droplet_w2)
       call mpi_irecv(droplet_w2,Arrive(iwest)*npvals,MPI_REAL, &
                      mywest,tag(iwest,2),MPI_COMM_WORLD,reqs(n_recv),ierr)
       !$acc end host_data
       index_w = n_recv
     end if
     if ( Arrive(ieast) .ne. 0 ) then
       n_recv = n_recv + 1
       tag(ieast,2)  = 1000 + Arrive(ieast)
       !$acc host_data use_device(droplet_e2)
       call mpi_irecv(droplet_e2,Arrive(ieast)*npvals,MPI_REAL, &
                      myeast,tag(ieast,2),MPI_COMM_WORLD,reqs(n_recv),ierr)
       !$acc end host_data
       index_e = n_recv
     end if
     if ( Arrive(inw) .ne. 0 ) then
       n_recv = n_recv + 1
       tag(inw,2) = 1000 + Arrive(inw)
       !$acc host_data use_device(droplet_nw2)
       call mpi_irecv(droplet_nw2,Arrive(inw)*npvals,MPI_REAL, &
                      mynw,tag(inw,2),MPI_COMM_WORLD,reqs(n_recv),ierr)
       !$acc end host_data
       index_nw = n_recv
     end if
     if ( Arrive(ine) .ne. 0 ) then
       n_recv = n_recv + 1
       tag(ine,2) = 1000 + Arrive(ine)
       !$acc host_data use_device(droplet_ne2)
       call mpi_irecv(droplet_ne2,Arrive(ine)*npvals,MPI_REAL, &
                      myne,tag(ine,2),MPI_COMM_WORLD,reqs(n_recv),ierr)
       !$acc end host_data
       index_ne = n_recv
     end if
     if ( Arrive(isw) .ne. 0 ) then
       n_recv = n_recv + 1
       tag(isw,2) = 1000 + Arrive(isw)
       !$acc host_data use_device(droplet_sw2)
       call mpi_irecv(droplet_sw2,Arrive(isw)*npvals,MPI_REAL, &
                      mysw,tag(isw,2),MPI_COMM_WORLD,reqs(n_recv),ierr)
       !$acc end host_data
       index_sw = n_recv
     end if
     if ( Arrive(ise) .ne. 0 ) then
       n_recv = n_recv + 1
       tag(ise,2) = 1000 + Arrive(ise)
       !$acc host_data use_device(droplet_se2)
       call mpi_irecv(droplet_se2,Arrive(ise)*npvals,MPI_REAL, &
                      myse,tag(ise,2),MPI_COMM_WORLD,reqs(n_recv),ierr)
       !$acc end host_data
       index_se = n_recv
     end if

     ! initiate the MPI non-blocking send interface to send the information of
     ! droplets that leave the current MPI region to the nearest neighbors 

     n_send = 8
     if ( Depart(inorth) .ne. 0 ) then
       n_send = n_send + 1
       tag(isouth,1) = 1000 + Depart(inorth)
       !$acc host_data use_device(droplet_n1)
       call mpi_isend(droplet_n1,Depart(inorth)*npvals,MPI_REAL, &
                      mynorth,tag(isouth,1),MPI_COMM_WORLD,reqs(n_send),ierr)
       !$acc end host_data
     end if
     if ( Depart(isouth) .ne. 0 ) then
       n_send = n_send + 1
       tag(inorth,1) = 1000 + Depart(isouth)
       !$acc host_data use_device(droplet_s1)
       call mpi_isend(droplet_s1,Depart(isouth)*npvals,MPI_REAL, &
                      mysouth,tag(inorth,1),MPI_COMM_WORLD,reqs(n_send),ierr)
       !$acc end host_data
     end if
     if ( Depart(iwest) .ne. 0 ) then
       n_send = n_send + 1
       tag(ieast,1) = 1000 + Depart(iwest)
       !$acc host_data use_device(droplet_w1)
       call mpi_isend(droplet_w1,Depart(iwest)*npvals,MPI_REAL, &
                      mywest,tag(ieast,1),MPI_COMM_WORLD,reqs(n_send),ierr)
       !$acc end host_data
     end if
     if ( Depart(ieast) .ne. 0 ) then
       n_send = n_send + 1
       tag(iwest,1) = 1000 + Depart(ieast)
       !$acc host_data use_device(droplet_e1)
       call mpi_isend(droplet_e1,Depart(ieast)*npvals,MPI_REAL, &
                      myeast,tag(iwest,1),MPI_COMM_WORLD,reqs(n_send),ierr)
       !$acc end host_data
     end if
     if ( Depart(inw) .ne. 0 ) then
       n_send = n_send + 1
       tag(ise,1) = 1000 + Depart(inw)
       !$acc host_data use_device(droplet_nw1)
       call mpi_isend(droplet_nw1,Depart(inw)*npvals,MPI_REAL, &
                      mynw,tag(ise,1),MPI_COMM_WORLD,reqs(n_send),ierr)
       !$acc end host_data
     end if
     if ( Depart(ine) .ne. 0 ) then
       n_send = n_send + 1
       tag(isw,1) = 1000 + Depart(ine)
       !$acc host_data use_device(droplet_ne1)
       call mpi_isend(droplet_ne1,Depart(ine)*npvals,MPI_REAL, &
                      myne,tag(isw,1),MPI_COMM_WORLD,reqs(n_send),ierr)
       !$acc end host_data
     end if
     if ( Depart(isw) .ne. 0 ) then
       n_send = n_send + 1
       tag(ine,1) = 1000 + Depart(isw)
       !$acc host_data use_device(droplet_sw1)
       call mpi_isend(droplet_sw1,Depart(isw)*npvals,MPI_REAL, &
                      mysw,tag(ine,1),MPI_COMM_WORLD,reqs(n_send),ierr)
       !$acc end host_data
     end if
     if ( Depart(ise) .ne. 0 ) then
       n_send = n_send + 1
       tag(inw,1) = 1000 + Depart(ise)
       !$acc host_data use_device(droplet_se1)
       call mpi_isend(droplet_se1,Depart(ise)*npvals,MPI_REAL, &
                      myse,tag(inw,1),MPI_COMM_WORLD,reqs(n_send),ierr)
       !$acc end host_data
     end if

     ! update the "holes" in "pdata" with the new droplets that
     ! enter the current MPI region from the nearest neighbors
     n = 1
     k = 1
     do while( n .le. n_recv )
       call mpi_waitany(n_recv,reqs(1:n_recv),indx,MPI_STATUS_IGNORE,ierr)
       if ( indx .eq. index_n ) then
          call update_new_droplet ( inorth, Arrive, ptrArrive, holes_ind, pdata, droplet_n2)
       else if ( indx .eq. index_s ) then
          call update_new_droplet ( isouth, Arrive, ptrArrive, holes_ind, pdata, droplet_s2)
       else if ( indx .eq. index_w ) then
          call update_new_droplet ( iwest, Arrive, ptrArrive, holes_ind, pdata, droplet_w2)
       else if ( indx .eq. index_e ) then
          call update_new_droplet ( ieast, Arrive, ptrArrive, holes_ind, pdata, droplet_e2)
       else if ( indx .eq. index_nw ) then
          call update_new_droplet ( inw, Arrive, ptrArrive, holes_ind, pdata, droplet_nw2)
       else if ( indx .eq. index_ne ) then
          call update_new_droplet ( ine, Arrive, ptrArrive, holes_ind, pdata, droplet_ne2)
       else if ( indx .eq. index_sw ) then
          call update_new_droplet ( isw, Arrive, ptrArrive, holes_ind, pdata, droplet_sw2)
       else if ( indx .eq. index_se ) then
          call update_new_droplet ( ise, Arrive, ptrArrive, holes_ind, pdata, droplet_se2)
       else
          write(*,*) indx," is not a nearest neighbor for myid = ",myid
          write(*,*) "Stop the program..."
          call stopcm1
       end if
       n = n + 1
     end do

    ! make sure that all the non-blocking MPI send operations are complete

    n_send = n_send - 8
    if ( n_send .ge. 1 ) then
       call mpi_waitall(n_send,reqs(9:9+n_send-1),status1(1:mpi_status_size,9:9+n_send-1),ierr)
    end if


     ! free up the memory for temporary variables
     !$acc end data

     if (allocated(droplet_n1))  deallocate(droplet_n1)
     if (allocated(droplet_n2))  deallocate(droplet_n2)
     if (allocated(droplet_s1))  deallocate(droplet_s1)
     if (allocated(droplet_s2))  deallocate(droplet_s2)
     if (allocated(droplet_w1))  deallocate(droplet_w1)
     if (allocated(droplet_w2))  deallocate(droplet_w2)
     if (allocated(droplet_e1))  deallocate(droplet_e1)
     if (allocated(droplet_e2))  deallocate(droplet_e2)
     if (allocated(droplet_ne1)) deallocate(droplet_ne1)
     if (allocated(droplet_ne2)) deallocate(droplet_ne2)
     if (allocated(droplet_nw1)) deallocate(droplet_nw1)
     if (allocated(droplet_nw2)) deallocate(droplet_nw2)
     if (allocated(droplet_se1)) deallocate(droplet_se1)
     if (allocated(droplet_se2)) deallocate(droplet_se2)
     if (allocated(droplet_sw1)) deallocate(droplet_sw1)
     if (allocated(droplet_sw2)) deallocate(droplet_sw2)


  end subroutine comm_droplet_value

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine makeContiguous(Depart,Arrive,num_fallout,holes_ind,backfill_ind,pdata,pdata_locind)
      ! =================================================================
      ! The GPU-enabled makeCongiguous subroutine will fill in "holes" in
      ! the pdata array such that the all active droplets are
      ! contained in the first part of the "pdata" array. This
      ! is achieved by moving droplets from the end of the "pdata" 
      ! array to fill in any remaining holes.
      ! =================================================================

      use input, only : myid,nparcelsLocalActive,nparcelsLocal,npvals,pract
      use constants, only : num_nn,undefined_index,neg_huge

      integer, intent(in) :: Depart(num_nn)                              ! Number of droplets that will enter each
                                                                         ! nearest neighbor at different directions
      integer, intent(in) :: Arrive(num_nn)                              ! Number of droplets that will enter the
                                                                         ! current MPI region from each nearest
      integer, intent(in), dimension(:) :: holes_ind                     ! location index in "pdata" that can 
                                                                         ! add a new droplet;
                                                                         ! not a fixed-size array but we do 
                                                                         ! not need to know its size
                                                                         ! neighbor at different directions
      integer, intent(in), dimension(:) :: backfill_ind
      integer, intent(in) :: num_fallout
      real, intent(inout), dimension(nparcelsLocal,npvals) :: pdata      ! droplet information
      integer, intent(inout), dimension(nparcelsLocal,3) :: pdata_locind ! x/y/z location index of each droplet

      integer :: num_tmp,numArrive,numDepart,num_holes,indx,k,n
      integer :: iv
      logical, parameter :: Debug=.false.

      !simplified code that fills the holes
      numDepart = SUM(Depart)
      numArrive = SUM(Arrive)
      num_tmp = numDepart+num_fallout
      if ( numArrive .lt. num_tmp ) then
        num_holes = num_tmp - numArrive             ! number of remaining "holes" we need to fill
                                                    ! to make the active dropelts in
                                                    ! the "pdata" contiguous

        !$acc parallel loop gang vector default(present)
        do n=1,num_holes
           indx = holes_ind(numArrive+n)            ! index of hole to be filled
           k    = backfill_ind(n)                   ! index of droplet to move into hole
           if(indx .lt. k) then
             if (Debug) then
               write(*,*) ' backfilled droplet # ',indx,' with droplet # ',k
               call flush(6)
             endif
             !$acc loop seq
             do iv=1,npvals
               pdata(indx,iv) = pdata(k,iv)
               pdata(k,iv)    = neg_huge
             enddo
             pdata_locind(indx,1) = pdata_locind(k,1)
             pdata_locind(indx,2) = pdata_locind(k,2)
             pdata_locind(indx,3) = pdata_locind(k,3)
             pdata_locind(k,1) = undefined_index
             pdata_locind(k,2) = undefined_index
             pdata_locind(k,3) = undefined_index
           endif
        enddo
      end if


      end subroutine makeContiguous

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  subroutine packDroplets( dir_idx, Depart, ptrDepart, Depart_ind, &
                                      pdata, droplet_buffer )

  ! this subroutine is used to:
  !      1. save the information of a droplet in "pdata" that will 
  !         leave the current MPI region to a buffer
  !      2. reset the corresponding slot in "pdata" to a negative value
  !      3. record this "hole" index in "pdata"

  use input, only : npvals,nparcelsLocal
  use constants, only: undefined_index, neg_huge

  implicit none

  integer, intent(in)    :: dir_idx                              ! neighbor's directional index 
  integer, intent(in), dimension(:) :: Depart                    ! number of droplets departing 
  integer, intent(in), dimension(:) :: ptrDepart                 ! Pointer into the departure index "Depart_ind"
  integer, intent(in), dimension(:) :: Depart_ind                !  Departure indices
  real, intent(inout), dimension(nparcelsLocal,npvals) :: pdata  ! the "pdata" with droplet information
  real, intent(out), dimension(:,:) :: droplet_buffer            ! array to save the information of leaving droplet
  integer :: i,iv,numDepart,np,iptr

     iptr = ptrDepart(dir_idx)
     numDepart   = Depart(dir_idx)
     !$acc parallel loop gang vector collapse(2) default(present)
     do i=1,numDepart
       do iv=1,npvals
         np = Depart_ind(iptr+i-1)
         droplet_buffer(i,iv) = pdata(np,iv)
         pdata(np,iv) = neg_huge
       enddo
     enddo

  end subroutine packDroplets
#endif

    subroutine setupIndexPointers(Cnt,ptrCnt)

      use constants, only: num_nn,inorth,isouth, &
                         iwest,ieast,inw,ine,isw,ise
      implicit none

      integer, intent(in) :: Cnt(num_nn) 
      integer, intent(inout) :: ptrCnt(num_nn)

      ! local variables
      integer :: n_idx,s_idx,e_idx,w_idx,nw_idx,ne_idx,sw_idx,se_idx


      ptrCnt(inorth) = 1
      ptrCnt(isouth) = ptrCnt(inorth) + Cnt(inorth)
      ptrCnt(iwest)  = ptrCnt(isouth) + Cnt(isouth)
      ptrCnt(ieast)  = ptrCnt(iwest)  + Cnt(iwest)
      ptrCnt(inw)    = ptrCnt(ieast)  + Cnt(ieast)
      ptrCnt(ine)    = ptrCnt(inw)    + Cnt(inw)
      ptrCnt(isw)    = ptrCnt(ine)    + Cnt(ine)
      ptrCnt(ise)    = ptrCnt(isw)    + Cnt(isw)

    end subroutine setupIndexPointers
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    subroutine setupDepartDroplet(Depart_ind,ptrDepart,Depart,pdata_neighbor)

      use input, only : mynw,mysw,myne,myse,myeast,mywest,mynorth,mysouth, &
        nparcelsLocal,nparcelsLocalActive,myid
      use constants, only: num_nn,inorth,isouth, &
                         iwest,ieast,inw,ine,isw,ise,undefined_index
#ifdef _OPENACC
      use cutensorex
#endif

      integer, intent(inout) :: ptrDepart(num_nn)
      integer, intent(in) :: Depart(num_nn)

      integer, intent(inout) :: Depart_ind(:)
      integer, intent(in) :: pdata_neighbor(:)

      ! local variables
      integer :: n_idx,s_idx,e_idx,w_idx,nw_idx,ne_idx,sw_idx,se_idx
      integer :: np

      integer :: numDepart,numerrors
      logical, parameter :: verbose = .false.

      numDepart=SUM(Depart)

      n_idx  = ptrDepart(inorth)
      s_idx  = ptrDepart(isouth)
      w_idx  = ptrDepart(iwest)
      e_idx  = ptrDepart(ieast)
      nw_idx = ptrDepart(inw)
      ne_idx = ptrDepart(ine)
      sw_idx = ptrDepart(isw)
      se_idx = ptrDepart(ise)
      if(verbose .and. (numDepart.gt.0)) then
         write(*,110) myid,'Depart:',Depart
         write(*,110) myid,'ptrDepart:',ptrDepart
      endif

#if (defined(_OPENACC) && defined(_COUNTPREFIX))

      if(Depart(inorth).gt.0) then
         call GatherIndices(inorth,ptrDepart,pdata_neighbor,Depart_ind)
      endif

      if(Depart(isouth).gt.0) then
         call GatherIndices(isouth,ptrDepart,pdata_neighbor,Depart_ind)
      endif

      if(Depart(iwest).gt.0) then
         call GatherIndices(iwest,ptrDepart,pdata_neighbor,Depart_ind)
      endif

      if(Depart(ieast).gt.0) then
         call GatherIndices(ieast,ptrDepart,pdata_neighbor,Depart_ind)
      endif

      if(Depart(inw).gt.0) then
         call GatherIndices(inw,ptrDepart,pdata_neighbor,Depart_ind)
      endif

      if(Depart(ine).gt.0) then
         call GatherIndices(ine,ptrDepart,pdata_neighbor,Depart_ind)
      endif

      if(Depart(isw).gt.0) then
         call GatherIndices(isw,ptrDepart,pdata_neighbor,Depart_ind)
      endif

      if(Depart(ise).gt.0) then
         call GatherIndices(ise,ptrDepart,pdata_neighbor,Depart_ind)
      endif

#else

      !$acc update host(pdata_neighbor)
      n_idx  = ptrDepart(inorth)
      s_idx  = ptrDepart(isouth)
      e_idx  = ptrDepart(ieast)
      w_idx  = ptrDepart(iwest)
      nw_idx = ptrDepart(inw)
      ne_idx = ptrDepart(ine)
      sw_idx = ptrDepart(isw)
      se_idx = ptrDepart(ise)
      do np=1, nparcelsLocalActive
        if(pdata_neighbor(np) .eq. inorth) then
          Depart_ind(n_idx) = np
          n_idx=n_idx+1
        else if(pdata_neighbor(np) .eq. isouth) then
          Depart_ind(s_idx) = np
          s_idx=s_idx+1
        else if(pdata_neighbor(np) .eq. iwest) then
          Depart_ind(w_idx) = np
          w_idx=w_idx+1
        else if(pdata_neighbor(np) .eq. ieast) then
          Depart_ind(e_idx) = np
          e_idx=e_idx+1
        else if(pdata_neighbor(np) .eq. inw) then
          Depart_ind(nw_idx) = np
          nw_idx=nw_idx+1
        else if(pdata_neighbor(np) .eq. ine) then
          Depart_ind(ne_idx) = np
          ne_idx=ne_idx+1
        else if(pdata_neighbor(np) .eq. isw) then
          Depart_ind(sw_idx) = np
          sw_idx=sw_idx+1
        else if(pdata_neighbor(np) .eq. ise) then
          Depart_ind(se_idx) = np
          se_idx=se_idx+1
        endif
      enddo
      !$acc update device(Depart_ind)
#endif

 110     format('myid: ',i4,' ',a10,' ',8(i8,1x))

     end subroutine setupDepartDroplet

#if (defined(_COUNTPREFIX) && defined(_OPENACC))
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

     subroutine GatherIndices(dir,ptrDepart,pdata_neighbor,Depart_ind)

        use constants, only : num_nn
        use input, only : nparcelsLocal, nparcelsLocalActive
#ifdef _OPENACC
        use cutensorex
#endif

        integer, intent(in)                             :: dir
        integer, intent(in), dimension(num_nn)          :: ptrDepart
        integer, intent(in), dimension(nparcelsLocal)   :: pdata_neighbor
        integer, intent(inout), dimension(:)            :: Depart_ind

        ! local variables
        integer :: np
        integer :: ptrdir
        integer, dimension(nparcelsLocal) :: Cidx

        !$acc data create(Cidx) present(pdata_neighbor)

        ptrdir = ptrDepart(dir)
        Cidx = count_prefix(mask=(pdata_neighbor.eq.dir),exclusive=.true.)

        !$acc parallel loop gang vector default(present)
        do np=1,nparcelsLocalActive-1
           if(Cidx(np) .ne. Cidx(np+1)) then
             Depart_ind(ptrdir+Cidx(np))=np
           endif
        enddo
        !$acc end parallel

        ! special treament for the last value in the array
        !$acc kernels default(present) 
        if(pdata_neighbor(nparcelsLocalActive).eq.dir) then
           Depart_ind(ptrdir+Cidx(nparcelsLocalActive)) = nparcelsLocalActive
        endif
        !$acc end kernels

        !$acc end data

      end subroutine GatherIndices
#endif

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

#ifdef MPI
     subroutine update_new_droplet ( dir_idx, Arrive, ptrArrive, holes_ind, pdata, droplet_buffer)

     ! this subroutine is used to update the new droplet information
     ! to the "pdata" array on a particular MPI rank after the MPI communication

     use input, only : npvals,nparcelsLocal

     implicit none

     integer, intent(in), dimension(:) :: holes_ind                  ! array to store the location index of "holes" in "pdata"
     integer, intent(in) :: dir_idx                               !  neighbor's directional index
     integer, intent(in), dimension(:) :: Arrive                    ! number of droplets departing
     integer, intent(in), dimension(:) :: ptrArrive                 ! Pointer into the departure index "Depart_ind"

     real, intent(in), dimension(:,:) :: droplet_buffer              ! array to save the information of new droplets
     real, intent(inout), dimension(nparcelsLocal,npvals) :: pdata   ! the "pdata" with droplet information

     integer :: numArrive,i,iptr,iv,ptr

     ptr = ptrArrive(dir_idx)
     numArrive = Arrive(dir_idx)
     !$acc parallel loop gang vector collapse(2) default(present)
     do i = 1, numArrive
       do iv=1,npvals
         iptr = holes_ind(ptr+i-1)
         pdata(iptr,iv) = droplet_buffer(i,iv)
       end do
     end do

     end subroutine update_new_droplet

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

     subroutine CollectDropCount ( nparcels_per_mpi, nparcelstot )

     ! This subroutine will collect the active droplet count from each
     ! MPI rank and assemble them into a single array for later usage.
     ! It works for > 2 billions droplets that can not be represented 
     ! by a normal integer type of length

     use constants, only : i8
     use input, only : myid,numprocs,ierr,nparcelsLocalActive
     use mpi

     implicit none

     integer, intent(out) :: nparcels_per_mpi(numprocs)   ! number of active droplets per MPI rank     
     integer(i8), intent(out) :: nparcelstot              ! total number of active droplets
                                                          ! over the whole domain

     ! Local variable
     integer :: i,tag

     ! Collect the number of active droplets from each MPI rank
     if ( myid .eq. 0 ) then
        nparcels_per_mpi(1) = nparcelsLocalActive
        nparcelstot = nparcels_per_mpi(1)
        do i = 2, numprocs
           tag = 98765 + i - 1
           call mpi_recv(nparcels_per_mpi(i),1,MPI_INTEGER,i-1, &
                         tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
           nparcelstot = nparcelstot + nparcels_per_mpi(i)
        end do
     else
        tag = 98765 + myid
        call mpi_send(nparcelsLocalActive,1,MPI_INTEGER, &
                      0,tag,MPI_COMM_WORLD,ierr)
     end if

     end subroutine CollectDropCount

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

     subroutine CollectDropInfo ( pdata, rbuf, nparcels_per_mpi, nparcelsout )

     ! This subroutine will collect the "pdata" information from each
     ! MPI rank and assemble them into a single big array "rbuf" for 
     ! writing out later.
     ! It does not work for > 2 billions droplets since we can not write
     ! out restart / pdata files with the current I/O interface

     use input, only : myid,numprocs,ierr,nparcelsLocalActive, &
                       nparcelsLocal,npvals,i8
     use mpi

     implicit none

     real, intent(in), dimension(nparcelsLocal,npvals) :: pdata
     integer, intent(in) :: nparcels_per_mpi(numprocs)   ! number of active droplets per MPI rank     
     integer, intent(in) :: nparcelsout                  ! total number of active droplets
                                                         ! to be written out to a binary file
     real, intent(out) :: rbuf(:,:)                      ! buffer to store all the active droplet 
                                                         ! information over the domain

     ! Local variable
     integer :: i,n,np,iv,tag
     integer :: beg_loc,end_loc
     integer(i8) :: tmp
     real, dimension(:,:), allocatable :: tmp_buf

     ! Collect the information of active droplets from each MPI rank
     if ( myid .eq. 0 ) then
        end_loc = min(nparcelsout,nparcelsLocalActive)
        do iv = 1, npvals
           do np = 1, end_loc 
              rbuf(np,iv) = pdata(np,iv)
           end do
        end do
        tmp = 1 
        do i = 2, numprocs
           allocate( tmp_buf(nparcels_per_mpi(i),npvals) )
           tag = 98765 + i - 1
           call mpi_recv(tmp_buf,nparcels_per_mpi(i)*npvals, &
                         MPI_REAL,i-1,tag,MPI_COMM_WORLD, &
                         MPI_STATUS_IGNORE,ierr)
           tmp = tmp + nparcels_per_mpi(i-1)
           ! We only save first "nparcelsout" droplets for I/O
           if ( tmp <= nparcelsout ) then
              beg_loc = tmp
              end_loc = min( beg_loc+nparcels_per_mpi(i)-1, nparcelsout )
              do iv = 1, npvals
                 do np = beg_loc, end_loc 
                    n = np - beg_loc + 1
                    rbuf(np,iv) = tmp_buf(n,iv)
                 end do
              end do
           end if
           deallocate(tmp_buf)
        end do
     else
        tag = 98765 + myid
        allocate( tmp_buf(nparcelsLocalActive,npvals) )
        do iv = 1, npvals
           do np = 1, nparcelsLocalActive
              tmp_buf(np,iv) = pdata(np,iv)
           end do
        end do
        call mpi_send(tmp_buf,nparcelsLocalActive*npvals, &
                      MPI_REAL,0,tag,MPI_COMM_WORLD,ierr)
        deallocate(tmp_buf)
     end if

     end subroutine CollectDropInfo

#endif

end module comm_droplet_module
