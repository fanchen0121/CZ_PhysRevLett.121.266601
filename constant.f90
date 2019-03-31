Module constant
implicit none
integer :: N=12,v1,v2,yita1,yita2
real(kind=8),allocatable :: Es(:,:,:),Ome_ks_reduce(:,:,:),Ome_ks(:,:,:)
real(kind=8) :: t1,t2,m1,m2,gamma,E1,E2,K1,K2,vx,vy,dkx,dky
real(kind=8),parameter :: pi=3.1415926
complex(kind=8) :: sx(2,2),sz(2,2),s0(2,2),I0(2,2)
complex(kind=8) :: sy(2,2)
integer(kind=4) :: myid,numprocs
end module constant

