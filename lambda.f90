program lambda
use constant
implicit none

integer :: nx,ny,i,j
real(kind=8) ::kx,ky,dEfs,dT,T,Efs
real(kind=8),allocatable :: obd(:,:)
allocate(Es(N,N,8),Ome_ks(N,N,8))
v1=2
v2=v1
t1=1.5
t2=t1
yita1=-1
yita2=1
m1=0.1
m2=m1
gamma=0.05
E1=0.02
E2=-0.08
K1=0.1*pi
K2=0.15*pi
vx=0.2
vy=0.0

sx(1,1)=(0,0)
sx(1,2)=(1,0)
sx(2,1)=(1,0)
sx(2,2)=(0,0)

sz(1,1)=(1,0)
sz(1,2)=(0,0)
sz(2,1)=(0,0)
sz(2,2)=(-1,0)

sy(1,1)=(0,0)
sy(1,2)=(0,-1)
sy(2,1)=(0,1)
sy(2,2)=(0,0)

s0(1,1)=(1,0)
s0(1,2)=(0,0)
s0(2,1)=(0,0)
s0(2,2)=(1,0)
I0=(0,0)

!call sth to get Es_reduce and Ome_reduce  just need to return a single value
!parallel for kx or ky?
!Then call mpi_reduce to combine them into Es(N,N,8) and Ome_ks(N,N,8)
!Es,Ome=pp(kx,ky)
!Es=np.array(np.reshape(Es,[N,N,8]))
!Ome_ks=np.array(np.reshape(Ome,[N,N,8]))

!kx=np.linspace(-0.3*np.pi,0*np.pi-0.001,N) # (-np.pi,0)
!ky=np.linspace(-1*np.pi,0*np.pi,N)

!do mm=mpime+1,N,nproc
!   kx=-0.3*pi+(mm-1.d0)/(N-1)*(0-(-0.3*pi))
!   ky=-1*pi
!   while (ky .lt. 0) do
!      call pp
!   enddo
!enddo

dkx=pi/(N-1)
dky=pi/(N-1)
do nx=1,N
   kx=-0.5*pi+(nx-1)*dkx
   !ky=-1*pi
   do ny=1,N
      ky=-0.5*pi+(ny-1)*dky
      call calc(kx,ky,Es(nx,ny,:),Ome_ks(nx,ny,:)) 
   enddo
enddo
write(*,*) Es(1,1,:)
!kx=-0.5*pi
!open(1,file="Evskx",status="replace")
!do nx=1,N
!write(*,"(8F20.10)") Ome_ks(nx,2,:)
!enddo
!close(1)

! call MPI_REDUCE

!dT=(100-0.5)/150  ! dT=(Tmax-Tmin)/(NT-1)
!dEfs=(0.5-(-0.5))/100 ! dEfs=(Efsmax-Efsmin)/(NEfs-1)
allocate(obd(1,1))
!allocate(obd(151,101))

!do i=1,151 ! for T
!   T=0.5+(i-1)*dT  ! Tmin+(i-1)*dT
    T=300
!   do j=1,101 ! for Efs
!      Efs=-0.5+(j-1)*dEfs ! 
    Efs=0
      call calc_2(T,Efs,obd(1,1)) 
!   enddo
!enddo
!write(*,*) obd(1,1)
end program

subroutine BerryCur_k(evac,eval,Ome_kx)
use constant
implicit none
real(kind=8)    :: eval(8),Ome_kx(8),temp,ome_kx1
complex(kind=8) :: V(8,8),evac(8,8) !,eval
complex(kind=8) :: dHd1_kx(2,2),dHd1_ky(2,2),dHd2_kx(2,2), dHd2_ky(2,2),dP_kx(2,2),dP_ky(2,2)
complex(kind=8) :: dH_kx(8,8),dH_ky(8,8)
complex(kind=8) :: nomi1,nomi2,nomi3,nomi4
integer :: i,j

    V=evac

    dHd1_kx=t1*s0+v1*yita1*sy
    dHd1_ky=v1*sx

    dHd2_kx=t2*s0+v2*yita2*sy
    dHd2_ky=v2*sx

    dP_kx=vx*sz
    dP_ky=(0,-1)*vy*s0

    dH_kx(1,1)=dHd1_kx(1,1)
    dH_kx(1,2)=dHd1_kx(1,2)
    dH_kx(1,3)=dP_kx(1,1)
    dH_kx(1,4)=dP_kx(1,2)
    dH_kx(1,5)=I0(1,1)
    dH_kx(1,6)=I0(1,2)
    dH_kx(1,7)=I0(1,1)
    dH_kx(1,8)=I0(1,2)

    dH_kx(2,1)=dHd1_kx(2,1)
    dH_kx(2,2)=dHd1_kx(2,2)
    dH_kx(2,3)=dP_kx(2,1)
    dH_kx(2,4)=dP_kx(2,2)
    dH_kx(2,5)=I0(2,1)
    dH_kx(2,6)=I0(2,2)
    dH_kx(2,7)=I0(2,1)
    dH_kx(2,8)=I0(2,2)


    dH_kx(3,1)=dP_kx(1,1)
    dH_kx(3,2)=dP_kx(1,2)
    dH_kx(3,3)=dHd1_kx(1,1)
    dH_kx(3,4)=dHd1_kx(1,2)
    dH_kx(3,5)=I0(1,1)
    dH_kx(3,6)=I0(1,2)
    dH_kx(3,7)=I0(1,1)
    dH_kx(3,8)=I0(1,2)

    dH_kx(4,1)=dP_kx(2,1)
    dH_kx(4,2)=dP_kx(2,2)
    dH_kx(4,3)=dHd1_kx(2,1)
    dH_kx(4,4)=dHd1_kx(2,2)
    dH_kx(4,5)=I0(2,1)
    dH_kx(4,6)=I0(2,2)
    dH_kx(4,7)=I0(2,1)
    dH_kx(4,8)=I0(2,2)

    dH_kx(5,1)=I0(1,1)
    dH_kx(5,2)=I0(1,2)
    dH_kx(5,3)=I0(1,1)
    dH_kx(5,4)=I0(1,2)
    dH_kx(5,5)=dHd2_kx(1,1)
    dH_kx(5,6)=dHd2_kx(1,2)
    dH_kx(5,7)=dP_kx(1,1)
    dH_kx(5,8)=dP_kx(1,2)

    dH_kx(6,1)=I0(2,1)
    dH_kx(6,2)=I0(2,2)
    dH_kx(6,3)=I0(2,1)
    dH_kx(6,4)=I0(2,2)
    dH_kx(6,5)=dHd2_kx(2,1)
    dH_kx(6,6)=dHd2_kx(2,2)
    dH_kx(6,7)=dP_kx(2,1)
    dH_kx(6,8)=dP_kx(2,2)

    dH_kx(7,1)=I0(1,1)
    dH_kx(7,2)=I0(1,2)
    dH_kx(7,3)=I0(1,1)
    dH_kx(7,4)=I0(1,2)
    dH_kx(7,5)=dP_kx(1,1)
    dH_kx(7,6)=dP_kx(1,2)
    dH_kx(7,7)=dHd2_kx(1,1)
    dH_kx(7,8)=dHd2_kx(1,2)

    dH_kx(8,1)=I0(2,1)
    dH_kx(8,2)=I0(2,2)
    dH_kx(8,3)=I0(2,1)
    dH_kx(8,4)=I0(2,2)
    dH_kx(8,5)=dP_kx(2,1)
    dH_kx(8,6)=dP_kx(2,2)
    dH_kx(8,7)=dHd2_kx(2,1)
    dH_kx(8,8)=dHd2_kx(2,2)


    dH_ky(1,1)=dHd1_ky(1,1)
    dH_ky(1,2)=dHd1_ky(1,2)
    dH_ky(1,3)=dP_ky(1,1)
    dH_ky(1,4)=dP_ky(1,2)
    dH_ky(1,5)=I0(1,1)
    dH_ky(1,6)=I0(1,2)
    dH_ky(1,7)=I0(1,1)
    dH_ky(1,8)=I0(1,2)

    dH_ky(2,1)=dHd1_ky(2,1)
    dH_ky(2,2)=dHd1_ky(2,2)
    dH_ky(2,3)=dP_ky(2,1)
    dH_ky(2,4)=dP_ky(2,2)
    dH_ky(2,5)=I0(2,1)
    dH_ky(2,6)=I0(2,2)
    dH_ky(2,7)=I0(2,1)
    dH_ky(2,8)=I0(2,2)

    dH_ky(3,1)=dP_ky(1,1)
    dH_ky(3,2)=dP_ky(1,2)
    dH_ky(3,3)=dHd1_ky(1,1)
    dH_ky(3,4)=dHd1_ky(1,2)
    dH_ky(3,5)=I0(1,1)
    dH_ky(3,6)=I0(1,2)
    dH_ky(3,7)=I0(1,1)
    dH_ky(3,8)=I0(1,2)


    dH_ky(4,1)=dP_ky(2,1)
    dH_ky(4,2)=dP_ky(2,2)
    dH_ky(4,3)=dHd1_ky(2,1)
    dH_ky(4,4)=dHd1_ky(2,2)
    dH_ky(4,5)=I0(2,1)
    dH_ky(4,6)=I0(2,2)
    dH_ky(4,7)=I0(2,1)
    dH_ky(4,8)=I0(2,2)

    dH_ky(5,1)=I0(1,1)
    dH_ky(5,2)=I0(1,2)
    dH_ky(5,3)=I0(1,1)
    dH_ky(5,4)=I0(1,2)
    dH_ky(5,5)=dHd2_ky(1,1)
    dH_ky(5,6)=dHd2_ky(1,2)
    dH_ky(5,7)=dP_ky(1,1)
    dH_ky(5,8)=dP_ky(1,2)

    dH_ky(6,1)=I0(2,1)
    dH_ky(6,2)=I0(2,2)
    dH_ky(6,3)=I0(2,1)
    dH_ky(6,4)=I0(2,2)
    dH_ky(6,5)=dHd2_ky(2,1)
    dH_ky(6,6)=dHd2_ky(2,2)
    dH_ky(6,7)=dP_ky(2,1)
    dH_ky(6,8)=dP_ky(2,2)

    dH_ky(7,1)=I0(1,1)
    dH_ky(7,2)=I0(1,2)
    dH_ky(7,3)=I0(1,1)
    dH_ky(7,4)=I0(1,2)
    dH_ky(7,5)=dP_ky(1,1)
    dH_ky(7,6)=dP_ky(1,2)
    dH_ky(7,7)=dHd2_ky(1,1)
    dH_ky(7,8)=dHd2_ky(1,2)

    dH_ky(8,1)=I0(2,1)
    dH_ky(8,2)=I0(2,2)
    dH_ky(8,3)=I0(2,1)
    dH_ky(8,4)=I0(2,2)
    dH_ky(8,5)=dP_ky(2,1)
    dH_ky(8,6)=dP_ky(2,2)
    dH_ky(8,7)=dHd2_ky(2,1)
    dH_ky(8,8)=dHd2_ky(2,2)

!    allocate(Ome_kx(size(eval)))
!    allocate(Ome_kx(8))
    do i =1,8 !size(eval)
        Ome_kx1=0
        do j=1,8 !size(eval) !do j in range(len(eval)):
            if (i .ne.j ) then
                nomi1=sum(matmul(conjg(V(:,i)),dH_kx)*V(:,j))
                nomi2=sum(matmul(conjg(V(:,j)),dH_ky)*V(:,i))
                nomi3=sum(matmul(conjg(V(:,i)),dH_ky)*V(:,j))
                nomi4=sum(matmul(conjg(V(:,j)),dH_kx)*V(:,i))
                temp=(aimag(nomi1*nomi2)-aimag(nomi3*nomi4))/(eval(i)-eval(j))**2
                !nomi_1=np.matmul(np.matmul(np.matrix.getH(V[:,i]),dH_kx),V[:,j]) ! getH complex conjugate
                !nomi_2=np.matmul(np.matmul(np.matrix.getH(V[:,j]),dH_ky),V[:,i])
                !nomi_3=np.matmul(np.matmul(np.matrix.getH(V[:,i]),dH_ky),V[:,j])
                !nomi_4=np.matmul(np.matmul(np.matrix.getH(V[:,j]),dH_kx),V[:,i])
                !temp=(np.imag(nomi_1[0,0]*nomi_2[0,0])-np.imag(nomi_3[0,0]*nomi_4[0,0]))/(eval[i]-eval[j])**2;
                Ome_kx1=Ome_kx1+temp  !! 
             endif
        enddo
        Ome_kx(i)=-1*Ome_kx1
    enddo   

end subroutine BerryCur_k

!function factor(T,E,Ef)
!  implicit none
!  real(kind=8) :: beta,T,E,Ef
!  integer :: ds(2),i,j
!  real(kind=8),allocatable :: n(:,:)

!  beta=1/(T*25.7*10**(-3)/298);
!  ds(1)=size(E,1)
!  ds(2)=size(E,2)

!  allocate(n(ds(1),ds(2)))

!  do i=1,ds(1)
!     do j=1,ds(2)
!        n(i,j)=1/(exp(beta*(E(i,j)-Ef))+1)
!     enddo
!  enddo
!return n
!end function factor

subroutine dfactor(T,E,Ef,n0)
  use constant
  implicit none

  real(kind=8) :: beta
  real(kind=8) :: t,e(N,N), ef, f
  integer :: ds(2),i,j
  real(kind=8) :: n0(N,N)

  beta=1/(T/298*25.7E-3)
!  write(*,*) beta
!  ds(1)=size(E,1)
!  ds(2)=size(E,2)

!  allocate(n(ds(1),ds(2)))
  do i=1,N !size(E,1)!ds(1)
     do j=1,N!size(E,2)!ds(2)
        f=1/(exp(beta*(E(i,j)-Ef))+1)
!        write(*,*) f
        n0(i,j)=f*(1-f)*(-1*beta)
     enddo
  enddo


!    beta=1/(T*25.7*10**(-3)/298);
!    ds=E.shape;
!    n=np.zeros(ds)
!   for i in range(ds[0]):
!        for j in range(ds[1]):
!            f=1/(np.exp(beta*(E[i,j]-Ef))+1);
!            n[i,j]=f*(1-f)*(-beta);

end subroutine dfactor

subroutine calc(kx,ky,es0,Ome_ks0)
  use constant
  implicit none
  real(kind=8) :: ky,kyp,kx,es0(8),Ome_ks0(8),kxp1,kxp2
  complex(kind=8) :: Hd1(2,2),Hd2(2,2),P(2,2),H(8,8),GM(2,2)
  integer :: Lwmax=1000
  integer :: info,Lwork
    real(kind=8), allocatable :: rwork(:)
    complex(kind=8), allocatable :: work(:)
    integer(kind=4) :: nwork=1
allocate(rwork(max(1,3*N-2)))
allocate(work(nwork))
        kyp=ky
        !kx=kx

        kxp1=kx+K1
        kxp2=kx+K2

        Hd1=(E1+t1*kxp1)*s0+v1*(kyp*sx + yita1*kxp1*sy)+m1/2*sz
        Hd2=(E2+t2*kxp2)*s0+v2*(kyp*sx + yita2*kxp2*sy)+m2/2*sz
        P=vx*kx*sz+((0,-1)*vy*kyp)*s0

        GM=gamma*s0

        H(1,1)=Hd1(1,1)
        H(1,2)=Hd1(1,2)
        H(1,3)=P(1,1)
        H(1,4)=P(1,2)
        H(1,5)=I0(1,1)
        H(1,6)=I0(1,2)
        H(1,7)=GM(1,1)
        H(1,8)=GM(1,2)

        H(2,1)=Hd1(2,1)
        H(2,2)=Hd1(2,2)
        H(2,3)=P(2,1)
        H(2,4)=P(2,2)
        H(2,5)=I0(2,1)
        H(2,6)=I0(2,2)
        H(2,7)=GM(2,1)
        H(2,8)=GM(2,2)

        H(3,1)=P(1,1)
        H(3,2)=P(1,2)
        H(3,3)=Hd1(1,1)
        H(3,4)=Hd1(1,2)
        H(3,5)=GM(1,1)
        H(3,6)=GM(1,2)
        H(3,7)=I0(1,1)
        H(3,8)=I0(1,2)

        H(4,1)=P(2,1)
        H(4,2)=P(2,2)
        H(4,3)=Hd1(2,1)
        H(4,4)=Hd1(2,2)
        H(4,5)=GM(2,1)
        H(4,6)=GM(2,2)
        H(4,7)=I0(2,1)
        H(4,8)=I0(2,2)

        H(5,1)=I0(1,1)
        H(5,2)=I0(1,2)
        H(5,3)=GM(1,1)
        H(5,4)=GM(1,2)
        H(5,5)=Hd2(1,1)
        H(5,6)=Hd2(1,2)
        H(5,7)=P(1,1)
        H(5,8)=P(1,2)

        H(6,1)=I0(2,1)
        H(6,2)=I0(2,2)
        H(6,3)=GM(2,1)
        H(6,4)=GM(2,2)
        H(6,5)=Hd2(2,1)
        H(6,6)=Hd2(2,2)
        H(6,7)=P(2,1)
        H(6,8)=P(2,2)

        H(7,1)=GM(1,1)
        H(7,2)=GM(1,2)
        H(7,3)=I0(1,1)
        H(7,4)=I0(1,2)
        H(7,5)=P(1,1)
        H(7,6)=P(1,2)
        H(7,7)=Hd2(1,1)
        H(7,8)=Hd2(1,2)

        H(8,1)=GM(2,1)
        H(8,2)=GM(2,2)
        H(8,3)=I0(2,1)
        H(8,4)=I0(2,2)
        H(8,5)=P(2,1)
        H(8,6)=P(2,2)
        H(8,7)=Hd2(2,1)
        H(8,8)=Hd2(2,2)
       call zheev("V","U",8,H,8,es0,work,-1,rwork,info)

       if(real(work(1)).gt.nwork) then
          nwork=nint(2*real(work(1)))
          deallocate(work)
          allocate(work(nwork))
       end if
       call zheev("V","U",8,H,8,es0,work,nwork,rwork,info)

       call BerryCur_k(H,es0,Ome_ks0)
        !H1=np.hstack((Hd1,P,I0,Gamma))
        !H2=np.hstack((P,Hd1,Gamma,I0))
        !H3=np.hstack((I0,Gamma,Hd2,P))
        !H4=np.hstack((Gamma,I0,P,Hd2))
        !H=np.vstack((H1,H2,H3,H4))

        !!!#### need to call sth. to find eigenvalue and eigenvectors of H
        !!!#### then es0 equals to the eigenvalue
        !!!#### call BerryCur_k(evec,eval,Ome_ks0)

        !eval,evec=la.eigh(H)
        !Es0.append(eval)
        !Ome_ks0.append(BerryCur_k(evec,eval));

        !return Es0,Ome_ks0
end subroutine calc

subroutine calc_2(T0,Ef0,Obd_temp)
use constant
implicit none
integer :: index_bands,i
real(kind=8) :: integraly(N),integralx,integral,temp(N,N),obd_temp,dEs(N,N),Bcur(N,N)
real(kind=8) :: T0,Ef0,df_dE(N,N)
    Obd_temp=0
    do index_bands=1,8 !for index_bands in range(min(Es.shape)):
        !!! dEs call sth to find the gradient of Es(:,:)
        !dEs=np.gradient(Es[:,:,index_bands],hx,axis=0); hx is the dkx
        call velocity(Es(:,:,index_bands),dEs)
!        if (index_bands .eq.1)   then
!        do i =1,N
!        write(*,*) "E"
!        write(*,*) Es(:,1,index_bands)
!        write(*,*) "dEs"
!        write(*,*) dEs(:,1)
!        enddo
!        endif
        
        Bcur=Ome_ks(:,:,index_bands)
        call dfactor(T0,Es(:,:,index_bands),Ef0,df_dE)
        temp=(dEs*Bcur)*df_dE*(Es(:,:,index_bands)-Ef0)**2
!        if (index_bands .eq.1) then
!        do i =1,6
!        write(*,*) temp(i,:)
!        enddo
!        endif

        temp=((dEs*Bcur)*df_dE)*(Es(:,:,index_bands)-Ef0)**2
!        if (index_bands .eq.1) then
!           do i =1,8
!           write(*,*) temp(1,:)
!           enddo
!        endif

        integraly=0.0
        do i =1, N-1
           integraly=integraly+(temp(:,i)+temp(:,i+1))*dky/2 !! used dky or dkx since they are equally spaced
        enddo
        integralx=0.0
        do i = 1,N-1
           integralx=integralx+(integraly(i)+integraly(i+1))*dkx/2
!        temp=-1*np.trapz(np.trapz(temp,ky,axis=1),kx,axis=0)/(max(kx)-min(kx))**2/(T0*25.7*10**(-3)/298)**2
        enddo
        integral=-1*integralx/pi**2/(T0/298*25.7E-3)**2
!        write(*,*) integral
!(max(kx)-min(kx))**2/(T0*25.7*10**(-3)/298)**2
        Obd_temp=Obd_temp+integral
   enddo
end subroutine calc_2

subroutine velocity(E,dE)
use constant
implicit none

real(kind=8) :: E(N,N),dE(N,N),Etmp(N)
integer :: i,j

do i=1,N
   Etmp=E(:,i)
   do j=1,N
     if (j.eq.1) then 
        dE(j,i)=(Etmp(j+1)-Etmp(j))/dkx
     elseif (j.eq.N) then 
         dE(j,i)=(Etmp(j)-Etmp(j-1))/dkx
     else
         dE(j,i)=(Etmp(j+1)-Etmp(j-1))/2/dkx
     endif
   enddo
enddo

end subroutine
