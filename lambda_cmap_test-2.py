#import kwant
import numpy as np
#import tinyarray as ty
import matplotlib.pyplot as plt
import scipy.sparse.linalg as sla
import scipy.linalg as la
from concurrent import futures
from functools import partial

N=12; ##N=1200

kx=np.linspace(-0.5*np.pi,0.5*np.pi,N) # (-np.pi,0)
ky=np.linspace(-0.5*np.pi,0.5*np.pi,N);
#print(kx)
Es=np.zeros((N,N,8))
Ome_ks=np.zeros((N,N,8))
#FC
#global v1 v2 t1 t2 yita1 yita2 m1 m2 gamma E1 E2 vx vy 
#global s0 sx sy sz I0

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
K1=0.1*np.pi
K2=0.15*np.pi
vx=0.2 
vy=0.0
#FC
sx=np.array([[0,1],[1,0]])
sy=np.array([[0,-1j],[1j,0]])
sz=np.array([[1,0],[0,-1]])
s0=np.array([[1,0],[0,1]]) 
I0=np.array([[0,0],[0,0]])
#FC
# parpool(3) 

# FC
def BerryCur_k(evac,eval):
    
    V=evac
    
    dHd1_kx=t1*s0+v1*yita1*sy;
    dHd1_ky=v1*sx;

    dHd2_kx=t2*s0+v2*yita2*sy;
    dHd2_ky=v2*sx;

    dP_kx=vx*sz;
    dP_ky=-1j*vy*s0;

    dH1=np.hstack((dHd1_kx,dP_kx,I0,I0))
    dH2=np.hstack((dP_kx,dHd1_kx,I0,I0))
    dH3=np.hstack((I0,I0,dHd2_kx,dP_kx))
    dH4=np.hstack((I0,I0,dP_kx,dHd2_kx))
    ####
    dH_kx=np.vstack((dH1,dH2,dH3,dH4))
    dH_kx=np.matrix(dH_kx)
    
    dH1=np.hstack((dHd1_ky,dP_ky,I0,I0))
    dH2=np.hstack((dP_ky,dHd1_ky,I0,I0))
    dH3=np.hstack((I0,I0,dHd2_ky,dP_ky))
    dH4=np.hstack((I0,I0,dP_ky,dHd2_ky))
    ###
    dH_ky=np.vstack((dH1,dH2,dH3,dH4))
    dH_ky=np.matrix(dH_ky)
#FC    
    Ome_kx=[]
    #print(len(eval))
    for i in range(len(eval)): 
        Ome_kx1=0
        for j in range(len(eval)): 
            #print(j)
            if i!=j:
                nomi_1=np.matmul(np.matmul(np.matrix.getH(V[:,i]),dH_kx),V[:,j]) # getH complex conjugate
                nomi_2=np.matmul(np.matmul(np.matrix.getH(V[:,j]),dH_ky),V[:,i])
                nomi_3=np.matmul(np.matmul(np.matrix.getH(V[:,i]),dH_ky),V[:,j])
                nomi_4=np.matmul(np.matmul(np.matrix.getH(V[:,j]),dH_kx),V[:,i])
                temp=(np.imag(nomi_1[0,0]*nomi_2[0,0])-np.imag(nomi_3[0,0]*nomi_4[0,0]))/(eval[i]-eval[j])**2;
                #print(temp)
                Ome_kx1=Ome_kx1+temp 
                
        Ome_kx.append(-1*Ome_kx1)
        
    return Ome_kx
#FC
def factor(T,E,Ef):
    beta=1/(T*25.7*10**(-3)/298); 
    ds=E.shape;
    for i in range(ds[0]):
        for j in range(ds[1]):
            f=1/(np.exp(beta*(E[i,j]-Ef))+1);
           # n[i,j]=f*(1-f)*(-beta);
            n[i,j]=f;
    return n
def dfactor(T,E,Ef):

    beta=1/(T*25.7*10**(-3)/298); 
    ds=E.shape;
    n=np.zeros(ds)
    for i in range(ds[0]):
        for j in range(ds[1]):
            f=1/(np.exp(beta*(E[i,j]-Ef))+1);
            n[i,j]=f*(1-f)*(-beta);
    return n
#FC

def calc(ky,kx):
        kyp=ky
        kx=kx
        Es0=[]
        Ome_ks0=[]
        
        kxp1=kx+K1
        kxp2=kx+K2
        Hd1=(E1+t1*kxp1)*s0+v1*(kyp*sx + yita1*kxp1*sy)+m1/2*sz
        Hd2=(E2+t2*kxp2)*s0+v2*(kyp*sx + yita2*kxp2*sy)+m2/2*sz
        P=vx*kx*sz+(-1j*vy*kyp)*s0
        Gamma=gamma*s0
            
        H1=np.hstack((Hd1,P,I0,Gamma))
        H2=np.hstack((P,Hd1,Gamma,I0))
        H3=np.hstack((I0,Gamma,Hd2,P))
        H4=np.hstack((Gamma,I0,P,Hd2))
        H=np.vstack((H1,H2,H3,H4))
            
        eval,evec=la.eigh(H)
           
        Es0.append(eval)  
        
        Ome_ks0.append(BerryCur_k(evec,eval));
        
        return Es0,Ome_ks0

def pp(kx,ky):
    Ome=[]
    Es=[]
    for kx0 in kx:
        with futures.ProcessPoolExecutor(max_workers=32) as executor:
            for result in executor.map(partial(calc, kx=kx0), ky):
                Es.append(result[0])
                Ome.append(result[1])
    return Es,Ome

def calc_2(T0,Ef0):
    Obd_temp=0
    Obd0=[]
    for index_bands in range(8): #range(min(Es.shape)):
        #print(hx)
        dEs=np.gradient(Es[:,:,index_bands],hx,axis=0);
        #print(Es[:,:,index_bands])
#        if index_bands == 0:
#            for i in range(N):
#                print(Es[:,0,index_bands])
#                print(dEs[:,0])
#        tmp1=np.zeros([N])
        #print(tmp1)
        Bcur=Ome_ks[:,:,index_bands]
        df_dE=dfactor(T0,Es[:,:,index_bands],Ef0)
#        if index_bands == 0:
#           for i in range(N):
#               print(df_dE[0,:])

#        if index_bands == 0:
#            for i in range(N):
#                print(((dEs*Bcur)*df_dE*(Es[:,:,index_bands]-Ef0)**2)[i,:])
        temp=((dEs*Bcur)*df_dE)*(Es[:,:,index_bands]-Ef0)**2;
        #print(temp)
        #print(temp.shape)
#        if index_bands == 0:
#           for i in range(N):
#               print(temp[i,:])

        temp=-1*np.trapz(np.trapz(temp,ky,axis=1),kx,axis=0)/(max(kx)-min(kx))**2/(T0*25.7*10**(-3)/298)**2; # already in T/eV
 

       #print(np.trapz(temp,ky,axis=1))
        #tmp1=(temp[N-1,]+temp[0,:])*np.pi/2
        #print(tmp1)
        #print(np.trapz(temp,ky,axis=1))
        #tmp1=0
        #for i in range(N-1):
        #        tmp1=tmp1+(temp[:,i]+temp[:,i+1])*(ky[i+1]-ky[i])/2
        #print(tmp1)
        #temp=0
        #for i in range(N-1):
        #    temp=temp+(tmp1[i]+tmp1[i+1])*(kx[i+1]-kx[i])/2 #(np.max(kx)-np.min(kx))/2
        #temp=-1*temp/(max(kx)-min(kx))**2/(T0*25.7*10**(-3)/298)**2
        Obd_temp=Obd_temp+temp;
        print(temp)
    Obd0.append(Obd_temp)
    
    return Obd0


def pp_2(Efs,Ts):
    Obd=[]
    Es=[]
    for Ef in Efs:
        with futures.ProcessPoolExecutor(max_workers=32) as executor:
            for result in executor.map(partial(calc_2, Ef0=Ef), Ts):
                #Es.append(result[0])
                Obd.append(result)
    return Obd


Es,Ome=pp(kx,ky)
Es=np.array(np.reshape(Es,[N,N,8]))
Ome_ks=np.array(np.reshape(Ome,[N,N,8]))

#for i in range(N):
#	print(Ome_ks[i,1,:])


hx=(np.max(kx)-np.min(kx))/(len(kx)-1);
###Efs=np.linspace(-0.5,0.5,101);   # 501

###Ts=np.linspace(0.5,100,151)      #551
#Ts=[0.005,20,50,100,300]
#Ts=[50]
Ts=np.linspace(300,400,1)
Efs=np.linspace(0,0.5,1)
Obd=pp_2(Efs,Ts)
print(Obd)

###Esky0=Es[:,N-1,:]
###Ome_ky0=Ome_ks[:,N-1,:]

#Esky0.tofile(file="zSpectrum_vx_0p40.csv", sep=",", format="%s")
#Ome_ky0.tofile(file="zBerryCur_vx_0p40.csv", sep=",", format="%s")

## THIS TWO SHOULD BE THE SAME WITH Dx

###np.savetxt('zvx0p2_cmap_test.csv',Obd,delimiter=',')   

