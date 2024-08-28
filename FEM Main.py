import numpy as np
import matplotlib.pyplot as plt

def shapefunc(p,co,cod,choice):
    if(choice==1):
        if (p==1):
            co=[[0.5,-0.5],
                [0.5,0.5]]
            cod=[[-0.5],
                [0.5]]
        elif (p==2):
            co=[[0,-0.5,0.5],
                [1,0,-1],
                [0,0.5,0.5]]
            cod=[[-0.5,1],
                [0,-2],
                [0.5,1]]
        elif (p==3):
            co=[[-0.0625,0.0625,0.5625,-0.5625],
                [0.5625,-1.6875,-0.5625,1.6875],
                [0.5625,1.6875,-0.5625,-1.6875],
                [-0.0625,-0.0625,0.5625,0.5625]]
            cod=[[0.0625,1.125,-1.6875],
                [-1.6875,-1.125,5.0625],
                [1.6875,-1.125,-5.0625],
                [-0.0625,1.125,1.6875]]
        elif (p==4):
            co=[[0,0.1667,-0.1667,-0.6667,0.6667],
                [0,-1.3333,2.6667,1.3333,-2.6667],
                [1,0,-5,0,4],
                [0,1.3333,2.6667,-1.3333,-2.6667],
                [0,-0.1667,-0.1667,0.6667,0.6667]]
            cod=[[0.1667,-0.3333,-2,2.6667],
                 [-1.3333,5.3333,4,-10.6667],
                 [0,-10,0,16],
                 [1.3333,5.3333,-4,-10.6667],
                 [-0.1667,-0.3333,2,2.6667]]
    else:
        if (p==1):
            co=[[0.5,-0.5],
                [0.5,0.5]]
            cod=[[-0.5],
                [0.5]]
        elif (p==2):
            co=[[0.5,-0.5,0,0,0],
                [-0.6124,0,0.6124,0,0],
                [0.5,0.5,0,0,0]]
            cod=[[-0.5,0],
                [0,1.2247],
                [0.5,0]]
        elif (p==3):
            co=[[0.5,-0.5,0,0,0],
                [0,-0.7906,0,0.7906,0],
                [-0.6124,0,0.6124,0,0],
                [0.5,0.5,0,0,0]]
            cod=[[-0.5,0,0,0],
                [-0.7906,0,2.3717,0],
                [0,1.2247,0,0],
                [0.5,0,0,0]]
        elif (p==4):
            co=[[0.5,-0.5,0,0,0],
                [0.2338,0,-1.4031,0,1.1693],
                [0,-0.7906,0,0.7906,0],
                [-0.6124,0,0.6124,0,0],
                [0.5,0.5,0,0,0]]
            cod=[[-0.5,0,0,0],
                [0,-2.8062,0,4.6771],
                [-0.7906,0,2.3717,0],
                [0,1.2247,0,0],
                [0.5,0,0,0]]
    return co,cod
def intg(re):
    coi=np.array([[0.23862,0.46791],[-0.23862,0.46791],[0.66121,0.36076],[-0.66121,0.36076],[0.93247,0.17132],[-0.93247,0.17132]])
    s=0
    for i in range(0,6):
        su=0
        for j in range(0,int(re.shape[0])):
            su=su+(re[j]*coi[i][0]**j)
        s=s+su*coi[i][1]
    return s
def elemat(ae,c,f,p):
    for i in range(0, p+1):
        for j in range(0, p+1):
            ke[i][j]=(2/h)*intg(np.polynomial.polynomial.polymul(np.polynomial.polynomial.polymul(cod[i],cod[j]),ae[0]))
            ge[i][j]=(h/2)*intg(np.polynomial.polynomial.polymul(np.polynomial.polynomial.polymul(co[i],co[j]),c[0]))
        fe[i][0]=(h/2)*intg(np.polynomial.polynomial.polymul(co[i],f[0]))
    return ke,ge,fe
def getval(co,x,p,d):
    if(d==0):
        s=0
        for i in range(0,p+1):
            s=s+co[i]*x**i
    else:
        s=0
        for i in range(0,p):
            s=s+co[i]*x**i
    return s
n=int(input("Enter number of elements "))
p=int(input("Enter order no. "))
choice=int(input("Enter type of Shape function (Lagrange-1/Heirarchic-2)"))
mptbc=0
if (n%2==0):
    mptbc+=int(input("Enter force at midpoint"))
co=np.zeros((5,5))
cod=np.zeros((5,5))
co,cod=shapefunc(p,co,cod,choice)
aeco=np.zeros((1,3))
cco=np.zeros((1,3))
fco=np.zeros((1,3))
bc1=int(input("Type of BC at x=0 (1 or 2 or 3)= "))
bc2=int(input("Type of BC at x=1 (1 or 2 or 3)= "))
if(bc1==2):
    force1=float(input("Enter force= "))
elif (bc1==1):
    dis1=float(input("Enter displacement= "))
else:
    spr1=float(input("Enter spring constant= "))
    dev1=float(input("Enter spring deviation= "))
if(bc2==2):
    force2=float(input("Enter force= "))
elif(bc2==1):
    dis2=float(input("Enter displacement= "))
else:
    spr2=float(input("Enter spring constant= "))
    dev2=float(input("Enter spring deviation= "))

aep=int(input("Enter order of ae "))
for i in range (0,aep+1):
    aeco[0][i]=int(input("Eneter coeff. "))
cp=int(input("Enter order of c "))
for i in range (0,cp+1):
    cco[0][i]=int(input("Eneter coeff. "))
fp=int(input("Enter order of t "))
for i in range (0,fp+1):
    fco[0][i]=int(input("Eneter coeff. "))
h=1/n
nodloc=np.zeros((n,2))
for i in range(0,n):
    for j in range(0,2):
        if (i==j==0):
            continue
        elif (j%2==0):
            nodloc[i][j]=nodloc[i-1][j+1]
        else:
            nodloc[i][j]=nodloc[i][j-1]+h
K=np.zeros((n*p+1,n*p+1))
G=np.zeros((n*p+1,n*p+1))
F=np.zeros((n*p+1,1))
Q=np.zeros((n*p+1,1))
re=np.zeros((1,5))
for i in range(0,n):
    ke=np.zeros((p+1,p+1))
    ge=np.zeros((p+1,p+1))
    fe=np.zeros((p+1,1))
    sunod=(nodloc[i][0]+nodloc[i][1])/2
    aecof=np.array([[aeco[0][0]+aeco[0][1]*sunod+aeco[0][2]*sunod**2,aeco[0][1]*(h/2)+aeco[0][2]*h*sunod,aeco[0][2]*(h**2/4)]])
    ccof=np.array([[cco[0][0]+cco[0][1]*sunod+cco[0][2]*sunod**2,cco[0][1]*(h/2)+cco[0][2]*h*sunod,cco[0][2]*(h**2/4)]])
    fcof=np.array([[fco[0][0]+fco[0][1]*sunod+fco[0][2]*sunod**2,fco[0][1]*(h/2)+fco[0][2]*h*sunod,fco[0][2]*(h**2/4)]])
    ke,ge,fe=elemat(aecof,ccof,fcof,p)
    if (i==0):
        K[:p+1,:p+1]+=ke
        G[:p+1,:p+1]+=ge
        F[:p+1,0:1]+=fe
    else:
        K[i*p:i*p+p+1,i*p:i*p+p+1]+=ke
        G[i*p:i*p+p+1,i*p:i*p+p+1]+=ge
        F[i*p:i*p+p+1,0:1]+=fe
KF=K+G
if (bc1==1):
    for i in range (1,n*p+1):
        F[i]=F[i]-dis1*KF[i][0]
    for i in range(0,n*p+1):
        for j in range (0,n*p+1):
            if(i==0 or j==0):
                KF[i][j]=0
    KF[0][0]=1
    F[0][0]=dis1
    Q[0][0]=0
if (bc2==1):
    for i in range (0,n*p):
        F[i]=F[i]-dis2*KF[i][n*p]
    for i in range(0,n*p+1):
        for j in range (0,n*p+1):
            if(i==n*p or j==n*p):
                KF[i][j]=0
    KF[n*p][n*p]=1
    F[n*p][0]=dis2
    Q[n*p][0]=0
if (bc1==2):
    Q[0][0]=force1
if (bc2==2):
    Q[n*p][0]=force2
if (bc1==3):
    KF[0][0]+=spr1
    Q[0][0]=spr1*dev1
if (bc2==3):
    KF[n*p][n*p]+=spr2
    Q[n*p][0]=spr2*dev2
if(mptbc!=0):
    v=int((n*p)/2)
    Q[v][0]=mptbc
U=np.linalg.inv(KF)@(F+Q)
print("U= ",U)
x=np.linspace(0,1,1000)
yh=[]
yhd=[]
for i in range(0,n):
    for j in range(0,1000):
        if(x[j]>=nodloc[i][0] and x[j]<=nodloc[i][1]):
            aux=(2*x[j]-(nodloc[i][0]+nodloc[i][1]))/h
            fr=0
            frd=0
            b=0
            for m in range(i*p,i*p+p+1):
                fr=fr+(U[m][0]*getval(co[b],aux,p,0))
                frd=frd+(U[m][0]*getval(cod[b],aux,p,1))*(2/h)
                b=b+1
            yh.append(round(fr,6))
            yhd.append(round(frd,6))
x=np.linspace(0,1,np.size(yh))
line1,=plt.plot(x,yh,label="FEM")
plt.legend(handles=[line1])
plt.title("Plot of elongation with x")
plt.ylabel("u")
plt.xlabel("x")
plt.show()
line1,=plt.plot(x,yhd,label="FEM")
plt.legend(handles=[line1])
plt.title("Plot of du/dx with x")
plt.ylabel("du/dx")
plt.xlabel("x")
plt.show()
