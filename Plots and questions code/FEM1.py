import numpy as np
import matplotlib.pyplot as plt
def FEMsolver(n,p,choice,bc1,bc2,aeco,cco,fco):
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
    co=np.zeros((5,5))
    cod=np.zeros((5,5))
    co,cod=shapefunc(p,co,cod,choice)
    if(bc1==2):
        force1=float(input("Enter force= "))
    elif (bc1==1):
        dis1=0
    else:
        spr1=float(input("Enter spring constant= "))
        dev1=float(input("Enter spring deviation= "))

    if(bc2==2):
        force2=0
    elif(bc2==1):
        dis2=float(input("Enter displacement= "))
    else:
        spr2=float(input("Enter spring constant= "))
        dev2=float(input("Enter spring deviation= "))
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
    # print(nodloc)
    K=np.zeros((n*p+1,n*p+1))
    G=np.zeros((n*p+1,n*p+1))
    F=np.zeros((n*p+1,1))
    Q=np.zeros((n*p+1,1))
    for i in range(0,n):
        ke=np.zeros((p+1,p+1))
        ge=np.zeros((p+1,p+1))
        fe=np.zeros((p+1,1))
        sunod=(nodloc[i][0]+nodloc[i][1])/2
        aecof=np.array([[aeco[0][0]+aeco[0][1]*sunod+aeco[0][2]*sunod**2,aeco[0][1]*(h/2)+aeco[0][2]*h*sunod,aeco[0][2]*(h/4)]])
        ccof=np.array([[cco[0][0]+cco[0][1]*sunod+cco[0][2]*sunod**2,cco[0][1]*(h/2)+cco[0][2]*h*sunod,cco[0][2]*(h/4)]])
        fcof=np.array([[fco[0][0]+fco[0][1]*sunod+fco[0][2]*sunod**2,fco[0][1]*(h/2)+fco[0][2]*h*sunod,fco[0][2]*(h/4)]])
        ke,ge,fe=elemat(aecof,ccof,fcof,p)
        # print(ke+ge,"\n",fe)
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
    # print(KF,"\n",F)
    U=np.linalg.inv(KF)@(F+Q)
    # print(U)
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
                yh.append(fr)
                yhd.append(frd)
    ye=[]
    yed=[]
    i=0
    while (i<1):
        ye.append(-i*i/2+(force2+1)*i+dis1)
        # ye.append(((np.e**(-i))*(np.e**(i) - 1)*(-np.e**(i) + np.e**(2)))/(np.e**(2) + 1))
        yed.append(-i+(force2+1))
        i=i+0.001
    # print(np.size(yh),np.size(ye))
    if (np.size(yh)<np.size(ye)):
        for i in range(0,(np.size(ye)-np.size(yh))):
            ye.pop()
            x=x[:-1]
    elif (np.size(ye)<np.size(yh)):
        for i in range(0,(np.size(yh)-np.size(ye))):
            yh.pop()
    return yh,yhd,ye,yed,x
    # print(np.size(yh))

nelem=[1,2,5,10,100]
order=[1,2]
choice=1
aeco=np.zeros((1,3))
cco=np.zeros((1,3))
fco=np.zeros((1,3))
aep=0
for i in range (0,aep+1):
    aeco[0][i]=1
cp=0
for i in range (0,cp+1):
    cco[0][i]=0
fp=0
for i in range (0,fp+1):
    fco[0][i]=1
bc1=1
bc2=2
for p in order:
    yh,yhd,ye,yed,x=FEMsolver(1,p,choice,bc1,bc2,aeco,cco,fco)
    line1, =plt.plot(x,yh, label="1 element")
    yh,yhd,ye,yed,x=FEMsolver(2,p,choice,bc1,bc2,aeco,cco,fco)
    line2, =plt.plot(x,yh, label="2 element")
    yh,yhd,ye,yed,x=FEMsolver(5,p,choice,bc1,bc2,aeco,cco,fco)
    line3, =plt.plot(x,yh, label="5 element")
    yh,yhd,ye,yed,x=FEMsolver(10,p,choice,bc1,bc2,aeco,cco,fco)
    line4, =plt.plot(x,yh, label="10 element")
    yh,yhd,ye,yed,x=FEMsolver(100,p,choice,bc1,bc2,aeco,cco,fco)
    line5, =plt.plot(x,yh, label="100 element")
    line6, =plt.plot(x,ye, label="exact")
    plt.legend(handles=[line1, line2, line3, line4, line5, line6])
    plt.title("Plot of FEM with exact solution")
    plt.xlabel("x")
    plt.ylabel("u")
    plt.show()
for p in order:
    yh,yhd,ye,yed,x=FEMsolver(1,p,choice,bc1,bc2,aeco,cco,fco)
    line1, =plt.plot(x,yhd, label="1 element")
    yh,yhd,ye,yed,x=FEMsolver(2,p,choice,bc1,bc2,aeco,cco,fco)
    line2, =plt.plot(x,yhd, label="2 element")
    yh,yhd,ye,yed,x=FEMsolver(5,p,choice,bc1,bc2,aeco,cco,fco)
    line3, =plt.plot(x,yhd, label="5 element")
    yh,yhd,ye,yed,x=FEMsolver(10,p,choice,bc1,bc2,aeco,cco,fco)
    line4, =plt.plot(x,yhd, label="10 element")
    yh,yhd,ye,yed,x=FEMsolver(100,p,choice,bc1,bc2,aeco,cco,fco)
    line5, =plt.plot(x,yhd, label="100 element")
    line6, =plt.plot(x,yed, label="exact")
    plt.legend(handles=[line1, line2, line3, line4, line5, line6])
    plt.title("Plot of derivatives")
    plt.xlabel("x")
    plt.ylabel("F")
    plt.show()
error=[]
lnn=[]
energy=[]
for n in nelem: 
    yh,yhd,ye,yed,x=FEMsolver(n,1,choice,bc1,bc2,aeco,cco,fco)
    er=np.log(np.linalg.norm(np.array(ye)-np.array(yh)))
    error.append(er)
    lnn.append(np.log(n))
    s=0
    for i in range(0,998):
        fx1=yhd[i]**2
        fx2=yhd[i+1]**2
        s=s+0.5*(fx1+fx2)*(x[i+1]-x[i])
    energy.append(s)
print(lnn,error)
line1,=plt.plot(nelem,energy, label="linear",marker=".")
error=[]
lnn=[]
energy=[]
eenergy=[]
for n in nelem: 
    yh,yhd,ye,yed,x=FEMsolver(n,2,choice,bc1,bc2,aeco,cco,fco)
    er=np.log(np.linalg.norm(np.array(ye)-np.array(yh)))
    error.append(er)
    lnn.append(np.log(n))
    s=0
    for i in range(0,998):
        fx1=yhd[i]**2
        fx2=yhd[i+1]**2
        s=s+0.5*(fx1+fx2)*(x[i+1]-x[i])
    energy.append(s)
    s=0
    for i in range(0,998):
        fx1=yed[i]**2
        fx2=yed[i+1]**2
        s=s+0.5*(fx1+fx2)*(x[i+1]-x[i])
    eenergy.append(s)
print(lnn,error)
line2,=plt.plot(nelem,energy, label="quadratic",marker=".")
line3,=plt.plot(nelem,eenergy, label="exact",marker=".")
plt.legend(handles=[line1, line2, line3])
plt.title("Plot of energy with increasing number of elements")
plt.ylabel("energy")
plt.xlabel("No. of elements")
plt.show()
# er=0
# for i in range(0,999):
#     er=er+((ye[i]-yh[i])**2)
# e=np.sqrt(er/1000)*100
# print ("RMSE error%= ",e,"%")
# plt.plot(yh)
# plt.plot(ye)
# plt.show()
# plt.plot(yhd)
# plt.plot(yed)
# plt.show()
