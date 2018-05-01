import numpy as np
import sympy as sp
"""
Dinicial no funciona
los sp.var() aún no funcionan -- cambiarlos por symbol
Faltan dudas
Falta iterar r2 
"""
#angulo para distancia Tierra-Sol
def d(dia):
    return(365/360)*(31+25+ dia)

def decimales(a,b,c):
    return(a+b/60+c/3600)
def squ(x):
    return(x**2)

#constantes
c= 173 #UA/dia~normal (velocidad de la luz) 
G=6.67408*10**-11 #N*m**2/kg**2
mu=G*(1.989*10**30 + 5.972*10**24)#kg
k=0.01720209895
E_=[[23,26,07.249],[23,26,07.293],[23,26,07.335]]#grados, minutos,segundos
E_=[np.sum(decimales(x[0],x[1],x[2])for x in E_)/len(E_)]#angulo ecliptica-eceleste
R=1#UA
Exct= 0.016
P=0#~luego se cambia ~creo que nunca se usa

#datos observaciones
Ar= np.array([[8,42,13.71520391],[8,42,3.6776478],[8,41,57.61776475],[8,41,45.35651401]])# reales
Dec=np.array([[31,57,48.37108478],[31,55,36.28691431],[31,53,43.2407565],[31,40,4.33621619]])# reales
datos=np.array([[15*decimales(x[0],x[1],x[2]) for x in Ar],[decimales(x[0],x[1],x[2]) for x in Dec]])#AR,DC
datosTransformados=[datos[:,0],datos[:,1],datos[:,2]]
datosExtra=[datos[:,3]]

#vector observacion i,j,k
def vector(ar,d):
    return np.cos(ar)*np.cos(d),np.sin(ar)*np.cos(d),np.sin(d)

vectores=[list(vector(i[0],i[1]))for i in datosTransformados]
p1=[vectores[0]]#observaciones
p2=[vectores[1]]#observaciones
p3=[vectores[2]]#observaciones
#print(vectores)
#norma vectores p
normaVectores=[np.linalg.norm(i)for i in vectores]#not sure ~no es necesario aún ¿?

#R (vector posición  entre sol y tierra)
def rTS (angulo):
    return(R*(1-squ(Exct))/1+(Exct*np.cos(angulo)))
R1=rTS(d(14))
R2=rTS(d(15))
R3=rTS(d(16))

#tiempo en kt (gaussian)
def dias(h,m,s):
    return((h+(s/60+m)/60)/24)
T1=dias(9,27,39)
T2=dias(11,21,33)
T3=dias(7,44,46)   
t3=k*(T3-T2)
t1=k*(T1-T2)
t=k*(T3-T1)
tiempos=[t,t1,t3]

#valor inicial de r2 (vector posicion de la observación central)
"""
Tal vez los siguientes se puedan optimizar en uno solo.
"""
Dinicial=np.dot(p1,np.cross(p2,p3))
print(Dinicial) ##-----error
def D1x (Rx):
    return np.dot(np.cross(Rx,p2),p3)#change
def D2x (Rx):
    return np.dot(np.cross(p1,Rx),p3)
def D3x (Rx):
    return np.dot(p1,np.cross(p2,Rx))

A1=t3/t
B1=(A1/6)*(squ(t)-squ(t3))
A3=-t1/t
B3=(A3/6)*(squ(t)-squ(t1))

A=(A1*D2x(R1)-D2x(R2)+A3*D2x(R3))/-Dinicial
B=(B1*D2x(R1)+B3*D2x(R3))/-Dinicial
Ee=-2*(np.dot(p2,R2))
F=squ(R2)
aa=-(squ(A)+A*Ee+F)
b=-mu(2*A*B+B*E)
cc=-squ(mu)*squ(B)

r_2=sp.Symbol('r')
equacion=r2**8+aa*(r2**6)+b*(r2**3)+cc
r2=sp.solve(equacion,r2)
print(r2)

# Iterate
"""
Desde aquí
"""
#truncated f & g
def f (r2,dr2,kt):
    return 1-mu*kt/2*(r2**3)+ mu*(np.dot(r2,dr2))/2*(r2**5)
f=[f(r2,dr2,kt)for kt in tiempos]
f1=f[0]
f3=f[2]

def g(r2,kt):
    return kt-(kt**3)*mu/6*(r2**3)
g=[g(r2,kt) for kt in tiempos]
g1=g[0]
g3=g[2]

#dr2
r1=f1*r2+g1*dr2
r3=f3*r2+g3*dr2
d1=-f3/(f1*g3-f3*g1)
d3=f1/(f1*g3-f3*g1)
dr2=d1*r1+d3*r3

#hallar c1 y c3
u=mu/r2**3
c1= g3/f1*g3-g1*f3#~~ A1+u*B1 
c2=-1
c3= -g1/f1*g3-g1*f3 #~~A3+u*B3

#hallar vectores posicion entre tierra y asteroide
P1=(c1*D1x(R1)+c2*D1(R2)+c3*D1(R3))/c1*Dinicial
P2=(c1*D2x(R1)+c2*D2(R2)+c3*D2(R3))/c2*Dinicial#~A+u*B ~ya entendí why not
P3=(c1*D3x(R1)+c2*D3(R2)+c3*D3(R3))/c3*Dinicial
print(P2)

#hallar posicion sol asteroide (2)
r1=P1-R1
r2=P2-R2
r3=P3-R3
r=[r1,r2,r3]
#hallar r y r.
#Corrección tiempo de la luz ~~en teoría desde aquí ya empieza
T1= T1-P1/c
T2= T2-P2/c
T3= T3-P3/c

t3=k*(T3-T2) #esto ya se hizo
t1=k*(T1-T2)
t=k*(T3-T1)
tiempos=[t,t1,t3]
"""
Jusqu'ici
"""
#rotar vectores r
def Rotar (E_,vector):
    return np.dot(np.matrix([[1,0,0],[0,np.cos(E_),np.sin(E_)],[0,-np.sin(E_),np.cos(E_)]]),vector)
r=[Rotar(E_,vector) for vector in r]
  
"""
ELEMENTOS ORBITALES
"""
def com(x1,x2):
    comun=[]
    for x in x1:
        if x in x2:
           comun.append(x)
           return (comun)

rm=r[1]-r[0]
rp= rm/(T2-T1)#¿?velocidad----hallar ~not sure
Norma_rm= np.linalg.norm(rm)
h=np.cross(rm,rp)

#Semieje mayor

a=sp.Symbol('a')
#eq_a=((2/Norma_rm)-(mu*(2/Norma_rm- 1/a))/mu)**-1
a=(2/Norma_rm)-(np.dot(rp,rp)/mu)**-1
#a=sp.solve(eq_a,a)

P=2*np.pi*a**(3/2) #~~Verificar donde se usa
#Excentricidad

e=np.sqrt((1-(np.linalg.norm(np.cross(rm,rp))**2))/mu*a)

#Inclinación

hz= squ(h)-squ(hx)+squ(hy)#hx hy ¿?--h[0] y h[1]-- esto no necesario
i=np.arccos(hz/np.linalg.norm(h))

#Longitud del nodo ascendente

hx=h*np.sin(i)*np.sin(O)
hy=-h*np.sin(i)*np.cos(O)
O1=np.arcsin(hx/h*np.sin(i))
O2=np.arccos(-hy/h*np.sin(i))
O= com(O1,O2)

#Perihelio
#hallar v True anomaly

v1=np.arccos(((a*(1-squ(e))/Norma_rm)-1)/e)
v2=np.arcsin(np.dot(rm,rp)*a*(1-squ(e))/np.linalg.norm(h)*Norma_rm)
v= np.rad2deg(com(v1,v2))

#hallar U

U1=np.arccos(np.dot(rm,np.cos(O)+np.sin(O))/Norma_rm) # = xcos(O) + ycos(O)/r
z=Norma_rm*np.cross(n,(rm/Norma_rm))
U2=np.arcsin(z[0]/Norma_rm*np.sin(i)*np.sin(O))
U=np.rad2deg(com(U1,U2))

W=(U-v)
w=[]
for x in W:
    if (0<=x<360):#x>=0 & x<360 ¿?
        w.append(x)
#Mean anomaly

E=np.arccos((1/e)*(1-Norma_rm/a)) 
if (0<=np.rad2deg(E)<=180):#~puede que no funcione, entonces for¿?
    M= E-e*np.sin(E) 
#Hay otra manera de M ~revisar
