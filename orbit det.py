import numpy as np
import sympy as sp
"""
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
G=1.48814*10**-34#AU**3/kg*día**2~~ 6.67408*10**-11 #N*m**2/kg**2
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
p1=vectores[0]#observaciones
p2=vectores[1]#observaciones
p3=vectores[2]#observaciones
#print(vectores)---[[],[]]
#norma vectores p
normaVectores=[np.linalg.norm(i)for i in vectores] #~no es necesario aún 

#R (vector posición  entre sol y tierra)
def rTS (angulo):
    return(R*(1-squ(Exct))/1+(Exct*np.cos(angulo)))
R1=rTS(d(14))#--magnitud 
R1=[R1*np.cos(d(14)),R1*np.sin(d(14)),0]#¿? No sé bien si aquí sí tenga el vector. En caso de que "sí"estaría en dos dimensiones
R2=rTS(d(15))
R2=[R2*np.cos(d(15)),R2*np.sin(d(15)),0]
R3=rTS(d(16))
R3=[R3*np.cos(d(16)),R3*np.sin(d(16)),0]
#print(R3,np.linalg.norm(R3))----[],#
#tiempo en kt (gaussian)

def dias(h,m,s):
    return((h+(s/60+m)/60)/24)
T1=dias(9,27,39)
T2=dias(11,21,33)
T3=dias(7,44,46)

def tiempos(T1,T2,T3):
    time=[]#t-t1-t3 
    time.append(k*(T3-T1))
    time.append(k*(T1-T2))
    time.append(k*(T3-T2))
    return(time)
t=tiempos(T1,T2,T3)
#print(t)---[]
#valor inicial de r2 (vector posicion de la observación central)
#print(p1)
#print(p2)--[]
#print(p3)
#print(R1)
Dinicial=np.dot(p1,np.cross(p2,p3))
#print(Dinicial) ---#
def DXx (Rx):
    if(Rx==R1):
        return np.dot(np.cross(Rx,p2),p3)
    elif (Rx==R2):
        return np.dot(np.cross(p1,Rx),p3)
    elif (Rx==R3):
        return np.dot(p1,np.cross(p2,Rx))
#print(DXx(R1))---#
A1=t[2]/t[0]#---#
B1=(A1/6)*(squ(t[0])-squ(t[2]))#---#
A3=-t[1]/t[0]#---#
B3=(A3/6)*(squ(t[0])-squ(t[1]))#---#

A=(A1*DXx(R1)-DXx(R2)+A3*DXx(R3))/-Dinicial#---#
B=(B1*DXx(R1)+B3*DXx(R3))/-Dinicial#---#

Ee=-2*(np.dot(p2,R2))#---#
F=squ(np.linalg.norm(R2))#---#
#print(Ee,F)
aa=-(squ(A)+A*Ee+F)#---#
b=-mu*(2*(A*B+B*Ee))#---#
cc=-squ(mu)*squ(B)#---#
#print(aa,b,cc)

r_2=sp.Symbol('r')
equacion=r_2**8+aa*(r_2**6)+b*(r_2**3)+cc
r2=sp.solve(equacion,r_2)
r2 = np.array(r2)
#r2=[x>0 for x in r2] 
print(r2[0])#---#¿? algunos son números y otros al parecer complejos;
"""
1. Está bien?
2. Cómo hago en ese caso para quitarlos
"""
r=[]#¿?en el for sí estoy reemplazandola?
# Iterate
"""
Desde aquí
"""
finale= False
while finale==False: #¿? cómo haría esto con for
    otro_r2=r2#¿?no estoy segura si estoy va dentro o fuera del for por lo que puede que al hacer la verificación al ser ambos 'r2', de 0
#truncated f & g
    dr2=sp.Symbol('d')
    def f (r2,dr2,kt):
        return 1-mu*kt/2*(r2**3)+ mu*(np.dot(r2,dr2))/2*(r2**5)
    f=[f(r2,dr2,kt)for kt in t]
    f1=f[0]
    f3=f[2]
    
    def g(r2,kt):
        return kt-(kt**3)*mu/6*(r2**3)
    g=[g(r2,kt) for kt in t]
    g1=g[0]
    g3=g[2]
    
    #dr2
    r1=f1*r2+g1*dr2
    r3=f3*r2+g3*dr2
    d1=-f3/(f1*g3-f3*g1)
    d3=f1/(f1*g3-f3*g1)
    dr2=d1*r1+d3*r3
    #dr2=sp.solve()
    #hallar c1 y c3
    
    c1= g3/f1*g3-g1*f3
    c2=-1
    c3= -g1/f1*g3-g1*f3 
    
    #hallar vectores posicion entre tierra y asteroide
    P1=(c1*DXx(R1)+c2*DXx(R2)+c3*DXx(R3))/c1*Dinicial
    P2=(c1*DXx(R1)+c2*DXx(R2)+c3*DXx(R3))/c2*Dinicial
    P3=(c1*DXx(R1)+c2*DXx(R2)+c3*DXx(R3))/c3*Dinicial
    print(P2)
    
    #hallar posicion sol asteroide 
    r1=P1-R1
    r2=P2-R2
    r3=P3-R3
    r=[r1,r2,r3]
    #hallar r y r.
    #Corrección tiempo de la luz ~~en teoría desde aquí ya empieza
    T1= T1-P1/c
    T2= T2-P2/c
    T3= T3-P3/c
    
    t=tiempos(T1,T2,T3)
    if(np.abs(otro_r2-r2)<=0.001):#¿? restando vectores ~ debería sacar norma o abs ya lo hace?
        finale=True
"""
Jusqu'ici
"""

#rotar vectores r
def Rotar (E_,vector):
    return np.dot(np.matrix([[1,0,0],[0,np.cos(E_),np.sin(E_)],[0,-np.sin(E_),np.cos(E_)]]),vector)
r=[Rotar(E_,vector) for vector in r]
r2=r[1] 
"""
ELEMENTOS ORBITALES
"""
def com(x1,x2):
    comun=[]
    for x in x1 and x2:
           comun.append(x)
           return (comun)
R2_punto=(R3-R2)/k*(T3-T2)#--¿?Método de Lagrange
p2_punto=(squ(t[2])*(p1-p2)-squ(t[0])*(p3-p2))/(t[0]*t[1]*t[2])#---manera de hacer esto más corta ¿? 
p2_2punto=-2((t[2]*(p1-p2)-t[0]*(p3-p2))/t[0]*t[1]*t[2])
P2_punto=-1/2*((1/np.linalg.norm(r2)**3)-((1+1/328900.5)/np.linalg.norm(R2)**3))*((np.dot(np.cross(p2,p2_2punto),R2)/(np.dot(np.cross(p2,p2_punto)),p2_2punto)))#--¿?supongo que es la norma la que se eleva 
r2_punto=p2_punto*p2+p2*P2_punto-R2_punto
rm=r2
rp=r2_punto
Norma_rm= np.linalg.norm(rm)
h=np.cross(rm,rp)

#Semieje mayor

a=sp.Symbol('a')
a=(2/Norma_rm)-(np.dot(rp,rp)/mu)**-1

P=2*np.pi*a**(3/2) #~~Verificar donde se usa
#Excentricidad

e=np.sqrt((1-(np.linalg.norm(np.cross(rm,rp))**2))/mu*a)

#Inclinación

hz= h[3]
i=np.arccos(hz/np.linalg.norm(h))

#Longitud del nodo ascendente

hx=h[0]
hy=-h[1]
O1=np.arcsin(hx/h*np.sin(i))
O2=np.arccos(-hy/h*np.sin(i))
O= com(O1,O2)

#Perihelio
#hallar v True anomaly

v1=np.arccos(((a*(1-squ(e))/Norma_rm)-1)/e)
v2=np.arcsin(np.dot(rm,rp)*a*(1-squ(e))/np.linalg.norm(h)*Norma_rm)
v= np.rad2deg(com(v1,v2))

#hallar U

U1=np.arccos(np.dot(rm,np.cos(O)+np.sin(O))/Norma_rm) 
n=[np.cos(O),np.sin(O),0]
z=Norma_rm*np.cross(n,(rm/Norma_rm)) 
U2=np.arcsin(z[0]/Norma_rm*np.sin(i)*np.sin(O))
U=np.rad2deg(com(U1,U2))

W=(U-v)
w=[]
for x in W:
    if (0<=x<360):
        w.append(x)
#Mean anomaly

E=np.arccos((1/e)*(1-Norma_rm/a)) 
M= E-e*np.sin(E) 
#Hay otra manera de M ~revisar
    
"""
    ¿?
"""
#f y g otra vez
n=np.sqrt(mu/a**3)
def fi(r2,DEi):#DEi= delta de la anomalía excéntrica (Ei-E.observación central) 
    return(1-((a/r2)*(1-np.cos(DEi))))
def gi(ti,DEi):
    return((ti-T2)+(1/n)*np.sin(DEi)-DEi)
