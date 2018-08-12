import numpy as np
"""
Falta iterar r2 
p=vector unitario tierra ceres
P=vector tierra ceres
r=vector sol ceres
R=vector tierra sol
T=tiempo observación
t=kt tiempo gauss

"""
#angulo para distancia Tierra-Sol
def d(dia):
    return(365/360)*(31+25+ dia)

def decimales(a,b,c):
    return(a+b/60+c/3600)
def squ(x):
    return(x**2)

#constantes
c= 173.1446 #UA/dia~normal (velocidad de la luz) 
G=1.48814*(10**-34)#AU**3/kg*día**2~~ 6.67408*10**-11 #N*m**2/kg**2
mu=G*((1.989*(10**30)) + (5.972*(10**24)))#kg
k= 0.01720209895 
#Epsilon=[[23,26,07.249],[23,26,07.293],[23,26,07.335]]#grados, minutos,segundos
#Epsilon=[np.deg2rad(np.sum(decimales(x[0],x[1],x[2])for x in E_)/len(E_))]#angulo ecliptica-eceleste =23.435358981481482
Epsilon=np.deg2rad(23.4352570810)
R=1#UA
Exct= 0.016
P=0#~luego se cambia ~creo que nunca se usa

#datos observaciones
Ar= np.array([[8, 42, 13.71520391],[8, 42, 3.6776478],[8, 41, 57.61776475],[8, 41, 45.35651401]])# reales
Dec=np.array([[31, 57, 48.37108478],[31, 55, 36.28691431],[31, 53, 43.2407565],[31, 40, 4.33621619]])# reales
datos=np.array([[np.deg2rad(15*decimales(x[0],x[1],x[2])) for x in Ar],[np.deg2rad(decimales(x[0],x[1],x[2])) for x in Dec]])#AR,DC
datosTransformados=[datos[:,0],datos[:,1],datos[:,2]]
datosExtra=[datos[:,3]]

#vector observacion i,j,k ~~vectores unitarios
def vector(ar,d):
    return np.cos(ar)*np.cos(d),np.sin(ar)*np.cos(d),np.sin(d)   

vectores=[list(vector(i[0],i[1]))for i in datosTransformados]
p1=np.array(vectores[0])#observaciones 
p2=np.array(vectores[1])#observaciones
p3=np.array(vectores[2])#observaciones
#print(p1)##---[]
#norma vectores p
normaVectores=[np.linalg.norm(i)for i in vectores] #~no es necesario aún 

#R (vector posición  entre sol y tierra)

R1=np.array([-9.878320430172510E-01,   1.032889954741025E-01,   4.478063517421658E-02])
R2=np.array([-9.900665016009336E-01,   8.627815975548332E-02,   3.740712265198442E-02])
R3=np.array([-9.913438623140280E-01,   7.509484630391650E-02,   3.255958030237219E-02])
#print(R2)

#tiempo en kt (gaussian)
T=[2458191.8942,2458192.9733,2458193.68163]

T1=np.array(T[0])
T2=np.array(T[1])
T3=np.array(T[2])
#print(T3)

def tiempos(T1,T2,T3):
    time=[]#t-t1-t3 
    time.append(k*(T3-T1))
    time.append(k*(T1-T2))
    time.append(k*(T3-T2))
    return(time)
t=tiempos(T1,T2,T3)
#print(t)#---[]
#valor inicial de r2 (vector posicion de la observación central)
#print(p1)
#print(p2)--np.array
#print(p3)
#print(R1)--np.array
Dinicial=np.dot(p1,np.cross(p2,p3))
#print(Dinicial)# ---#

def D1x (j):
    return np.dot(np.cross(j,p2),p3)
def D2x (j):
    return np.dot(np.cross(p1,j),p3)
def D3x (j):
    return np.dot(p1,np.cross(p2,j))
#print(D1x(R1))#---#

a1=t[2]/t[0]#---#
a2=-1
a3=-t[1]/t[0]#---#

P1=((a1*D1x(R1))+(a2*D1x(R2))+(a3*D1x(R3)))/(a1*Dinicial)
P2=((a1*D2x(R1))+(a2*D2x(R2))+(a3*D2x(R3)))/(a2*Dinicial)
P3=((a1*D3x(R1))+(a2*D3x(R2))+(a3*D3x(R3)))/(a3*Dinicial)
#print(P1)#---# -1.4860910079700003

r1o=np.array((p1*P1)-R1)
r2o=np.array((p2*P2)-R2)
r3o=np.array((p3*P3)-R3)
#print(r2o)#---vector [ 1.73409729 -0.95695466 -0.75102124]
r2o_punto=(r3o-r1o)/t[0]
#print(r2o_punto)#  [-363.0428774815191, -1320.8112113396828, 0.0027111788399157956] - con np.array--esto[-3.63042877e+02 -1.32081121e+03  2.71117884e-03]
rs=[]

def f (r2,r2_punto,kt):
    return 1-(mu*kt/2*(np.linalg.norm(r2)**3))+((kt**3)*mu*(np.dot(r2,r2_punto))/2*(np.linalg.norm(r2)**5))+(((kt**4)/24*(np.linalg.norm(r2)**3))*((3*((np.dot(r2_punto,r2_punto)/(np.linalg.norm(r2)**2))-(1/(np.linalg.norm(r2)**3))))-(15*((np.dot(r2,r2_punto)/(np.linalg.norm(r2)**2))**2))+(1/(np.linalg.norm(r2)**3))))
def g(r2,r2_punto,kt):
    return kt-(kt**3)*mu/6*(np.linalg.norm(r2)**3)-((kt**4)*(np.dot(r2,r2_punto)/4*(np.linalg.norm(r2)**5)))
# Iterate


#print(r2o)
#print(r2o_punto)
finale=False
r2_ciclo=r2o
r2_punto_ciclo= r2o_punto
while finale==False: 
    
    #Corrección tiempo de la luz ~~en teoría desde aquí ya empieza
    T1= T1-(P1/c)
    T2= T2-(P2/c)
    T3= T3-(P3/c)
    
    t=tiempos(T1,T2,T3)
    #truncated f & g        
    """
    f=[f(r2_ciclo,r2_punto_ciclo,kt)for kt in t]
    f1=f[0]
    f3=f[2]
    #print(f)
    """
    f1=f(r2_ciclo,r2_punto_ciclo,T[1])
    f3=f(r2_ciclo,r2_punto_ciclo,T[2])
    print(f1)
    print(f3)
    """
    g=[g(r2_ciclo,r2_punto_ciclo,kt) for kt in t]
    g1=g[0]
    g3=g[2]
    #print(g)
    """
    g1=g(r2_ciclo,r2_punto_ciclo,T[0])
    g3=g(r2_ciclo,r2_punto_ciclo,T[2])
   
    #hallar c1 y c3
    
    c1= g3/((f1*g3)-(g1*f3))
    c2=-1
    c3= -g1/((f1*g3)-(g1*f3)) 
    
    #hallar vectores posicion entre tierra y asteroide
    P1=(c1*D1x(R1)+c2*D1x(R2)+c3*D1x(R3))/c1*Dinicial
    P2=(c1*D2x(R1)+c2*D2x(R2)+c3*D2x(R3))/c2*Dinicial
    P3=(c1*D3x(R1)+c2*D3x(R2)+c3*D3x(R3))/c3*Dinicial
   
    #hallar posicion sol asteroide 
    r1=np.array((p1*P1)-R1)
    r2=np.array((p2*P2)-R2)
    r3=np.array((p3*P3)-R3)
    #print(r2_punto_ciclo)
    r2_punto=(r3-r1)/t[0]

    r=[r1,r2,r3]
    #print(r2)
    #hallar r y r.
    #print(r2_ciclo,r2)
    if(np.abs(np.linalg.norm(r2_ciclo)-np.linalg.norm(r2))<=0.00001 and np.abs(np.linalg.norm(r2_punto_ciclo)-np.linalg.norm(r2_punto))<=0.00001):
        if(np.abs(np.linalg.norm(r2_ciclo)-np.linalg.norm(r2))>= -0.00001 and np.abs(np.linalg.norm(r2_punto_ciclo)-np.linalg.norm(r2_punto))>= -0.00001):
            finale=True
            print("aaa")
        #print(r2_ciclo,r2)
    else:
        r2_ciclo = r2
        r2_punto_ciclo = r2_punto
        finale=False
        #print("bbb")
r_new = r
r2s_punto = r2_punto
#print(r2)
#print(r2s_punto)
#print(rs,r2s_punto)
"""
Jusqu'ici
"""


#rotar vectores r
def Rotar (Ep,vector):
    #return np.dot(np.matrix([[1,0,0],[0,np.cos(Ep),np.sin(Ep)],[0,-np.sin(Ep),np.cos(Ep)]]),vector)
    return np.dot(np.array([[1,0,0],[0,np.cos(Ep),np.sin(Ep)],[0,-np.sin(Ep),np.cos(Ep)]]),vector)
r= [Rotar(Epsilon,vector) for vector in r_new]
r2n_punto= Rotar(Epsilon,r2s_punto)
r2=r[1]
#print(r2n_punto)
 
"""
ELEMENTOS ORBITALES
"""
def com(x1,x2):
    comun=[]
    for x in x1 and x2:
           comun.append(x)
           return comun
rm=r2
rp=r2n_punto
#print (r2,rp, np.shape(rp))
Norma_rm= np.linalg.norm(rm)
h=np.cross(rm,rp)
v2=np.dot(rp,rp)
#Semieje mayor(a)

a= 1/((2/Norma_rm)-v2)
#print(h[2])

#Excentricidad(e)

#e=np.sqrt((1-(np.linalg.norm(h)**2))/mu*a)
e=np.sqrt((1-((np.linalg.norm(h)**2))/a))
#print(e)
#Inclinación(i)

hz= h[2]
i=np.arccos(hz/np.linalg.norm(h))#entre 0 y 90°

#Longitud del nodo ascendente(O-omega)

hx=h[0]
hy=-h[1]
O1=np.arcsin(hx/(np.linalg.norm(h)*np.sin(i)))
O2=np.arccos(-hy/(np.linalg.norm(h)*np.sin(i)))
#print(O1,O2)
O= com(O1,O2)

#Perihelio(w-omega)
#hallar v True anomaly

v1=np.arccos(((a*(1-squ(e))/Norma_rm)-1)*1/e)
#v2=np.arcsin((np.dot(rm,rp)*a*(1-squ(e)))/np.linalg.norm(h)*Norma_rm)
v2=np.arcsin(((np.dot(rm,rp)*a*(1-squ(e)))/np.linalg.norm(h)*Norma_rm)*1/e)
v= np.rad2deg(com(v1,v2))

#hallar U--revisar x,y,z

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
#Mean anomaly(M)

E=np.arccos((1/e)*(1-(Norma_rm/a)))#---verificar mismo cuadrante que v 
M= E-(e*np.sin(E)) 
#Hay otra manera de M ~revisar

