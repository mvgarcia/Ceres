import numpy as np
import sympy as sp
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
G=1.48814*10**-34#AU**3/kg*día**2~~ 6.67408*10**-11 #N*m**2/kg**2
mu=G*(1.989*10**30 + 5.972*10**24)#kg
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

RA=np.array([[11, 36, 6.30],[11, 40, 3.62],[11, 43, 10.25]])
DEC=np.array([[2, 35, 0.5],[2, 9, 27.1],[1, 49, 19]])
data=np.array([[np.deg2rad(15*decimales(x[0],x[1],x[2])) for x in RA],[np.deg2rad(decimales(x[0],x[1],x[2])) for x in DEC]])#AR,DC
transformed_data=[data[:,0],data[:,1],data[:,2]]

vectores2=[list(vector(i[0],i[1]))for i in transformed_data]
R1=np.array(vectores2[0])
R2=np.array(vectores2[1])
R3=np.array(vectores2[2])
#print(R2)
'''
def rTS (angulo):
    return(R*(1-squ(Exct))/1+(Exct*np.cos(angulo)))
R1=rTS(d(14))#--magnitud 
R1=np.array([R1*np.cos(d(14)),R1*np.sin(d(14)),0])
R2=rTS(d(15))
R2=np.array([R2*np.cos(d(15)),R2*np.sin(d(15)),0])
R3=rTS(d(16))
R3=np.array([R3*np.cos(d(16)),R3*np.sin(d(16)),0])
#print(R3,np.linalg.norm(R3))----[],#
'''
#tiempo en kt (gaussian)

T=[2458191.8942,2458192.9733,2458193.82275]

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
#print(t)---[]
#valor inicial de r2 (vector posicion de la observación central)
#print(p1)
#print(p2)--np.array
#print(p3)
#print(R1)--np.array
Dinicial=np.dot(p1,np.cross(p2,p3))
#print(Dinicial) ---#
'''def DXx (Rx):
    if(np.all(Rx==R1)):
        return np.dot(np.cross(R1,p2),p3)
    elif (np.all(Rx==R2)):
        return np.dot(np.cross(p1,R2),p3)
    elif (np.all(Rx==R3)):
        return np.dot(p1,np.cross(p2,R3))
    '''
def D1x (j):
    return np.dot(np.cross(j,p2),p3)
def D2x (j):
    return np.dot(np.cross(p1,j),p3)
def D3x (j):
    return np.dot(p1,np.cross(p2,j))
#print(DXx(R1))---#
"""A1=t[2]/t[0]#---#
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

r_2=sp.Symbol('r', real=True)
equacion=r_2**8+aa*(r_2**6)+b*(r_2**3)+cc
r2=sp.solve(equacion,r_2)
r2 = np.array(r2)
r2_inicial = r2[r2>0]
"""
a1=t[2]/t[0]#---#
a2=-1
a3=-t[1]/t[0]#---#

P1=(a1*D1x(R1)+a2*D1x(R2)+a3*D1x(R3))/a1*Dinicial
P2=(a1*D2x(R1)+a2*D2x(R2)+a3*D2x(R3))/a2*Dinicial
P3=(a1*D3x(R1)+a2*D3x(R2)+a3*D3x(R3))/a3*Dinicial
#print(P1)#---# 2.9212279935443062e-06

r1o=np.array((p1*P1)-R1)
r2o=np.array((p2*P2)-R2)
r3o=np.array((p3*P3)-R3)
#print(r2o)#---vector [ 2.81140060e-01 -9.54691768e-01  1.52068239e-06]
r2o_punto=(r3o-r1o)/t[0]
#print(r2o_punto)#  [-363.0428774815191, -1320.8112113396828, 0.0027111788399157956] - con np.array--esto[-3.63042877e+02 -1.32081121e+03  2.71117884e-03]
rs=[]

def f (r2,r2_punto,kt):
    return 1-(mu*kt/2*(np.linalg.norm(r2)**3))+((kt**3)*mu*(np.dot(r2,r2_punto))/2*(np.linalg.norm(r2)**5))+(((kt**4)/24*(np.linalg.norm(r2)**3))*((3*((np.dot(r2_punto,r2_punto)/(np.linalg.norm(r2)**2))-(1/(np.linalg.norm(r2)**3))))-(15*((np.dot(r2,r2_punto)/(np.linalg.norm(r2)**2))**2))+(1/(np.linalg.norm(r2)**3))))
def g(r2,kt):
    return kt-(kt**3)*mu/6*(np.linalg.norm(r2)**3)-((kt**4)*(np.dot(r2,r2_punto)/4*(np.linalg.norm(r2)**5)))
# Iterate
    '''
rs=[[],[],[]]
  
for i in range(len(r2_inicial)):
    x = r2_inicial[i]
    #r2=x    
    R2_punto=(R3-R2)/k*(t[2]-t[1])#--¿?Método de Lagrange
    p2_punto=(squ(t[2])*(p1-p2)-squ(t[0])*(p3-p2))/(t[0]*t[1]*t[2]) 
    p2_2punto=-2*((t[2]*(p1-p2)-t[0]*(p3-p2))/t[0]*t[1]*t[2])
    P2_punto=-(1/2)*((1/x**3)-((1+(1/328900.5))/np.linalg.norm(R2)**3))*((np.dot(np.cross(p2,p2_2punto),R2))/(np.dot(np.cross(p2,p2_punto),p2_2punto))) 
    AA=np.dot(np.cross(p2,p2_punto),R2)/np.dot(np.cross(p2,p2_punto),p2_2punto)
    BB=((1+(1/328900.5))/np.linalg.norm(R2)**3)*AA
    P2= (AA/x**3)-BB
    #print(P2)--#
    r2o=(P2*p2)-R2
    print(r2o)
    dr2=(P2_punto*p2)+(p2*p2_punto)-R2_punto
    print(dr2)
    #print(P2_punto)---#
    finale= False
    r = []
    #print (np.shape(dr2))
    while finale==False: 
        otro_r2=r2o
    #truncated f & g        
        
        f=[f(otro_r2,dr2,kt)for kt in t]
        f1=f[0]
        f3=f[2]
        #print(np.shape(f),f)
        
        g=[g(otro_r2,kt) for kt in t]
        g1=g[0]
        g3=g[2]
        #print(g)
        
        #dr2
        r1=f1*otro_r2+g1*dr2
        r3=f3*otro_r2+g3*dr2
        d1=-f3/(f1*g3-f3*g1)
        d3=f1/(f1*g3-f3*g1)
        dr2=d1*r1+d3*r3
        #dr2=sp.solve()
        #hallar c1 y c3
        
        c1= g3/((f1*g3)-(g1*f3))
        c2=-1
        c3= -g1/((f1*g3)-(g1*f3)) 
        
        #hallar vectores posicion entre tierra y asteroide
        P1=(c1*D1x(R1)+c2*D1x(R2)+c3*D1x(R3))/c1*Dinicial
        P2=(c1*D2x(R1)+c2*D2x(R2)+c3*D2x(R3))/c2*Dinicial
        P3=(c1*D3x(R1)+c2*D3x(R2)+c3*D3x(R3))/c3*Dinicial
        #print(P2)
        
        #hallar posicion sol asteroide 
        r1=P1-R1
        r2=P2-R2
        r3=P3-R3
        r=[r1,r2,r3]
        #print(r2)
        #hallar r y r.
        #Corrección tiempo de la luz ~~en teoría desde aquí ya empieza
        T1= T1-P1/c
        T2= T2-P2/c
        T3= T3-P3/c
        
        t=tiempos(T1,T2,T3)
        if(np.abs(np.linalg.norm(otro_r2)-np.linalg.norm(r2))<=0.001):
            finale=True
        else:
            otro_r2=r2
    rs[i] = r 
    
'''

#print(r2o)
#print(r2o_punto)
finale=False
while finale==False: 
    r2_ciclo=r2o
    r2_punto_ciclo= r2o_punto
    
    #Corrección tiempo de la luz ~~en teoría desde aquí ya empieza
    T1= T1-(P1/c)
    T2= T2-(P2/c)
    T3= T3-(P3/c)
    
    t=tiempos(T1,T2,T3)
    #truncated f & g        
    
    f=[f(r2_ciclo,r2_punto_ciclo,kt)for kt in t]
    f1=f[0]
    f3=f[2]
    #print(f)
    
    g=[g(r2_ciclo,kt) for kt in t]
    g1=g[0]
    g3=g[2]
    #print(g)
    """ 
    #dr2
    r1=f1*r2_ciclo+g1*r2_punto
    r3=f3*r2_ciclo+g3*r2_punto
    d1=-f3/(f1*g3-f3*g1)
    d3=f1/(f1*g3-f3*g1)
    dr2=d1*r1+d3*r3
    #dr2=sp.solve()
    """
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

    """
    r1=P1-R1
    r2=P2-R2
    r3=P3-R3
    """
    r=[r1,r2,r3]
    #print(r2)
    #hallar r y r.
    #print(r2_ciclo,r2)
    if(np.abs(np.linalg.norm(r2_ciclo)-np.linalg.norm(r2))<=0.0001 and np.abs(np.linalg.norm(r2_punto_ciclo)-np.linalg.norm(r2_punto))<=0.0001):
        if(np.abs(np.linalg.norm(r2_ciclo)-np.linalg.norm(r2))>= -0.0001 and np.abs(np.linalg.norm(r2_punto_ciclo)-np.linalg.norm(r2_punto))>= -0.0001):
            finale=True
            print("aaa")
        #print(r2_ciclo,r2)
    else:
        r2_ciclo = r2
        r2_punto_ciclo = r2_punto
        finale=False
        print("bbb")
rs = r
r2s_punto = r2_punto
print(r2)
print(r2_punto)
#print(rs,r2s_punto)
"""
Jusqu'ici
"""


#rotar vectores r
def Rotar (Ep,vector):
    return np.dot(np.matrix([[1,0,0],[0,np.cos(Ep),np.sin(Ep)],[0,-np.sin(Ep),np.cos(Ep)]]),vector)
r=[Rotar(Epsilon,vector) for vector in rs]
r2n_punto=[Rotar(Epsilon,r2s_punto)]
r2=r[1]
 
"""
ELEMENTOS ORBITALES
"""
def com(x1,x2):
    comun=[]
    for x in x1 and x2:
           comun.append(x)
           return (comun)
rm=r2
rp=r2n_punto
print (r2,rp[0], np.shape(rp))
Norma_rm= np.linalg.norm(rm)
h=np.cross(rm,rp)
v2=np.dot(rp,rp)
#Semieje mayor(a)
"""
a=sp.Symbol('a')
a=(2/Norma_rm)-(np.dot(rp,rp)/mu)**-1
"""
a=1/((2/rm)-v2)
#P=2*np.pi*a**(3/2) #~~Verificar donde se usa

#Excentricidad(e)

#e=np.sqrt((1-(np.linalg.norm(h)**2))/mu*a)
e=np.sqrt((1-(np.linalg.norm(h)**2))/a)

#Inclinación(i)

hz= h[3]
i=np.arccos(hz/np.linalg.norm(h))#entre 0 y 90°

#Longitud del nodo ascendente(O-omega)

hx=h[0]
hy=-h[1]
O1=np.arcsin(hx/np.linalg.norm(h)*np.sin(i))
O2=np.arccos(-hy/np.linalg.norm(h)*np.sin(i))
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
    
'''
    ¿?

#f y g otra vez
n=np.sqrt(mu/a**3)
def fi(r2,DEi):#DEi= delta de la anomalía excéntrica (Ei-E.observación central) 
    return(1-((a/r2)*(1-np.cos(DEi))))
def gi(ti,DEi):
    return((ti-T2)+(1/n)*np.sin(DEi)-DEi)
'''