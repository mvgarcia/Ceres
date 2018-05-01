import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
from scipy.optimize import curve_fit

def dec(a,b,c):
    return(a+b/60+c/3600)

def inverso(d):
    g=np.floor(d)
    md=(d-g)*60
    m=np.floor(md)
    sd=(md-m)*60
    return np.array([g,m,sd])

elsecoord=np.genfromtxt('campo.dat')
coordenadas=[]
for i in elsecoord:
    ar=15*dec(i[0],i[1],i[2])
    dc=dec(i[3],i[4],i[5])
    coordenadas.append(np.array([ar,dc]))
elsepx=np.genfromtxt('pixel.dat')
cerespx=np.genfromtxt('ceres.dat')[:,1:]

'''cerespx=[]
for i in elsecoord:
    x=i[1]
    y=i[2]
    cerespx.append(np.array([x,y]))
'''
distgrados=[]
for i in combinations(range(len(coordenadas)),2):
    e1=coordenadas[i[0]]
    e2=coordenadas[i[1]]
    distgrados.append(np.linalg.norm(e1-e2))
    
distpx=[]
for i in combinations(range(len(elsepx)),2):
    e1=elsepx[i[0]]
    e2=elsepx[i[1]]
    distpx.append(np.linalg.norm(e1-e2))

#plt.scatter(distpx,distgrados)
#plt.xlabel('pixeles')
#plt.ylabel('coordenadas')
#plt.show()

#def ajuste(a,t):
    #return a*t
'''
ajuste
'''
ajuste=lambda t,a:a*t
m,cov=curve_fit(ajuste,distpx,distgrados)
m=m[0] #pendiente
px=elsepx*m 
px=-px

'''
Matriz de rotación
'''
def RotationMatrix(angle):
    return np.matrix([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])

def Rotate(vector,angle):
    return np.array(np.dot(RotationMatrix(angle),vector))[0]
'''
centros de masa
'''
#cmpx=np.sum(elsepx,0)/len(elsepx)
cmpx=np.sum(px,0)/len(px)
cmcoord=np.sum(coordenadas,0)/len(coordenadas)
#print(cmcoord)
'''
angulos de coordenadas con su centro de masa
'''
anglec=[]
for i in coordenadas:
    vector=i-cmcoord
    angles=np.arctan2(vector[1],vector[0])
    anglec.append(angles)
#print(anglec)
'''
angulos de pixeles con su centro de masa
'''
anglepx=[]
#for i in elsepx:
for i in px:
    vector=np.subtract(i,cmpx)
    angles=np.arctan2(vector[1],vector[0])
    anglepx.append(angles)
#print(anglepx)
'''
restan
'''
a1=[np.subtract(anglec,anglepx)]
a1=a1[0]#####
#print(a1)

'''
Arreglar ángulos
'''
a11=[]
for i in a1:
    a11.append(np.deg2rad(np.rad2deg(i)%360))

#print(a11)
#plt.figure()
#plt.plot(np.rad2deg(a11)-360)

a1=a11


'''
promedio
'''
angle=np.mean(a1)
'''
angulo
'''
#angle=angle[0]
#print (np.rad2deg(angle)-360)
'''
rotar px con angulo
'''
px_rotado=[]
for i in px:
    px_rotado.append(Rotate(i,angle))
px=px_rotado
'''
resta de coordenadas con los pixeles rotados
'''
c=[np.subtract(coordenadas,px_rotado)]
c=c[0]
#print(c)
'''
promedio de la diferencia de coordenadas a pixeles

h es la diferencia en AR
d es la diferencia en DEC
'''
h=[np.sum(c[:,0])/len(c[:,0])][0]
d=[np.sum(c[:,1])/len(c[:,1])][0]

'''
funcion transformar general
'''
def Transform(pixeles):
    R=pixeles
    R=-m*R
    Ro=np.array([Rotate(i,angle)for i in R])
    R=Ro
    RR=[i+np.array([h,d]) for i in R]
    return RR

T=Transform(elsepx)
Errores=abs(np.array([T[i]-coordenadas[i] for i in range(len(T))]))

ErrAR=np.mean(Errores[:,0])
ErrDEC=np.mean(Errores[:,1])
#print(inverso(ErrAR/15),inverso(ErrDEC))

ce=Transform(cerespx)

for i in ce:
    print(inverso(i[0]/15),inverso(i[1]))

