import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import pandas as pd
def dydt(y,t):
    x=y[0]
    h=y[1]
    W=y[2]
    V=y[3]
    G=y[4]
    g=9.81 #cte de gravedad
    cp=1.33E-06 #consumo especifico   
    n=0.7 #rendimiento de la helice
    alpha=0.026 #angulo de ataque
    CL=0.56 #cl max crucero
    S=33.93 #Superficie alar
    r=1.225*(1-h*22.557E-6)**4.256 #variacion de la densidad en funcion de la altura
    P0=652400*(r/1.225) #potencia maxima en el eje
    L=0.5*S*CL*r*V**2 #Sustentacion
    Cd=0.0290-0.00593*CL+0.04464*CL**2 #Coef de resistencia crucero
    D=0.5*r*S*Cd*V**2 #Resistencia
    Nm=2 #numero de motores
    ###############
    #Ecuaciones#
    dxdt=V*np.cos(G)
    dhdt=V*np.sin(G)
    dWdt=-cp*Nm*P0/n
    dvdt=(g/W)*((Nm*n*P0/V)*np.cos(alpha)-D-W*np.sin(G)-V*dWdt)
    dGdt=(g/(W*V))*((Nm*n*P0/V)*np.sin(alpha)+L-W*np.cos(G))
   
    return np.array([dxdt,dhdt,dWdt,dvdt,dGdt])

    
    
#vector tiempo
t=np.linspace(0,2400,100000)
#condiciones iniciales
y0=np.array([0,500,5627*9.81,70,0.052])

#integracion de las ecuaciones
sol=odeint(dydt,y0,t)
X=sol[:,0]
H=sol[:,1]
W=sol[:,2]
v=sol[:,3]
G=sol[:,4]

#grafico derivadas
S=33.93 #superficie
r=1.225*(1-H*22.557E-6)**4.256 #densidad
cp=1.33E-06#consumo especifico
P0=652400*r/1.225 #potencia
g=9.81 #cte de gravedad
n=0.7 #rendimiento de la helice
alpha=0.026 #angulo de ataque
CL=0.56 #cl max crucero
Cd=0.0290-0.00593*CL+0.04464*CL**2 #Coef de resistencia crucero
D=0.5*r*S*Cd*v**2 #resistencia
L=0.5*S*CL*r*v**2 #sustentacion

dxdt=v*np.cos(G)
dhdt=v*np.sin(G)
dWdt=-cp*P0/n
dvdt=(g/W)*((n*P0/v)*np.cos(alpha)-D-W*np.sin(G)-v*dWdt)
dGdt=(g/(W*v))*((n*P0/v)*np.sin(alpha)+L-W*np.cos(G))

#factor de carga
N0=L/W
#consumo de combustible
Peso=5627*9.81
Wmax=18491.85
Pw=((Peso-W)/W)
#evolvente
datos=pd.read_csv('evolvente.csv',names=['h','Vmax','Vmin','Vs'])
x=datos.iloc[:,1]
y=datos.iloc[:,0]
x1=datos.iloc[:,2]
x3=datos.iloc[:,3]

plt.xlabel('V(m/s)')
plt.ylabel('H(m)')
plt.plot(x,y,label='Vmin')
plt.plot(x1,y,label='Vmax')
plt.plot(v,H,label='Equilibrio')
plt.plot(x3,y,label='Vs')
plt.legend()
plt.savefig('evolvente.png')
plt.show()

#graficos de las derivadas
fig1, axs=plt.subplots(3,2,figsize=(15,15))


ax=axs[0,0]
ax.plot(t,X)
ax.set_ylabel('X(m)')
ax.set_xlabel('t(s)')

ax=axs[0,1]
ax.plot(t,H)
ax.set_ylabel('H(m)')
ax.set_xlabel('t(s)')

ax=axs[1,0]
ax.plot(t,W)
ax.set_ylabel('W(N)')
ax.set_xlabel('t(s)')

ax=axs[1,1]
ax.plot(t,v)
ax.set_ylabel('V(m/s)')
ax.set_xlabel('t(s)')

ax=axs[2,0]
ax.plot(t,G)
ax.set_ylabel('gamma(deg)')
ax.set_xlabel('t(s)')


fig1.tight_layout()


plt.savefig('simulacion.png')
#grafico  las variables de estado en funcion del tiempo
fig2, axs=plt.subplots(3,2,figsize=(15,15))
ax=axs[0,0]
ax.plot(t,dxdt)
ax.set_ylabel('X(m/s)')
ax.set_xlabel('t(s)')

ax=axs[0,1]
ax.plot(t,dhdt)
ax.set_ylabel('H(m/s)')
ax.set_xlabel('t(s)')

ax=axs[1,0]
ax.plot(t,dWdt)
ax.set_ylabel('W(N/s)')
ax.set_xlabel('t(s)')

ax=axs[1,1]
ax.plot(t,dvdt)
ax.set_ylabel('V(m/s^2)')
ax.set_xlabel('t(s)')

ax=axs[2,0]
ax.plot(t,dGdt)
ax.set_ylabel('gamma(deg)')
ax.set_xlabel('t(s)')

plt.savefig('simulacion2.png')

#grafico factor de carga y %consumido de combustible
f, axarr = plt.subplots(2, sharex=True,figsize=(15,15))
axarr[0].plot(t, N0)
axarr[0].set_ylabel('Fc')
axarr[0].set_xlabel('t(s)')
axarr[1].plot(t, Pw)
axarr[1].set_ylabel('W%')
axarr[1].set_xlabel('t(s)')
plt.savefig('simulacion3.png')


plt.show()
