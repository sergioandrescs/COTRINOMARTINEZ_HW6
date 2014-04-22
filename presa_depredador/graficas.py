import numpy as np 
import pylab as plt

datos = np.loadtxt('datos.txt')
for i in range(30):
    x=datos[i*1000:((i*1000)+999),0]
    y=datos[i*1000:((i*1000)+999),1]
    plt.plot(x,y)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('X(0)='+str(30-i)+'\\Y(0)='+str(datos[0,1]))
    plt.savefig('presa_depredador'+str(30-i)+'.png',dpi=200)


