import numpy as np 
import pylab as plt

datos = np.loadtxt('datos.txt')
for i in range(30):
    x=datos[i*1000:((i*1000)+999),0]
    y=datos[i*1000:((i*1000)+999),1]
    plt.scatter(x,y)
    plt.xlabel('Presas')
    plt.ylabel('Cazadores')
    plt.title('Presas iniciales:'+str(30-i)+', Cazadores iniciales:'+str(int(datos[0,1])))
    plt.savefig('presa_depredador'+str(30-i)+'.png',dpi=200)
    plt.show()


