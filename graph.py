import numpy as np
import matplotlib.pyplot as plt
data = np.loadtxt('data.dat')
print np.shape (data)

fig = plt.figure(figsize=(7.5,7.0))     
plt.plot(data[:,0], data[:,1], 'ko')

filename = 'Tres_Cuerpos' 
plt.savefig(filename + '.pdf',format = 'pdf', transparent=True)
plt.show()
plt.close()
