import numpy as np

fil = open('test.txt','w')
x = np.zeros((8))
for i in range(8):
    x[i] = np.random.uniform(0,1)

fil.write(str(x[0])+' '+str(x[1])+' '+str(x[2])+' '+str(x[3])+'\n')
fil.write(str(x[0])+' '+str(x[1])+' '+str(x[2])+' '+str(x[3]))
fil.close()
