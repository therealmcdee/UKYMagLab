import numpy as np
import matplotlib.pyplot as plt


## DISTANCES ARE IN CENTIMETERS

N = 500              # number of points
ny = 20              # number of points in y-direction
nx = 5               # number of points in one direction on the plane (squared gives you number in plane)
Lytot = 2*15.24      # target length to map 
dx = 1               # furthest point on plane to go to  

sigfig = 2           # the mapper's lowest resolution is 1/750 cm^(-1) ~ 0.001 cm


dy = round(ny/Lytot, sigfig)

x = np.linspace(-dx, dx, nx)
y = np.linspace(0, Lytot, ny)


mesh = np.zeros((N, 3))
cnt = 0
for i in range(ny):
    for j in range(nx):
        for k in range(nx):
            if (i%2)==0:
                if (j%2)==0:
                    mesh[cnt][0] = x[j]
                    mesh[cnt][1] = y[i]
                    mesh[cnt][2] = x[k]
                else:
                    mesh[cnt][0] = x[j]
                    mesh[cnt][1] = y[i]
                    mesh[cnt][2] = -x[k] 
            else:
                if (j%2)==0:
                    mesh[cnt][0] = -x[j]
                    mesh[cnt][1] = y[i]
                    mesh[cnt][2] = -x[k]
                else:
                    mesh[cnt][0] = -x[j]
                    mesh[cnt][1] = y[i]
                    mesh[cnt][2] = x[k]
            cnt += 1

with open('stc_coilb_map_positions.txt', 'w+') as outf:
    for i in range(len(mesh)):
        outf.write(f'{round(mesh[i][0],sigfig)} {round(mesh[i][1],sigfig)} {round(mesh[i][2],sigfig)}\n')


fig = plt.figure()
ax = fig.add_subplot(projection='3d')

ax.scatter(mesh[:,0], mesh[:,1], mesh[:,2])

fig2 = plt.figure()
plane0 = mesh[(mesh[:,1]==0)]
plane1 = mesh[(mesh[:,1]==mesh[1,1])]
plt.plot(plane0[:,0], plane0[:,2])
print(plane0[-1])
print(plane1[0])


plt.show()
