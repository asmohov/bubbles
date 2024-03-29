#view slice at y=0 of B field spheromak initial condition

import sys
import time
sys.path.append('/home/asmohov/athena/vis/python/')
sys.path.append('~/.local/lib/python3.8/site-packages/')
sys.path.append('~/working')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import athena_read

#movie making section
sys.path.append('/home/asmohov/athena/vis/python/')
sys.path.append('~/.local/lib/python3.8/site-packages/')
sys.path.append('~/working')
import numpy as np
import time
from matplotlib import pyplot as plt
from celluloid import Camera as cam
import ffmpeg
import matplotlib as mpl
import athena_read
import IPython
from base64 import b64encode
start = time.time()
z=[]
fig, ax = plt.subplots(figsize=(10,20),dpi=50)
Cam = cam(fig)
outrange = range(7501)
file_path = '../mag_bubble'
for i in outrange:
        if z != []:
            z.remove()
        if i<10:
            fname = file_path+'/rt.block0.out2.0000'+str(i)+'.vtk'
        elif 9<i<100:
            fname = file_path+'/rt.block0.out2.000'+str(i)+'.vtk'
        elif 99<i<1000:
            fname = file_path+'/rt.block0.out2.00'+str(i)+'.vtk'
        else:
            fname = file_path+'/rt.block0.out2.0'+str(i)+'.vtk'
        #print(fname)
       #section to make plot
        if i%50== 0:
            data = athena_read.vtk(fname)
            #print(data)
            x_coords = data[0][:-1]
            z_coords = data[2][:-1]
            y_coords = data[1][:-1]#+.03125
            #print(x_coords)
            splice = np.array(data[3]['Bcc'])
            #print(splice.shape)
            splice = splice[:,16,:,:]

            rho_splice = np.array(data[3]['rho'])
            rho_splice = rho_splice[:,:,16]
            #print(rho_splice.shape)
            #print(len(x_coords))
            #print(len(z_coords))
            Bx = splice[:,:,0]
            By = splice[:,:,1]
            Bz = splice[:,:,2]
            #print('len bx ',len(Bx))
            #print('len bz ',len(Bz))
            #print('len by ',len(By))
            Bmag = np.sqrt(Bx*Bx+By*By+Bz*Bz)

            #fig, ax = plt.subplots(figsize=(8,8))
            ax.set_xlim(-1,1)
            ax.set_ylim(0,4)
            ax.pcolormesh(x_coords,z_coords,rho_splice,shading='gouraud',cmap='binary_r')
            q = ax.quiver(x_coords,z_coords,Bx,Bz,Bmag,cmap='plasma')

            #print('snapping')
            Cam.snap()
print('Beginning animation')
anim = Cam.animate(interval=100)
#writervideo = anim.FFMpegWriter(fps=60)
print('done animating')
anim.save('29feb_test.mp4')
#Video("bubble_test.mp4",embed=True)
print('Run time is ',(time.time()-start),' seconds' )
#play video inline
print('done')
