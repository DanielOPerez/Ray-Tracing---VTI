import numpy as np
import scipy.optimize as optimize

from sys import argv
from ray_mod import straight_ray, refracted_ray

#=============================================================================
# Datos de entrada
#=============================================================================
z=np.loadtxt(argv[1], usecols=[0])
vp0=np.loadtxt(argv[1], usecols=[1])
vs0=np.loadtxt(argv[1], usecols=[2])
ani=np.loadtxt(argv[1], usecols=range(3,6))

s=np.loadtxt(argv[2], ndmin=2)
r=np.loadtxt(argv[3], ndmin=2)

n_s=np.size(s,0)
n_r=np.size(r,0)

times=np.zeros([n_s,n_r,3])
ray_name=['p','sh', 'sv']

for i in np.arange(n_s):
    for j in np.arange(n_r):
        #indexo las disc entre la fuente y el receptor
        #revisar, debe haber uan forma mas facil de hacerlo
        ind_min=np.min(np.where(z>np.minimum(s[i,2],r[j,2])))
        ind_max=np.max(np.where(z<np.maximum(s[i,2],r[j,2])))

        ind1=np.arange(ind_min,ind_max+1)
        ind2=np.arange(ind_min,ind_max+2)

        if np.size(ind1)==0:
            #llamar a straigth ray
            for k in np.arange(len(ray_name)):
                times[i,j,k]=straight_ray(vp0[ind2],
                                          vs0[ind2],
                                          ani[ind2,:],
                                          s,r,
                                          ray_name[k])                                                                       
        else:
            #llamar a refracted ray

            if s[i,2] > r[j,2]:
                ind1=np.flipud(ind1)
                ind2=np.flipud(ind2)

            esp=np.abs(np.diff(np.insert(z[ind1],
                                         (0,np.size(ind1)),
                                         [s[i,2],r[j,2]])))    

            x0=np.linspace(0,np.sqrt((r[0,0]-s[0,0])**2+(r[0,1]-s[0,1])**2),
                           np.size(esp)-1)

            args_c=(esp,vp0[ind2],vs0[ind2],ani[ind2,:],s,r)
            for k in np.arange(len(ray_name)):
                args=args_c+(ray_name[k],)
                x_opt=optimize.fmin_cg(refracted_ray,
                                       x0,
                                       fprime=None,
                                       args=args,
                                       gtol=1.e-4,
                                       disp=False)
                times[i,j,k]=refracted_ray(x_opt,*args)
                
                print(times[i,j,k])
