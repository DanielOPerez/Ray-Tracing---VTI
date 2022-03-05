import numpy as np

#=======================================================================================
# Straight ray
#=======================================================================================

def straight_ray(vp0,vs0,ani,s,r,raytype):

    d=np.sqrt((r[0,0]-s[0,0])**2+
              (r[0,1]-s[0,1])**2+
              (r[0,2]-s[0,2])**2)
 
    if (r[0,2]-s[0,2])==0.0:
        theta=np.pi/2.0
    else:
        theta=np.arctan(((r[0,0]-s[0,0])**2+(r[0.1]-s[0,1])**2)/np.abs(r[0,2]-s[0,2]))
    
    if raytype == 'p':
        v=vp0*(1.0+ani[:,1]*(np.sin(theta)*np.cos(theta))**2+ani[:,0]*np.sin(theta)**4)
    elif raytype == 'sh':
        v=vs0*(1.0+ani[:,2]*np.sin(theta)**2)
    elif raytype == 'sv':
        v=vs0*(1.0+(vp0/vs0)**2*(ani[:,0]-ani[:,1])*(np.sin(theta)*np.cos(theta))**2)
        
    time=np.divide(d,v)
    
    return time

#=======================================================================================
# Refracted ray
#=======================================================================================

def refracted_ray(x, *args):

    esp,vp0,vs0,ani,s,r,raytype=args
    
    dx=np.diff(np.insert(x,(0,np.size(x)),
                         [0,
                          np.sqrt((r[0,0]-s[0,0])**2+
                                  (r[0,1]-s[0,1])**2)]))

    theta=np.arctan(np.divide(dx,esp))
    d=np.sqrt(esp**2+dx**2)

    if raytype == 'p':
        v=vp0*(np.ones(np.size(dx))+
               ani[:,1]*(np.sin(theta)*np.cos(theta))**2+
               ani[:,0]*np.sin(theta)**4)
    elif raytype == 'sh':
        v=vs0*(np.ones(np.size(dx))+ani[:,2]*np.sin(theta)**2)
    elif raytype == 'sv':
        v=vs0*(np.ones(np.size(dx))+
            (vp0/vs0)**2*(ani[:,0]-ani[:,1])*(np.sin(theta)*np.cos(theta))**2)

    time=np.sum(np.divide(d,v))

    return time
