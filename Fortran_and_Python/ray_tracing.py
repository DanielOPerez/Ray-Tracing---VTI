import numpy as np

from cffi import FFI
from sys import argv


#=============================================================================
# FFI  python <--> C <--> fortran
#=============================================================================
ffi = FFI()
lib = ffi.dlopen("./libray_tracing.so")
#For each function/subroutine you want to call from Python,
#you should provide a C-signature to CFFI:
ffi.cdef("void ray_tracing_wrap(double*, double*, double*, double*, \
double*, double*, int, int, int, double*, double*, double*);")


#=============================================================================
# Datos de entrada
#=============================================================================
z=np.loadtxt(argv[1], usecols=[0])
vp0=np.loadtxt(argv[1], usecols=[1])
vs0=np.loadtxt(argv[1], usecols=[2])
ani=np.loadtxt(argv[1], usecols=range(3,6))

s=np.loadtxt(argv[2], ndmin=2)
r=np.loadtxt(argv[3], ndmin=2)

n_d=np.size(z,0)
n_s=np.size(s,0)
n_r=np.size(r,0)

tp_est=np.zeros([n_r,n_s])
tsh_est=np.zeros([n_r,n_s])
tsv_est=np.zeros([n_r,n_s])


#np.array --> pointer
ffi_z=ffi.cast("double*", z.ctypes.data)
ffi_vp0=ffi.cast("double*", vp0.ctypes.data)
ffi_vs0=ffi.cast("double*", vs0.ctypes.data)
ffi_ani=ffi.cast("double*", ani.ctypes.data)
ffi_r=ffi.cast("double*", r.ctypes.data)
ffi_s=ffi.cast("double*", s.ctypes.data)
ffi_tp_est=ffi.cast("double*", tp_est.ctypes.data)
ffi_tsh_est=ffi.cast("double*", tsh_est.ctypes.data)
ffi_tsv_est=ffi.cast("double*", tsv_est.ctypes.data)


lib.ray_tracing_wrap(ffi_z, \
                     ffi_vp0, \
                     ffi_vs0, \
                     ffi_ani, \
                     ffi_s, \
                     ffi_r, \
                     n_d, \
                     n_s, \
                     n_r, \
                     ffi_tp_est, \
                     ffi_tsh_est, \
                     ffi_tsv_est)


np.savetxt('tp_est',tp_est)
np.savetxt('tsh_est',tsh_est)
np.savetxt('tsv_est',tsv_est)




