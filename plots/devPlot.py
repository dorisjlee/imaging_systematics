import  matplotlib.pyplot as plt
import numpy as np
V = np.array([[0.,0.],[0.,0.]])
def N(xgrid,ygrid,m,Vinv):
  """
  2d normal (Gaussian) distribution.
  #http://en.wikipedia.org/wiki/Gaussian_function#Two-dimensional_Gaussian_function
  xgrid,ygrid is the output of a meshgrid.
  #m is the mean (x_0,y_0)
  Vinv is a 2x2 np.array, you already took the inverse!
  """
  assert Vinv[0,1] == Vinv[0,1]
  exparg = (xgrid - m[0])*Vinv[0,0]*(xgrid - m[0]) + \
           2*(xgrid - m[0])*Vinv[1,0]*(ygrid -m[1]) + \
           (ygrid - m[1])*Vinv[1,1]*(ygrid - m[1])
  
  return (np.linalg.det(Vinv))**0.5/(2.*np.pi)*np.exp(-0.5*exparg)
nx = 2048
ny = 1361
psize = 0.396 # arcseconds per pixel.

# coordinates of center of pixels
x = (np.arange(0,nx,1)+0.5)*psize
y = (np.arange(0,ny,1)+0.5)*psize

xM, yM = np.meshgrid(x,y,indexing='xy')
def Iconv(mdev):
    I_conv = np.zeros((1361, 2048))
    for a in mdev.keys():
        v=mdev[a]
        m = np.array([100,300])
        v_inv = np.linalg.inv(np.array([[v**2,0],[0,v**2]]))
        I_conv += a*N(xM,yM,m,v_inv)
    return I_conv
mdev8 = {0.00262:0.00113, 0.02500:0.00475,0.13413:0.01462,0.51326:0.03930,1.52005:0.09926,3.56204:0.24699,6.44845:0.63883,8.10105:1.92560}
I_conv= Iconv(mdev8)
# plotI(I_conv8)
n=0
rlist = []
Iconvlst = []
m = np.array([100,300]) #center
# too large can not run on ipynb
for i in np.arange(I_conv.shape[0]-1):
    for j in np.arange(I_conv.shape[1]-1):
        print(i,j)
        r = np.sqrt(abs(i-m[0])**2+abs(j-m[1])**2) #treating center as origin
#         print I_conv[i][j]
        rlist.append(r)
        Iconvlst.append(I_conv[i][j])
print rlist
print Iconvlst
plt.figure()
plt.plot(rlist,Iconvlst)
plt.show()
plt.figure()
plt.loglog(rlist,Iconvlst)
plt.show()