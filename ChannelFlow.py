from pylab import *

close('all')

Re2000 = loadtxt('profiles/Re2000.prof', comments='%')

yt = Re2000[:,0]
yp = Re2000[:,1]
Up = Re2000[:,2]
up = Re2000[:,3]
vp = Re2000[:,4]
wp = Re2000[:,5]
uvp = Re2000[:,10]
uwp = Re2000[:,11]
vwp = Re2000[:,12]

dUdy = zeros(size(yt))

figure()
plot(yp, up)
plot(yp, vp)
plot(yp, wp)

figure()
plot(yp, uvp)

show()
