from pylab import *

close('all')

Re = 2003
nu = 1/0.485e+5
u_tau = 0.41302030e-1
delta = Re*nu/u_tau
kappa =0.4

prof = loadtxt('profiles/Re2000.prof', comments='%')
bal = loadtxt('balances/Re2000.bal.kbal', comments='%')

yt = prof[:,0]
yp = prof[:,1]
Up = prof[:,2]
up = prof[:,3]
vp = prof[:,4]
wp = prof[:,5]
uvp = prof[:,10]
uwp = prof[:,11]
vwp = prof[:,12]

epsilon = -bal[:,2]
p_diff = bal[:,5]
t_diff = bal[:,6]

k = 0.5*(up*up+vp*vp+wp*wp)*u_tau**2

N = size(yp)

# figure()
# plot(yp, log(yp)/kappa+5)
# plot(yp, Up)

dUdyp = zeros(N)
dkdyp = zeros(N)
for i in range(0,N):
  if i > 0 and i < N-1:
    dym = yp[i+1]-yp[i]
    dyp = yp[i]-yp[i-1]
    dUdyp[i] = Up[i+1]*dym/((dym+dyp)*dyp) \
        + Up[i]*(dyp-dym)/(dyp*dym) \
        - Up[i-1]*dyp/((dym+dyp)*dym)
    dkdyp[i] = k[i+1]*dym/((dym+dyp)*dyp) \
        + k[i]*(dyp-dym)/(dyp*dym) \
        - k[i-1]*dyp/((dym+dyp)*dym)
  elif i > 0:
    dy = yp[i]-yp[i-1]
    dUdyp[i] = (Up[i]-Up[i-1])/dy
    dkdyp[i] = (k[i]-k[i-1])/dy
  elif i < N-1:
    dy = yp[i+1]-yp[i]
    dUdyp[i] = (Up[i+1]-Up[i])/dy
    dkdyp[i] = (k[i+1]-k[i])/dy

dUdy = u_tau*Re/delta*dUdyp
dkdy = Re/delta*dkdyp
l_m = kappa*delta*yt
l_exact = sqrt(abs(uvp*u_tau**2/dUdy**2))

figure()
plot(yp, l_m)
plot(yp, l_exact)

uu_exact = up*up*u_tau**2-2/3*k
vv_exact = vp*vp*u_tau**2-2/3*k
ww_exact = wp*wp*u_tau**2-2/3*k
uv_exact = uvp*u_tau**2
uw_exact = uwp*u_tau**2
vw_exact = vwp*u_tau**2

nu_t_ke = 0.09*k**2/epsilon
uu_ke = zeros(N)
vv_ke = zeros(N)
ww_ke = zeros(N)
uv_ke = -nu_t_ke*dUdy
uw_ke = zeros(N)
vw_ke = zeros(N)

figure()
plot(yp, uu_exact, label='uu', color='b')
plot(yp, vv_exact, label='vv', color='r')
plot(yp, ww_exact, label='ww', color='g')
plot(yp, uv_exact, label='uv', color='m')
plot(yp, uw_exact, label='uw', color='c')
plot(yp, vw_exact, label='vw', color='y')

plot(yp, uu_ke, linestyle='--', label='uu', color='b')
plot(yp, vv_ke, linestyle='--', label='vv', color='r')
plot(yp, ww_ke, linestyle='--', label='ww', color='g')
plot(yp, uv_ke, linestyle='--', label='uv', color='m')
plot(yp, uw_ke, linestyle='--', label='uw', color='c')
plot(yp, vw_ke, linestyle='--', label='vw', color='y')
legend()

figure()
transport_exact = (0.5*t_diff + p_diff)*u_tau**3
transport_ke = nu_t_ke*dkdy
plot(yp, transport_exact, label='exact')
plot(yp, transport_ke, label='ke')
legend()

show()
