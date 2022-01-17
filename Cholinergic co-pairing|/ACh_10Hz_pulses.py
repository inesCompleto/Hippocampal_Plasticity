import sk_dsp_comm.sigsys as ss
import scipy.signal as signal
import numpy as np
import matplotlib.pyplot as plt


def periodic_rect(t,init,tau,T):
    """
    The period is set by T and tau is the pulse width
    """
    x = np.zeros(len(t))
    x += ss.rect(np.mod(t - init,T),tau*2)
    return x

t = np.arange(0,2000,.02)
x = periodic_rect(t,0,5,100)

plt.figure()
plt.plot(t,x)
plt.plot(t,x)
plt.grid()
#xlim([-5,5])
plt.xlabel(r'Time (s)')
plt.ylabel(r'$x_2(t)$ and $x_2^2(t)$');



x2 = periodic_rect(t,0,405,700)
x1 = x* x2
#ss.rect(t-405/2,405)

plt.figure()
plt.subplot(211)
plt.plot(t,x1)
plt.xlabel(r'Time (s)')
plt.ylabel(r'$x_1(t)$')



plt.subplot(212)
plt.plot(t,x2)
plt.xlabel(r'Time (s)')
plt.ylabel(r'$x_2(t)$')
plt.grid()
plt.tight_layout()