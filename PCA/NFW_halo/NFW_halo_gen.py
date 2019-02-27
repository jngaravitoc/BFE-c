# Hacked version of NFW halo generation code, original from:
# https://github.com/poveda-ruiz/HaloGenerator/blob/master/src/new/main.py

# usage python3 NFW_halo_gen.py 10000 20 300 NFW_n1E4.txt

# to -d0 
# Implement henrquist halo

import numpy as np, pylab, emcee
import sys

n = int(sys.argv[1])
c = int(sys.argv[2])
rcut = int(sys.argv[3])
M = float(sys.argv[4])
filename = sys.argv[5]

def lnprob(r):
    if 1 > r > 0:
        return np.log(r/((1+r*c)*(1+r*c)))
    return -np.inf


def lnprob_hern(r):
    if 1 > r > 0:
        return np.log(r/((1+r*c)*(1+r*c)))
    return -np.inf

ndim, nwalkers = 1, 2
p0 = [[0.01],[0.09]]

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob)
sampler.run_mcmc(p0, n/2)

r = np.concatenate((sampler.chain[0],sampler.chain[1]))
r = np.sort(r)
r = [x[0] for x in r]


theta = np.arccos(2*np.random.random(n)-1)
phi =  2.0 * np.pi * np.random.random(n)

x = r * np.sin(theta) * np.cos(phi) * rcut
y = r * np.sin(theta) * np.sin(phi) * rcut
z = r * np.cos(theta) * rcut

m_part = np.ones(len(x))*M/len(x)

f = open(filename, 'w')
for i in range(len(x)):
    f.write('%f \t %f \t %f \t %f \n'%(x[i], y[i], z[i], m_part[i]))
f.close()

#pylab.plot(x,y,'.k')
#pylab.show()
