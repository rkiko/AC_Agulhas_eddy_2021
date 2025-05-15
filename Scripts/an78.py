import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
from pathlib import Path
home = str(Path.home())

x0,xf=-1000,1001
x=np.r_[-1000:1001]
mean=48
sigma=75
y=1/(sigma*np.sqrt(2*np.pi))*np.e**(-(x-mean)**2/(2*sigma**2))
prob_greater_than_0 = (1 - norm.cdf(0, mean, sigma))*100

plt.figure(figsize=(8, 5))
plt.plot(x, y, color='blue')
plt.vlines(0,ymin=-1,ymax=1)
plt.xlabel('X')
plt.ylabel('Density')
plt.xlim(x0,xf)
plt.ylim(0,1.05/(sigma*np.sqrt(2*np.pi)))
plt.title('%d$\pm$%d, %d%% of values are > 0' % (mean,sigma,prob_greater_than_0))
plt.grid()
plt.savefig("%s/GIT/AC_Agulhas_eddy_2021/Plots/an78/Gauss_%dpm%d.pdf" % (home,mean,sigma), dpi=200)
plt.close()


mean=71
sigma=354
y=1/(sigma*np.sqrt(2*np.pi))*np.e**(-(x-mean)**2/(2*sigma**2))
prob_greater_than_0 = (1 - norm.cdf(0, mean, sigma))*100

plt.figure(figsize=(8, 5))
plt.plot(x, y, color='blue')
plt.vlines(0,ymin=-1,ymax=1)
plt.xlabel('X')
plt.ylabel('Density')
plt.xlim(x0,xf)
plt.ylim(0,1.05/(sigma*np.sqrt(2*np.pi)))
plt.title('%d$\pm$%d, %d%% of values are > 0' % (mean,sigma,prob_greater_than_0))
plt.grid()
plt.savefig("%s/GIT/AC_Agulhas_eddy_2021/Plots/an78/Gauss_%dpm%d.pdf" % (home,mean,sigma), dpi=200)
plt.close()



