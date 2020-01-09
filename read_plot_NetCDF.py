# Author: JingBang Liu

import netCDF4 as NC
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

# This function transfer coordinate of a matrix to normal real life coordinate
def transfer_axes(M):
	s = M.shape
	N = np.zeros((s[1],s[0]))
	for i in range(s[1]):
		for j in range(s[0]):
			N[i,j] = M[j,s[1]-1-i]
	return N

# Import data from netCDF file
dat = NC.Dataset("electrostatics_data.nc","r",format="NETCDF4")

# Extract data that we are going to use
Ex = dat.variables['Ex_field_intensity']
Ey = dat.variables['Ey_field_intensity']
pos = dat.variables['particle_position']

# create axis for pcolor
x_axis = np.linspace(-1,1,100)
y_axis = np.linspace(-1,1,100)

# show two plots in one window
fig, (ax0,ax1) = plt.subplots(1,2)

# plot
c = ax0.pcolor(x_axis,y_axis,transfer_axes(Ex),cmap='RdBu',)
fig.colorbar(c,ax=ax0)
ax0.set_title('Ex')
ax0.set_xlabel('x')
ax0.set_ylabel('y')
c = ax1.scatter(pos[:,0],pos[:,1])
ax1.set_title('Trajectory')
ax1.set_xlabel('x')
ax1.set_ylabel('y')
plt.show()
