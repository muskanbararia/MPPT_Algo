#Jupyter Notebook- Python PV-Lib methods to simulate results of Solar Tracking Algorithms

%matplotlib inline
import matplotlib.pyplot as plt
try:
    import seaborn as sns
    sns.set(rc={"figure.figsize": (12, 6)})
except ImportError:
    print('We suggest you install seaborn using conda or pip and rerun this cell')

# built in python modules
import datetime
import logging
import os
import inspect

# python add-ons
import numpy as np
import pandas as pd

#Importiing pvLib - PVLIB Python is a community supported tool that provides a set of functions and classes for simulating the performance of photovoltaic energy systems. 
import pvlib
#Import Location from pvLib to calculate location based parameters
from pvlib.location import Location
#Cosd/sinD returns cos and sine values for given angles
from pvlib.tools import cosd, sind

#Setting variable location as BBSR with Lat = 20.2961 & Longitude as 85.8245, Time zone Asia/Kolkata, Altitude 45
tus = Location(20.2961, 85.8245, 'Asia/Kolkata', 45, 'BBSR')

#Defining the Y-axis for Fig. 1 with Date/Time and plotting frequency as 1 min.
times = pd.date_range(start=datetime.datetime(2015,6,24), end=datetime.datetime(2015,6,25), freq='1Min', tz=tus.tz)

#Generating result for plotting X-axis based on input for "times" variable
#Using "get_solarposition" to obtain apparent elevation, apparent azimuth, elevation, azimuth, apparent and zenith
pyephem_ephem = pvlib.solarposition.get_solarposition(times, tus.latitude, tus.longitude, method='pyephem')

#Generating heads to plot
print(pyephem_ephem.head())
#Printing the plot
pyephem_ephem.plot()

#Storing $pyephem_ephem in $emphemout
ephemout = pyephem_ephem

#Allocating variables to use in Fig 2.
#Solar azimuth angle in decimal degrees
azimuth = ephemout['azimuth']
#Solar apparent azimuth in decimal degrees
apparent_azimuth = ephemout['azimuth']
#Solar apparent zenith angle in decimal degrees
apparent_zenith = ephemout['apparent_zenith']
#The tilt of axis of rotation with respect to horizontal
axis_tilt = 10
#Value donating the compas direction along which the axis of rotation lies
axis_azimuth = 170

latitude = 32
max_angle = 65 
backtrack = True
# A value denotaing the ground coverage ratio of a tracker
gcr = 2.0/7.0

#Acquiring Azimuth values and setting it on $times
times = azimuth.index

#Setting Variables
az = apparent_azimuth - 180
apparent_elevation = 90 - apparent_zenith
x = cosd(apparent_elevation) * sind(az)
y = cosd(apparent_elevation) * cosd(az)
z = sind(apparent_elevation)

#Definging Plot for Fig 2.
earth_coords = pd.DataFrame({'x':x,'y':y,'z':z})

#Printing Fig. 2 plot
earth_coords.plot()
#Giving the plot A header
plt.title('sun position in Earth coordinate system')

#To Transform to panel coordinate syste
axis_azimuth_south = axis_azimuth - 180

print('cos(axis_azimuth_south)={}, sin(axis_azimuth_south)={}'
      .format(cosd(axis_azimuth_south), sind(axis_azimuth_south)))
print('cos(axis_tilt)={}, sin(axis_tilt)={}'
      .format(cosd(axis_tilt), sind(axis_tilt)))

#Changing the output for transformed table
xp = x*cosd(axis_azimuth_south) - y*sind(axis_azimuth_south);
yp = (x*cosd(axis_tilt)*sind(axis_azimuth_south) +
      y*cosd(axis_tilt)*cosd(axis_azimuth_south) -
      z*sind(axis_tilt))
zp = (x*sind(axis_tilt)*sind(axis_azimuth_south) +
      y*sind(axis_tilt)*cosd(axis_azimuth_south) +
      z*cosd(axis_tilt))

#Storing result in multidimensional array/ Dataframe
panel_coords = pd.DataFrame({'x':xp,'y':yp,'z':zp})

#Printing The plot
panel_coords.plot()
plt.title('sun position in panel coordinate system')

# Calculate angle from X-Y plane to calculate sun vector onto X-Z plane nad then obtain $wid
wid = pd.Series(90 - np.degrees(np.arctan2(zp, xp)), index=times)

# filter for sun above panel horizon
wid[zp <= 0] = np.nan

wid.plot(label='tracking angle')
ephemout['apparent_elevation'].plot(label='apparent elevation')
plt.legend()
plt.title('Ideal tracking angle without backtracking')

# angle from x-y plane to projection of sun vector onto x-z plane
tmp = np.degrees(np.arctan(zp/xp))  
# Obtain wid by translating tmp to convention for rotation angles.
# Have to account for which quadrant of the x-z plane in which the sun 
# vector lies.  Complete solution here but probably not necessary to 
# consider QIII and QIV.
wid = pd.Series(index=times)
wid[(xp>=0) & (zp>=0)] =  90 - tmp[(xp>=0) & (zp>=0)];  # QI
wid[(xp<0)  & (zp>=0)] = -90 - tmp[(xp<0)  & (zp>=0)];  # QII
wid[(xp<0)  & (zp<0)]  = -90 - tmp[(xp<0)  & (zp<0)];   # QIII
wid[(xp>=0) & (zp<0)]  =  90 - tmp[(xp>=0) & (zp<0)];   # QIV

# filter for sun above panel horizon
wid[zp <= 0] = np.nan

wid.plot(label='tracking angle')
ephemout['apparent_elevation'].plot(label='apparent elevation')
plt.legend()
plt.title('Ideal tracking angle without backtracking')

if backtrack:
    axes_distance = 1/gcr
    temp = np.minimum(axes_distance*cosd(wid), 1)

    # backtrack angle
    # (always positive b/c acosd returns values between 0 and 180)
    wc = np.degrees(np.arccos(temp))

    v = wid < 0
    widc = pd.Series(index=times)
    widc[~v] = wid[~v] - wc[~v]; # Eq 4 applied when wid in QI
    widc[v] = wid[v] + wc[v];    # Eq 4 applied when wid in QIV
else:
    widc = wid

widc.plot(label='tracking angle')
#pyephemout['apparent_elevation'].plot(label='apparent elevation')
plt.legend(loc=2)
plt.title('Ideal tracking angle with backtracking')

tracking_angles = pd.DataFrame({'with backtracking':widc,'without backtracking':wid})
tracking_angles.plot()
plt.legend()

#Applying restriction to maximum angle 
tracker_theta = widc.copy()
tracker_theta[tracker_theta > max_angle] = max_angle
tracker_theta[tracker_theta < -max_angle] = -max_angle

tracking_angles['with restriction'] = tracker_theta
tracking_angles.plot()

#Calculating vector for panel-normal
panel_norm = np.array([sind(tracker_theta), tracker_theta*0,cosd(tracker_theta)])

panel_norm_df = pd.DataFrame(panel_norm.T, columns=('x','y','z'), index=times)
panel_norm_df.plot()
plt.title('panel normal vector components in panel coordinate system')
plt.legend()

un_vec = np.array([xp, yp, zp])

panel_coords = pd.DataFrame(sun_vec.T, columns=('x','y','z'), index=times)

panel_coords.plot()
plt.title('sun position in panel coordinate system')

#Plotting angle of incidence
aoi = np.degrees(np.arccos(np.abs(np.sum(sun_vec*panel_norm, axis=0))))
aoi = pd.Series(aoi, index=times)

aoi.plot()
plt.title('angle of incidence')

#Plotting cos of angle of incidence
cosd(aoi).plot()

# Calculate standard rotation matrix
print('cos(axis_azimuth_south)={}, sin(axis_azimuth_south)={}'
      .format(cosd(axis_azimuth_south), sind(axis_azimuth_south)))
print('cos(axis_tilt)={}, sin(axis_tilt)={}'
      .format(cosd(axis_tilt), sind(axis_tilt)))

#rotate about x-axis by angle -axis_tilt so that y-axis is also parallel to earth surface
rot_x = np.array([[1, 0, 0], [0, cosd(-axis_tilt), -sind(-axis_tilt)], [0, sind(-axis_tilt), cosd(-axis_tilt)]])

# panel_norm_earth contains the normal vector expressed in earth-surface coordinates
# (z normal to surface, y aligned with tracker axis parallel to earth)
panel_norm_earth = np.dot(rot_x, panel_norm).T

# projection to plane tangent to earth surface,
# in earth surface coordinates
projected_normal = np.array([panel_norm_earth[:,0], panel_norm_earth[:,1], panel_norm_earth[:,2]*0]).T

# calculate magnitudes
panel_norm_earth_mag = np.sqrt(np.nansum(panel_norm_earth**2, axis=1))
projected_normal_mag = np.sqrt(np.nansum(projected_normal**2, axis=1))
#print('panel_norm_earth_mag={}, projected_normal_mag={}'.format(panel_norm_earth_mag, projected_normal_mag))

#First rotate about x-axis by angle -axis_tilt so that y-axis is also parallel to earth surface, then project.
projected_normal = (projected_normal.T / projected_normal_mag).T

panel_norm_earth_df = pd.DataFrame(panel_norm_earth, columns=('x','y','z'), index=times)
panel_norm_earth_df.plot()
plt.title('panel normal vector components in Earth coordinate system')

projected_normal_df = pd.DataFrame(projected_normal, columns=('x','y','z'), index=times)
projected_normal_df.plot()

plt.title('panel normal vector projected to surface in Earth coordinate system')

#importing pvsystem to study Sandia National Lab's Database
from pvlib import pvsystem

#pvlib can import TMY2 and TMY3 data. Here, we import the example files. tmy3 is for Sandpoint and tmy2 is for Miami
pvlib_abspath = os.path.dirname(os.path.abspath(inspect.getfile(pvlib)))
tmy3_data, tmy3_metadata = pvlib.tmy.readtmy3(os.path.join(pvlib_abspath, 'data', '703165TY.csv'))
tmy2_data, tmy2_metadata = pvlib.tmy.readtmy2(os.path.join(pvlib_abspath, 'data', '12839.tm2'))
pvlib.pvsystem.systemdef(tmy3_metadata, 0, 0, .1, 5, 5)
pvlib.pvsystem.systemdef(tmy2_metadata, 0, 0, .1, 5, 5)

#Angle of incidence modifiers. Here we are demonstarting Ashrae modifier using .ashraeiam
angles = np.linspace(-180,180,3601)
ashraeiam = pd.Series(pvsystem.ashraeiam(angles, .05), index=angles)
ashraeiam.plot()
plt.ylabel('ASHRAE modifier')
plt.xlabel('input angle (deg)')

#Angle of incidence modifiers. Here we are demonstarting physical modifier using .physicaliam
angles = np.linspace(-180,180,3601)
physicaliam = pd.Series(pvsystem.physicaliam(angles), index=angles)
physicaliam.plot()
plt.ylabel('physical modifier')
plt.xlabel('input index')

#Here we are comparing the above two modifiers
plt.figure()
ashraeiam.plot(label='ASHRAE')
physicaliam.plot(label='physical')
plt.ylabel('modifier')
plt.xlabel('input angle (deg)')
plt.legend()

#PV system efficiency can vary by up to 0.5% per degree C, so it's important to accurately model cell and module temperature. The sapm_celltemp function uses plane of array irradiance, ambient temperature, wind speed, and module and racking type to calculate cell and module temperatures. The default parameter set is open_rack_cell_glassback.

# scalar inputs
pvsystem.sapm_celltemp(900, 5, 20) # irrad, wind, temp

# vector inputs
times = pd.DatetimeIndex(start='2015-01-01', end='2015-01-02', freq='12H')
temps = pd.Series([0, 10, 5], index=times)
irrads = pd.Series([0, 500, 0], index=times)
winds = pd.Series([10, 5, 0], index=times)
#Plotting temp_cell vs temp_module
pvtemps = pvsystem.sapm_celltemp(irrads, winds, temps)
pvtemps.plot()


