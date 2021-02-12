#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pyproj
import timeit
import numpy as np
import matplotlib.pyplot as plt

pyproj.datadir.set_data_dir('/usr/local/share/proj')
wgs84 = {'proj': 'longlat', 'datum': 'WGS84', 'no_defs': None, 'type': 'crs'}
geod = pyproj.Geod(ellps='WGS84')
control_s = {
    'Canada_Atlas': [(-98-13/60, 61+39/60), (-135, 40), (-55, 40)],
    'Canada_Wall_Map': [(-150, 60), (-97.5, 50), (-45, 60)],#!
     'NW_South_America': [(-69, -25), (-55, 10), (-85, 10)],
     'Australia': [(134, -8), (110, -32), (158, -32)],#!
     'S_South_America': [(-43, -18), (-72, -18), (-72, -56)],
     'E_South_America': [(-63-33/60, 8+8/60),
                         (-58-33/60, -34-35/60), (-35-13/60, -5-47/60)],
     'Europe_Wall_Map': [(15, 72), (-8, 33), (38, 33)],#~
     'Strebe_Africa': [(0,22), (22.5, -22), (45, 22)],
     'North_America_Wall_Map': [(-150, 55), (-92.5, 10), (-35, 55)],#!
     'Africa_Wall_Map': [(-19-3/60, 24+25/60), (20, -35), (59+3/60, 24+25/60)],
     #  'Hemisphere': [(-180,0), (0, -60), (0, 60)],#!
     'South_America_Wall_Map': [(-80, 9), (-71, -53), (-35, -6)],
}

for name, ctrlpts in control_s.items():
    centerx = np.mean([ctrlpts[0][0], ctrlpts[1][0], ctrlpts[2][0]])
    #centery = np.mean([ctrlpts[0][1], ctrlpts[1][1], ctrlpts[2][1]])
    xlow = np.clip(centerx - 89.5, -180, 180)
    xhigh = np.clip(centerx + 89.5, -180, 180)
    x, y = np.array(np.meshgrid(#np.linspace(xlow, xhigh, 180),
                                np.linspace(-179.5, 179.5, 360),
                                np.linspace(-89.5, 89.5, 180)))
    mattristring = {'proj': 'mattri', #counterclockwise
                  'lon_1': ctrlpts[0][0],
                  'lon_2': ctrlpts[1][0],
                  'lon_3': ctrlpts[2][0],
                  'lat_1': ctrlpts[0][1],
                  'lat_2': ctrlpts[1][1],
                  'lat_3': ctrlpts[2][1]}
    tr_mattri = pyproj.transformer.Transformer.from_crs(wgs84, mattristring)
    fwd = tr_mattri.transform(x, y)
    inv = tr_mattri.transform(*fwd,
                              direction=pyproj.enums.TransformDirection.INVERSE)
    lt = pyproj.Proj(mattristring)    
    dist = geod.inv(x,y,inv[0],inv[1])[2]
    print(name, ", ", np.nanmax(dist))
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.set_title(name)
    cs = ax.contourf(x, y, dist, levels=[0.0001,0.001,0.01,0.1,1])
    ax.scatter(*zip(*ctrlpts))
    ax.axis('equal')
    fig.colorbar(cs)

#%%
for name, ctrlpts in control_s.items():
    mattristring = {'proj': 'mattri', #counterclockwise
                  'lon_1': ctrlpts[0][0],
                  'lon_2': ctrlpts[1][0],
                  'lon_3': ctrlpts[2][0],
                  'lat_1': ctrlpts[0][1],
                  'lat_2': ctrlpts[1][1],
                  'lat_3': ctrlpts[2][1]}    
    tr_mattri = pyproj.transformer.Transformer.from_crs(wgs84, mattristring)
    for xy in ctrlpts:
        xy2 = tr_mattri.transform(*tr_mattri.transform(*xy),
                            direction=pyproj.enums.TransformDirection.INVERSE)
        print(name, ": ", np.round(xy[0]-xy2[0]), np.round(xy[1]-xy2[1]))
