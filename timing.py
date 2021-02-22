#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pyproj
import timeit
import numpy as np

pyproj.datadir.set_data_dir('/usr/local/share/proj')
wgs84 = {'proj': 'longlat', 'datum': 'WGS84', 'no_defs': None, 'type': 'crs'}
control_s = {
    'Canada_Atlas': [(-98-13/60, 61+39/60), (-135, 40), (-55, 40)],
    'Canada_Wall_Map': [(-150, 60), (-97.5, 50), (-45, 60)],
    'NW_South_America': [(-69, -25), (-55, 10), (-85, 10)],
    'Australia': [(134, -8), (110, -32), (158, -32)],
    'S_South_America': [(-43, -18), (-72, -18), (-72, -56)],
    'E_South_America': [(-63-33/60, 8+8/60),
                        (-58-33/60, -34-35/60), (-35-13/60, -5-47/60)],
    'Europe_Wall_Map': [(15, 72), (-8, 33), (38, 33)],
    'Strebe_Africa': [(0,22), (22.5, -22), (45, 22)],
    'North_America_Wall_Map': [(-150, 55), (-92.5, 10), (-35, 55)],
    'Africa_Wall_Map': [(-19-3/60, 24+25/60), (20, -35), (59+3/60, 24+25/60)],
    'Hemisphere': [(-180,0), (0, -60), (0, 60)],
    'South_America_Wall_Map': [(-80, 9), (-71, -53), (-35, -6)],
}
x, y = np.array(np.meshgrid(np.linspace(-179.5, 179.5, 360),
                               np.linspace(-89.5, 89.5, 180))).reshape(2,-1)

pnames = ['No-op', 'Matrix', 'Chamberlin', 'TPEQD']
allresults = {name: np.inf for name in pnames}
print('Control_points, No-op, Matrix, Chamberlin, TPEQD')
for name, ctrlpts in control_s.items():
    inits = {'Chamberlin': {'proj': 'chamb', #clockwise
                      'lon_1': ctrlpts[0][0],
                      'lon_2': ctrlpts[2][0],
                      'lon_3': ctrlpts[1][0],
                      'lat_1': ctrlpts[0][1],
                      'lat_2': ctrlpts[2][1],
                      'lat_3': ctrlpts[1][1]},
             'Matrix': {'proj': 'mattri', #counterclockwise
                      'lon_1': ctrlpts[0][0],
                      'lon_2': ctrlpts[1][0],
                      'lon_3': ctrlpts[2][0],
                      'lat_1': ctrlpts[0][1],
                      'lat_2': ctrlpts[1][1],
                      'lat_3': ctrlpts[2][1]},
             'TPEQD': {'proj': 'tpeqd',
                      'lon_1': ctrlpts[0][0],
                      'lon_2': ctrlpts[1][0],
                      'lat_1': ctrlpts[0][1],
                      'lat_2': ctrlpts[1][1]},
             'No-op': wgs84
             }
    transformers = {}
    for pn, p in inits.items():
        transformers[pn] = pyproj.transformer.Transformer.from_crs(wgs84, p)
    
    results = {}
    for trn, tr in transformers.items():
        xport = {'tr': tr, 'x': x, 'y': y}
        results[trn] = min(timeit.repeat(stmt='tr.transform(x, y)', 
                              globals=xport, number=100)) 
        allresults[trn] = min(allresults[trn], results[trn])
    print(name, ", ", end="")
    for name in pnames:
        print(results[name], ", ", end="")
    print()
print('Overall, ', end="")
for name in pnames:
    print(allresults[name], ", ", end="")
print()
