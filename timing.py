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
res_chamb_a = np.inf
res_mattri_a = np.inf
res_noop_a = np.inf
print('Control_points, No-op, Chamberlin, Matrix')
for name, ctrlpts in control_s.items():
    chamstring = {'proj': 'chamb', #clockwise
                  'lon_1': ctrlpts[0][0],
                  'lon_2': ctrlpts[2][0],
                  'lon_3': ctrlpts[1][0],
                  'lat_1': ctrlpts[0][1],
                  'lat_2': ctrlpts[2][1],
                  'lat_3': ctrlpts[1][1]}
    mattristring = {'proj': 'mattri', #counterclockwise
                  'lon_1': ctrlpts[0][0],
                  'lon_2': ctrlpts[1][0],
                  'lon_3': ctrlpts[2][0],
                  'lat_1': ctrlpts[0][1],
                  'lat_2': ctrlpts[1][1],
                  'lat_3': ctrlpts[2][1]}
    tr_chamb = pyproj.transformer.Transformer.from_crs(wgs84, chamstring)
    tr_mattri = pyproj.transformer.Transformer.from_crs(wgs84, mattristring)
    tr_noop = pyproj.transformer.Transformer.from_crs(wgs84, wgs84) 
    res_chamb = min(timeit.repeat(stmt='tr_chamb.transform(x, y)', 
                              globals=globals(), number=100))
    res_mattri = min(timeit.repeat(stmt='tr_mattri.transform(x, y)', 
                              globals=globals(), number=100))
    res_noop = min(timeit.repeat(stmt='tr_noop.transform(x, y)', 
                              globals=globals(), number=100))
    res_chamb_a = min(res_chamb_a, res_chamb)
    res_mattri_a = min(res_mattri_a, res_mattri)
    res_noop_a = min(res_noop_a, res_noop)    
    print(name, ", ", res_noop, ", ", res_mattri, ", ", res_chamb)
print("Total, ", res_noop_a, ", ", res_mattri_a, ", ", res_chamb_a)