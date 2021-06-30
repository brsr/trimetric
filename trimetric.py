#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import geopandas
import pandas as pd
#import shapely
from shapely.geometry import Point, LineString#, MultiPolygon, Polygon
import pyproj
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
#import mapproj
import functools
#import os
# os.chdir('Code/mapproj')
np.seterr(divide='ignore')
pyproj.datadir.set_data_dir('/usr/local/share/proj')
world = geopandas.read_file(geopandas.datasets.get_path('naturalearth_lowres'))

def geodesics(lon, lat, geod, n=100, includepts=False):
    """Draw geodesics between each adjacent pair of points given by
    lon and lat.
    """
    lon2 = np.roll(lon, -1, axis=0)
    lat2 = np.roll(lat, -1, axis=0)
    result = []
    for l, t, l2, t2 in zip(lon, lat, lon2, lat2):
        g = geod.npts(l, t, l2, t2, n)
        g.insert(0, (l, t))
        g.append((l2, t2))
        result.append(LineString(g))
    ctrlboundary = geopandas.GeoSeries(result)
    if includepts:
        controlpts = arraytoptseries(np.array([lon, lat]))
        ctrlpoly = geopandas.GeoSeries(pd.concat([ctrlboundary, controlpts],
                                            ignore_index=True))
        return ctrlpoly
    else:
        return ctrlboundary

def graticule(spacing1=15, spacing2=1,
              lonrange = [-180, 180], latrange = [-90, 90],
              crs=world.crs):
    """
    Create a graticule (or another square grid)
    """
    a = int((lonrange[1] - lonrange[0])//spacing2)
    b = int((latrange[1] - latrange[0])//spacing1)
    c = int((lonrange[1] - lonrange[0])//spacing1)
    d = int((latrange[1] - latrange[0])//spacing2)
    plx = np.linspace(lonrange[0], lonrange[1], num=a + 1)
    ply = np.linspace(latrange[0], latrange[1], num=b + 1)
    mex = np.linspace(lonrange[0], lonrange[1], num=c + 1)
    mey = np.linspace(latrange[0], latrange[1], num=d + 1)
    parallels = np.stack(np.meshgrid(plx, ply), axis=-1).transpose((1,0,2))
    meridians = np.stack(np.meshgrid(mex, mey), axis=-1)
    gratlist = [parallels[:, i] for i in range(parallels.shape[1])]
    gratlist += [meridians[:, i] for i in range(meridians.shape[1])]
    gratl2 = [LineString(line) for line in gratlist]
    grat = geopandas.GeoSeries(gratl2)
    grat.crs = crs
    return grat

def ptseriestoarray(ser):
    """Convert a geopandas GeoSeries containing shapely Points
    (or LineStrings of all the same length) to an array of
    shape (2, n) or (3, n).
    """
    return np.stack([x.coords for x in ser], axis=-1).squeeze()

def arraytoptseries(arr, crs=world.crs):
    """Convert an array of shape (2, ...) or (3, ...) to a
    geopandas GeoSeries containing shapely Point objects.
    """
    if arr.shape[0] == 2:
        result = geopandas.GeoSeries([Point(x[0], x[1])
                                    for x in arr.reshape(2, -1).T])
    else:
        result = geopandas.GeoSeries([Point(x[0], x[1], x[2])
                                    for x in arr.reshape(3, -1).T])
    result.crs = crs
    return result

def latlon_vec(ll, scale=np.pi/180):
    """Convert longitude and latitude to 3-vector

    >>> ll = np.arange(6).reshape(2,3)*18
    >>> UnitVector.transform_v(ll)
    array([[5.87785252e-01, 2.93892626e-01, 4.95380036e-17],
           [0.00000000e+00, 9.54915028e-02, 3.59914664e-17],
           [8.09016994e-01, 9.51056516e-01, 1.00000000e+00]])
    """
    lon, lat = ll[0]*scale, ll[1]*scale
    x = np.cos(lat)*np.cos(lon)
    y = np.cos(lat)*np.sin(lon)
    z = np.sin(lat)
    return np.stack([x, y, z], axis=0)

def circgen(lon, lat, r, cross=False):
    az = 0
    circs = []
    for x, y, ri in zip(lon, lat, r):
        circ = []
        for az in range(361):
            out = geod.fwd(x, y, az, ri)
            circ.append(out[:2])
        circs.append(LineString(circ))
        if cross:
            for ax in [0, 180], [90, 270]:
                p1 = geod.fwd(x, y, ax[0], ri)[:2]
                p2 = geod.fwd(x, y, ax[1], ri)[:2]
                crox = geod.npts(p1[0], p1[1], p2[0], p2[1], 181)
                circs.append(LineString(crox))
    return geopandas.GeoSeries(circs)

r = 6371
geod = pyproj.Geod(a=r)
grat = graticule()

cyclic = [0, 1, 2, 0]
cyclic2 = [1,2,0,1]
pnames = ['Chamberlin', 'Matrix']
adegpts = np.array(np.meshgrid(np.linspace(-179.5, 179.5, 360),
                               np.linspace(-89.5, 89.5, 180)))
degpts = arraytoptseries(adegpts)
degpts.crs = world.crs
vdegpts = latlon_vec(adegpts)
wght = np.cos(adegpts[1]*np.pi/180)

north = geopandas.GeoSeries(Point([0,90]))
north.crs = world.crs

control_points = {
    'Canada_Atlas': geopandas.GeoSeries([Point(-98-13/60, 61+39/60), Point(-135, 40), Point(-55, 40)]),
    'Canada_Wall_Map': geopandas.GeoSeries([Point(-150, 60), Point(-97.5, 50), Point(-45, 60)]),
    'NW_South_America': geopandas.GeoSeries([Point(-69, -25), Point(-55, 10), Point(-85, 10)]),
    'Australia': geopandas.GeoSeries([Point(134, -8), Point(110, -32), Point(158, -32)]),
    'S_South_America': geopandas.GeoSeries([Point(-43, -18), Point(-72, -18), Point(-72, -56)]),
    'E_South_America': geopandas.GeoSeries([Point(-63-33/60, 8+8/60),
                                  Point(-58-33/60, -34-35/60), Point(-35-13/60, -5-47/60)]),
    'Europe_Wall_Map': geopandas.GeoSeries([Point(15, 72), Point(-8, 33), Point(38, 33)]),
    #'Strebe_Africa': geopandas.GeoSeries([Point(0,22), Point(22.5, -22), Point(45, 22)]),
    'North_America_Wall_Map': geopandas.GeoSeries([Point(-150, 55), Point(-92.5, 10), Point(-35, 55)]),
    'Africa_Wall_Map': geopandas.GeoSeries([Point(-19-3/60, 24+25/60), Point(20, -35), Point(59+3/60, 24+25/60)]),
    #'Hemisphere': geopandas.GeoSeries([Point(-180,0), Point(0, -60), Point(0, 60)]),
    'South_America_Wall_Map': geopandas.GeoSeries([Point(-80, 9), Point(-71, -53), Point(-35, -6)]),    
}#south america last so we can use it in the construction figure
#control_points = {'South_America_Wall_Map': 
#                  control_points['South_America_Wall_Map']}
focus = {
    'Canada_Atlas': world.index[world.name == 'Canada'],
    'Canada_Wall_Map': world.index[world.name == 'Canada'],
    'NW_South_America': world.index[world.name.isin(['Peru', 'Ecuador',
                                                     'Colombia', 'Suriname', 'Venezuela', 'Guyana'])],
    'Australia': world.index[world.name == 'Australia'],
    'S_South_America': world.index[world.name.isin(['Bolivia', 'Chile',
                                                    'Argentina', 'Uruguay', 'Paraguay'])],
    'E_South_America': world.index[world.name == 'Brazil'],
    # don't focus on France because of French Guiana
    # and Russia because half of it is in Asia
    'Europe_Wall_Map': world.index[(world.continent == 'Europe') &
                                   (world.name != 'France') &
                                   (world.name != 'Russia')],
    'Strebe_Africa': world.index[world.continent == 'Africa'],
    'South_America_Wall_Map': world.index[world.continent == 'South America'],
    'North_America_Wall_Map': world.index[world.continent == 'North America'],
    'Africa_Wall_Map': world.index[world.continent == 'Africa'],
}

exclude = {
    'Canada_Atlas': world.index[world.continent.isin(['Antarctica', 'Oceania'])],
    'Canada_Wall_Map': world.index[world.continent.isin(['Antarctica', 'Oceania'])],
    'NW_South_America': world.index[world.continent.isin(['Asia', 'Oceania'])],
    'Australia': world.index[world.continent.isin(['North America',
                                                   'South America',
                                                   'Africa', 'Europe'])],
    'S_South_America': world.index[world.continent.isin(['Asia'])],
    'E_South_America': world.index[world.continent.isin(['Asia', 'Oceania'])],
    'Europe_Wall_Map': world.index[world.continent.isin(['Antarctica', 'Oceania'])],
    'Strebe_Africa': world.index[world.continent.isin(['North America'])],
    'South_America_Wall_Map': world.index[world.continent.isin(['Asia', 'Oceania'])],
    'North_America_Wall_Map': world.index[world.continent.isin(['Antarctica', 'Oceania'])],
    'Africa_Wall_Map': world.index[world.continent.isin(['North America'])],
    'Hemisphere': world.index[~world.continent.isin(['North America',
                                                     'South America'])],
}

extents = {'South_America_Wall_Map': ([-3500, 4500], [-4000, 4000])}
scalelevels = {'South_America_Wall_Map': np.linspace(0.96, 1.14, 10)}
omegalevels = {'South_America_Wall_Map': np.linspace(0, 10, 6)}
dlevels = {'South_America_Wall_Map': np.linspace(0, 150, 7)}

for controlpts in control_points.values():
    controlpts.crs = world.crs

pdindex = pd.MultiIndex.from_product(
    [control_points.keys(), pnames], names=['ctrlpts', 'proj'])
# %%
figsize=(7.2, 3.5)
cptable = pd.DataFrame(dtype=float,
                       columns=['pt1_lon', 'pt1_lat',
                                'pt2_lon', 'pt2_lat',
                                'pt3_lon', 'pt3_lat',
                                'len23', 'len31', 'len12',
                                'asymmetry', 'area'])

pdtable = pd.DataFrame(index=pdindex, dtype=float,
                       columns=['avgomega', 'maxomega',
                                'avgld', 'maxld', 'minscale', 'maxscale', 'scalerat'])

for name, controlpts in control_points.items():
    print(name)
    actrlpts = ptseriestoarray(controlpts)
    # stats about the control triangle
    linelengths = geod.line_lengths(actrlpts[0, cyclic2],
                                    actrlpts[1, cyclic2])
    asymmetry = 3*max(linelengths)/sum(linelengths) - 1
    cptable.loc[name] = np.concatenate([actrlpts.T.flatten(),
                                        linelengths,
                                        [asymmetry],
                                        [geod.polygon_area_perimeter(actrlpts[0], actrlpts[1])[0]]])

    # determine index for interior of triangle
    ar = np.linalg.inv(latlon_vec(actrlpts))
    adegpts_ar = np.tensordot(ar, vdegpts, axes=(1,0))
    zone_ar = np.signbit(adegpts_ar).sum(axis=0)
    index_ar = zone_ar == 0
    wght_m = wght.copy()
    wght_m[~index_ar] = np.nan
    # initialize projections
    chamstring = {'proj': 'chamb',
                  'lon_1': actrlpts[0, 0],
                  'lon_2': actrlpts[0, 2],
                  'lon_3': actrlpts[0, 1],
                  'lat_1': actrlpts[1, 0],
                  'lat_2': actrlpts[1, 2],
                  'lat_3': actrlpts[1, 1],
                  'R': geod.a}
    mattristring = {'proj': 'mattri',
                  'lon_1': actrlpts[0, 0],
                  'lon_2': actrlpts[0, 1],
                  'lon_3': actrlpts[0, 2],
                  'lat_1': actrlpts[1, 0],
                  'lat_2': actrlpts[1, 1],
                  'lat_3': actrlpts[1, 2],
                  'R': geod.a}
    ct = pyproj.Proj(chamstring)
    lt = pyproj.Proj(mattristring)
    # various quantities to map
    gd = geodesics(actrlpts[0], actrlpts[1], geod, n=100)
    gd.crs = world.crs
    antipodes = np.array([actrlpts[0] - 180, -actrlpts[1]])
    antipodes[0] = antipodes[0] + np.where(antipodes[0] < -180, 360, 0)
    center = actrlpts.mean(axis=1)
    gratrange = np.array([-90, 90]) + center[0]
    # to avoid the antipodal part
    grat2 = graticule(lonrange=gratrange)
    cpts = np.array(np.meshgrid(np.linspace(-60, 60, 9) + center[0],
                                np.linspace(-60, 60, 9) +
                                np.clip(center[1], -25, 25))).reshape(2, -1)
    tissotr = np.ones(cpts.shape[1])*500
    tissot = circgen(cpts[0], cpts[1], tissotr, cross=True)
    tissot.crs = world.crs
    # perform the projections
    lt_factors = lt.get_factors(adegpts[0], adegpts[1])
    degpts_lt = degpts.to_crs(mattristring)
    adegpts_lt = ptseriestoarray(degpts_lt).reshape(adegpts.shape)
    controlpts_lt = controlpts.to_crs(mattristring)
    world_lt = world.to_crs(mattristring)
    grat_lt = grat.to_crs(mattristring)
    grat2_lt = grat2.to_crs(mattristring)
    gd_lt = gd.to_crs(mattristring)
    tissot_lt = tissot.to_crs(mattristring)
    north_lt = north.to_crs(mattristring)

    ct_factors = ct.get_factors(adegpts[0], adegpts[1])
    adegpts_ct = ptseriestoarray(degpts.to_crs(chamstring)).reshape(adegpts.shape)
    controlpts_ct = controlpts.to_crs(chamstring)
    world_ct = world.to_crs(chamstring)
    grat_ct = grat.to_crs(chamstring)
    grat2_ct = grat2.to_crs(chamstring)
    gd_ct = gd.to_crs(chamstring)
    tissot_ct = tissot.to_crs(chamstring)

    #rotate lt so north is up
    angle = 90 - np.arctan2(north_lt[0].y, north_lt[0].x)*180/np.pi
    degpts_lt = degpts_lt.rotate(angle, origin=(0, 0))
    adegpts_lt = ptseriestoarray(degpts_lt).reshape(adegpts.shape)
    controlpts_lt = controlpts_lt.rotate(angle, origin=(0, 0))
    world_lt = world_lt.rotate(angle, origin=(0, 0))
    grat_lt = grat_lt.rotate(angle, origin=(0, 0))
    grat2_lt = grat2_lt.rotate(angle, origin=(0, 0))
    gd_lt = gd_lt.rotate(angle, origin=(0, 0))
    tissot_lt = tissot_lt.rotate(angle, origin=(0, 0))

    # affine transform ct to line up with lt
    uno = np.ones((1, 3))
    ltz = np.concatenate([ptseriestoarray(controlpts_lt), uno])
    ctz = np.concatenate([ptseriestoarray(controlpts_ct), uno])
    # want matrix such that ltz = M @ ctz
    #M = ltz @ inv(ctz)
    # this isn't necessarily numerically stable but oh well
    matrix = ltz @ np.linalg.inv(ctz)
    affine = [matrix[0, 0], matrix[0, 1], matrix[1, 0], matrix[1, 1],
              matrix[0, 2], matrix[1, 2]]
    adegpts_ct = (np.tensordot(matrix[:2, :2], adegpts_ct, axes=(1, 0))
                  + matrix[:2, 2][:, np.newaxis, np.newaxis])
    controlpts_ct = controlpts_ct.affine_transform(affine)
    world_ct = world_ct.affine_transform(affine)
    grat_ct = grat_ct.affine_transform(affine)
    grat2_ct = grat2_ct.affine_transform(affine)
    gd_ct = gd_ct.affine_transform(affine)
    tissot_ct = tissot_ct.affine_transform(affine)

    tgtpts = ptseriestoarray(controlpts_lt)

    # distortion calculations
    omega_lt = lt_factors.angular_distortion
    scale_lt = lt_factors.areal_scale
    omega_ct = ct_factors.angular_distortion
    scale_ct = ct_factors.areal_scale

    omegas = [omega_ct, omega_lt]
    scales = [scale_ct, scale_lt]

    unprojlengths = []
    for i in range(3):
        f, b, l = geod.inv(adegpts[0], adegpts[1],
                           actrlpts[0, i]*np.ones(adegpts.shape[1:]),
                           actrlpts[1, i]*np.ones(adegpts.shape[1:]))
        unprojlengths.append(l)
    unprojlengths = np.array(unprojlengths)

    lengthdistorts = []
    for dp in [adegpts_ct, adegpts_lt]:
        x = dp[:, np.newaxis] - tgtpts[..., np.newaxis, np.newaxis]
        y = np.sqrt(np.sum(x**2, axis=0))
        result = np.nansum(abs(y - unprojlengths), axis=0)
        lengthdistorts.append(result)

    for omega, scale, ld, pname in zip(omegas, scales, lengthdistorts, pnames):
        avgomega = np.nansum(omega*wght_m)/np.nansum(wght_m)
        maxomega = np.nanmax(omega[index_ar])
        avgld = np.nansum(ld*wght_m)/np.nansum(wght_m)
        maxld = np.nanmax(ld[index_ar])
        avgscale = np.nansum(scale*wght_m)/np.nansum(wght_m)
        minscale = np.nanmin(scale[index_ar])
        maxscale = np.nanmax(scale[index_ar])
        scalerat = maxscale/minscale - 1
        pdtable.loc[name, pname] = [avgomega, maxomega, avgld, maxld,
                                    minscale, maxscale, scalerat*100]
    # %
    excl = exclude[name]
    world_ctv = world_ct.drop(excl)
    world_ltv = world_lt.drop(excl)
    # %#figure: comparison of projections (omit, use the next one)
    fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True)
    ax = axes[0]
    world_ctv.plot(ax=ax, color='k')
    grat_ct.plot(ax=ax, color='lightgrey', linewidth=1)
    controlpts_ct.plot(ax=ax, color='r')
    #gd_ct.plot(ax=ax, color='green', linestyle=':')
    ax.set_title('Chamberlin')
    ax = axes[1]
    world_ltv.plot(ax=ax, color='k')
    grat_lt.plot(ax=ax, color='lightgrey', linewidth=1)
    controlpts_lt.plot(ax=ax, color='r')
    #gd_lt.plot(ax=ax, color='green', linestyle=':')
    ax.set_title('Matrix')
    fig.savefig(name + '_whole.png', bbox_inches = 'tight')
    # bounds for plots
    try:
        xbounds, ybounds = extents[name]
    except KeyError:
        bd = pd.concat([gd_ct.bounds, gd_lt.bounds])
        try:
            bd = pd.concat([bd,
                            world_ct.loc[focus[name]].bounds,
                            world_lt.loc[focus[name]].bounds])
        except KeyError:
            pass
        bd[~np.isfinite(bd)] = np.nan
        xbounds = [bd.minx.min(), bd.maxx.max()]
        xbounds[0] = xbounds[0]*1.1 if xbounds[0] < 0 else xbounds[0]*0.9
        xbounds[1] = xbounds[1]*1.1 if xbounds[1] > 0 else xbounds[1]*0.9
        ybounds = [bd.miny.min(), bd.maxy.max()]
        ybounds[0] = ybounds[0]*1.1 if ybounds[0] < 0 else ybounds[0]*0.9
        ybounds[1] = ybounds[1]*1.1 if ybounds[1] > 0 else ybounds[1]*0.9
        xsize = xbounds[1] - xbounds[0]
        ysize = ybounds[1] - ybounds[0]
        if xsize > ysize:
            ymean = tgtpts[1].mean()
            ybounds = ymean - xsize/2, ymean + xsize/2,
    # %#figure: zoom comparison of projections
    fig, axes = plt.subplots(1, 2, figsize=figsize, sharex=True, sharey=True)
    ax = axes[0]
    world_ctv.plot(ax=ax, color='k')
    grat2_ct.plot(ax=ax, color='lightgrey', linewidth=1)#, alpha=0.5)
    #gd_ct.plot(ax=ax, color='green', linestyle=':')
    controlpts_ct.plot(ax=ax, color='r')#, marker='x')
    ax.set_title('Chamberlin')
    ax = axes[1]
    world_ltv.plot(ax=ax, color='k')
    grat2_lt.plot(ax=ax, color='lightgrey', linewidth=1)#, alpha=0.5)
    #gd_lt.plot(ax=ax, color='green', linestyle=':')
    controlpts_lt.plot(ax=ax, color='r')#, marker='x')
    ax.set_xlim(*xbounds)
    ax.set_ylim(*ybounds)
    ax.set_title('Matrix')
    fig.savefig(name + '_zoom.eps', bbox_inches = 'tight')
    # %#figure: tissot
    fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True)
    axes[0].set_title('Chamberlin')
    axes[1].set_title('Matrix')
    for gdx, worldx, grat2x, tissots, ax in zip(
        [gd_ct, gd_lt],
        [world_ctv, world_ltv],
        [grat2_ct, grat2_lt],
            [tissot_ct, tissot_lt], axes):
        worldx.plot(ax=ax, color='#B4B4B4')
        grat2x.plot(ax=ax, color='lightgrey', linewidth=1)#, alpha=0.5)
        tissots.plot(ax=ax, color='k')
        ax.scatter(tgtpts[0], tgtpts[1], color='r')
        ax.axis('equal')
        #gdx.plot(ax=ax, color='g', linestyle=':')
        #ax.tick_params(axis='y', labelrotation=90)
    ax.set_xlim(*xbounds)
    ax.set_ylim(*ybounds)
    plt.draw()
    axes[0].set_yticklabels(axes[0].get_yticklabels(), va='center')
    fig.savefig(name + '_tissot.eps', bbox_inches = 'tight')
    # %#figure: scale
    try:
        levels = scalelevels[name]
    except KeyError:
        px = pdtable.loc[name]
        pmax = px.maxscale.max()
        pmin = np.clip(px.minscale.min(), 0, pmax)
        levels = np.linspace(pmin, pmax, 6)

    fig, axes = plt.subplots(1, 3, gridspec_kw={"width_ratios": [1, 1, 0.1]},
                             figsize=figsize)
    axes[0].set_title('Chamberlin')
    axes[1].set_title('Matrix')
    axes[0].get_shared_y_axes().join(axes[0], axes[1])
    axes[1].set_yticklabels([])
    rindex = np.isfinite(scales[0]) & np.isfinite(scales[1])
    for scale, dp, gdx, worldx, grat2x, ax in zip(scales,
                                                  [adegpts_ct, adegpts_lt],
                                                  [gd_ct, gd_lt],
                                                  [world_ctv, world_ltv],
                                                  [grat2_ct, grat2_lt], axes):
        #worldx.plot(ax=ax, color='#B4B4B4')
        #grat2x.plot(ax=ax, color='lightgrey', linewidth=1)#, alpha=0.5)
        #cs = ax.contour(dp[0], dp[1],
        #                np.where(rindex, scale, np.nan),
        #                colors='k', levels=levels)
        cs = ax.contourf(dp[0], dp[1], scale,
                        #np.where(rindex, scale, np.nan),
                        cmap='gray',                        
                        levels=levels)        
        #ax.clabel(cs, fmt='%1.2f', inline_spacing=-8)
        thing = scale.copy()
        thing[~index_ar] = np.nan
        ah = np.unravel_index(np.nanargmax(thing), thing.shape)
        al = np.unravel_index(np.nanargmin(thing), thing.shape)
        ax.scatter(dp[0][al], dp[1][al], color='w', marker='+')
        ax.scatter(tgtpts[0], tgtpts[1], color='r')
        ax.axis('equal')
        #gdx.plot(ax=ax, color='k', linestyle=':')
        #ax.tick_params(axis='y', labelrotation=90)
    ax.set_xlim(*xbounds)
    ax.set_ylim(*ybounds)
    plt.draw()
    axes[0].set_yticklabels(axes[0].get_yticklabels(), va='center')
    fig.colorbar(cs, axes[2])
    fig.savefig(name + '_scale.eps', bbox_inches = 'tight')
    # %#figure: omega
    try:
        levels = omegalevels[name]
    except KeyError:
        levels = np.linspace(0, np.floor(px.maxomega.max()), 6)
    fig, axes = plt.subplots(1, 3, gridspec_kw={"width_ratios": [1, 1, 0.1]},
                             figsize=figsize)
    axes[0].set_title('Chamberlin')
    axes[1].set_title('Matrix')
    axes[0].get_shared_y_axes().join(axes[0], axes[1])
    axes[1].set_yticklabels([])
    param = omegas[0] - omegas[1]
    for omega, dp, gdx, worldx, grat2x, ax in zip(omegas,
                                                  [adegpts_ct, adegpts_lt],
                                                  [gd_ct, gd_lt],
                                                  [world_ctv, world_ltv],
                                                  [grat2_ct, grat2_lt], axes):
        #worldx.plot(ax=ax, color='#B4B4B4')
        #grat2x.plot(ax=ax, color='lightgrey', linewidth=1)#, alpha=0.5)
        #cs = ax.contour(dp[0], dp[1],
        #                np.where(rindex, omega, np.nan), colors='k',
        #                levels=levels)
        cs = ax.contourf(dp[0], dp[1], omega,
                        #np.where(rindex, omega, np.nan),
                        cmap='gray',                        
                        levels=levels)
        
        #ax.clabel(cs, fmt='%1.0f', inline_spacing=0)
        thing = omega.copy()
        thing[~index_ar] = np.nan
        ah = np.unravel_index(np.nanargmax(thing), thing.shape)
        al = np.unravel_index(np.nanargmin(thing), thing.shape)
        ax.scatter(dp[0][al], dp[1][al], color='w', marker='+')
        ax.scatter(tgtpts[0], tgtpts[1], color='r')
        ax.axis('equal')
        #gdx.plot(ax=ax, color='k', linestyle=':')
        #ax.tick_params(axis='y', labelrotation=90)
        ax.contour(dp[0], dp[1],
                   np.where(rindex, param, np.nan),
                   levels=[0], colors='b', linestyles='--')
    ax.set_xlim(*xbounds)
    ax.set_ylim(*ybounds)
    plt.draw()
    axes[0].set_yticklabels(axes[0].get_yticklabels(), va='center')
    fig.colorbar(cs, axes[2])
    fig.savefig(name + '_omega.eps', bbox_inches = 'tight')
    # %#figure: deviation in distance
    try:
        levels = dlevels[name]
    except KeyError:
        levels = np.linspace(0, np.ceil(px.maxld.max()), 6)
    fig, axes = plt.subplots(1, 3, gridspec_kw={"width_ratios": [1, 1, 0.1]},
                             figsize=figsize)
    axes[0].set_title('Chamberlin')
    axes[1].set_title('Matrix')
    axes[0].get_shared_y_axes().join(axes[0], axes[1])
    axes[1].set_yticklabels([])
    param = lengthdistorts[0] - lengthdistorts[1]
    #levels = [2, 5, 10, 20, 50, 100, 200, 500]
    for ld, dp, gdx, worldx, grat2x, ax in zip(lengthdistorts,
                                               [adegpts_ct, adegpts_lt],
                                               [gd_ct, gd_lt],
                                               [world_ctv, world_ltv],
                                               [grat2_ct, grat2_lt], axes):
        #worldx.plot(ax=ax, color='#B4B4B4')
        #grat2x.plot(ax=ax, color='lightgrey', linewidth=1)#, alpha=0.5)
        #cs = ax.contour(dp[0], dp[1],
        #                np.where(rindex, ld, np.nan),
        #                colors='k', levels=levels)
        cs = ax.contourf(dp[0], dp[1], ld,
                        #np.where(rindex, ld, np.nan),
                        cmap='gray', #norm=matplotlib.colors.LogNorm(),
                        levels=levels)
        
        #ax.clabel(cs, fmt='%1.0f', inline_spacing=-2)
        thing = ld.copy()
        thing[~index_ar] = np.nan
        ah = np.unravel_index(np.nanargmax(thing), thing.shape)
        #al = np.unravel_index(np.nanargmin(thing), thing.shape)
        ax.scatter(dp[0][ah], dp[1][ah], color='k', marker='+')
        ax.scatter(tgtpts[0], tgtpts[1], color='r', zorder=8)
        #gdx.plot(ax=ax, color='k', linestyle=':')
        ax.axis('equal')
        #ax.tick_params(axis='y', labelrotation=90)
        ax.contour(dp[0], dp[1],
                   np.where(rindex, param, np.nan),
                   levels=[0], colors='b', linestyles='--')
    ax.set_xlim(*xbounds)
    ax.set_ylim(*ybounds)
    plt.draw()
    axes[0].set_yticklabels(axes[0].get_yticklabels(), va='center')
    fig.colorbar(cs, axes[2])
    fig.savefig(name + '_distance.eps', bbox_inches = 'tight')
# %%
#pdtable.loc['Hemisphere', 'scalerat'] = np.inf
#pdtable.loc['Hemisphere', 'maxomega'] = 180
pdtable.to_csv('cham_matr_stats.csv', index_label='name')
cptable.to_csv('control_triangles.csv', index_label='name')
pdtablenoh = pdtable#.drop('Hemisphere')
area = cptable.area#.drop('Hemisphere')
#sl = cptable[['len12','len23','len31']].drop('Hemisphere')
#aspect = sl.max(axis=1) - sl.min(axis=1)
aorder = cptable.sort_values('area')#'symmetry')
cham_maxo = aorder.join(pdtablenoh.xs('Chamberlin', level=1))['maxomega']
matrix_maxo = aorder.join(pdtablenoh.xs('Matrix', level=1))['maxomega']
cham_avgo = aorder.join(pdtablenoh.xs('Chamberlin', level=1))['avgomega']
matrix_avgo = aorder.join(pdtablenoh.xs('Matrix', level=1))['avgomega']
cham_maxd = aorder.join(pdtablenoh.xs('Chamberlin', level=1))['maxld']
matrix_maxd = aorder.join(pdtablenoh.xs('Matrix', level=1))['maxld']
cham_avgd = aorder.join(pdtablenoh.xs('Chamberlin', level=1))['avgld']
matrix_avgd = aorder.join(pdtablenoh.xs('Matrix', level=1))['avgld']
cham_sr = aorder.join(pdtablenoh.xs('Chamberlin', level=1))['scalerat']
matrix_sr = aorder.join(pdtablenoh.xs('Matrix', level=1))['scalerat']
labels = aorder.index
ticks = np.arange(len(labels))
#%%
cptablefmt = cptable.copy()
lens = ['len23', 'len31', 'len12']
cptablefmt[lens] = cptable[lens].round()
cptablefmt.area = (cptable.area/1E6).round(2)

def decdeg2dm(dd, suffix=['W','E']):
    if dd == 0:
        return '0\\degree'
    elif dd >= 0:
        suf = suffix[1]
    else:
        suf = suffix[0]
    dd = abs(dd)
    degrees,minutes = divmod(dd*60,60)
    dstring = "{}\\degree".format(degrees)
    mstring = "{}\'".format(minutes) if minutes != 0 else ''
    return dstring + mstring + suf
lons = ['pt1_lon', 'pt2_lon', 'pt3_lon']
lats = ['pt1_lat', 'pt2_lat', 'pt3_lat']
cptablefmt[lons] = np.vectorize(decdeg2dm)(cptable[lons].values)
dd2 = functools.partial(decdeg2dm, suffix=['S','N'])
cptablefmt[lats] = np.vectorize(dd2)(cptable[lats].values)
cptablefmt.to_csv('control_triangles.csv', index_label='name')
# %% scatter plots of the previous: omega
fig, ax = plt.subplots(figsize=(7.8, 3))
ax.grid(b=True, which='major', color='#666666', linestyle='-', axis='x')
ax.scatter(cham_maxo, ticks, edgecolor='k', zorder=6,
           color='lightblue', label='Chamberlin max')
ax.scatter(cham_avgo, ticks, edgecolor='k', zorder=6,
           color='dodgerblue', label='Chamberlin average')
ax.scatter(matrix_maxo, ticks, edgecolor='k',  marker='s', zorder=5,
           color='orange', label='Matrix max')
ax.scatter(matrix_avgo, ticks, edgecolor='k',  marker='s', zorder=5,
           color='orangered', label='Matrix average')
ax.set_xlabel('$\omega$, degrees')
ax.set_yticks(ticks)
ax.set_yticklabels(labels)
ax.legend(loc='best')
fig.subplots_adjust(left=0.3)
fig.savefig('omegaplot.eps', bbox_inches = 'tight')
# %%D
fig, ax = plt.subplots(figsize=(7.8, 3))
ax.grid(b=True, which='major', color='#666666', linestyle='-', axis='x')
#ax.scatter(cham_maxd, ticks, edgecolor='k', zorder=6,
#           color='lightblue', label='Chamberlin max')
ax.scatter(np.fmax(cham_maxd, matrix_maxd), ticks, edgecolor='k',
           marker='s', zorder=5, color='lightgrey', label='Both max')
ax.scatter(cham_avgd, ticks, edgecolor='k', zorder=6,
           color='dodgerblue', label='Chamberlin average')
ax.scatter(matrix_avgd, ticks, edgecolor='k',  marker='s', zorder=5,
           color='orangered', label='Matrix average')
ax.set_xlabel('$D$, km')
ax.set_yticks(ticks)
ax.set_yticklabels(labels)
ax.legend(loc='best')
fig.subplots_adjust(left=0.3)
fig.savefig('distanceplot.eps', bbox_inches = 'tight')
# %% scale
fig, ax = plt.subplots(figsize=(7.8, 3))
ax.grid(b=True, which='major', color='#666666', linestyle='-', axis='x')
ax.scatter(cham_sr, ticks, edgecolor='k', color='dodgerblue',
           label='Chamberlin', zorder=6)
ax.scatter(matrix_sr, ticks, edgecolor='k', color='orange', marker='s',
           label='Matrix', zorder=5)
ax.set_xlabel('scale ratio, %')
ax.set_yticks(ticks)
ax.set_yticklabels(labels)
ax.legend(loc='best')
ax.set_xlim(left=0)
fig.subplots_adjust(left=0.3)
fig.savefig('scaleplot.eps', bbox_inches = 'tight')
# %%figure: show construction
testpt = [-70, -5]
f, b, r = geod.inv(testpt[0]*np.ones(actrlpts.shape[1]),
                   testpt[1]*np.ones(actrlpts.shape[1]),
                   actrlpts[0], actrlpts[1])

az = 0
circs = circgen(actrlpts[0], actrlpts[1], r)
stestpt = geopandas.GeoSeries(Point(testpt))
stestpt.crs = world.crs
sctpt = stestpt.to_crs(chamstring).affine_transform(affine)
ctpt = np.array(sctpt[0].xy)
#ctpt = ct.transform(testpt[0], testpt[1])
sltpt = stestpt.to_crs(mattristring).rotate(angle, origin=(0, 0))
ltpt = np.array(sltpt[0].xy)

fig, axes = plt.subplots(1, 3, figsize=(8, 2))
ax = axes[0]
#world.plot(ax=ax, color='k')
grat.plot(ax=ax, color='lightgrey', linewidth=1)#, alpha=0.5)
controlpts.plot(ax=ax, color='green', marker='x')
gd.plot(ax=ax, color='green', linestyle=':')
circs.plot(ax=ax, color='b')
ax.scatter(*testpt, color='b')
#antigd.plot(ax=ax, color='yellow')
ax.set_aspect('equal', 'box')
ax.set_xlim(-105, -15)
ax.set_ylim(-60, 30)
xticks = ['$90\\degree$W','$60\\degree$W','$30\\degree$W',]
yticks = ['$60\\degree$S', '$30\\degree$S', '$0\\degree$', '$30\\degree$N']
ax.set_xticks(np.linspace(-90, -30, 3))
ax.set_xticklabels(xticks)
ax.set_yticks(np.linspace(-60, 30, 4))
ax.set_yticklabels(yticks)
ax.annotate("$v_1$", actrlpts[:, 0], xytext=(0, 5), ha='right',
            textcoords="offset points")
ax.annotate("$v_2$", actrlpts[:, 1], xytext=(5, -5), ha='left',
            textcoords="offset points")
ax.annotate("$v_3$", actrlpts[:, 2], xytext=(5, -5), ha='left',
            textcoords="offset points")
ax.annotate("$v$", testpt, xytext=(8, -10), ha='center',
            textcoords="offset points")
ax.set_title('a')
ax = axes[1]
#world_lt.plot(ax=ax, color='k')
#grat2_lt.plot(ax=ax, color='lightgrey', linewidth=1)
#gd_lt.plot(ax=ax, linestyle=':', color='green')
ax.plot(tgtpts[0, cyclic], tgtpts[1, cyclic], color='green',
        marker='x', linestyle=':')
for xi, yi, ri in zip(tgtpts[0], tgtpts[1], r):
    c = plt.Circle((xi, yi), ri, color='b', fill=False)
    ax.add_artist(c)

ri2 = np.roll(r, 1)**2
rj2 = np.roll(r, -1)**2
zi = np.roll(tgtpts, 1, axis=-1)
zj = np.roll(tgtpts, -1, axis=-1)
# should work as long as xi != xj
y = np.array([-5E6, 5E6])[..., np.newaxis]
x = ((ri2 - rj2)/2 + y*(zi[1]-zj[1]))/(zj[0]-zi[0])
ax.plot(x, y, color='r', linestyle='--')
ax.scatter(ctpt[0], ctpt[1], color='b')
ax.scatter(ltpt[0], ltpt[1], color='r')
ax.set_aspect('equal', 'box')
ax.set_xlim(-5E3, 5E3)
ax.set_ylim(-5E3, 5E3)
ax.annotate("$p_1$", tgtpts[:, 0], xytext=(0, 5), ha='right',
            textcoords="offset points")
ax.annotate("$p_2$", tgtpts[:, 1], xytext=(5, -5), ha='left',
            textcoords="offset points")
ax.annotate("$p_3$", tgtpts[:, 2], xytext=(5, -5), ha='left',
            textcoords="offset points")
ax.annotate("$p$", (ctpt[0], ctpt[1]), xytext=(0, 10),
            ha='left', textcoords="offset points")
ax.set_title('b')

# figure out how to inset this later? nah
ax = axes[2]
#gd_lt.plot(linestyle=':', ax=ax, color='green')
#ax.plot(tgtpts[0, index], tgtpts[1, index], color='green', marker='o')
for xi, yi, ri in zip(tgtpts[0], tgtpts[1], r):
    c = plt.Circle((xi, yi), ri, color='b', fill=False)
    ax.add_artist(c)
ax.scatter(ctpt[0], ctpt[1], color='b')
ax.scatter(ltpt[0], ltpt[1], color='r')
ax.plot(x, y, color='r', linestyle='--')
ax.set_aspect('equal', 'box')
#ax.set_xlim(ctpt[0] -500, ctpt[0]+500)
#ax.set_ylim(ctpt[1] -500, ctpt[1]+500)
ax.set_xlim(-750, -580)
ax.set_ylim(1550, 1700)
ax.annotate("$p_\ell$", (ltpt[0], ltpt[1]), xytext=(0, 8),
            ha='right', textcoords="offset points")
ax.annotate("$p_c$", (ctpt[0], ctpt[1]), xytext=(5, 0),
            ha='left', textcoords="offset points")
ax.set_title('c')
fig.savefig('construction.eps', bbox_inches = 'tight')
