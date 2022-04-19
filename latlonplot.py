#!/usr/bin/.venv python
# -------------------------------------------------------------------
#   Filename:  latlonplot.py
#   Purpose:   
#   Author:    Will-Sturgeon
#   Comments:  
#   Version:   0.1
# -------------------------------------------------------------------

# -----------------------------------------------------------------------
# ---------------Import required Modules (Python and ...) ---------------
# -----------------------------------------------------------------------

import pyssht
import numpy as np
from scipy import sparse
import pygmt

import pandas as pd
import xarray as xr
from scipy.interpolate import griddata

# -----------------------------------------------------------------------
# ------------------------------ FUNCTIONS ------------------------------
# -----------------------------------------------------------------------


# -----------------------------------------------------------------------
# ---------------------------- LOAD .NPZ file ---------------------------
# -----------------------------------------------------------------------

L=28
pathmatrix = sparse.load_npz("/Users/willsturgeon/Documents/SGLOBEv2/npz_files/0S018.npz")

pathmatrix = pathmatrix.toarray()
pathx=(pathmatrix)
pathsum = (pathmatrix.sum(axis=0)).reshape((L, 2*L-1))

theta, phi = pyssht.sample_positions(L, Grid=True)
lat = np.degrees(np.pi / 2 - theta).flatten()
lon = np.degrees(phi).flatten()

data = pathsum.flatten()

# -----------------------------------------------------------------------
# --------------------------- PLOT .NPZ FILE ----------------------------
# -----------------------------------------------------------------------

coordinates0 = np.column_stack((lon,lat))

# ## Create structured data for plotting
minlon, maxlon = 0,360
minlat, maxlat = -90,90
step = 0.5

lons = np.arange(minlon, maxlon, step)
lats = np.arange(minlat, maxlat, step)

## interpolate data on spatial grid
xintrp, yintrp = np.meshgrid(lons, lats)
z1 = griddata(coordinates0, data, (xintrp, yintrp), method='cubic') #cubic interpolation
xintrp = np.array(xintrp, dtype=np.float32)
yintrp = np.array(yintrp, dtype=np.float32)

## xarray dataarray for plotting using pygmt
da = xr.DataArray(z1,dims=("lat", "long"),
coords={"long": lons, "lat": lats},)

frame =  ["WSen"]

# Visualization
fig = pygmt.Figure()

# make color pallets
lim=abs(max(data.min(),data.max()))
# print(f'{data.min():.2f}/{data.max():.2f}')
pygmt.makecpt(
    #cmap='red,white,blue',
    cmap='plasma',
    series=f'{data.min()}/{data.max()}/0.01',
    #series=f'-{lim}/{lim}/0.01',
    continuous=True
)

#plot tomography
fig.grdimage(
    region="g",
    grid=da,
    projection='W360/10c',
    interpolation='l'
    )

# plot coastlines
fig.coast(
    region="g", 
    shorelines=True,
    frame=frame,
    area_thresh=1000
    )

## Plot colorbar
# Default is horizontal colorbar
fig.colorbar(
    frame='+l"Ray coverage"'
    )

# save figure
save_fig = 1
if not save_fig:
    fig.show() 
    #fig.show(method='external') #open with the default pdf reader
else:
#    fig.savefig("0S018_raydensity.png", crop=True, dpi=300, transparent=True)
    fig.savefig("0S018_raydensity.pdf", crop=True, dpi=720)
    print('Figure saved!')