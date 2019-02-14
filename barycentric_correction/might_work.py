# Distilled from:
# https://stackoverflow.com/questions/52900678/coordinates-transformation-in-astropy

from astropy import coordinates as coord
from astropy import units as uu
import astropy.time as astt
#from astropy.coordinates import EarthLocation as 


ipath = './sqa0m801-en03-20151019-0007-e00.fits.fz'
ivals, hkeys = pf.getdata(ipath, header=True)

# For LCO FITS data (WGS84, meters):
geo_keys = ['OBSGEO-X', 'OBSGEO-Y', 'OBSGEO-Z']
geo_xyz  = [hkeys.get(kk) for kk in geo_keys]
#cartrep = coord.CartesianRepresentation(*geo_xyz, unit=uu.m)
obs_loc  = coord.EarthLocation.from_geocentric(*geo_xyz, unit=uu.m)
obstime  = astt.Time( .... , scale='utc', location=obs_loc)

ip_peg = coord.SkyCoord("23:23:08.55", "+18:24:59.3",
               unit=(uu.hourangle, uu.deg), frame='icrs')

targ_coo = ip_peg
#targ_coo = coord.SkyCoord( ... )

ltt_bary = obstime.light_travel_time(targ_coo)
ltt_bary_jpl = obstime.light_travel_time(targ_coo, ephemeris='jpl')

time_barycen = obstime.tdb + ltt_bary
#loc = EarthLocation(lon=

