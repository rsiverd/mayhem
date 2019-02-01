import numpy as np
import os, sys, time
import brute_astrometry
bfa = brute_astrometry.BruteForceAstrometry()

rot_deg = 13.091
rot_deg = -13.091
#rot_deg = 0

nres_rmatrix = bfa.rotation_matrix(np.radians(rot_deg))
sys.stderr.write("\nnres_rmatrix (%.4f deg):\n%s\n\n" 
        % (rot_deg, str(nres_rmatrix)))

x0, y0 = 0.0, 0.0
x0 = 500
y0 = -500
#origin = np.array([x0, y0])[:, None]
#poly_xform = np.hstack((origin, nres_rmatrix))
origin = np.array([x0, y0])
poly_xform = np.vstack((origin, nres_rmatrix))
sys.stderr.write("\npolynomial transformation:\n")
sys.stderr.write("%s\n\n" % str(poly_xform))

## Save ASCII polynomial transformation for is3_interp:
#isis_xform = poly_xform[::-1, :].T
isis_xform = poly_xform[:, ::-1]
ncoeff = isis_xform.shape[0]
save_xform = 'nres_xform.txt'
sys.stderr.write("Saving transform to '%s' ... " % save_xform)

with open(save_xform, 'w') as f:
    f.write("%d\n" % ncoeff)
    for row in isis_xform:
        f.write("%s\n" % ' '.join(['%6f'%x for x in row]))
sys.stderr.write("done.\n")

