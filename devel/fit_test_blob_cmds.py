import numpy as np
import fcheck
reload(fcheck)

trnum = 7
trcoo = trace_pixel_pos[trnum]

spec_blob = spec_data[trcoo]
lamp_blob = lamp_data[trcoo]

# ----------------------------------------------------------------------- 
# Start by comparing the results of a fit for the first column:
lamp_vals_i = lamp_blob[:, 0]
spec_vals_i = spec_blob[:, 0]

slope_np_lsq = fcheck.fit_yxratio_numpy_lstsq(lamp_vals_i, spec_vals_i)
slope_direct = fcheck.fit_yxratio_direct_sums(lamp_vals_i, spec_vals_i)
sys.stderr.write("slope_np_lsq: %15.8f\n" % slope_np_lsq)
sys.stderr.write("slope_direct: %15.8f\n" % slope_direct)

# ----------------------------------------------------------------------- 
# Try the whole spectrum (direct method):
ratio_allpix = fcheck.blob_fit_yxratio_direct(lamp_blob, spec_blob)

# Check with plot ...


