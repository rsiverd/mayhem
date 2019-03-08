
#def nearest_xy(x1, y1, x2, y2, xmax=110):
#    tdiff = []
#    for tx,ty in zip(x1, y1):
#        which = (np.abs(tx - x2) <= xmax)
#        if (np.sum(which) > 0):
#            tdiff.append(np.min(np.hypot(tx - x2[which], ty - y2[which])))
#    return np.min(tdiff)
#
### Closest separations:
#def get_minseps(ridges):
#    ndiffs = len(ridges) - 1
#    separations = []
#    for ii in range(ndiffs):
#        x1, y1 = ridges[ii]
#        x2, y2 = ridges[ii+1]
#        tsep = [np.min(np.hypot(x2 - tx, y2 - ty)) for tx,ty in zip(x1,y1)]
#        separations.append(np.min(tsep))
#    return np.array(separations)
#
### Closest separations:
#def get_minseps2(ridges):
#    ndiffs = len(ridges) - 1
#    separations = []
#    for ii in range(ndiffs):
#        x1, y1 = ridges[ii]
#        x2, y2 = ridges[ii+1]
#        sys.stderr.write("matching orders %d,%d\n" % (ii, ii+1))
#        separations.append(nearest_xy(x1, y1, x2, y2))
#    return np.array(separations)
#
##sys.stderr.write("Timing approach 1 ... ")
##tik = time.time()
##yay_seps_1 = get_minseps(fib0_ridges)
##tok = time.time()
##sys.stderr.write("took %.3f seconds.\n" % (tok-tik))
## TAKES 20 SECONDS
#
#sys.stderr.write("Timing approach 2 ... ")
#tik = time.time()
#yay_seps_2 = get_minseps2(fib0_ridges)
#tok = time.time()
#sys.stderr.write("took %.3f seconds.\n" % (tok-tik))
## TAKES 14 SECONDS


