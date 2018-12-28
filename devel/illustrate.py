
def lcheck(thar_data, line_centers, tnum):
    plt.clf()

    odata = thar_data[tnum]
    peaks = line_centers[tnum]

    plot(odata['xpix'], odata['spec'], c='b')
    for ii in peaks:
        axvline(odata['xpix'][ii], ls='--', c='r')

    return

