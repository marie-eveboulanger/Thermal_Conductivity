"""
This subsubmodule contains detailled versions of useful functions
"""


def compute_kxx(I, dTx, width, thickness, length, Resistance=5000):
    """
    Computes the longtitudinal thermal grandient aka kxx

    Parameters:
        ----------------------------------------------------------------------------
        I:      1d array
        The values of the current that flows through the heater.
        dTx:    1d array
        The values of the longitudinal temperature gradient in the sample.
        width:  float
        The sample's width in meters
        thickness:  1d array
        The sample's thickness in meters
        length: 1d array
        The sample's length in meters
        Resistance: float or int
        The heater's resistance in ohm
        """

    Q = Resistance*I*I
    alpha = (width*thickness)/length  # The geometric factor
    kxx = Q/(dTx*alpha)

    return kxx


def compute_kxy(kxx, dTx, dTy, width, length, Resistance=5000):
    """
    Computes the transverse thermal grandient aka kxy

    Parameters:
    ----------------------------------------------------------------------------
    kxx:      1d array
    The values of the longitudinal thermal gradient.
    dTx:    1d array
    The values of the longitudinal temperature gradient in the sample.
    dTy:    1d array
    The values of the transverse temperature gradient in the sample.
    width:  float
    The sample's width in meters
    length: 1d array
    The sample's length in meters
    """

    geo_factor = length/width
    delta_ratio = dTy/dTx
    kxy = kxx*delta_ratio*geo_factor

    return kxy
