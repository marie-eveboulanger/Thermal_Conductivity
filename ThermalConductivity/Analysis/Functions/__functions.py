"""
This subsubmodule contains detailled versions of useful functions
"""
import numpy as np
import numpy.polynomial as npp
from ThermalConductivity.Thermometry import seebeck_thermometry


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


def tallahassee_temp(Resistance_heat_off, Resistance_heat_on, T0, order=8):
    """
    Computes either T+ or T- since it's the same calculations

    Parameters:
    ----------------------------------------------------------------------------
    Resistance_heat_off:    1d array
                        The values of the resistance without heat in volts
    Resistance_heat_on:     1d array
                        The values of the resistance with heat in volts
    T0:                     1d array
                        The values of the probe's base temperature in Kelvin
    order:                  int
                        The order of the polynomial fit
    """

    # Polynomial fit on the heat off values to "calibrate"
    # Log is used to increase the precision of the fit
    Coefficients = npp.polyfit(np.log(Resistance_heat_off), np.log(T0), order)

    # Use the fit to get the calibrated heat on temperature
    Temperature = np.exp(npp.polyval(np.log(Resistance_heat_on), Coefficients))

    return Temperature


def compute_thermocouple(dT_heat_off, dT_heat_on, T_reference, gain=1000):
    """
    Computes the calibrated difference in temperature from a thermocouple

    Parameters:
    ----------------------------------------------------------------------------
    dT_heat_off:    1d array
                    The difference in temperatures when the heater is off
    dT_heat_on:     1d array
                    The difference in temperatures when the heater is on
    T_ref:          The temperature of the point at which the thermocouple
                    is connected
    gain:           float or int
                    The preamp gain
    """

    seebeck_coefficients = seebeck_thermometry(T_reference)
    delta_T = abs(dT_heat_on-dT_heat_off)/seebeck_coefficients/gain

    return delta_T


def vti_calibration_loop(dT_abs_off, dT_abs_on, dTx_off, dTx_on, T0):
    """
    Computes a calibration loop for the VTI's thermocouple

    Parameters:
    ----------------------------------------------------------------------------
    dT_abs_off:     1d array
                    The readings of the absolute thermocouple with no heat
    dT_abs_on:      1d array
                    The readings of the absolute thermocouple with heat
    dTx_off:        1d array
                    The reading of the thermocouple along x with no heat
    dTx_on:         1d array
                    The reading of the thermocouple along x with heat
    T0:             1d array
                    The temperature of the probe
    """

    # First iteration
    dT_abs = 0*T0
    dTx = 0*T0
    T_minus = 0*T0
    T_plus = 0*T0
    T_av = 0*T0
    previous_T_av = T_av+1000

    while abs(previous_T_av.sum()-T_av.sum()) > 1e-10:
        previous_T_av = T_av
        T_ref_1 = T0+dT_abs/2
        T_ref_2 = T_minus+dTx/2
        dT_abs = compute_thermocouple(dT_abs_off, dT_abs_on, T_ref_1)
        dTx = compute_thermocouple(dTx_off, dTx_on, T_ref_2)
        T_minus = T0+dT_abs
        T_plus = T_minus+dTx
        T_av = T_minus+dTx/2

    return T_av, dTx, T_plus, T_minus
