"""
Functions used for tracer injection discharge calculations.
"""

def calc_Q_const(C_stream, pump_rate, C_injection):
    """
    Calculates discharge for constant rate injection.

    Parameters
    ----------
    C_stream : float
        The concentration values from the stream.
    
    pump_rate : float
        The rate at which the injection solution is being pumped into the stream. Discharge will be calculated in units of this pump rate.

    C_injection : float
        The concentration of the injection solution (must be in same units as stream concentration).

    Returns
    -------
    Discharge : float
         Discharge in units of the provided pump rate.

    """
    return C_injection*pump_rate/C_stream

def calc_Q_pulse(C_stream, V_dye, C_dye, dt = 5., C_background=0.):
    """
    Calculates discharge for a pulse injection.

    Parameters
    ----------
    C_stream : float
        The concentration values from the stream, trimmed to the desired time range.
    
    V_dye : float
        The volume of concentrated dye solution (usually 1 mL).

    C_dye : float
        The concentration of the concentrated dye solution (must be in same units as stream concentration).
    
    dt : float
         The time difference between concentration measurements. Default is 5 seconds.

    C_background : float
        The background dye concentration. Default is 0.

    Returns
    -------
    Discharge : float
         Discharge in units of the provided volume/dt.
    """  

    return C_dye*V_dye/(sum(C_stream.values - C_background)*dt)




