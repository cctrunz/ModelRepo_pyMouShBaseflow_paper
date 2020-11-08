from typing import Optional, List, Dict, Tuple, Union, Any

import pandas as pd
import numpy as np


from .constants import ICE_DENSITY, WATER_DENSITY

IntOrFloat = Union[float, int]


def range_from_errors(value: IntOrFloat,
                      error: IntOrFloat,
                      value_in_range=False
                      ) -> Union[Tuple[float, float, float], Tuple[float, float]]:
    return ((value - error, value, value + error) if value_in_range
            else (value - error, value + error))


def overburden_pressure_water_level(ice_thickness: IntOrFloat
                                    ) -> float:
    return (ICE_DENSITY / WATER_DENSITY) * ice_thickness


def mH2O_to_gwl(water_level: IntOrFloat,
                surface_elev: IntOrFloat,
                ice_thickness: IntOrFloat
                ) -> IntOrFloat:
    return water_level + bed_elev(surface_elev, ice_thickness)


def gwl_to_mH2O(gwl: IntOrFloat,
                surface_elev: IntOrFloat,
                ice_thickness: IntOrFloat
                ) -> IntOrFloat:
    return gwl - bed_elev(surface_elev, ice_thickness)


def bed_elev(surface_elev: IntOrFloat,
             ice_thickness: IntOrFloat
             ) -> IntOrFloat:
    return surface_elev - ice_thickness


def overburden_fraction(water_level: IntOrFloat,
                        ice_thickness: IntOrFloat
                        ) -> float:
    return water_level / overburden_pressure_water_level(ice_thickness)


def gwl_overburden_fraction(gwl: IntOrFloat,
                            surface_elev: IntOrFloat,
                            ice_thickness: IntOrFloat
                            ) -> float:
    return overburden_fraction(gwl_to_mH2O(gwl, surface_elev, ice_thickness),
                               ice_thickness)


def percent_overburden(water_level: IntOrFloat,
                       ice_thickness: IntOrFloat
                       ) -> float:
    return overburden_fraction(water_level, ice_thickness) * 100


def wl_to_overburden_fraction(water_level: IntOrFloat,
                              ice_thickness: IntOrFloat,
                              error: Optional[IntOrFloat] = None,
                              value_in_range: Optional[bool] = False,
                              ) -> Any:
    if error is not None:
        return tuple(map(lambda z:
                         overburden_fraction(water_level, z),
                         range_from_errors(ice_thickness, error,
                                           value_in_range=value_in_range)))
    return overburden_fraction(water_level, ice_thickness)


def gwl_to_overburden_fraction(gwl: IntOrFloat,
                               surface_elev: IntOrFloat,
                               ice_thickness: IntOrFloat,
                               error: Optional[IntOrFloat] = None,
                               value_in_range: Optional[bool] = False,
                               ) -> Any:
    return wl_to_overburden_fraction(
        gwl_to_mH2O(gwl, surface_elev, ice_thickness),
        ice_thickness, error=error, value_in_range=value_in_range)


# print(gwl_to_overburden_fraction(550, 769, 500, error=100))

# if water_level_reference == 'wlb':
#     overburden_equ = lambda z: water_level / ((DENSITY_RATIO)*z)

# return (lambda z: (water_level-(surface_elev-z))/(DENSITY_RATIO*z) if error == 0 else
#         tuple(map(overburden_equ, range_from_errors(ic65e_thickness, error,
#                                                     keep_value=add_midpoint))))

# return (tuple(map(lambda z: wlb / ((ICE_DENSITY/WATER_DENSITY) * z),
#                  range_from_errors(ice_thickness, error,
#                                    keep_value=add_midpoint)))
#         if error != 0
#         else wlb / ((ICE_DENSITY/WATER_DENSITY) * ice_thickness)))

# def lc_overburden_fraction(water_level,
#                            water_level_type='gwl',
#                            add_midpoint=False):
#     """overburden fraction for low camp"""
#     return (overburden_fraction(water_level_type,
#                                 765.8, 503, 100,
#                                 add_midpoint=add_midpoint)
#             if water_level_type == 'gwl'
#             else wlb_overburden_fraction=(water_level,
#                                           765.8, 503, 100,
#                                           add_midpoint=add_midpoint))

# def add_overburden_fraction_axis(ax_gwl: matplotlib.axes._subplots.AxesSubplot,
#                                  ice_thickness=503.0,
#                                  surface_elev=765.8):
#     """twinx axis with ground water level altitude in overburden pressure

#     Arguments:
#         ax_gwl {matplotlib.axes._subplots.AxesSubplot} -- gwl axis

#     Keyword Arguments:
#         ice_thickness {float} -- [description] (default: {503.0})
#         surface_elev {float} -- [description] (default: {765.8})

#     creates axis named ax_fob with limits matching gwl for inputs
#     """
#     import matplotlib
#     y1, y2 = ax_gwl.get_ylim()
#     ax_fob.set_ylim(overburden_fraction(y1, surface_elev, ice_thickness,))
