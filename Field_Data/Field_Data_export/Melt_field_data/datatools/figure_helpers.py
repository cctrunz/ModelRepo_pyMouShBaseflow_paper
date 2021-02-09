# -*- coding: utf-8 -*-
"""A collection of plotting functions
    """

from datatools.units import convert


def fig_size_mm(w, h, dec=2):
    """return figure size tuple in inches with precision dec {default=2})"""
    return tuple(map(lambda x: round(convert('mm', 'in', x), dec), (w, h)))
