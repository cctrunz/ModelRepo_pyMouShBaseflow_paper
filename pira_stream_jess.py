
from pathlib import Path
import pickle

import pandas as pd


def c_rolling(data: FrameOrSeries,
              window: str,
              window_func='mean',
              min_periods: Optional[int] = None,
              win_type: Optional[str] = None
              ) -> FrameOrSeries:
    """apply and center datetime index of pandas rolling window function

    Args:
        data (pd.Series, pd.DataFrame): pandas datetime index obj.
        window (str): pandas recognized time string
        window_func (str): rolling window function to apply.
            Defaults to mean (rolling average)
        min_periods (int): Minimum number of observations in window 
            required to have a value. Defaults to None.
        win_type (str): Window type for rolling window calculation.
            If None, all points are evenly weighted. Defaults to None.

    Returns:
        Series or DataFrame: rolling window function with centered index
    """
    rolled = getattr(data.rolling(
        window, min_periods=min_periods, win_type=win_type), window_func)()
    rolled.index = rolled.index - (pd.Timedelta(window) / 2)
    return rolled


def pira_stream(for_diurnal_calc=False):
    """generate and return a dataframe with all pira stream measurements

    returns
        pira_stream (pd.DataFrame):

    """
    with open('./data/ADAT18_HYD_1.txt', 'rb') as file:
        pira = pickle.load(file)

    stage_m = (2.26-convert('in', 'm', pira.Stage[:'2018-9-27'])).dropna()

    pira_dye = pd.read_csv(
        Path.cwd()/'data_csv/PIRA_Qdye.csv', index_col=0, parse_dates=True)
    pira_dye.rename(
        columns={'Discharge (m^3/s)': 'dye_discharge'}, inplace=True)

    pira_discharge = pd.read_csv(script_dir.parent/'C_DISCHARGE-MEASUREMENT' /
                                 'IceRay_LowCamp_2018'/'Low-Camp-Discharge.csv',
                                 header=30, index_col=0, parse_dates=True,
                                 names=['date', 'flow_rate'])
    pira_level = pd.read_csv(script_dir.parent/'C_DISCHARGE-MEASUREMENT' /
                             'IceRay_LowCamp_2018'/'Low-Camp-Level.csv',
                             header=30, names=['date', 'ir_stage'],
                             index_col=0, parse_dates=True)

    r_trunc = c_rolling(pira_level[:'2018-7-12 23:10'], '2H').append(
        c_rolling(pira_level['2018-7-12 23:10':], '2H', min_periods=100))

    pira_stream = stage_m.dropna()[:'2018-09-27'].to_frame().join(
        r_trunc.resample('15T').interpolate(
            method='linear', limit_direction='forward'), how='outer')
    pira_stream['combined_stage'] = pira_stream['ir_stage'].dropna().append(
        pira_stream['Stage']['2018-07-19 19':] - 0.1)
    # add discharge measurments (yeilds 15minute sampling rate)
    pira_stream['ir_discharge'] = c_rolling(pira_discharge[1:], '2H')
    pira_stream['dye_discharge'] = c_rolling(pira_dye, '2H').resample(
        '1T').interpolate(method='linear', limit_direction='forward', limit=5)

    return pira_stream if for_diurnal_calc is False else c_rolling(
        pira_stream['combined_stage'][:'2018-08-20'], '6H')
