import pandas as pd
from datetime import timedelta, datetime
# import numpy as np
import scipy
from scipy import stats, signal
# ------------------------------------------------------------------------------

from pandas import read_csv, DataFrame
from pathlib import Path


def JAR1_df(file=Path.cwd().joinpath('data_csv', 'JAR117_WEA_2.csv')):
    avgT = []
    df = pd.read_csv(file, parse_dates=True, index_col=0)
    for idx, data in df[['Temp_1G', 'Temp_2H', 'Temp_1I', 'Temp_2J']].iterrows():
        avgT.append(data.mean())
    df['Temp_avg'] = avgT
    return df


# gcnet_stndata = {
    # 'JAR1': }


# class gcnet:

#     def __init__(self, stn='JAR1'):


def read_hobo_csv(csv_file, all_columns=False):
    """
    Reads data from a csv file exported from HOBOware.

    Parameters
    ----------
    csv_file : string
        A string containing the file name of the csv file to be read.
    all_columns : boolean (optional)
        Determines whether to read in all columns or just ones that we
        search for and relabel
        (RH, DewPt, Abs Pres, Temp, Attached, Stopped, Connected, EOF,
        Cond High Range, 
        Cond Low Range). 

        Default = False

    Returns
    -------
    df : pandas.DataFrame
        DataFrame containing data from HOBO csv file.
    """
    skiprows = 1
    index_col = 1
    parse_dates = True
    df = read_csv(csv_file, skiprows=skiprows, index_col=index_col,
                  parse_dates=parse_dates, na_values=['-888.88', '-888.9'])
    # Convert column names into something nicer
    columns = df.columns
    old_columns = columns
    rename_dict = {}
    # cond_count = 0
    solar_count = 1
    """
    If you want to search for a term and replace a column name with something different 
    use python tuples instead, example:

    new_names = (('find_this', 'replace_with_this'),) 

    use it like this:

    for old, new in  new_names:
        if old in label:
            new_name = new
            wantcol = True
    """
    new_names = ['RH', 'Gust', 'Wind Speed',
                 'Wind Direction', 'DewPt', 'Abs Pres', 'Rain', 'Temp']
    # cols = (old_columns, new_names)

    del df['#']
    for label in columns:
        new_name = label

        for name in new_names:
            if name in label:
                new_name = name

        # Account for multiple Solar Radiation Sensors and name them differently.
        # NOTE: current code only allows for up to two solar rad sensors <JZM>
        if 'Solar' in label:
            if solar_count == 1:
                new_name = 'Solar1'
                solar_count = 2
            elif solar_count == 2:
                new_name = 'Solar2'
                solar_count = 3
            else:
                print(">2 Solar Rad Sensors Detecteds")

        rename_dict[label] = new_name
    df = df.rename(columns=rename_dict)

    return df

# def import_hobo(file):
#     '''import hobo weather data from input file and convert to local time'''
#     data = pd.read_csv(file, index_col=1, skiprows=2, parse_dates=True)
#     data.index = data.index - timedelta(hours=3)
#     data['SolarAdj'] = scipy.ndimage.filters.median_filter(data['Solar'], 10)
#     return data


def read_and_rename_hobo(file, tz_info='UTC'):
    # print("reading in data from: {}".format(file))
    df = read_hobo_csv(file)
    # Determine and rename incoming and reflected solar data
    if 'Solar2' in df:
        if df['Solar1'].mean() < df['Solar2'].mean():
            df = df.rename(columns={'Solar1': 'Reflected',
                                    'Solar2': 'Solar'})
        elif df['Solar1'].mean() > df['Solar2'].mean():
            df = df.rename(columns={'Solar1': 'Solar',
                                    'Solar2': 'Reflected'})
        df['Solar_corrected'] = scipy.ndimage.filters.median_filter(
            df['Solar'], 10)
    elif df.index[0].year == 2018:
        if 'LOWC' in file:
            df = df.rename(columns={'Solar1': 'Solar'})
            df['Solar_corrected'] = scipy.ndimage.filters.median_filter(
                df['Solar'], 10)
        elif 'HIGH' in file:
            df = df.rename(columns={'Solar1': 'Reflected'})
    df = df.rename(columns={'Wind Speed': 'Wind_speed'})
    df.tz_localize(tz_info)
    return df


def read_JAR1_data(file, num_cols=20):
    """Read in GC-NET weather data

    Arguments:
        file {str} -- file path

    Keyword Arguments:
        num_cols {int} -- number of data columns (default: {20})
                            acceptable values {22}

    Returns:
        jar1 {dataframe} -- weather data for station
                            index: date

    """

    if num_cols == 22:
        col_names = ['stnnum', 'year', 'doy', 'SW_down', 'SW_up',
                     'Net_Radiation', 'Temp_1G', 'Temp_2H', 'Temp_1I',
                     'Temp_2J', 'RH_1K', 'RH_2L',
                     'Wind_Speed', 'Wind_Dir',
                     'Atm_Pressure', 'Snow_Height_1R', 'Snow_Height_2S',
                     'RAD1', 'RAD2', 'NetRadMax', 'Albedo', 'Zenith_Angle']
        skip = 23
    elif num_cols == 20:
        col_names = ['stnnum', 'year', 'doy', 'SW_down', 'SW_up',
                     'Net_Radiation', 'Temp_1G', 'Temp_2H', 'Temp_1I',
                     'Temp_2J', 'RH_1K', 'RH_2L',
                     'Atm_Pressure', 'Snow_Height_1R', 'Snow_Height_2S',
                     'RAD1', 'RAD2', 'NetRadMax', 'Albedo', 'Zenith_Angle']
        skip = 21

    jar1 = pd.read_csv(file, sep=' ', skiprows=skip,
                       names=col_names, na_values=['999.0000', '999.0'])

    df = []
    for idx, val in jar1.iterrows():
        if val.year == 2017:
            epoch = pd.to_datetime('2017-01-01')
        elif val.year == 2018:
            epoch = pd.to_datetime('2018-01-01')
        elif val.year == 2019:
            epoch = pd.to_datetime('2019-01-01')
        else:
            print(f'{idx}, {val.year}')
        df.append({'Date': epoch+timedelta(days=(val.doy-1))})
    df = pd.DataFrame(df)
    jar1['Date'] = df
    jar1 = jar1.set_index('Date')

#     P_atm_mbar = jar1.Atm_Pressure.resample('15T').pad()
#     df2 = pd.DataFrame(P_atm_mbar)
#     #* Convert mbar to ft H2O and m H2O
#     df2['P_atm_ftH2O'] = df2.Atm_Pressure * mbar2ftH2O
#     df2['P_atm_mH2O'] = df2.Atm_Pressure * mbar2mH2O
    return jar1


class weather_data:
    def __init__(self, file: str):
        self.file = file
        self.data = read_and_rename_hobo(self.file)
        self.Temp = self.data['Temp']
        self.relative_humidity = self.data['RH']
        if 'Solar' in self.data:
            self.solar = self.data['Solar']
            self.solar_corrected = self.data['Solar_corrected']
        if 'Reflected' in self.data:
            self.reflected = self.data['Reflected']
        if 'Wind Direction' in self.data:
            self.wind_direction = self.data['Wind Direction']
            self.data = self.data.rename(columns={'Wind Direction':
                                                  'Wind_direction'})
        self.wind_speed = self.data['Wind_speed']
        self.gust = self.data['Gust']
        if 'Rain' in self.data:
            self.rain = self.data['Rain']
        self.units = 'Temp - Deg C\n Melt: mm mw h^-1'

    def __repr__(self):
        return 'weather data object'

    def __str__(self):
        return 'Weather data from {} to {}'.format(
            self.Temp.index[0], self.Temp.index[-1])

    def calc_melt(self, incoming_shortwave_radiation=None):
        '''
        Calculate hourly melt rates (mm m.w. equivalent h^-1)
        using the enhanced temperature-index glacier melt model


        M = {   TF * T + SRF * (1 - alpha) * G  if T >  TT
                0                               if T <= TT
        where
            M - Melt Rate (mm per hour)
            TF - Temp Factor (mm h^-1 C^-1)
            T - mean air Temp of each time step (C)
            SRF - Shortwave radiation factor (m^2 mm W^-1 h^-1)
            alpha - albedo
            G - incoming shortwave radiation (W m^-2)
            TT - threshold Temp (0 deg C)

        Ref:
        ----
        Pellicciotti et al. (2005). An enhanced temperature - index 
            glacier melt model including the shortwave radiation 
            balance: development and testing for Haut Glacier 
            d'Arolla, Switzerland. J. Glac. 51(175), 573-587.

        Parameters
        ---------
        temperature :
            Air Temperature
            Units: degree C

        incoming_shortwave_radiation :
            Incoming shortwave solar radiation
            Units: W m^-2
        albedo :

        Other Parameters
        ----------------
        temperature_threshold : 
            temperature at or below no melt will occur
            deg C
            default value = 0


        Output
        ------
            df : pd.dataframe
                data frame (time-indexed)
                melt rates: mm w.e. h^-1

            temp_threshold : optional, int
                determines the cut off for melt to occur
                default = 0 deg C
            ice_snow_transition : optional, string
                date string in a format readable by pd.to_datetime()
                if ice_snow_transition=2017, '2017-08-17 00:00:00' is used
                if ice_snow_transition=None, no transition is given
                    bare ice albedo used for entire calc period

        '''
        from math import isnan
        # Define constants:
        TF = 0.05
        SRF = 0.0094
        threshold_temp = 0
        melt_rate = []

        if incoming_shortwave_radiation is not None:
            self.data['Solar'] = incoming_shortwave_radiation
        # Determine if there is both incoming & reflected rad
        try:
            print("Calculating melt rate with daily albedo values\n\
                    melt rate units: mm w.e. h^-1")
            wea_df = pd.DataFrame({'Temp': self.Temp,
                                   'Solar': self.solar_corrected})
            # wea_df = dropna(wea_df)
            daily_peak_solar = self.solar_corrected.between_time('15:00', '15:10',
                                                                 include_start=True)
            daily_peak_reflected = self.reflected.between_time('15:00',
                                                               '15:10', include_start=True)

            alpha = daily_peak_reflected / daily_peak_solar
            albedo_df = pd.DataFrame(alpha, columns=['Albedo'])

            df = []
            for j in albedo_df.index.to_period("D"):
                for i, k in wea_df[str(j)].iterrows():
                    df.append({'Albedo': albedo_df[str(j)].get_values()[0][0],
                               'Temp': k[0], 'Solar': k[1], 'Date': i})
            df = pd.DataFrame(df).set_index('Date')
            melt_rate = []
            albedo_corrected = []
            for index, j in df.iterrows():
                if isnan(j['Temp']) == True or j['Temp'] <= threshold_temp:
                    melt_rate.append({'Date': index, 'melt_rate': 0})
                # elif j['Temp'] <= threshold_temp:
                #     melt_rate.append({'Date': index, 'melt_rate': 0})
                elif j['Temp'] > threshold_temp:
                    if j.Albedo <= 1:
                        melt_rate.append({'Date': index, 'melt_rate':
                                          (TF * j['Temp'])
                                          + (SRF * (1 - j['Albedo']) * j['Solar'])})
                        albedo_corrected.append({'Date': index,
                                                 'Albedo': j['Albedo']})
                        last_useable_albedo = j.Albedo
                        # last_useable_albedo_date = index

                    else:
                        try:
                            melt_rate.append({'Date': index, 'melt_rate':
                                              (TF * j['Temp'])
                                              + (SRF * (1 - last_useable_albedo) * j['Solar'])})
                            albedo_corrected.append({'Date': index,
                                                     'Albedo': last_useable_albedo})
                            # print(f'')
                            # print(f'{index} {j.Albedo:0.4f} --> {last_useable_albedo:0.4f} from {last_useable_albedo_date} ')
                        except:
                            continue
            albedo_corrected = pd.DataFrame(albedo_corrected).set_index('Date')

        except AttributeError:
            try:
                print("Calculating melt rate with albedo=0.7\n\
                        melt rate units: mm w.e. h^-1")
                for index, j in self.data.iterrows():
                    if isnan(j['Temp']) == True or j['Temp'] <= threshold_temp:
                        melt_rate.append({'Date': index, 'melt_rate': 0})
                    elif j['Temp'] > threshold_temp:
                        melt_rate.append({'Date': index, 'melt_rate':
                                          (TF * j['Temp'])
                                          + (SRF * (1 - 0.7) * j['Solar'])
                                          })

                melt_rateDF = pd.DataFrame(melt_rate)
                melt_rateDF = melt_rateDF.set_index(['Date'])

            except:
                print("ERROR: Data set has no incoming solar radiation data\n\
                        can not preform melt rate calculation")

        self.melt_rate = melt_rateDF
        self.data['Melt_rate'] = self.melt_rate

        try:
            self.data['Albedo'] = albedo_corrected['Albedo']
        except:
            print('No Albedo Calculated')
        return melt_rateDF
