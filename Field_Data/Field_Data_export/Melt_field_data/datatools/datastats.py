# Perform daily calculations on data
import pandas as pd
from datetime import timedelta

def calcDailyMinMax(data, col_name='val'):
    """
    Determine daily min and max values with indexes
    
    Arguments:
        data {series, float64} -- input data, index - datetime
    
    Returns:
        df_min {dataframe, float64} -- daily minimum values
        df_max {dataframe, float64} -- daily maximum values
    """
    #check if the input data's dtype is float64
    if data.dtype != 'float64':
        try:
            data = data.astype(float)
        except:
            print(f'ERROR:input data dtype is not float64')

    issame = 1
    df_max = []
    df_min = []

    for day in data.index.to_period('D'):
        #Only perform calculations once per unique day in index
        if day != issame:
            issame = day
            try:
                #Return maximum/minimum values and index
                valMin = data[str(day)].min()
                idxMin = data[str(day)].idxmin()
                # df_min.append({'Date': idxMin, 'Velocity': valMin})


                valMax = data[str(day)].max()
                idxMax = data[str(day)].idxmax()

                #Make sure maximum values occur after minimums
                if idxMax < idxMin:
                    # next_day = str(pd.to_datetime(str(day))+timedelta(days=1))
                    valMax = data[str(idxMin):str(day)].max()
                    idxMax = data[str(idxMin):str(day)].idxmax()
                
                #if daily maximum is after 23:00 allow to go into next day
                if idxMax.hour >= 23:
                    t_start = str(idxMax - timedelta(hours=6))
                    t_end = str(idxMax + timedelta(hours=6))
                    valMax = data[t_start:t_end].max()
                    idxMax = data[t_start:t_end].idxmax()


                if idxMin.hour >= 23:
                    t_start = str(idxMax - timedelta(hours=idxMax.hours))
                    t_end = str(idxMax)
                    valMin = data[t_start:t_end].min()
                    idxMin = data[t_start:t_end].idxmin()
                    
                #append data to series
                df_min.append({'date': idxMin, col_name: valMin})
                df_max.append({'date': idxMax, col_name: valMax})
                

            except:
                continue

    #convert to dataframes
    df_max = pd.DataFrame(df_max)
    df_max['date'] = pd.to_datetime(df_max['date'])
    df_max = df_max.set_index('date')
    df_min = pd.DataFrame(df_min)
    df_min['date'] = pd.to_datetime(df_min['date'])
    df_min = df_min.set_index('date')
    return df_min, df_max



