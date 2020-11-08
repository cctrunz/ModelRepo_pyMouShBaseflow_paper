
from math import atan, sin, cos, sqrt
from datetime import timedelta
import pandas as pd

#Function to calculate strain rates between stations
def calcStrainRates(stn1_obj, stn2_obj):
    '''Calculate longitudinal, lateral, and vertical strain
    rates between two gps stations following Hoffman et al., 2011
    Input:
        stn1 - class: gps
        stn2 - class: gps
            
    Output:
        strain_rate pd.df
            strain_rate.longitudinal
            strain_rate.lateral
            strain_rate.vertical

    Equations:
        strain_rate.longitudinal = (1 / initial_distance.xflow) 
                       * (change_in_distance.xflow / 24)
        strain_rate.lateral = (1 / initial_distance.xtran)
                               * (change_in_distance.xtran / 24)
        strain_rate.vertical = -1 * (strain_rate.longitudinal + 
                                     strain_rate.lateral)
    '''
    import pickle
    #? Determine which stations are being input to pull bg vals:
    if 'usf1' in stn1_obj.file or 'usf1' in stn2_obj.file:
        if 'usf3' in stn1_obj.file or 'usf3' in stn2_obj.file:
            with open('/Users/jzm/WorkingDir/GrIS/X_PYTHON-CODES/data/bg_vals/lmid_jnih.txt','rb') as file:
                bg = pickle.load(file)
        elif 'usf2' in stn1_obj.file or 'usf2' in stn2_obj.file:
            with open('/Users/jzm/WorkingDir/GrIS/X_PYTHON-CODES/data/bg_vals/lmid_jeme.txt','rb') as file:
                bg = pickle.load(file)
        elif 'jeme' in stn1_obj.file or 'jeme' in stn2_obj.file:
            with open('/Users/jzm/WorkingDir/GrIS/X_PYTHON-CODES/data/bg_vals/lmid_jeme.txt','rb') as file:
                bg = pickle.load(file)
            
    omega = bg.omega
    dist_bg = bg.dist
    xflow_bg = bg.deltaXflow
    xtran_bg = bg.deltaXtran


    from datetime import timedelta        
    #* Convert station object to dataframe for manipulation, resmple to 30 min
    stn1 = stn1_obj.data.resample('30T', loffset=timedelta(minutes=15)).mean().dropna()
    stn2 = stn2_obj.data.resample('30T', loffset=timedelta(minutes=15)).mean().dropna()


    #! Needs to be tested: may pull errors if there
    #!  is no data at those times. 
    strain_rate = []
    stn_separation = []
    deltaXX0 = deltaYY0 = 0
    counter = 0
    delta_t = timedelta(days=1)
    for idx, val in stn1.iterrows():
        


        try:
            # print(f'data available for station 1 on {day}')
            stn2N_vals = stn2_NEcomp.dNorth[str(day)]
            stn2E_vals = stn2_NEcomp.dEast[str(day)]
            if len(stn2N_vals) > 1:
                stn2N_vals = stn2N_vals[0]
                stn2E_vals = stn2E_vals[0]
            elif len(stn2N_vals) == 0:
                deltaXX0 = deltaYY0 = 0
                # print('No data for station 2')
                continue
            # print(f'Data available for both stationos on day {day}')
            # Calculate distance between stations for
            # Normal coordinate system: North, East, Direct
            deltaN = abs(stn1N_vals) - abs(stn2N_vals)
            deltaE = abs(stn1E_vals) - abs(stn2E_vals)
            direct_dist = sqrt(deltaN**2 + deltaE**2)

            #Transform to transverse/along flow directions
            deltaXX = calc_dist_longitudinal(direct_dist, omega)
            deltaYY = calc_dist_lateral(direct_dist, omega)
            stn_separation.append({'Date': pd.to_datetime(str(day)),
                                    'xflow': deltaXX,
                                    'xtran': deltaYY})

                #calculate strain rate
            if deltaXX0 != 0:
                strainXX = calc_strain_rate(bg.deltaXflow, deltaXX0, deltaXX)
                strainYY = calc_strain_rate(bg.deltaXtran, deltaYY0, deltaYY)
                strainZZ = -1 * (strainXX + strainYY)
                strain_rate.append({'longitudinal':strainXX,
                                    'lateral': strainYY,
                                    'vertical': strainZZ,
                                    'Date': pd.to_datetime(str(day))})
                counter =+ 1  
            deltaXX0 = deltaXX
            deltaYY0 = deltaYY        
        except:
            # deltaXX0 = deltaYY0 = 0
            print(f'Exception Triggered')


    print(f'Number of strain calcs = {counter}')
    strain_df = pd.DataFrame(strain_rate)
    strain_df.index = strain_df['Date']
    return strain_df, stn_separation
    



def calc_dist(stn1, stn2, time):
    delta_north = abs(stn1.dnorth[time]) - abs(stn2.dnorth[time])
    delta_east = abs(stn1.deast[time]) - abs(stn2.deast[time])
    distance=sqrt(delta_north**2 + delta_east**2)
    return delta_north, delta_east, distance

def calc_dist_longitudinal(dist, omega):
    dist_xx = dist * cos(omega)
    return dist_xx

def calc_dist_lateral(dist, omega):
    dist_yy = dist * sin(omega)
    return dist_yy

def calc_strain_rate(dist0, L0, L1):
    strain_rate = (L1-L0)/(dist0*24)
    return strain_rate





# 	def ldetrend(self, v0, t0):        # Linearly detrend v1
# 	    slope, intercept, r_value, p_value, std_err = stats.linregress(t0,v0)
# 	    v1 = v0 - intercept - slope*t0
# 	    return v1


'''

        # Determine start end end dates
        start_date = self.stn1_date[0].replace(hour=0, minute=0, second=0)
        end_date = self.stn1_date[-1].replace(hour=0, minute=0, second=0)
        t0 = start_date
        t1 = t0+timedelta(days=1)
        strain_rate = []
        timevec = []
        for i in range((end_date - start_date).days):
            # only calculate strainrate if there is data
            if (t1 in self.stn1.index) == True &\
                   (t1 in self.stn2.index) == True &\
                   (t0 in self.stn1.index) == True &\
                   (t0 in self.stn2.index) == True:
                # dist_start = self.calc_dist_UTM(t0)
                # dist_end = self.calc_dist_UTM(t1)
                # dL = dist_end - dist_start
                strain = dL / (dist_start * calc_interval)
                strain_rate.append(strain)
                timevec.append(t1)
                t0 = t0 + timedelta(days=1)
                t1 = t1 + timedelta(days=1)
            else:
                t0 = t0 + timedelta(days=1)
                t1 = t1 + timedelta(days=1)
        strain_rate_vec = pd.DataFrame(
            {'Date': timevec, 'strain_rate': strain_rate})
        strain_rate_vec = strain_rate_vec.set_index(['Date'])
        return strain_rate_vec
'''
