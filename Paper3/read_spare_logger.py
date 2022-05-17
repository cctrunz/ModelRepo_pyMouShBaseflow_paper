import pandas as pd
import matplotlib.pyplot as plt

pira = pd.read_csv('SpareLogger_Table1.dat',delimiter=',',skiprows=[0,2,3],index_col=0,parse_dates=True)

plt.figure()
plt.plot(pira.Lvl*(-1))