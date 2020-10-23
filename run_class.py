from PyMouSh import MoulinShape
import numpy as np







Q_mean = 1
dQ = 0.1


#Initial paramters -- moulin shape
initial_moulin_radius = moulin.initial_moulin_radius()
Qin = moulin.sinusoidal(Q_mean,dQ)
moulin = MoulinShape(Qin,initial_moulin_radius)


sim = moulin.run1step()
