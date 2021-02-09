import matplotlib.pyplot as plt
import numpy as np

def generate_semicircle(center_x, center_y, radius, stepsize=0.1):
    """
    generates coordinates for a semicircle, centered at center_x, center_y
    """        
    y = np.arange(center_y, center_y+radius+stepsize, stepsize)
    x = np.sqrt(radius**2 - y**2)
    # since each x value has two corresponding y-values, duplicate x-axis.
    # [::-1] is required to have the correct order of elements for plt.plot. 
    y = np.concatenate([y,y[::-1]])
    # concatenate y and flipped y. 
    x = np.concatenate([x,-x[::-1]])
    return x, y + center_y

x,y = generate_semicircle(0,0,0.1, 0.01)
plt.plot(x, y)