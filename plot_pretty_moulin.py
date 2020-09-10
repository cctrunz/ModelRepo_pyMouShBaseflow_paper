import matplotlib.pyplot as plt

'''Moulin graph alone'''
fig,ax = plt.subplots()
def live_plot(hw,Mx_upstream,Mx_downstream,z):
    ax.clear() #clear figure content -- especially important for ploting in the loop
    ax.axhspan(0, hw, facecolor ='lightblue', alpha = 1,zorder=1)
    ax.axhspan(-100, 0, facecolor ='peru', alpha = 1,zorder=1)
    ax.plot(Mx_upstream,z,color='black') #plot major axis on the left
    ax.plot(Mx_downstream,z,color='black')  #plot minor axis on the right
    ax.set_xlim([-10,10])    
    ax.fill_betweenx(z,-10,Mx_upstream,color='aliceblue',zorder=2)
    ax.fill_betweenx(z,Mx_downstream,10,color='aliceblue',zorder=2)   