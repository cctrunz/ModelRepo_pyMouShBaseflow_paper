import matplotlib.pyplot as plt
import numpy as np

fig1 = plt.figure(figsize=(25,8))
grid = plt.GridSpec(2,4)#, wspace=-0.7)
ax1 = fig1.add_subplot(grid[0:3, 0])
ax7 = fig1.add_subplot(grid[1, 1:4])#, sharey=ax1) 
ax8 = ax7.twinx()
ax9 = fig1.add_subplot(grid[0, 1:4], sharex=ax7)   


dlim = 10
dbound = 10
dtick = [-10,5,0,5,10]
    

def live_plot(results,dt,z,Qin,time,H,t_real,h_real,idx=0):
    Mx_upstream = results['Mx_upstream']
    Mx_downstream = results['Mx_downstream']

    hw = results['hw']
    SCs = results['SCs']
    Qout = results['Qout']

    #MOULIN SHAPE
    ax1.clear() #clear figure content -- especially important for ploting in the loop
    ax1.plot(Mx_upstream[idx],z,color='black') #plot major axis on the left
    ax1.plot(Mx_downstream[idx],z,color='black')  #plot minor axis on the right
    
    ax7.plot(t_real/3600,h_real,'-',color='black')      
    ax7.plot(time[0:idx+1]/3600,hw[0:idx+1],'-',color='blue')     
    ax8.plot(time[0:idx+1]/3600,SCs[0:idx+1],'-',color='red')  
    
    ax9.plot(time/3600,Qin,'--',color='blue')#,label='Qin')
    ax9.plot(time[0:idx+1]/3600,Qout[0:idx+1],color='red')#,'-',color='red',label='Qout')
    ax9.legend(['Qin','Qout'],loc="upper right")
    

    ax1.set_ylabel('z(m)')
    
    ax1.set_xlabel('(m)')

    ax7.set_ylabel('hw',color='blue')
    ax8.set_ylabel('Scs',color='red')
    ax9.set_ylabel('Qin')
    #ax9.set_ylabel('Qout')

    
    ax1.set_xlim([-5,5]) 

    ax7.set_xlim([0,max(time)/3600])
    ax8.set_xlim([0,max(time)/3600])
    ax9.set_xlim([0,max(time)/3600])
    ax9.set_xlim([0,max(time/3600)])
    
    ax7.set_ylim([0,H])
    ax8.set_ylim([min(SCs),max(SCs)])
    ax9.set_ylim([0,max(Qin)])
    ax9.set_ylim([0,max(Qin)])
    
    ax7.yaxis.tick_left()
    ax7.yaxis.set_label_position("left")
    ax7.tick_params(axis='y', labelcolor='blue')
    
    ax8.yaxis.tick_right()
    ax8.yaxis.set_label_position("right")
    ax8.tick_params(axis='y', labelcolor='red')
    
    
    
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    
    
    ax7.spines['top'].set_visible(False)
    ax7.spines['bottom'].set_visible(False)
    
    ax8.spines['top'].set_visible(False)
    ax8.spines['bottom'].set_visible(False)
    
    ax9.spines['top'].set_visible(False)
    ax9.spines['right'].set_visible(False)
    
    ax9.spines['top'].set_visible(False)
    ax9.spines['right'].set_visible(False)
    #ax9.spines['bottom'].set_visible(False)
    

    
    ax7.axes.xaxis.set_visible(False)
    
    
    ax1.spines['bottom'].set_position(('zero'))
 
    
    ax1.spines['left'].set_bounds(0,H)
    ax1.spines['bottom'].set_bounds(-5,5)
 
    
    ax1.set_xticks([-4,-3,-2,-1,0,1,2,3,4]) 

    ax9.set_xticks(np.round(np.linspace(1,max(time)/3600,20)))
    ax1.set_xticklabels([-4,-3,-2,-1,0,1,2,3,4]) 

    ax9.set_xticklabels(np.round(np.linspace(1,max(time)/3600,20)))
    
 
    
    ax1.axhspan(0, hw[idx], facecolor ='lightblue', alpha = 1,zorder=1)
    ax1.axhspan(-140, 0, facecolor ='peru', alpha = 1,zorder=1)
    ax1.fill_betweenx(z,-10,Mx_upstream[idx],color='aliceblue',zorder=2)
    ax1.fill_betweenx(z,Mx_downstream[idx],10,color='aliceblue',zorder=2)        



    
