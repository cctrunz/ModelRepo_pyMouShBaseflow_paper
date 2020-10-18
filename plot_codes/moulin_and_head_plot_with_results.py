import matplotlib.pyplot as plt
import numpy as np

fig1 = plt.figure(figsize=(25,8))
grid = plt.GridSpec(4,4)#, wspace=-0.7)
ax1 = fig1.add_subplot(grid[0:4, 0])
ax2 = fig1.add_subplot(grid[0:4, 1:4])#, sharey=ax1)  #hw
ax3 = fig1.add_subplot(grid[2, 1:4])#ax2.twinx() #SCs
ax4 = fig1.add_subplot(grid[3, 1:4], sharex=ax2)    #Qin-Qout
ax5 = ax4.twinx()#fig1.add_subplot(grid[2, 1:4], sharex=ax2)    #Qin-Qout

dlim = 10
dbound = 10
dtick = [-10,5,0,5,10]
    

def live_plot(results,dt,z,Qin,time,H,t_real,h_real,idx=0,
              hw_min = 200
              ):
    

    
    


    Mx_upstream = results['Mx_upstream']
    Mx_downstream = results['Mx_downstream']

    hw = results['hw']
    SCs = results['SCs']
    Qout = results['Qout']

    #MOULIN SHAPE
    ax1.clear() #clear figure content -- especially important for ploting in the loop
    ax1.plot(Mx_upstream[idx],z,color='black') #plot major axis on the left
    ax1.plot(Mx_downstream[idx],z,color='black')  #plot minor axis on the right
    
    ax2.clear()
    ax2.plot(t_real/3600,h_real,'-',color='black')      
    ax2.plot(time[0:idx+1]/3600,hw[0:idx+1],'-',color='blue')   
    ax2.legend(['measured','simulated'],loc="upper left")
    
    ax3.clear()
    ax3.plot(time[0:idx+1]/3600,SCs[0:idx+1],'-',color='red')  
      
    ax4.clear()
    ax5.clear()
    ax4.plot(time/3600,Qin,'--',color='blue')#,label='Qin')    
    ax5.plot(time[0:idx+1]/3600,Qout[0:idx+1],color='red')#,'-',color='red',label='Qout')
    
    ax1.axhspan(0, hw[idx], facecolor ='lightblue', alpha = 1,zorder=1)
    ax1.axhspan(-140, 0, facecolor ='peru', alpha = 1,zorder=1)
    ax1.fill_betweenx(z,-10,Mx_upstream[idx],color='aliceblue',zorder=2)
    ax1.fill_betweenx(z,Mx_downstream[idx],10,color='aliceblue',zorder=2)     
    
    '''Moulin'''
    ax1.set_ylabel('z(m)')    
    ax1.set_xlabel('(m)')
    ax1.set_xlim([-5,5]) 
    ax1.set_ylim([-50,H]) 
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['bottom'].set_position(('zero')) 
    ax1.spines['left'].set_bounds(0,H)
    ax1.spines['bottom'].set_bounds(-5,5)
    ax1.set_xticks([-4,-3,-2,-1,0,1,2,3,4]) 
    ax1.set_xticklabels([-4,-3,-2,-1,0,1,2,3,4]) 
    ax1.set_yticks(np.arange(0,H+1,100)) 
    ax1.set_yticklabels(np.arange(0,H+1,100))

    '''Head'''
    ax2.set_ylabel('hw',color='blue')
    ax2.set_xlim([0,max(time)/3600])
    ax2.set_ylim([-50,H]) 
    ax2.yaxis.tick_left()
    ax2.yaxis.set_label_position("left")
    ax2.tick_params(axis='y', labelcolor='blue')
    ax2.spines['top'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_color('blue')
    ax2.axes.xaxis.set_visible(False)
    ax2.spines['left'].set_bounds(hw_min,H)
    ax2.set_yticks(np.arange(hw_min,H+1,100)) 
    ax2.set_yticklabels(np.arange(hw_min,H+1,100))

    
    '''Subglacial channel'''
    ax3.set_ylabel('Scs',color='red')
    ax3.set_xlim([0,max(time)/3600])
    ax3.set_ylim([min(SCs),max(SCs)])
    ax3.yaxis.tick_right()
    ax3.yaxis.set_label_position("right")
    ax3.tick_params(axis='y', labelcolor='red')
    ax3.spines['top'].set_visible(False)
    ax3.spines['bottom'].set_visible(False)
    ax3.spines['left'].set_visible(False)
    ax3.spines['right'].set_color('red')
    ax3.axes.xaxis.set_visible(False)
    
    
    '''Qin'''
    ax4.set_ylabel('Qin',color='blue')
    ax4.set_xlim([0,max(time)/3600])
    ax4.set_ylim([min(Qin),max(Qin)])
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)
    ax4.spines['left'].set_color('blue')
    ax4.tick_params(axis='y', labelcolor='blue')
    ax4.set_xticks(np.round(np.linspace(1,max(time)/3600,20)))
    ax4.set_xticklabels(np.round(np.linspace(1,max(time)/3600,20)))
    
    '''Qout'''
    ax5.set_ylabel('Qout',color='red')
    ax5.set_xlim([0,max(time)/3600])
    ax5.set_ylim([min(Qout),max(Qout)])
    ax5.yaxis.tick_right()
    ax5.yaxis.set_label_position("right")
    ax5.tick_params(axis='y', labelcolor='red')
    ax5.spines['top'].set_visible(False)
    ax5.spines['left'].set_visible(False)
    ax5.spines['bottom'].set_visible(False)
    ax5.spines['right'].set_color('red')

    
    
    #make backgroud transparent    
    ax1.patch.set_alpha(0)
    ax2.patch.set_alpha(0)
    ax3.patch.set_alpha(0)
    ax4.patch.set_alpha(0)
    ax5.patch.set_alpha(0)
    
 
    
   



    
