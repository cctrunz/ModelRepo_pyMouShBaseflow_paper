import matplotlib.pyplot as plt
import numpy as np

fig1 = plt.figure(figsize=(20,15))
grid = plt.GridSpec(4,7)#, wspace=-0.7)
ax1 = fig1.add_subplot(grid[0:2, 0:3])
ax2 = fig1.add_subplot(grid[0:2, 3], sharey=ax1)
ax3 = fig1.add_subplot(grid[0:2, 4], sharey=ax1)
ax4 = fig1.add_subplot(grid[0:2, 5], sharey=ax1)
ax5 = fig1.add_subplot(grid[0:2, 6], sharey=ax1)
ax7 = fig1.add_subplot(grid[3, 0:7]) 
ax8 = ax7.twinx()
ax9 = fig1.add_subplot(grid[2, 0:7])  


dlim = 10
dbound = 10
dtick = [-10,5,0,5,10]
    

def live_plot(results,dt,z,Qin,time,H,plot_interval):
    Mx_upstream = results['Mx_upstream']
    Mx_downstream = results['Mx_downstream']
    dTM = results['dTM']
    dC_major = results['dC_major']
    dC_minor = results['dC_minor']
    dE_major = results['dE_major']
    dE_minor = results['dE_minor']
    dr_major = results['dr_major']
    dr_minor = results['dr_minor']
    #dGlen = results['dGlen']
    #t = results['t']
    hw = results['hw']
    SCs = results['SCs']
    #z = results['z']
    #Qin = results['Qin']
    Qout = results['Qout']


    #time = results['time']
    #H = results['H']
    mts_to_cmh = 100*60*60/dt
    
    idx_plot = np.arange(0,len(time),plot_interval)
    idx_save = 0
    
    for idx in idx_plot:
        idx_save = idx_save+1
        # if idx_save/len(idx_plot)*100 %10 == 0:
        #     print('percentage plot:',idx_save/len(idx_plot))

        #MOULIN SHAPE
        ax1.clear() #clear figure content -- especially important for ploting in the loop
        ax1.plot(Mx_upstream[idx],z,color='black') #plot major axis on the left
        ax1.plot(Mx_downstream[idx],z,color='black')  #plot minor axis on the right
        
        #MELT
        ax2.clear()
        ax2.plot(dTM[idx]*mts_to_cmh,z,color='black') 
        #ax2.plot(dTM[wet]*mts_to_cmh,z[wet],color='black') #plot major axis on the left    
        #ax2.plot(dPD[~wet]*mts_to_cmh,z[~wet],color='black') #plot major axis on the left
        
        #CREEP
        ax3.clear()
        ax3.plot(dC_major[idx]*mts_to_cmh,z,color='grey') #plot major axis on the left
        ax3.plot(dC_minor[idx]*mts_to_cmh,z,color='black')  #plot minor axis on the right
        
        ax4.clear()
        ax4.plot(dE_major[idx]*mts_to_cmh,z,color='grey') #plot major axis on the left
        ax4.plot(dE_minor[idx]*mts_to_cmh,z,color='black')  #plot minor axis on the right  
        
        ax5.clear()
        ax5.plot(dr_major[idx]*mts_to_cmh,z,color='grey') #plot major axis on the left
        ax5.plot(dr_minor[idx]*mts_to_cmh,z,color='black')  #plot minor axis on the right  
    
        
        ax7.plot(time[0:idx+1]/3600,hw[0:idx+1],'-',color='blue') 
        ax8.plot(time[0:idx+1]/3600,SCs[0:idx+1],'-',color='red')  
    
         
        ax9.plot(time/3600,Qin,'--',color='blue')#,label='Qin')
        ax9.plot(time[0:idx+1]/3600,Qout[0:idx+1],color='red')#,'-',color='red',label='Qout')
        ax9.legend(['Qin','Qout'],loc="upper right")
        
        
        ax2.set_title('Turbulent Melting')
        ax3.set_title('Creep Deformation')
        ax4.set_title('Elastic Deformation')
        ax5.set_title('Change in Radius')
    
        ax1.set_ylabel('z(m)')
        
        ax1.set_xlabel('(m)')
        ax2.set_xlabel('(cm/h)')
        ax3.set_xlabel('(cm/h)')
        ax4.set_xlabel('(cm/h)')
        ax5.set_xlabel('(cm/h)')
        ax7.set_ylabel('hw',color='blue')
        ax8.set_ylabel('Scs',color='red')
        ax9.set_ylabel('Qin')
        #ax9.set_ylabel('Qout')
    
        
        ax1.set_xlim([-5,5]) 
        ax2.set_xlim([-dlim,dlim])
        ax3.set_xlim([-dlim,dlim])
        ax4.set_xlim([-dlim,dlim])
        ax5.set_xlim([-dlim,dlim])
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
        
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        
        ax3.spines['top'].set_visible(False)
        ax3.spines['right'].set_visible(False)
        ax3.spines['left'].set_visible(False)
        
        ax4.spines['top'].set_visible(False)
        ax4.spines['right'].set_visible(False)
        ax4.spines['left'].set_visible(False)
        
        ax5.spines['top'].set_visible(False)
        ax5.spines['right'].set_visible(False)
        ax5.spines['left'].set_visible(False)
        
        
        ax7.spines['top'].set_visible(False)
        ax7.spines['bottom'].set_visible(False)
        
        ax8.spines['top'].set_visible(False)
        ax8.spines['bottom'].set_visible(False)
        
        ax9.spines['top'].set_visible(False)
        ax9.spines['right'].set_visible(False)
        
        ax9.spines['top'].set_visible(False)
        ax9.spines['right'].set_visible(False)
        #ax9.spines['bottom'].set_visible(False)
        
        
        
        ax2.axes.yaxis.set_visible(False)
        ax3.axes.yaxis.set_visible(False)
        ax4.axes.yaxis.set_visible(False)
        ax5.axes.yaxis.set_visible(False)
    
        
        ax7.axes.xaxis.set_visible(False)
        
        
        ax1.spines['bottom'].set_position(('zero'))
        ax2.spines['bottom'].set_position(('zero'))
        ax3.spines['bottom'].set_position(('zero'))
        ax4.spines['bottom'].set_position(('zero'))
        ax5.spines['bottom'].set_position(('zero'))
    
        
        ax1.spines['left'].set_bounds(0,H)
        ax1.spines['bottom'].set_bounds(-5,5)
        ax2.spines['bottom'].set_bounds(-dbound,dbound)
        ax3.spines['bottom'].set_bounds(-dbound,dbound)
        ax4.spines['bottom'].set_bounds(-dbound,dbound)
        ax5.spines['bottom'].set_bounds(-dbound,dbound)
    
        
        ax1.set_xticks([-4,-3,-2,-1,0,1,2,3,4]) 
        ax2.set_xticks(dtick) 
        ax3.set_xticks(dtick)
        ax4.set_xticks(dtick)
        ax5.set_xticks(dtick)
    
        ax9.set_xticks(np.round(np.linspace(1,max(time)/3600,20)))
        ax1.set_xticklabels([-4,-3,-2,-1,0,1,2,3,4]) 
        ax2.set_xticklabels(dtick)         
        ax3.set_xticklabels(dtick)         
        ax4.set_xticklabels(dtick)         
        ax5.set_xticklabels(dtick)         
    
        ax9.set_xticklabels(np.round(np.linspace(1,max(time)/3600,20)))
        
        #ax3.patch.set_facecolor('red')
        ax2.patch.set_alpha(0)
        ax2.patch.set_alpha(0)
        ax3.patch.set_alpha(0)
        ax4.patch.set_alpha(0)
        ax5.patch.set_alpha(0)
    
        
        ax1.axhspan(0, hw[idx], facecolor ='lightblue', alpha = 1,zorder=1)
        ax1.axhspan(-140, 0, facecolor ='peru', alpha = 1,zorder=1)
        ax1.fill_betweenx(z,-10,Mx_upstream[idx],color='aliceblue',zorder=2)
        ax1.fill_betweenx(z,Mx_downstream[idx],10,color='aliceblue',zorder=2)        
        ax2.fill_betweenx(z,0,dTM[idx]*mts_to_cmh,color='orangered')
    
        ax3.fill_betweenx(z,dC_minor[idx]*mts_to_cmh,0,where=dC_minor[idx]>0,facecolor=('orangered'))
        ax3.fill_betweenx(z,dC_minor[idx]*mts_to_cmh,0,where=dC_minor[idx]<0,facecolor=('lightgreen'))
        ax4.fill_betweenx(z,dE_minor[idx]*mts_to_cmh,0,where=dE_minor[idx]>0,facecolor=('orangered'))
        ax4.fill_betweenx(z,dE_minor[idx]*mts_to_cmh,0,where=dE_minor[idx]<0,facecolor=('lightgreen'))
        ax5.fill_betweenx(z,dr_minor[idx]*mts_to_cmh,0,where=dr_minor[idx]>0,facecolor=('orangered'))
        ax5.fill_betweenx(z,dr_minor[idx]*mts_to_cmh,0,where=dr_minor[idx]<0,facecolor=('lightgreen'))
        
        plt.savefig('Movies/Figure_movie_%s'%idx_save)

    
