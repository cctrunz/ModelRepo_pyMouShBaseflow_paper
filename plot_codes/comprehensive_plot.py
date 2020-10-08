import matplotlib.pyplot as plt
import numpy as np

fig1 = plt.figure(figsize=(20,15))
grid = plt.GridSpec(4, 16, wspace=-0.7)
ax1 = fig1.add_subplot(grid[0:2, 0:4])
ax2 = fig1.add_subplot(grid[0:2, 5:9], sharey=ax1)
ax3 = fig1.add_subplot(grid[0:2, 7:11], sharey=ax1)
ax4 = fig1.add_subplot(grid[0:2, 9:13], sharey=ax1)
ax5 = fig1.add_subplot(grid[0:2, 11:15], sharey=ax1)
ax6 = fig1.add_subplot(grid[0:2, 13:17], sharey=ax1)
ax7 = fig1.add_subplot(grid[2, 0:16]) 
ax8 = ax7.twinx()
ax9 = fig1.add_subplot(grid[3, 0:16])  
ax10 = fig1.add_subplot(grid[3, 0:16]) 
    

def live_plot(Mx_upstream,Mx_downstream,dTM,dPD,dC_major,dC_minor,dE_major,dE_minor,dr_major,dr_minor,dGlen,t,hw,SCs,z,Qin,Qout,idx,wet,mts_to_cmh,time,H):

    ax1.clear() #clear figure content -- especially important for ploting in the loop
    ax1.plot(Mx_upstream,z,color='black') #plot major axis on the left
    ax1.plot(Mx_downstream,z,color='black')  #plot minor axis on the right
              
    ax2.clear()
    ax2.plot(dTM[wet]*mts_to_cmh,z[wet],color='black') #plot major axis on the left
    ax2.plot(dPD[~wet]*mts_to_cmh,z[~wet],color='black') #plot major axis on the left
    
    ax3.clear()
    ax3.plot(dC_major*mts_to_cmh,z,color='grey') #plot major axis on the left
    ax3.plot(dC_minor*mts_to_cmh,z,color='black')  #plot minor axis on the right
    
    ax4.clear()
    ax4.plot(dE_major*mts_to_cmh,z,color='grey') #plot major axis on the left
    ax4.plot(dE_minor*mts_to_cmh,z,color='black')  #plot minor axis on the right  
    
    ax5.clear()
    ax5.plot(dr_major*mts_to_cmh,z,color='grey') #plot major axis on the left
    ax5.plot(dr_minor*mts_to_cmh,z,color='black')  #plot minor axis on the right  
    
    ax6.clear()
    ax6.plot(dGlen*mts_to_cmh,z,color='black') #plot major axis on the left
    
    ax7.plot(t/3600,hw,'.',color='blue')
    
    ax8.plot(t/3600,SCs,'.',color='red')       
    ax9.plot(t/3600,Qin[idx],'.',color='blue')
    ax10.plot(t/3600,Qout,'.',color='red')
    
    
    ax2.set_title('Turbulent Melting')
    ax3.set_title('Creep Deformation')
    ax4.set_title('Elastic Deformation')
    ax5.set_title('Change in Radius')
    ax6.set_title('Ice Motion')
    
    ax1.set_ylabel('z(m)')
    
    ax1.set_xlabel('(m)')
    ax2.set_xlabel('(cm/h)')
    ax3.set_xlabel('(cm/h)')
    ax4.set_xlabel('(cm/h)')
    ax5.set_xlabel('(cm/h)')
    ax6.set_xlabel('(cm/h)')
    
    ax1.set_xlim([-10,10]) 
    ax2.set_xlim([-5,5])
    ax3.set_xlim([-5,5])
    ax4.set_xlim([-5,5])
    ax5.set_xlim([-5,5])
    ax6.set_xlim([-5,5])
    ax7.set_xlim([0,max(time)/3600])
    ax8.set_xlim([0,max(time)/3600])
    ax9.set_xlim([0,max(time)/3600])
    ax10.set_xlim([0,max(time/3600)])
    
    ax7.set_ylim([0,H])
    #ax8.set_ylim([0,3])
    #ax9.set_ylim([0,max(Qin)])
    #ax10.set_ylim([0,max(Qin)])
    
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
    
    ax6.spines['top'].set_visible(False)
    ax6.spines['right'].set_visible(False)
    ax6.spines['left'].set_visible(False)
    
    ax7.spines['top'].set_visible(False)
    ax7.spines['bottom'].set_visible(False)
    
    ax8.spines['top'].set_visible(False)
    ax8.spines['bottom'].set_visible(False)
    
    ax9.spines['top'].set_visible(False)
    ax9.spines['right'].set_visible(False)
    
    ax10.spines['top'].set_visible(False)
    ax10.spines['right'].set_visible(False)
    #ax10.spines['bottom'].set_visible(False)
    
    
    
    ax2.axes.yaxis.set_visible(False)
    ax3.axes.yaxis.set_visible(False)
    ax4.axes.yaxis.set_visible(False)
    ax5.axes.yaxis.set_visible(False)
    ax6.axes.yaxis.set_visible(False)
    
    ax7.axes.xaxis.set_visible(False)
    
    
    ax1.spines['bottom'].set_position(('zero'))
    ax2.spines['bottom'].set_position(('zero'))
    ax3.spines['bottom'].set_position(('zero'))
    ax4.spines['bottom'].set_position(('zero'))
    ax5.spines['bottom'].set_position(('zero'))
    ax6.spines['bottom'].set_position(('zero'))
    
    ax1.spines['left'].set_bounds(0,H)
    ax1.spines['bottom'].set_bounds(-10,10)
    ax2.spines['bottom'].set_bounds(-1,1)
    ax3.spines['bottom'].set_bounds(-1,1)
    ax4.spines['bottom'].set_bounds(-1,1)
    ax5.spines['bottom'].set_bounds(-1,1)
    ax6.spines['bottom'].set_bounds(-1,1)
    
    
    ax1.set_xticks([-8,-6,-4,-2,0,2,4,6,8]) 
    ax2.set_xticks([-1,0,1]) 
    ax3.set_xticks([-1,0,1])
    ax4.set_xticks([-1,0,1])
    ax5.set_xticks([-1,0,1])
    ax6.set_xticks([-1,0,1])
    ax10.set_xticks(np.round(np.linspace(1,max(time)/3600,20)))
    ax1.set_xticklabels([8,6,4,2,0,2,4,6,8]) 
    ax2.set_xticklabels([-1,0,1])         
    ax3.set_xticklabels([-1,0,1])         
    ax4.set_xticklabels([-1,0,1])         
    ax5.set_xticklabels([-1,0,1])         
    ax6.set_xticklabels([-1,0,1])
    ax10.set_xticklabels(np.round(np.linspace(1,max(time)/3600,20)))
    
    #ax3.patch.set_facecolor('red')
    ax2.patch.set_alpha(0)
    ax2.patch.set_alpha(0)
    ax3.patch.set_alpha(0)
    ax4.patch.set_alpha(0)
    ax5.patch.set_alpha(0)
    ax6.patch.set_alpha(0)
    
    ax1.axhspan(0, hw, facecolor ='lightblue', alpha = 1,zorder=1)
    ax1.axhspan(-140, 0, facecolor ='peru', alpha = 1,zorder=1)
    ax1.fill_betweenx(z,-10,Mx_upstream,color='aliceblue',zorder=2)
    ax1.fill_betweenx(z,Mx_downstream,10,color='aliceblue',zorder=2)        
    ax2.fill_betweenx(z[wet],0,dTM[wet]*mts_to_cmh,color='orangered')
    ax2.fill_betweenx(z[~wet],0,dPD[~wet]*mts_to_cmh,color='orangered')
    ax3.fill_betweenx(z,dC_minor*mts_to_cmh,0,where=dC_minor>0,facecolor=('orangered'))
    ax3.fill_betweenx(z,dC_minor*mts_to_cmh,0,where=dC_minor<0,facecolor=('lightgreen'))
    ax4.fill_betweenx(z,dE_minor*mts_to_cmh,0,where=dE_minor>0,facecolor=('orangered'))
    ax4.fill_betweenx(z,dE_minor*mts_to_cmh,0,where=dE_minor<0,facecolor=('lightgreen'))
    ax5.fill_betweenx(z,dr_minor*mts_to_cmh,0,where=dr_minor>0,facecolor=('orangered'))
    ax5.fill_betweenx(z,dr_minor*mts_to_cmh,0,where=dr_minor<0,facecolor=('lightgreen'))
    ax6.fill_betweenx(z,0,dGlen*mts_to_cmh,color='grey')
    
    
