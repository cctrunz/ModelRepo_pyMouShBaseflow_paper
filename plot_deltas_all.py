import matplotlib.pyplot as plt
import numpy as np

fig1 = plt.figure(figsize=(20,10))
grid = plt.GridSpec(2, 7, wspace=0)
ax1 = fig1.add_subplot(grid[0, 0:3])
ax2 = fig1.add_subplot(grid[0, 1:4], sharey=ax1)
ax3 = fig1.add_subplot(grid[0, 2:5], sharey=ax1)
ax4 = fig1.add_subplot(grid[0, 3:6], sharey=ax1)
ax5 = fig1.add_subplot(grid[0, 4:7], sharey=ax1)

ax6 = fig1.add_subplot(grid[1, 0:3])
ax7 = fig1.add_subplot(grid[1, 1:4], sharey=ax6)
ax8 = fig1.add_subplot(grid[1, 2:5], sharey=ax6)
ax9 = fig1.add_subplot(grid[1, 3:6], sharey=ax6)
ax10 = fig1.add_subplot(grid[1, 4:7], sharey=ax6)




def live_plot(dTM,dPD,dOC,dC_major,dC_minor,dE_major,dE_minor,dr_major,dr_minor,dGlen,z,mts_to_cmh):

    ax1.clear()
    ax1.plot(dTM*mts_to_cmh,z,color='black') #plot major axis on the left
    ax6.clear()
    ax6.plot(dTM*mts_to_cmh,z,color='black') #plot major axis on the left

              
    ax2.clear()
    ax2.plot(dOC*mts_to_cmh,z,color='black') #plot major axis on the left
    ax7.clear()
    ax7.plot(dPD*mts_to_cmh,z,color='black') #plot major axis on the left

    
    ax3.clear()
    ax3.plot(dC_major*mts_to_cmh,z,color='black') #plot major axis on the left
    ax8.clear()
    ax8.plot(dC_minor*mts_to_cmh,z,color='black')  #plot minor axis on the right
    
    ax4.clear()
    ax4.plot(dE_major*mts_to_cmh,z,color='black') #plot major axis on the left
    ax9.clear()
    ax9.plot(dE_minor*mts_to_cmh,z,color='black')  #plot minor axis on the right  
    
    ax5.clear()
    ax5.plot(dr_major*mts_to_cmh,z,color='black') #plot major axis on the left
    ax10.clear()
    ax10.plot(dr_minor*mts_to_cmh,z,color='black')  #plot minor axis on the right  

    ax1.fill_betweenx(z,0,dTM*mts_to_cmh,color='orangered')
    ax6.fill_betweenx(z,0,dTM*mts_to_cmh,color='orangered')

    ax2.fill_betweenx(z,0,dOC*mts_to_cmh,color='orangered')
    ax7.fill_betweenx(z,0,dPD*mts_to_cmh,color='orangered')

    ax3.fill_betweenx(z,dC_major*mts_to_cmh,0,where=dC_major>0,facecolor=('orangered'))
    ax3.fill_betweenx(z,dC_major*mts_to_cmh,0,where=dC_major<0,facecolor=('lightgreen'))
    ax4.fill_betweenx(z,dE_major*mts_to_cmh,0,where=dE_major>0,facecolor=('orangered'))
    ax4.fill_betweenx(z,dE_major*mts_to_cmh,0,where=dE_major<0,facecolor=('lightgreen'))
    ax5.fill_betweenx(z,dr_major*mts_to_cmh,0,where=dr_major>0,facecolor=('orangered'))
    ax5.fill_betweenx(z,dr_major*mts_to_cmh,0,where=dr_major<0,facecolor=('lightgreen'))

    ax8.fill_betweenx(z,dC_minor*mts_to_cmh,0,where=dC_minor>0,facecolor=('orangered'))
    ax8.fill_betweenx(z,dC_minor*mts_to_cmh,0,where=dC_minor<0,facecolor=('lightgreen'))
    ax9.fill_betweenx(z,dE_minor*mts_to_cmh,0,where=dE_minor>0,facecolor=('orangered'))
    ax9.fill_betweenx(z,dE_minor*mts_to_cmh,0,where=dE_minor<0,facecolor=('lightgreen'))
    ax10.fill_betweenx(z,dr_minor*mts_to_cmh,0,where=dr_minor>0,facecolor=('orangered'))
    ax10.fill_betweenx(z,dr_minor*mts_to_cmh,0,where=dr_minor<0,facecolor=('lightgreen'))
    
    

    

    ax1.set_title('Turbulent Melting')
    ax2.set_title('Open Channel')
    ax3.set_title('Creep Deformation major')
    ax4.set_title('Elastic Deformation major')
    ax5.set_title('Change in Radius major')
    ax6.set_title('Turbulent Melting')
    ax7.set_title('Potential drop')
    ax8.set_title('Creep Deformation minor')
    ax9.set_title('Elastic Deformation minor')
    ax10.set_title('Change in Radius minor')
    
    

    
    ax1.set_xlabel('(cm/h)')
    ax2.set_xlabel('(cm/h)')
    ax3.set_xlabel('(cm/h)')
    ax4.set_xlabel('(cm/h)')
    ax5.set_xlabel('(cm/h)')
    ax6.set_xlabel('(cm/h)')
    ax7.set_xlabel('(cm/h)')
    ax8.set_xlabel('(cm/h)')
    ax9.set_xlabel('(cm/h)')
    ax10.set_xlabel('(cm/h)')
    
    
    ax1.set_xlim([-5,5])   
    ax2.set_xlim([-5,5])
    ax3.set_xlim([-5,5])
    ax4.set_xlim([-5,5])
    ax5.set_xlim([-5,5])
    ax6.set_xlim([-5,5])   
    ax7.set_xlim([-5,5])
    ax8.set_xlim([-5,5])
    ax9.set_xlim([-5,5])
    ax10.set_xlim([-5,5])
    
   

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
    
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_visible(False)

    ax6.spines['top'].set_visible(False)
    ax6.spines['right'].set_visible(False)
    ax6.spines['left'].set_visible(False)

    ax7.spines['top'].set_visible(False)
    ax7.spines['right'].set_visible(False)
    ax7.spines['left'].set_visible(False)

    ax8.spines['top'].set_visible(False)
    ax8.spines['right'].set_visible(False)
    ax8.spines['left'].set_visible(False)

    ax9.spines['top'].set_visible(False)
    ax9.spines['right'].set_visible(False)
    ax9.spines['left'].set_visible(False)

    ax10.spines['top'].set_visible(False)
    ax10.spines['right'].set_visible(False)
    ax10.spines['left'].set_visible(False)

    
    ax1.axes.yaxis.set_visible(False)
    ax2.axes.yaxis.set_visible(False)
    ax3.axes.yaxis.set_visible(False)
    ax4.axes.yaxis.set_visible(False)
    ax5.axes.yaxis.set_visible(False)
    ax6.axes.yaxis.set_visible(False)
    ax7.axes.yaxis.set_visible(False)
    ax8.axes.yaxis.set_visible(False)
    ax9.axes.yaxis.set_visible(False)
    ax10.axes.yaxis.set_visible(False)
    
    
    ax1.spines['bottom'].set_position(('zero'))
    ax2.spines['bottom'].set_position(('zero'))
    ax3.spines['bottom'].set_position(('zero'))
    ax4.spines['bottom'].set_position(('zero'))
    ax5.spines['bottom'].set_position(('zero'))
    ax7.spines['bottom'].set_position(('zero'))
    ax8.spines['bottom'].set_position(('zero'))
    ax9.spines['bottom'].set_position(('zero'))
    ax10.spines['bottom'].set_position(('zero'))

    
    ax1.spines['bottom'].set_bounds(-1,1)
    ax2.spines['bottom'].set_bounds(-1,1)
    ax3.spines['bottom'].set_bounds(-1,1)
    ax4.spines['bottom'].set_bounds(-1,1)
    ax5.spines['bottom'].set_bounds(-1,1)   
    ax6.spines['bottom'].set_bounds(-1,1)
    ax7.spines['bottom'].set_bounds(-1,1)
    ax8.spines['bottom'].set_bounds(-1,1)
    ax9.spines['bottom'].set_bounds(-1,1)
    ax10.spines['bottom'].set_bounds(-1,1)
    
    
    ax1.set_xticks([-1,0,1])
    ax2.set_xticks([-1,0,1]) 
    ax3.set_xticks([-1,0,1])
    ax4.set_xticks([-1,0,1])
    ax5.set_xticks([-1,0,1])
    ax6.set_xticks([-1,0,1])
    ax7.set_xticks([-1,0,1]) 
    ax8.set_xticks([-1,0,1])
    ax9.set_xticks([-1,0,1])
    ax10.set_xticks([-1,0,1])
    

    ax1.set_xticklabels([-1,0,1])
    ax2.set_xticklabels([-1,0,1])         
    ax3.set_xticklabels([-1,0,1])         
    ax4.set_xticklabels([-1,0,1])         
    ax5.set_xticklabels([-1,0,1])  
    ax6.set_xticklabels([-1,0,1])
    ax7.set_xticklabels([-1,0,1])         
    ax8.set_xticklabels([-1,0,1])         
    ax9.set_xticklabels([-1,0,1])         
    ax10.set_xticklabels([-1,0,1])        
    

    
    #ax3.patch.set_facecolor('red')
    ax1.patch.set_alpha(0)
    ax2.patch.set_alpha(0)
    ax3.patch.set_alpha(0)
    ax4.patch.set_alpha(0)
    ax5.patch.set_alpha(0)
    ax6.patch.set_alpha(0)
    ax7.patch.set_alpha(0)
    ax8.patch.set_alpha(0)
    ax9.patch.set_alpha(0)
    ax10.patch.set_alpha(0)

    
    

    
    
