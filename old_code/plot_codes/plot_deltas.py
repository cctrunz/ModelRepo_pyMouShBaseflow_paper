import matplotlib.pyplot as plt
import numpy as np

fig1 = plt.figure(figsize=(20,5))
grid = plt.GridSpec(1, 7, wspace=0)
ax1 = fig1.add_subplot(grid[0, 0:3])
ax2 = fig1.add_subplot(grid[0, 1:4], sharey=ax1)
ax3 = fig1.add_subplot(grid[0, 2:5], sharey=ax1)
ax4 = fig1.add_subplot(grid[0, 3:6], sharey=ax1)
ax5 = fig1.add_subplot(grid[0, 4:7], sharey=ax1)



def live_plot(dTM,dPD,dC_major,dC_minor,dE_major,dE_minor,dr_major,dr_minor,dGlen,z,wet,mts_to_cmh):

              
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
    
    ax1.clear()
    ax1.plot(dGlen*mts_to_cmh,z,color='black') #plot major axis on the left
    

    
    ax2.set_title('Turbulent Melting')
    ax3.set_title('Creep Deformation')
    ax4.set_title('Elastic Deformation')
    ax5.set_title('Change in Radius')
    ax1.set_title('Ice Motion')
    

    

    ax2.set_xlabel('(cm/h)')
    ax3.set_xlabel('(cm/h)')
    ax4.set_xlabel('(cm/h)')
    ax5.set_xlabel('(cm/h)')
    ax1.set_xlabel('(cm/h)')
    
    ax1.set_xlim([-10,10])     
    ax2.set_xlim([-5,5])
    ax3.set_xlim([-5,5])
    ax4.set_xlim([-5,5])
    ax5.set_xlim([-5,5])
    ax1.set_xlim([-5,5])
   

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

    
  
    ax2.axes.yaxis.set_visible(False)
    ax3.axes.yaxis.set_visible(False)
    ax4.axes.yaxis.set_visible(False)
    ax5.axes.yaxis.set_visible(False)
    ax1.axes.yaxis.set_visible(False)
    

    ax2.spines['bottom'].set_position(('zero'))
    ax3.spines['bottom'].set_position(('zero'))
    ax4.spines['bottom'].set_position(('zero'))
    ax5.spines['bottom'].set_position(('zero'))
    ax1.spines['bottom'].set_position(('zero'))
    

    ax2.spines['bottom'].set_bounds(-1,1)
    ax3.spines['bottom'].set_bounds(-1,1)
    ax4.spines['bottom'].set_bounds(-1,1)
    ax5.spines['bottom'].set_bounds(-1,1)
    ax1.spines['bottom'].set_bounds(-1,1)
    
    

    ax2.set_xticks([-1,0,1]) 
    ax3.set_xticks([-1,0,1])
    ax4.set_xticks([-1,0,1])
    ax5.set_xticks([-1,0,1])
    ax1.set_xticks([-1,0,1])


    ax2.set_xticklabels([-1,0,1])         
    ax3.set_xticklabels([-1,0,1])         
    ax4.set_xticklabels([-1,0,1])         
    ax5.set_xticklabels([-1,0,1])         
    ax1.set_xticklabels([-1,0,1])

    
    #ax3.patch.set_facecolor('red')
    ax2.patch.set_alpha(0)
    ax2.patch.set_alpha(0)
    ax3.patch.set_alpha(0)
    ax4.patch.set_alpha(0)
    ax5.patch.set_alpha(0)
    ax1.patch.set_alpha(0)
    
   
    ax2.fill_betweenx(z[wet],0,dTM[wet]*mts_to_cmh,color='orangered')
    ax2.fill_betweenx(z[~wet],0,dPD[~wet]*mts_to_cmh,color='orangered')
    ax3.fill_betweenx(z,dC_minor*mts_to_cmh,0,where=dC_minor>0,facecolor=('orangered'))
    ax3.fill_betweenx(z,dC_minor*mts_to_cmh,0,where=dC_minor<0,facecolor=('lightgreen'))
    ax4.fill_betweenx(z,dE_minor*mts_to_cmh,0,where=dE_minor>0,facecolor=('orangered'))
    ax4.fill_betweenx(z,dE_minor*mts_to_cmh,0,where=dE_minor<0,facecolor=('lightgreen'))
    ax5.fill_betweenx(z,dr_minor*mts_to_cmh,0,where=dr_minor>0,facecolor=('orangered'))
    ax5.fill_betweenx(z,dr_minor*mts_to_cmh,0,where=dr_minor<0,facecolor=('lightgreen'))
    ax1.fill_betweenx(z,0,dGlen*mts_to_cmh,color='grey')
    
    
