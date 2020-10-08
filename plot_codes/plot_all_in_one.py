import matplotlib.pyplot as plt
import numpy as np

def generate_semicircle(center_x, center_y, radius, stepsize=0.01):
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


#def set_up_figure():
# fig = plt.figure(figsize=(25,8))
# grid = plt.GridSpec(11, 41, wspace=-0.7)
# ax1 = fig.add_subplot(grid[1:9, 7:12]) #moulin plot
# ax2 = fig.add_subplot(grid[1:9, 15:26], sharey=ax1) #total
# ax3 = fig.add_subplot(grid[1:9, 19:30], sharey=ax1) #melt
# ax4 = fig.add_subplot(grid[1:9, 23:34], sharey=ax1) #creep
# ax5 = fig.add_subplot(grid[1:9, 27:38], sharey=ax1) #elastic
# ax6 = fig.add_subplot(grid[1:9, 31:42], sharey=ax1) #refreezing

# ax7 = fig.add_subplot(grid[1:9, 0:1], sharey=ax1) #glen ice motion
# ax8 = fig.add_subplot(grid[1:9, 2:3], sharey=ax1) #temperature
# ax9 = fig.add_subplot(grid[1:9, 15:17], sharey=ax1)  #hw
# ax10 = fig.add_subplot(grid[0:1, 15:17]) #Qin-Qout
# ax11 = fig.add_subplot(grid[6:8, 15:16]) #S
    

def live_plot(Mx_upstream,Mx_downstream,dTM,dPD,dC_major,dC_minor,dE_major,dE_minor,dr_major,dr_minor,dGlen,t,hw,SCs,z,Qin,Qout,idx,wet,mts_to_cmh,time,H,results,T_far):

    
    fig = plt.figure(figsize=(25,8))
    grid = plt.GridSpec(11, 41, wspace=-0.7)
    ax1 = fig.add_subplot(grid[1:9, 7:12]) #moulin plot
    ax2 = fig.add_subplot(grid[1:9, 15:26], sharey=ax1) #total
    ax3 = fig.add_subplot(grid[1:9, 19:30], sharey=ax1) #melt
    ax4 = fig.add_subplot(grid[1:9, 23:34], sharey=ax1) #creep
    ax5 = fig.add_subplot(grid[1:9, 27:38], sharey=ax1) #elastic
    ax6 = fig.add_subplot(grid[1:9, 31:42], sharey=ax1) #refreezing
    
    ax7 = fig.add_subplot(grid[1:9, 0:1], sharey=ax1) #glen ice motion
    ax8 = fig.add_subplot(grid[1:9, 2:3], sharey=ax1) #temperature
    ax9 = fig.add_subplot(grid[1:9, 15:17], sharey=ax1)  #hw
    ax10 = fig.add_subplot(grid[0:1, 15:17]) #Qin-Qout
    ax11 = fig.add_subplot(grid[6:8, 15:16]) #S
    
    '''clear figures'''    
    ax1.clear() #clear figure content -- especially important for ploting in the loop
    ax2.clear()
    ax3.clear()
    ax4.clear()
    ax5.clear()
    ax6.clear()
    ax7.clear()
    ax8.clear()
    ax9.clear()
    ax10.clear()
    ax11.clear()
    
    '''plot'''
    #moulin
    ax1.plot(Mx_upstream,z,color='black') #plot major axis on the left
    ax1.plot(Mx_downstream,z,color='black')  #plot minor axis on the right
    #total
    ax2.plot(dr_major*mts_to_cmh,z,color='grey') #plot major axis on the left
    ax2.plot(dr_minor*mts_to_cmh,z,color='black')  #plot minor axis on the right      
    #melt    
    ax3.plot(dTM[wet]*mts_to_cmh,z[wet],color='black') #plot major axis on the left
    ax3.plot(dPD[~wet]*mts_to_cmh,z[~wet],color='black') #plot major axis on the left
    #creep    
    ax4.plot(dC_major*mts_to_cmh,z,color='grey') #plot major axis on the left
    ax4.plot(dC_minor*mts_to_cmh,z,color='black')  #plot minor axis on the right
    #elastic 
    ax5.plot(dE_major*mts_to_cmh,z,color='grey') #plot major axis on the left
    ax5.plot(dE_minor*mts_to_cmh,z,color='black')  #plot minor axis on the right      
    #refreezing   
    ax6.plot(np.zeros(len(z)),z,color='black') #plot major axis on the left
    #glen-ice motion
    ax7.plot(dGlen*mts_to_cmh,z,color='black') #plot major axis on the left
    #temperature
    ax8.plot(T_far,z,'-',color='black')
    #hw    
    ax9.plot(time[0:idx]/3600,results['hw'][0:idx],'-',color='blue')
    #Qin-Qout     
    ax10.plot(time[0:idx]/3600,Qin[0:idx],'-',color='blue')
    ax10.plot(time[0:idx]/3600,results['Qout'][0:idx],'-',color='red')
    #
    # if SCs < 0.1:
    #     radius = 0.1
    # else:
    #     radius = np.sqrt(np.abs(SCs))
    
    
    radius = np.sqrt(SCs+0.01/np.pi)
    x,y = generate_semicircle(0, 0, radius, stepsize=radius*1e-3)
    ax11.plot(x,y,'-',color='black') 
    ax11.plot(x,np.zeros(len(x)),'-',color='black') 
    
    # if idx < 6:
    #     radius = np.sqrt(SCs+0.01/np.pi)
    #     x,y = generate_semicircle(0, 0, radius, stepsize=radius*1e-3)
    #     ax11.plot(x,y,'-',color='black') 
    #     ax11.plot(x,np.zeros(len(x)),'-',color='black') 
    # else:
    #     colors = [1,0.9,0.8,0.7,0.6]
    #     for i in [0,1,2,3,4,5]:
    #         radius = np.sqrt(results['SCs'][idx-i]+0.01/np.pi)
    #         x,y = generate_semicircle(0, 0, radius, stepsize=radius*1e-3)
    #         ax11.plot(x,y,'-',color='black') 
    #         ax11.plot(x,np.zeros(len(x)),'-',color='black') 
            

    
    '''fill'''
    ax1.axhspan(0, hw, facecolor ='lightblue', alpha = 1,zorder=1)
    ax1.axhspan(-140, 0, facecolor ='peru', alpha = 1,zorder=1)
    ax1.fill_betweenx(z,-10,Mx_upstream,color='aliceblue',zorder=2)
    ax1.fill_betweenx(z,Mx_downstream,10,color='aliceblue',zorder=2)     
    #subglacial stream
    ax1.fill_between([Mx_downstream[0],10],0,0.015*H, color='lightblue',zorder=3)
    ax1.plot([Mx_downstream[0],10],[0,0],color='black',zorder=3)
    ax1.plot([Mx_downstream[0],10],[0.015*H,0.015*H],color='black',zorder=3)
    #supraglacial stream
    ax1.fill_between([-10,Mx_upstream[-1]],H,0.985*H, color='lightblue',zorder=3)
    ax1.plot([-10,Mx_upstream[-1]],[0.985*H,0.985*H],color='black',zorder=3)
    
    ax2.fill_betweenx(z,dr_minor*mts_to_cmh,0,where=dr_minor>0,facecolor=('orangered'))
    ax2.fill_betweenx(z,dr_minor*mts_to_cmh,0,where=dr_minor<0,facecolor=('lightgreen'))       
    ax3.fill_betweenx(z[wet],0,dTM[wet]*mts_to_cmh,color='orangered')
    ax3.fill_betweenx(z[~wet],0,dPD[~wet]*mts_to_cmh,color='orangered')
    ax4.fill_betweenx(z,dC_minor*mts_to_cmh,0,where=dC_minor>0,facecolor=('orangered'))
    ax4.fill_betweenx(z,dC_minor*mts_to_cmh,0,where=dC_minor<0,facecolor=('lightgreen'))
    ax5.fill_betweenx(z,dE_minor*mts_to_cmh,0,where=dE_minor>0,facecolor=('orangered'))
    ax5.fill_betweenx(z,dE_minor*mts_to_cmh,0,where=dE_minor<0,facecolor=('lightgreen'))   
    
    #ax6.fill_betweenx(z,0,*mts_to_cmh,color='grey')
    #ax7.fill_betweenx(z,0,dGlen*mts_to_cmh,color='grey')
    
    '''Title'''
    ax2.set_title('Total')
    ax3.set_title('Turbulent Melting')
    ax4.set_title('Creep Deformation')
    ax5.set_title('Elastic Deformation')
    ax6.set_title('Refreezing')
    ax7.set_title('Ice Motion')
    ax8.set_title('Temperature')
    #ax9.set_title('Head')
    ax10.set_title('Qin and Qout')
    #ax11.set_title('SCs')$
    
    #ax8.text(0, 0, 'Temperature')
    
    '''label'''
    ax1.set_ylabel('z(m)')    
    ax1.set_xlabel('r(m)')
    ax2.set_xlabel('(cm/h)')
    ax3.set_xlabel('(cm/h)')
    ax4.set_xlabel('(cm/h)')
    ax5.set_xlabel('(cm/h)')
    ax6.set_xlabel('(cm/h)')
    ax7.set_xlabel('(cm/h)')
    ax8.set_xlabel('(K)')
    ax11.set_xlabel('r(m)')
    ax10.set_ylabel('($m^3/s$)')


    
    # ax7.yaxis.tick_left()
    # ax7.yaxis.set_label_position("left")
    # ax7.tick_params(axis='y', labelcolor='blue')
    
    # ax8.yaxis.tick_right()
    # ax8.yaxis.set_label_position("right")
    # ax8.tick_params(axis='y', labelcolor='red')
    
    
    '''spine visibility'''
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
    ax7.spines['right'].set_visible(False)
    ax7.spines['left'].set_visible(False)
    
    ax8.spines['top'].set_visible(False)
    ax8.spines['right'].set_visible(False)
    ax8.spines['left'].set_visible(False)
    
    ax9.spines['top'].set_visible(False)
    ax9.spines['right'].set_visible(False)
    ax9.spines['left'].set_visible(False)
    
    ax10.spines['top'].set_visible(False)
    ax10.spines['bottom'].set_visible(False)
    ax10.spines['right'].set_visible(False)
    
    ax11.spines['top'].set_visible(False)
    ax11.spines['right'].set_visible(False)
    ax11.spines['left'].set_visible(False)
    
    
    
    
    ax2.axes.yaxis.set_visible(False)
    ax3.axes.yaxis.set_visible(False)
    ax4.axes.yaxis.set_visible(False)
    ax5.axes.yaxis.set_visible(False)
    ax6.axes.yaxis.set_visible(False)
    ax7.axes.yaxis.set_visible(False)
    ax8.axes.yaxis.set_visible(False)
    ax9.axes.yaxis.set_visible(False)
    ax11.axes.yaxis.set_visible(False)
    
    ax10.axes.xaxis.set_visible(False)
    
    '''spine position'''
    ax1.spines['bottom'].set_position(('zero'))
    ax2.spines['bottom'].set_position(('zero'))
    ax3.spines['bottom'].set_position(('zero'))
    ax4.spines['bottom'].set_position(('zero'))
    ax5.spines['bottom'].set_position(('zero'))
    ax6.spines['bottom'].set_position(('zero'))
    ax7.spines['bottom'].set_position(('zero'))
    ax8.spines['bottom'].set_position(('zero'))
    ax9.spines['bottom'].set_position(('data', H))


    '''spine bounds'''
    ax1.spines['left'].set_bounds(0,H)
    ax11.spines['left'].set_bounds(0,5)
    ax1.spines['bottom'].set_bounds(-4,4)
    #ax10.spines['left'].set_bounds()
    ax1.spines['bottom'].set_bounds(-3,3)
    ax2.spines['bottom'].set_bounds(-1,1)
    ax3.spines['bottom'].set_bounds(-1,1)
    ax4.spines['bottom'].set_bounds(-1,1)
    ax5.spines['bottom'].set_bounds(-1,1)
    ax6.spines['bottom'].set_bounds(-1,1)
    ax7.spines['bottom'].set_bounds(-1,1)
    ax8.spines['bottom'].set_bounds(min(T_far)-2,max(T_far)+2)
    ax11.spines['bottom'].set_bounds(0,2)

    '''xlim''' 
    ax1.set_xlim([-10,10]) 
    ax2.set_xlim([-7,7])
    ax3.set_xlim([-7,7])
    ax4.set_xlim([-7,7])
    ax5.set_xlim([-7,7])
    ax6.set_xlim([-7,7])
    ax7.set_xlim([-2,2])
    ax8.set_xlim([min(T_far)-2,max(T_far)+2]) #Temperature
    ax9.set_xlim([0,max(time)/3600])#hw
    ax10.set_xlim([0,max(time)/3600])#Qin
    ax11.set_xlim([-2,2]) #SCs
    
    '''ylim'''
    ax11.set_ylim([0,2])
    ax10.set_ylim([min(Qin),max(Qin)])
    
    '''ticks position'''
    ax1.set_xticks([-3,-2,-1,0,1,2,3]) 
    ax2.set_xticks([-1,0,1]) 
    ax3.set_xticks([-1,0,1])
    ax4.set_xticks([-1,0,1])
    ax5.set_xticks([-1,0,1])
    ax6.set_xticks([-1,0,1])
    ax7.set_xticks([-1,0,1])
    ax8.set_xticks(np.round(np.linspace(min(T_far),max(T_far),3),0))
    ax10.set_xticks(np.round(np.linspace(1,max(time)/3600,20),0))
    ax11.set_xticks([0,1,2])

    '''ticks values'''
    ax1.set_xticklabels([-3,-2,-1,0,1,2,3]) 
    ax2.set_xticklabels([-1,0,1])         
    ax3.set_xticklabels([-1,0,1])         
    ax4.set_xticklabels([-1,0,1])         
    ax5.set_xticklabels([-1,0,1])         
    ax6.set_xticklabels([-1,0,1])
    ax6.set_xticklabels([-1,0,1])
    ax7.set_xticklabels([-1,0,1])
    ax8.set_xticklabels(np.round(np.linspace(min(T_far),max(T_far),3),0))
    ax10.set_xticklabels(np.round(np.linspace(1,max(time)/3600,20),0))
    ax11.set_xticklabels([0,1,2])
    
    '''Background transparency'''
    ax2.patch.set_alpha(0)
    ax3.patch.set_alpha(0)
    ax4.patch.set_alpha(0)
    ax5.patch.set_alpha(0)
    ax6.patch.set_alpha(0)
    ax7.patch.set_alpha(0)
    ax8.patch.set_alpha(0)
    ax9.patch.set_alpha(0)
    ax10.patch.set_alpha(0)
    ax11.patch.set_alpha(0)
    

    
    

