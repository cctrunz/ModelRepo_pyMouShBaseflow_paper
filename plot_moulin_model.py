import matplotlib.pyplot as plt
import numpy as np

def plot_geom(  results, time, z, set_xlim_moulin=False,
                ax2_varname=False,
                ax3_varname=False, 
                ax4_varname=False, 
                ax5_varname=False, 
                ax6_varname=False):#,
                #color_scale=rainbow): this is not working yet
    
    '''Create horizontal subplots with any vertical data (in the z direction) that we want to input. 
    x axis is either the moulin radius or the delta radius for each physical mechanism
    - inputs for varname:  ['varname_major','varname_minor'] or 'varname'
    -- results is the dictionnary containing all the saved results from run_moulin_model.py
    -- plot the moulin radius in the first graph on the left
    -- you can enter up to 5 additional variable names of z data to plot. Name must be presesnt in the results dictionnary. 
    - (say where the major axis is, .. or show it in the figure..)'''
    
    fig = plt.figure(figsize=(15,5))
  
    #Define grid specifications 
    if ax2_varname == False:
        grid = plt.GridSpec(1, 1, wspace=0, hspace=0)
    if ax2_varname != False and ax3_varname == False:
        grid = plt.GridSpec(1, 2, wspace=0, hspace=0)
    if ax3_varname != False and ax4_varname == False:
        grid = plt.GridSpec(1, 3, wspace=0, hspace=0)
    if ax4_varname != False and ax5_varname == False:
        grid = plt.GridSpec(1, 4, wspace=0, hspace=0)
    if ax5_varname != False and ax6_varname == False:
        grid = plt.GridSpec(1, 5, wspace=0, hspace=0)        
    if ax6_varname != False:
        grid = plt.GridSpec(1, 6, wspace=0, hspace=0)
              
    #  define position of axis in the grid
    ax1 = fig.add_subplot(grid[0, 0])
    if ax2_varname !=False:
        ax2 = fig.add_subplot(grid[0, 1], sharey=ax1)
        if ax3_varname !=False:
            ax3 = fig.add_subplot(grid[0, 2], sharey=ax1)
            if ax4_varname !=False:
                ax4 = fig.add_subplot(grid[0, 3], sharey=ax1)
                if ax5_varname !=False:
                    ax5 = fig.add_subplot(grid[0, 4], sharey=ax1)
                    if ax6_varname !=False:
                        ax6 = fig.add_subplot(grid[0, 5], sharey=ax1)

    colors = [plt.cm.rainbow(i) for i in np.linspace(0, 1, len(time))] 
        
    #plot 
    for i in np.arange(0,len(time),100):
        #Plot the moulin radius
        ax1.plot(results['Mx_upstream'][i],z,color=colors[i]) #plot major axis on the left
        ax1.plot(results['Mx_downstream'][i],z,color=colors[i])  #plot minor axis on the right
        ax1.title.set_text('Moulin')
        if set_xlim_moulin != False:
            ax1.set_xlim(set_xlim_moulin)
    
        if ax2_varname !=False:
            #identify if single array or if major and minor arrays
            ax2.title.set_text(ax2_varname)
            if len(ax2_varname) == 2:
                ax2.plot(-results[ax2_varname[0]][i],z,color=colors[i]) #plot major axis on the left
                ax2.plot(results[ax2_varname[1]][i],z,color=colors[i])   #plot minor axis on the right             
            else:
                ax2.plot(results[ax2_varname][i],z,color=colors[i])               
                
        if ax3_varname !=False:
            #identify if single array or if major and minor arrays
            ax3.title.set_text(ax3_varname)
            if len(ax3_varname) == 2:
                ax3.plot(-results[ax3_varname[0]][i],z,color=colors[i]) #plot major axis on the left
                ax3.plot(results[ax3_varname[1]][i],z,color=colors[i])   #plot minor axis on the right             
            else:
                ax3.plot(results[ax3_varname][i],z,color=colors[i])            

        if ax4_varname !=False:
            #identify if single array or if major and minor arrays
            ax4.title.set_text(ax4_varname)
            if len(ax4_varname) == 2:
                ax4.plot(-results[ax4_varname[0]][i],z,color=colors[i]) #plot major axis on the left
                ax4.plot(results[ax4_varname[1]][i],z,color=colors[i])   #plot minor axis on the right             
            else:
                ax4.plot(results[ax4_varname][i],z,color=colors[i])    

        if ax5_varname !=False:
            #identify if single array or if major and minor arrays
            ax5.title.set_text(ax5_varname)
            if len(ax5_varname) == 2:
                ax5.plot(-results[ax5_varname[0]][i],z,color=colors[i]) #plot major axis on the left
                ax5.plot(results[ax5_varname[1]][i],z,color=colors[i])   #plot minor axis on the right             
            else:
                ax5.plot(results[ax5_varname][i],z,color=colors[i])        

        if ax6_varname !=False:
            #identify if single array or if major and minor arrays
            ax6.title.set_text(ax6_varname)
            if len(ax6_varname) == 2:
                ax6.plot(-results[ax6_varname[0]][i],z,color=colors[i]) #plot major axis on the left
                ax6.plot(results[ax6_varname[1]][i],z,color=colors[i])   #plot minor axis on the right             
            else:
                ax6.plot(results[ax6_varname][i],z,color=colors[i])        
       
        
def plot_2Darray(results,array2d):
    ''' plots map of results with time as the x axis and z as the y axis
    input:
        results: dictionnary 
        'varname':  2d array '''
    plt.figure()
    #prepare for coordinates
    # real_x = time
    # real_y = z

    #plot image
    plt.imshow(np.rot90(array2d))#,extent=extent)#,origin='lower'
    plt.colorbar()
    #plt.title()
    #plt.clim(clim_min,clim_max)   
    
    # plt.gca().set_xticks(range(len(real_x)))
    # plt.gca().set_yticks(range(len(real_y)))
    # plt.gca().set_xticklabels(real_x)
    # plt.gca().set_yticklabels(real_y)
    

def plot_2Darray_with_1Darray(results, time,array2d,array1d):
    ''' plots map of results with time as the x axis and z as the y axis
    input:
        results: dictionnary 
        'varname':  2d array '''
        
    fig, ax = plt.subplots(2, 1, sharex=True)
    
    
    #plot array
    
    ax.plot(time,array1d)

    #plot image
    im = ax.imshow(np.rot90(array2d))#,extent=extent)#,origin='lower'
    plt.colorbar(im)
    #plt.clim(clim_min,clim_max)   
    # real_x = time
    # real_y = z
    # plt.gca().set_xticks(range(len(real_x)))
    # plt.gca().set_yticks(range(len(real_y)))
    # plt.gca().set_xticklabels(real_x)
    # plt.gca().set_yticklabels(real_y)

def plot_1Darray_timeserie(results, time, array1d):
    plt.figure()
    plt.plot(time,results[array1d],label = array1d)
    plt.legend()