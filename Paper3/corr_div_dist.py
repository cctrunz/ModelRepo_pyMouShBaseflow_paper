import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme(style="ticks")#, font_scale=1.25)
#palette=sns.color_palette("colorblind")


name=['A','AA','B','BB','X','Y','YY','Z']
div=[1,1,9,7,16,11,20,20]
dist=[8,50,28,12,85,50,100,30]

#plt.plot(div,dist,'o')

fig, ax = plt.subplots()
ax.scatter(div, dist,marker='+',color='black')

for i, txt in enumerate(name):
    ax.annotate(txt, (div[i]+0.2, dist[i]+1))
    
ax.set_ylabel('Distance from previous moulin (m)')
ax.set_xlabel('Angle between stream and ice motion ($°$)')
ax.set_xticks([0,2,4,6,8,10,12,14,16,18,20])
ax.set_yticks([0,20,40,60,80,100])

sns.despine(offset=5,trim=True)
#sns.despine(ax=ax[0],bottom=True)
plt.tight_layout()


plt.savefig('corr_dis_div.pdf')