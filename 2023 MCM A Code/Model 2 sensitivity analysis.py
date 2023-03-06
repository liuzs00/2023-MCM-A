import numpy as np
import matplotlib.pyplot as plt   


high1=[13,13,17]
mean1=[65/6,69/6,84/6]
low1=[9,10,12]

high2=[155.539,153.779,153.302]
mean2=[120.555,102.63,102.245]
low2=[119.830,128.375,144.753]

x=["90%N*","N*","110%N*"]






fig,axs = plt.subplots(1,2,figsize=(10,10))
axs[0].plot(x,high1,"navy", label="High")
axs[0].plot(x,mean1,"black", label="Mean")
axs[0].plot(x,low1,"blue", label="Low")
axs[0].fill_between(x, high1,mean1,color='darkorange',alpha=.5)
axs[0].fill_between(x, mean1,low1,color='orange',alpha=.5)
axs[0].set_ylabel('Total Species' )
axs[0].set_xlabel('N*' )
axs[0].set_title('N* VS The number of Species')
axs[0].legend()
axs[0].set_ylim(0,20)

axs[1].plot(x,high2,"navy", label="High")
axs[1].plot(x,mean2,"black", label="Mean")
axs[1].plot(x,low2,"blue", label="Low")
axs[1].fill_between(x, high2,mean2,color='darkorange',alpha=.5)
axs[1].fill_between(x, mean2,low2,color='orange',alpha=.5)
axs[1].set_ylabel('Total Utility' )
axs[1].set_xlabel('N*' )
axs[1].set_title('N* VS Optimized Z')
axs[1].legend()

plt.show()