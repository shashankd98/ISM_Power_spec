import numpy as np 
from scipy import optimize
from astropy.io import fits
from matplotlib import pyplot as plt

# This produces a plot of spatially averaged T_B as a function of v_lsr

image=fits.open('../Data/lab_data.fit')
b=20
l=20
numv=78

# Average T_B
def TBave(size):
    tb=np.zeros(numv)
    for i in range(-size,size+1):
        for j in range(-size,size+1):
            tb+=image[0].section[:,b+i,l+j]
    tb=tb/(2*size+1)**2
    return tb

def double_gaussian(v, a1,v01,sig1,a2,v02,sig2):
    return a1*np.exp(-1/2*(v-v01)**2/sig1**2)+a2*np.exp(-1/2*(v-v02)**2/sig2**2)

def gaussian(v,a,v0,sig):
    return a*np.exp(-1/2*(v-v0)**2/sig**2)

def main():
    # Plot average T_B as a function of v_lsr
    vlsr=np.linspace(-40,40,78)
    sizes=[5,10,15,20]
    for size in sizes:
        tb=TBave(size)
        plt.plot(vlsr,tb,label=str(size)+"x"+str(size))
        print("Running...")

    # Fit the 5x5 image to a double Gaussian function
    tb5=TBave(5)
    params= optimize.curve_fit(double_gaussian,vlsr,tb5)[0]
    fit1=gaussian(vlsr,params[0],params[1],params[2])
    fit2=gaussian(vlsr,params[3],params[4],params[5])

    plt.plot(vlsr,fit1,label='WNM fit')
    plt.plot(vlsr,fit2,label='CNM fit')
    plt.xlabel("$V_{lsr}$ (km/s)")
    plt.ylabel("$T_B$ (K)")
    plt.title("LAB Survey")
    plt.legend(title="Img size")
    plt.show()

if __name__ == "__main__":
    main()
