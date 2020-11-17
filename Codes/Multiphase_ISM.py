import numpy as np 
from astropy.io import fits
from matplotlib import pyplot as plt
from scipy import stats
from scipy.signal import correlate2d
from scipy import fftpack
import azimuthal_average

# Parameters of the FITS file, obtained from the FITS header file
v1=-40e3
delv= 1.030571969000E+03
numv=78

# The coordinates of the FITS file corresponding to the galactic coordinates (l,b)=(198,32)
b0=20
l0=20

# Size of the image (deg x deg)
size=20

image=fits.open('../Data/lab_data.fit')
data=image[0]

def TBave(vlsr):
    # This function returns the average T_B of the image at a particular v_lsr.
    tb=0
    nv=int((vlsr-v1)/delv)
    for i in range(-size,size+1):
        for j in range(-size,size+1):
            tb+=data.section[nv,b0+i,l0+j]
    tb=tb/(2*size+1)**2
    return tb

def power_spec(vlsr,boxsize):
    # This function calculates the power spectrum of the image
    nv=int((vlsr-v1)/delv)
    thetab=np.zeros((41,41))
    tba=TBave(vlsr)

    # Boxcar smoothing with width=boxsize
    for i in range(2*size+1):
        for j in range(2*size+1):
                for b in range(-int(boxsize/2),int(boxsize/2)+1):
                    thetab[i][j]+=(data.section[nv+b,b0+(i-size),l0+(j-size)]-tba)*np.sin(i*np.pi/(2*size))*np.sin(j*np.pi/(2*size))
                thetab[i][j]=thetab[i][j]/boxsize

    # Autocorrelation and Fourier transform
    autocorr=correlate2d(thetab,thetab)/(2*size+1)**4
    power2D=fftpack.fft2(autocorr)
    power2D=fftpack.fftshift(power2D)
    power2D=np.abs(power2D)

    n=len(power2D[0])

    # Create a grid of k values corresponding to the 2D power spectrum image
    kx=np.linspace(-1,1,n)
    ky=np.linspace(-1,1,n)
    kx=np.outer(np.ones(n),kx)
    ky=np.outer(ky,np.ones(n))
    k=np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            k[i][j]=np.sqrt(kx[i][j]**2+ky[i][j]**2)

    # Azimuthal averaging
    k1D,power1D=azimuthal_average.aave(k,power2D,40)

    # Calculating the positions of the array at k=0.3, k=0.9 
    n3=0
    n9=0
    minn3=10
    minn9=10
    for i in range(len(k1D)):
        if abs(k1D[i]-0.3)<minn3: 
            n3=i
            minn3=abs(k1D[i]-0.3)
        if abs(k1D[i]-0.9)<minn9: 
            n9=i
            minn9=abs(k1D[i]-0.9)

    k1D=np.log10(k1D[1:])
    power1D=np.log10(power1D[1:])

    # Fitting the log-log power spectrum to a linear function
    fit=stats.linregress(k1D[n3:n9],power1D[n3:n9])
    slope=fit[0]+0.2
    intercept=fit[1]
    error=fit[4]
    x=np.linspace(np.min(k1D),np.max(k1D),100)
    y=slope*x+intercept

    power1D=10**(power1D)
    k1D=10**(k1D)
    y=10**(y)
    x=10**(x)

    # Uncomment this block of code to produce the P(k) vs k plot
    # plt.loglog(k1D,power1D,label='Calculation')
    # plt.loglog(x,y,label="Fit")
    # plt.xlabel("k (deg^-1)")
    # plt.ylabel("P(k) (arbitrary units)")
    # plt.title("v="+str(vlsr/1e3)+" km/s"+"\n beta="+str('%.2f'%slope))
    # plt.legend()
    # plt.show()

    return slope,error

def main():
    box_width=1
    vlsr=np.linspace(-30,20,51)
    slope=np.zeros(51)
    error=np.zeros(51)

    tba=np.zeros(51)
    for i in range(51):
        tba[i]=TBave(10**3*vlsr[i])

    for i in range(51):
        if i%10==0: print("Running...")
        slope[i],error[i]=power_spec(10**3*vlsr[i],box_width)
        
    plt.scatter(vlsr,slope,s=20)
    
    plt.scatter(vlsr[int((0+15)/2)],slope[int((0+15)/2)],c='g',marker='v',s=100,label="IVC")
    plt.scatter(vlsr[int((16+25)/2)],slope[int((16+25)/2)],c='m',marker='v',s=100,label="WNM (-ve v)")
    plt.scatter(vlsr[int((26+45)/2)],slope[int((26+45)/2)],c='k',marker='v',s=100,label="CNM")
    plt.scatter(vlsr[int((46+50)/2)],slope[int((46+50)/2)],c='r',marker='v',s=100,label="WNM (+ve v)")
    plt.axvline(x=-15,c='k')
    plt.axvline(x=-5,c='k')
    plt.axvline(x=15,c='k')
    plt.errorbar(vlsr,slope,yerr=error)
    plt.legend(title="Median values")
    plt.title("LAB Survey\nMultiphase ISM power-law index values")
    plt.xlabel("v_lsr (km/s)")
    plt.ylabel("beta")
    plt.show()
    
    print("\n Median values:")
    print("IVC: "+str(slope[int((0+15)/2)]))
    print("WNM -ve: "+str(slope[int((16+25)/2)]))
    print("CNM: "+str(slope[int((26+45)/2)]))
    print("WNM +ve: "+str(slope[int((46+50)/2)]))
    print("Overall: "+str(slope[int((0+50)/2)]))
    print("\n Column density weighted average:")
    print("IVC: " +str(np.sum(np.multiply(slope[0:15],tba[0:15]))/np.sum(tba[0:15])))
    print("WNM -ve: " +str(np.sum(np.multiply(slope[16:25],tba[16:25]))/np.sum(tba[16:25])))
    print("CNM: " +str(np.sum(np.multiply(slope[26:45],tba[26:45]))/np.sum(tba[26:45])))
    print("WNM +vs: " +str(np.sum(np.multiply(slope[46:50],tba[46:50]))/np.sum(tba[46:50])))
    print("Overall: " +str(np.sum(np.multiply(slope[:],tba[:]))/np.sum(tba[:])))

if __name__ == "__main__":
    main()

# Uncomment this to produce a P(k) vs k plot at a particular v_lsr
# print(power_spec(9.28e3,1)[0])