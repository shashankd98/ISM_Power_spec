import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
from scipy import fftpack
from scipy import stats
from scipy.signal import correlate2d
import azimuthal_average

def main():
    hdulist = fits.open('../Data/planck_data.fit')  
    col_map=hdulist[0].data   

    plt.imshow(col_map,extent=(-7.5,7.5,-7.5,7.5))
    plt.title("Planck 857 GHz Survey \n (centered at (198,32) )")
    plt.xlabel("Galactic longitude l (deg)")
    plt.ylabel("Galactic latitude b (deg)")
    plt.show()

    nl=int(15*(60/5)+1)
    nb=int(15*(60/5)+1)

    # Power spectrum calculation 
    thetab=np.zeros((nl,nb))
    ave=np.average(col_map)
    for i in range(nl):
        for j in range(nb):
            thetab[i][j]=(col_map[i][j]-ave)*np.sin(i*np.pi/nl)*np.sin(j*np.pi/nb)

    # Autocorrelation and Fourier transform
    autocorr=correlate2d(thetab,thetab)/(nl*nb)**2
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
    k1D,power1D=azimuthal_average.aave(k,power2D,100)

    k1D=np.log10(k1D[1:])
    power1D=np.log10(power1D[1:])

    # Fit the log-log power spectrum to a linear function
    fit=stats.linregress(k1D[:],power1D[:])
    slope=fit[0]
    intercept=fit[1]
    error=fit[4]
    x=np.linspace(np.min(k1D),np.max(k1D),1000)
    y=slope*x+intercept

    power1D=10**(power1D)
    k1D=10**(k1D)
    x=10**(x)
    y=10**(y)

    plt.title("Planck Survey Power Spectrum \n beta="+str('%.2f'%slope))
    plt.loglog(k1D,power1D,label='Calculation')
    plt.loglog(x,y,label='Fit')
    plt.legend()
    plt.ylabel("P(k) (arbitrary units)")
    plt.xlabel("k (deg^-1)")
    plt.show()

    print("beta="+str(slope))
    print("error="+str(error))

if __name__ == "__main__":
    main()