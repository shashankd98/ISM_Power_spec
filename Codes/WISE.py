import numpy as np 
from astropy.io import fits
from matplotlib import pyplot as plt
from scipy import stats
from scipy import fftpack
from scipy.signal import correlate2d
import azimuthal_average

def main():

    # Read the FITS file
    hdulist=fits.open('../Data/WISE_data.fit')
    image=hdulist[0]
    data=image.section[:][:]

    # Uncomment this block to plot the column density map
    # plt.title("WISE 12 um Survey \n (centered at (198,32) )")
    # plt.xlabel("Galactic longitude l (deg)")
    # plt.ylabel("Galactic latitude b (deg)")
    # plt.imshow(data,extent=(-2.35,2.35,-2.35,2.35))
    # plt.show()

    n=len(data[0])

    # Power spectrum calculation
    tba=np.average(data)
    thetab=np.zeros((n,n),float)
    for i in range(n):
        for j in range(n):
            thetab[i][j]=(data[i][j]-tba)*np.sin(i*np.pi/n)*np.sin(j*np.pi/n)

    # Autocorrelation and Fourier transform
    autocorr=correlate2d(thetab,thetab,mode='full')/n**4
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
    k1D,power1D=azimuthal_average.aave(k,power2D,int(n/2))

    k1D=np.log10(k1D[1:])
    power1D=np.log10(power1D[1:])

    # Fit the log-log power spectrum to a linear function
    fit=stats.linregress(k1D[2:20],power1D[2:20])
    slope=fit[0]
    intercept=fit[1]
    error=fit[4]
    x=np.linspace(k1D[0],k1D[50],1000)
    y=slope*x+intercept

    power1D=10**(power1D)
    k1D=10**(k1D)
    x=10**(x)
    y=10**(y)

    plt.loglog(k1D,power1D,label="Calculation")
    plt.loglog(x,y,label="Fit")
    plt.title("WISE Survey Power Spectrum \n beta="+str('%.2f'%slope))
    plt.legend()
    plt.ylabel("P(k) (arbitrary units)")
    plt.xlabel("k (deg^-1)")
    plt.show()

    print("beta="+str(slope))
    print("error="+str(error))

if __name__ == "__main__":
    main()