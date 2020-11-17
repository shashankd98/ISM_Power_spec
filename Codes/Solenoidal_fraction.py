import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
from scipy import fftpack
from scipy.optimize import curve_fit
import azimuthal_average

# This calculates the solenoidal fraction in the Orion B cloud using IRAM-30 data. 
w0=fits.open('../Data/13co10-w0.fits')[0].data
w0=w0[0]
w1=fits.open('../Data/13co10-w1c.fits')[0].data
w1=w1[0]
w2=fits.open('../Data/13co10-w2c.fits')[0].data
w2=w2[0]

ny=len(w1[0])
nx=len(np.transpose(w1)[0])
nz=ny
kspace2D=np.zeros((nx,ny))
kspace3D=np.zeros((nx,ny,nz))
cosine=np.zeros((nx,ny,nz))
kx=np.linspace(-1,1,nx)
ky=np.linspace(-1,1,ny)
kz=np.linspace(-1,1,nz)

# Create a 2D and 3D array containing the k values
# cosine is a 3D array containing the values of (kx^2+ky^2)/k^2
for i in range(nx):
    for j in range(ny):
        kspace2D[i][j]=np.sqrt(kx[i]**2+ky[j]**2)
        if kspace2D[i][j]==0: kspace2D[i][j]=kspace2D[i][j-1]/2
    
for i in range(nx):
    for j in range(ny):
        kspace3D[i][j][:]=np.sqrt(kspace2D[i][j]**2+kz[:]**2)
        cosine[i][j][:]=kspace2D[i][j]**2/(kspace3D[i][j][:])**2

# Applying the window function
w0p=np.zeros((nx,ny))
w1p=np.zeros((nx,ny))
w2p=np.zeros((nx,ny))
for i in range(nx):
    for j in range(ny):
        w0p[i][j]=w0[i][j]*np.sin(np.pi*i/nx)*np.sin(np.pi*j/ny)
        w1p[i][j]=w1[i][j]*np.sin(np.pi*i/nx)*np.sin(np.pi*j/ny)
        w2p[i][j]=w2[i][j]*np.sin(np.pi*i/nx)*np.sin(np.pi*j/ny)

def power(k,a,b):
    return a*k**b

def power_spectrum():
    # Calculate the power spectrum of W0 and W1
    power2D=np.fft.fft2(w0p)/(nx*ny)
    power2D=fftpack.fftshift(power2D)
    power2D=np.abs(power2D)**2
    k1D,power1D=azimuthal_average.aave(kspace2D,power2D,100)

    # Fit to the power function
    f_fit=curve_fit(power,k1D[2:30],power1D[2:30],p0=[1e-5,-2])[0]
    pcov=curve_fit(power,k1D[2:30],power1D[2:30],p0=[1e-5,-2])[1]
    perr=np.sqrt(np.diag(pcov))
    lin=np.linspace(np.min(k1D),np.max(k1D),100)
    out=power(lin,f_fit[0],f_fit[1])
    plt.title("Power spectrum of W0")
    plt.xlabel("k")
    plt.ylabel("P(k) (arbitrary units)")
    plt.loglog(k1D,power1D)
    plt.loglog(lin,out)
    plt.show()

    a1=f_fit[0]
    b1=f_fit[1]
    print("\nFit to P(k)=a*k^b:")
    print("\na1="+str(a1))
    print("b1="+str(b1))
    print("Uncertainity in b1="+str(perr[1]))

    power2D=np.fft.fft2(w1p)/(nx*ny)
    ppower2D=fftpack.fftshift(power2D)
    ppower2D=np.abs(ppower2D)**2
    kp1D,ppower1D=azimuthal_average.aave(kspace2D,ppower2D,100)

    fp_fit=curve_fit(power,kp1D[3:30],ppower1D[3:30],p0=[1e-5,-2])[0]
    pcov2=curve_fit(power,kp1D[3:30],ppower1D[3:30],p0=[1e-5,-2])[1]
    perr2=np.sqrt(np.diag(pcov2))
    linp=np.linspace(np.min(kp1D),np.max(kp1D),100)
    outp=power(linp,fp_fit[0],fp_fit[1])

    plt.title("Power spectrum of W1")
    plt.xlabel("k")
    plt.ylabel("P(k) (arbitrary units)")
    plt.loglog(kp1D,ppower1D)
    plt.loglog(linp,outp)
    plt.show()
    a2=fp_fit[0]
    b2=fp_fit[1]
    print("\na2="+str(a2))
    print("b2="+str(b2))
    print("Uncertainity in b2="+str(perr2[1]))

    return a1,b1,a2,b2

def constants():
    # Calculate the A and B constants
    a1,b1,a2,b2=power_spectrum()
    sum1=np.sum(power(kspace3D,a1,b1)/kspace3D)
    sum2=np.sum(power(kspace2D,a1,b1)/kspace2D)
    sum3=np.sum(power(kspace3D,a2,b2)*cosine/kspace3D)
    sum4=np.sum(power(kspace2D,a2,b2)/kspace2D)
    
    k02=np.sqrt(kx[int(nx/2)]**2+ky[int(ny/2)]**2)
    k03=np.sqrt(kx[int(nx/2)]**2+ky[int(ny/2)]**2+kz[int(nz/2)]**2)
    f02=power(k02,a1,b1)
    f03=power(k03,a1,b1)

    A=(sum1-f03)/(sum2-f02)
    B=sum3/sum4
    print("\nA="+str(A))
    print("B="+str(B))
    return A,B

def g12():
    # Calculate the g12 constant
    psvs=np.average(w1**2)
    ps=np.average(w0**2)
    pvs=np.average(w2)
    p=np.average(w0)
    factor=(psvs/ps)/(pvs/p)
    print("\ng12="+str(factor))
    return factor

def main():
    # Put all the functions together and calculate R
    A,B=constants()
    g=g12()
    w0a=np.average(w0)
    w0sa=np.average(w0**2)
    w1sa=np.average(w1**2)
    w2a=np.average(w2)
    R=abs(w1sa/w0sa)*(w0sa/w0a**2/(1+A*(w0sa/w0a**2-1)))/(g*w2a/w0a)*B
    print("\nSolenoidal fraction R="+str(R))

if __name__ == "__main__":
    main()


