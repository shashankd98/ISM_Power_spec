import numpy as np 
from astropy.io import fits
from matplotlib import pyplot as plt

# This generates column density maps for each of the phases of the ISM

v1=-40e3
delv= 1.030571969000E+03
numv=78

b0=20
l0=20

size=20

def img(medium):
    vivc=[-30e3,-15e3]
    vwnmn=[-14e3,-5e3]
    vcnm=[-4e3,15e3]
    vwnmp=[16e3,20e3]
    vwhole=[-30e3,20e3]
    nv1=0
    nv2=0
    if medium=="ivc": 
        nv1=int((vivc[0]-v1)/delv)
        nv2=int((vivc[1]-v1)/delv)
    if medium=="wnmn": 
        nv1=int((vwnmn[0]-v1)/delv)
        nv2=int((vwnmn[1]-v1)/delv)
    if medium=="cnm": 
        nv1=int((vcnm[0]-v1)/delv)
        nv2=int((vcnm[1]-v1)/delv)
    if medium=="wnmp": 
        nv1=int((vwnmp[0]-v1)/delv)
        nv2=int((vwnmp[1]-v1)/delv)
    if medium=="whole":
        nv1=int((vwhole[0]-v1)/delv)
        nv2=int((vwhole[1]-v1)/delv)

    image=fits.open('../Data/lab_data.fit')
    tb=np.zeros((2*size+1,2*size+1))
    for i in range(0,2*size+1):
        for j in range(0,2*size+1):
            for v in range(nv1,nv2+1):
                tb[i][j]+=image[0].section[v,b0-i+size,l0+j-size]
            tb[i][j]/=(nv2-nv1+1)
    return tb

def main():
    choice="whole" #can be 'ivc', 'wnmn', 'cnm', 'wnmp', or 'whole'
    plt.imshow(img(choice),extent=(188,208,22,42)) 
    plt.xlabel("Galactic longitube l (deg)")
    plt.ylabel("Galactic latitude b (deg)")
    if choice=="ivc": plt.title("LAB Survey\nIntermediate velocity gas component")
    elif choice=="wnmn": plt.title("LAB Survey\nWarn neutral medium (-ve velocity)")
    elif choice=="cnm": plt.title("LAB Survey\nCold neutral medium")
    elif choice=="wnmp": plt.title("LAB Survey\nWarn neutral medium (+ve velocity)")
    elif choice=="whole": plt.title("LAB Survey\nWhole velocity range")
    plt.colorbar()
    plt.show()

if __name__ == "__main__":
    main()