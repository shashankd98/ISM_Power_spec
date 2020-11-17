import numpy as np

def aave(k,power2D,nbins):
    n=len(power2D[0])
    m=len(np.transpose(power2D)[0])

    # Reshape k and power2D into 1D arrays
    power=np.reshape(power2D,n*m)
    ka=np.reshape(k,n*m)

    # Sort power2D into bins based on their k values
    k1D=np.zeros(nbins-1)
    power1D=np.zeros(nbins-1)
    kbin=np.linspace(np.min(k),np.max(k),nbins)
    for b in range(nbins-1):
        num=0
        for i in range(n*m):
            if ka[i]>=kbin[b] and ka[i]<=kbin[b+1] and power[i]>1e-20:
                num+=1
                k1D[b]+=ka[i]
                power1D[b]+=power[i]
        if num!=0: 
            k1D[b]/=num
            power1D[b]/=num

    return k1D,power1D