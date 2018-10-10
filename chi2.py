execfile("utilities.py")
from histogram import histo

h_chi2=histo(20,0,40)


# chi^2-distribution for 10 degreees of freedom
# 11 bins in histogramkmet
# 1 parameter fit

h=histo(11,0,1) # 11 bins between 0 and 1

for j in range(0,1000): # 1000 chi^2 values
    if j%100==0: print j 
    for i in range(0,1000): # 1000 iterations of...
        h.fill(flat(0,1)) # ... throwing numbers in histogram from a flat distribution.

    h_chi2.fill(h.fit("p0")["chi2"]) # Fit to a constant.
                                     # Histogram chi^2 value
    h.clear() # Rinse, repeat

h_chi2.plot() # Plot chi^2-distribution
