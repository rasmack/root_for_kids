execfile("utilities.py")
from histogram import histo

h_chi2=histo(20,0,40)


# Her tjekker vi chi^2-fordelingen for 10 frihedsgrader
# 11 bins i histogrammet
# 1 parameter i regressionen

h=histo(11,0,1) # 11 bins mellem 0 og 1

for j in range(0,1000): # Generer 1000 vaerdier af chi^2
    if j%100==0: print j # Hold styr paa, hvor langt vi er.
    for i in range(0,1000): # 1000 iterationer af...
        h.fill(flat(0,1)) # ... at smide tal i histogrammet efter en
                          # flad fordeling

    h_chi2.fill(h.fit("p0")["chi2"]) # Fit til en konstant og smid
                                     # chi^2 i det rigtige histogram.
    h.clear() # Slet indholdet af arbejdshistogrammet og goer klar til
              # at proeve igen.

h_chi2.plot() # Plot chi^2-fordelingen
