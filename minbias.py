execfile("utilities.py")
from minbiasdata import list as events
from histogram import histo
import math as math



# Declaring histograms

h_pt=histo(100,0.5,2.5) # 100 intervals between 0.5 and 2.5
h_e=histo(100,1.5,3)
h_eta=histo(100,-3.5,3.5)

# Dictionary for counting particles
flavours={}

# Examples...


for event in events: # Looping over all events

    for [px,py,pz,name] in event: # Filling variables

        pt=math.sqrt(px*px+py*py) # Calculating transverse momentum
        h_pt.fill(pt) # Filling it to histogram

        e=math.sqrt(px*px+py*py+pz*pz+mass[name]*mass[name]) # Energy of particle
        h_e.fill(e) # Filling

        theta=np.arctan2(pt,pz) # Polar angle wrt to beam axis

        eta=-np.log(np.tan(theta/2.)) # Pseudo-rapidity
        h_eta.fill(eta) # Filling


        if name in flavours: 
            flavours[name]+=1 # Count particle by name
        else: # If this is the first instance, start counter
            flavours[name]=1
            
print flavours # Printing flavour

fitres=h_e.fit("pow",True) # Fitting energy distribution to power law

print fitres["prob"] # Printing fit probability

print h_eta.fit("p2",True)["prob"] # Fitting eta-distribution w
                                   # quadratic and plotting.
                                   # Printing fit probability


fitres=h_pt.fit("exp",True) # Fitting p_t to exponential and plotting

print fitres["prob"]


# Here's "the business end of things". Looking for
# Lambda0-particles. These decay to proton + pi-, and antiparticle
# decays to antiproton + pi+.

# Strategy
# - identify proton-pion pairs close in angle
# - add momenta
# - add energies
# - mass of the pair from E^2 = p^2 + m^2

# Histogram for masses
h_m=histo(10,1.11,1.12)

for event in events:
    for [px,py,pz,name] in event: # Loop over all particles

        if name!="p+" and name !="pbar-": continue # If not a p or pbar, move on.

        for [pxpi,pypi,pzpi,namepi] in event: # Loop over the other particles and assign values

            if not((name=="p+" and namepi=="pi-") or (name=="pbar-" and namepi=="pi+")): continue # If not correct pairing, move on.
            scalarproduct=px*pxpi+py*pypi+pz*pzpi
            pp=math.sqrt(px*px+py*py+pz*pz)
            ppi=math.sqrt(pxpi*pxpi+pypi*pypi+pzpi*pzpi)
            costheta=scalarproduct/(pp*ppi) # cos(v)
            if costheta<0.995: continue # Particles must be close
            Ep=math.sqrt(px**2+py**2+pz**2+mass[name]**2) # Energy of proton
            Epi=math.sqrt(pxpi**2+pypi**2+pzpi**2+mass[namepi]**2) # Energy of pion
            p=[px+pxpi,py+pypi,pz+pzpi] # Sum of momenta
            E=Ep+Epi # Sum of energies
            m=math.sqrt(E**2-(p[0]**2+p[1]**2+p[2]**2)) # Mass
            h_m.fill(m) # Fill mass to histogram

fitres=h_m.fit("norm",True) # Fit and plot
print fitres["prob"] # Print fit probability
print "mass: ",fitres["parameters"][1] # Print measured Lambda0-mass
                                       # Table value:
                                       # 1.115683 GeV
