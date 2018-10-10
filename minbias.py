execfile("utilities.py")
from minbiasdata import list as events
from histogram import histo
import math as math



# Deklaration af histogrammer

h_pt=histo(100,0.5,2.5) # 100 intervaller mellem 0.5 og 2.5
h_e=histo(100,1.5,3)
h_eta=histo(100,-3.5,3.5)

# Deklaration af dictionary til at taelle partikler
flavours={}

# Her er en bunke eksempler paa ting man kan goere. Vaer opmaerksom
# paa, at danske karakterer er en daarlig ide.


for event in events: # Vi looper over containeren af events og giver
                     # hver enkelt event navnet "event".

    for [px,py,pz,name] in event: # Hver event er jo en liste af
                                  # partikler repraesenteret ved en
                                  # liste.

        pt=math.sqrt(px*px+py*py) # Beregner den transverse impuls
        h_pt.fill(pt) # Fylder den i histogram

        e=math.sqrt(px*px+py*py+pz*pz+mass[name]*mass[name]) #Beregner energi af partikel ud fra E^2=p^2+m^2
        h_e.fill(e) # Fyler energi i histogram

        theta=np.arctan2(pt,pz) # Beregner vinkel til z-aksen
                                # (beam-aksen). Bemaerk, hvordan
                                # funktionen tager (y,x) som argument
                                # i stedet for y/x. Det betyder, at
                                # den kan vaelge vinklen i den rigtige
                                # kvadrant.

        eta=-np.log(np.tan(theta/2.)) # Beregner en stoerrelse kaldet pseudorapiditet
        h_eta.fill(eta) # Fylder den i histogram


        if name in flavours: # Hvis navnet allerede findes i listen
                             # over partikler, saa laeg 1 til
                             # taelletallet.
            flavours[name]+=1
        else: # Hvis ikke, saa opret den og sig, at vi har fundet en.
            flavours[name]=1
            
print flavours # Skriver fordelingen vlandt partikelflavours ud.

fitres=h_e.fit("pow",True) # Fitter energifordelingen til en
                           # potensfunktion og plotter histogram +
                           # fit. Gemmer informationer om regressionen
                           # i variablen fitres.

print fitres["prob"] # Udskriver sandsynligheden for at faa dette
                     # chi^2 eller et der er hoejere.

print h_eta.fit("p2",True)["prob"] # Fitter eta-fordelingen med et
                                   # andengradspolynomium og
                                   # plotter. Udskriver
                                   # sandsynligheden for at faa dette
                                   # chi^2 eller hoejere.

fitres=h_pt.fit("exp",True) # Fitter fordelingen af den transverse
                            # impuls til en eksponentiel
                            # udvikling. Plotter.
print fitres["prob"]


# Dette her er "the business end of things". Her skal vi kigge efter
# Lambda0-partikler. Disse henfalder ofte til proton + pi-, og
# antipartiklen henfalder til antiproton + pi+.

# Strategien er 
# - at identificere alle proton-pion par der ligger taet paa hinanden i vinkel.
# - udregne vektorsummen af deres impulser
# - udregne summen af deres energier.
# - udregne massen af parret ud fra E^2 = p^2 + m^2

# Histogram til de invariante masser
h_m=histo(10,1.11,1.12)

for event in events:
    for [px,py,pz,name] in event: # Loop over alle partikler

        if name!="p+" and name !="pbar-": continue # Hvis partiklen
                                                   # ikke er en proton
                                                   # eller antiproton,
                                                   # saa gaa videre
                                                   # til den naeste.
        for [pxpi,pypi,pzpi,namepi] in event: # For hver partikel,
                                              # loop over alle de
                                              # andre og tildel
                                              # vaerdier som vist til
                                              # venstre.
            if not((name=="p+" and namepi=="pi-") or (name=="pbar-" and namepi=="pi+")): continue # Hvis ikke en korrekt parring af partikler, saa gaa videre.
            prikprodukt=px*pxpi+py*pypi+pz*pzpi
            pp=math.sqrt(px*px+py*py+pz*pz)
            ppi=math.sqrt(pxpi*pxpi+pypi*pypi+pzpi*pzpi)
            costheta=prikprodukt/(pp*ppi) # Udregn cos(v), hvor vinklen er vinklen mellem partiklerne
            if costheta<0.995: continue # Partiklerne skal vaere taet paa hinanden.
            Ep=math.sqrt(px**2+py**2+pz**2+mass[name]**2) # Udregn energi af proton (Bemaerk, at vi traekker paa viden om massen)
            Epi=math.sqrt(pxpi**2+pypi**2+pzpi**2+mass[namepi]**2) # Energi af pion
            p=[px+pxpi,py+pypi,pz+pzpi] # Vektorsum af impulser
            E=Ep+Epi # Sum af energier
            m=math.sqrt(E**2-(p[0]**2+p[1]**2+p[2]**2)) # Udregn masse
            h_m.fill(m) # Fyld masse i histogram

fitres=h_m.fit("norm",True) # Fit og plot
print fitres["prob"] # Udskriv fitsandsynlighed
print "mass: ",fitres["parameters"][1] # Udskriv den maalte masse af
                                       # Lambda0. Tabelvaerdi er
                                       # 1.115683 GeV
