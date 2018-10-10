smeared=False

if smeared:
    from higgs_smeared import list as events
else:
    from higgsdata import list as events


from histogram import histo
import math as math

if smeared:
    h=histo(100,115,135)
else:
    h=histo(51,124.99,125.01)


for event in events:
    for i in range (0,len(event)):
        for j in range (0,i):
            [px1,py1,pz1,name1]=event[i]
            [px2,py2,pz2,name2]=event[j]
            p=[px1+px2,py1+py2,pz1+pz2]
            p2=p[0]**2+p[1]**2+p[2]**2
            E1=math.sqrt(px1*px1+py1*py1+pz1*pz1)
            E2=math.sqrt(px2*px2+py2*py2+pz2*pz2)
            E=E1+E2
            h.fill(math.sqrt(E*E-p2))

if smeared:
    # Fitting Gauss
    fitres=h.fit("norm",True)
else:
    # Fitting Breit-Wigner
    fitres=h.fit("resonans",True)

print fitres["prob"]
print fitres["parameters"]

