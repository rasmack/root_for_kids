execfile("utilities.py")
from jpsidata import list as events


h_m=histo(100,3.06,3.14)
h_spec=histo(100,0,100)
h_pt=histo(100,2,6)

for event in events:
    for electron in event:
        if electron[3]!="e-": continue
        for positron in event:
            if positron[3]!="e+": continue
            [px1,py1,pz1,pid1]=electron
            h_pt.fill( math.sqrt(px1**2+py1**2))
            E1=math.sqrt(px1**2+py1**2+pz1**2+mass[pid1]**2)
            [px2,py2,pz2,pid2]=positron
            E2=math.sqrt(px2**2+py2**2+pz2**2+mass[pid2]**2)
            p=[px1+px2,py1+py2,pz1+pz2]                                                                                                                                           
            E=E1+E2
            m=math.sqrt(E**2-(p[0]**2+p[1]**2+p[2]**2))
            h_m.fill(m)
            h_spec.fill(m)
print h_pt.fit("exp",True)["prob"]
h_spec.plot()
fitres=h_m.fit("norm",True)
print "Fit probability:", fitres["prob"]
print "Fitted mass: ", fitres["parameters"][1], " (Table value: 3.096916 GeV)"

