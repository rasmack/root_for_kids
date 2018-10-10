execfile("utilities.py")
from histogram import histo



# Fitting gaussian + linear background
h=histo(100,0,10)
for i in range(0,1000):
    h.fill(linbg(0,10,-1.1,15))

g=histo(100,0,10)
for i in range(0,1000):
    g.fill(norm(3,0.5))

hist=h+g
print hist.fit("norm+p1",True)



# Fitting exponential
h=histo(100,0,10)
data=[]
c=800
k=-0.23
for i in range (0, 100): data.append(flat(0.9,1.1)*c*np.exp(k*(float(i)+0.5)/10))
h.setcontent(data)
h.computeerrors()
print h.fit("exp",True)["prob"]


# Demonstrating central limit theorem
h=histo(10,2,7)
for i in range(0,1000):
    h.fill(
        flat(1,3)
        +flat(0.5,3)
        +flat(0,1.5)
        )
print h.fit("norm",True)["prob"]


"""
h=histo(10,0,1)
g=histo(10,0,40)

for j in range(0,1000):
    for i in range(0,1000):
        h.fill(flat(0,1))
    fitres=h.fit("p0")
#    print fitres["chi2"]
    g.fill(fitres["chi2"])
    h.clear()
g.plot()
#g.dump()
"""
"""
h=histo(10,0,1)

for i in range(0,1000):
    h.fill(flat(0,0.3)+flat(0,0.3)+flat(0,0.3))
fitres=h.fit("norm",True)
print "Chi^2: ",fitres["chi2"], " Sandsynlighed: ", fitres["prob"]
"""

"""
h=histo(10,0,1)
for i in range(0,1000):
    h.fill(flat(0,1))

print h.fit("p0",True)["prob"]
print h.fit("p1",True)["prob"]
print h.fit("p2",True)["prob"]
"""

"""
h=histo(10,0,1)
chi2sum=0
for i in range(0,1000):
    for j in range(0,10000):
        h.fill(flat(0,1))
    chi2sum+=h.fit("p0")["prob"]
    average=chi2sum/float(i)
    print average
    h.clear()
"""

"""
h=histo(10,0,1)
for i in range(0,10000):
    h.fill(norm(0.5,0.1))
print h.fit("norm")
"""

"""
h=histo(10,0,1)
for i in range(0,1000):
    for j in range(0,10000):
        h.fill(secret())
    print h.fit("p0")["chi2"], sum(h.getcontent())
"""
