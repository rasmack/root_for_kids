import matplotlib.pyplot as plt
from scipy.stats import chi2 as chi2dist
import numpy as np
import math as math
import scipy.optimize as optimization
import operator

class histo:
    def __init__(self,nbins,xlow,xhigh):
        self.data=[]
        self.errors=[]
        self.bincenter=[]
        self.binlow=[]
        self.binhigh=[]
        self.nbins=nbins
        self.xlow=float(xlow)
        self.xhigh=float(xhigh)
        for x in range (0,nbins):
            self.data.append(0)
            self.errors.append(0)
            binwidth=(1.*xhigh-xlow)/(1.*nbins)
            self.bincenter.append(xlow+binwidth*(x+0.5))
            self.binlow.append(xlow+binwidth*x)
            self.binhigh.append(xlow+binwidth*(x+1))
        self.manerr=False

    def __add__(self,x):
        c1=np.array(self.data)
        c2=np.array(x.getcontent())
        e1=np.array(self.geterrors())
        e2=np.array(x.geterrors())
        c=c1+c2
        e=np.sqrt(e1*e1+e2*e2)
        y=histo(x.getnbins(),x.getxlow(),x.getxhigh())
        y.setcontent(c.tolist())
        y.seterr(e.tolist())
        return y

    def clear(self):
        self.data=[]
        self.errors=[]
        self.manerr=False
        for x in range (0,self.nbins):
            self.data.append(0)
            self.errors.append(0)

    def setcontent(self,content):
        self.data=content

    def seterr(self,errors):
        self.errors=errors
        manerr=True

    def getcontent(self):
        return self.data

    def geterrors(self):
        self.computeerrors()
        return self.errors

    def seterrors(self,errors):
        self.errors=errors
        self.manerr=True

    def getnbins(self):
        return self.nbins

    def getxlow(self):
        return self.xlow

    def getxhigh(self):
        return self.xhigh


    def fill(self,x):
        if x<self.xlow: return
        bin=int((x-self.xlow)/(self.xhigh-self.xlow)*self.nbins)
        if(0<=bin<self.nbins):
            self.data[bin]+=1.
        if(x==self.binhigh[self.nbins-1]):
            self.data[self.nbins-1]+=1.

    def computeerrors(self):
        if not self.manerr:
            for x in range (0,self.nbins):
                if(self.data[x]==0):
                    self.errors[x]=math.sqrt(2.3)
                else:
                    self.errors[x]=math.sqrt(self.data[x])

    def mean(self):
        return sum(self.data*np.array(self.bincenter))/sum(self.data)

    def spread(self):
        mu=self.mean()
        offset=np.array(self.bincenter)-mu
        varians=sum(self.data*(offset*offset))/sum(self.data)
        return math.sqrt(varians)

    def dump(self):
        self.computeerrors()
        print self.data
        print self.errors

    def entries(self):
        return sum(self.data)

    def plot(self):
        plt.clf()
        self.computeerrors()
        plt.errorbar(self.bincenter, self.data, yerr=self.errors, fmt='o')
        if self.entries()!=0:
            plt.axis([
                    self.xlow,
                    self.xhigh,
                    min(0,min(np.array(self.data)-np.array(self.errors))),
                    max(np.array(self.data)+np.array(self.errors))
                    ])
        plt.show()
        plt.close()

    def fit(self,type="flat",plot=False):
        self.computeerrors()

        # Defining fit functions
        def gauss(input, a, b, c):
            return map(lambda x: a*np.exp(-0.5*((x-b)/c)**2),input)

        def gaussp0(input, a, b, c,k):
            return map(lambda x: a*np.exp(-0.5*((x-b)/c)**2)+k,input)

        def gaussp1(input, a, b, c,slope,intercept):
            return map(lambda x: a*np.exp(-0.5*((x-b)/c)**2)+slope*x+intercept,input)

        def bw(input,k,m,gamma):
            return map(lambda x: k/((x**2-m**2)**2+m**2*gamma**2),input)

        def bwlin(input,k,m,gamma,a,b):
            return map(lambda x: k/((x**2-m**2)**2+m**2*gamma**2+a*x+b),input)

        def p0(input,a):
            return map(lambda x: a,input)
        
        def p1(input, a,b):
            return map(lambda x: a*x+b,input)

        def p2(input, a,b,c):
            return map(lambda x: a*x*x+b*x+c,input)

        def exp(input, c,k):
            return map(lambda x: c*np.exp(k*x),input)
        
        def pow(input,b,a):
            return map(lambda x: b*x**a,input)

        # Setting up axes
        if(plot):
            plt.errorbar(self.bincenter, self.data, yerr=self.errors, fmt='o')
            if self.entries()!=0:
                plt.axis([
                        self.xlow,
                        self.xhigh,
                        min(0,min(np.array(self.data)-np.array(self.errors))),
                        max(np.array(self.data)+np.array(self.errors))
                        ])

        # Performing the actual fit based on the option given
        if(type=="norm"):
            initial=[self.entries(),self.mean(),self.spread()]
            fit=optimization.curve_fit(gauss, np.array(self.bincenter), np.array(self.data), initial, np.array(self.errors))
            scale=fit[0][0]
            mu=fit[0][1]
            sigma=fit[0][2]
            chi2=sum(((gauss(np.array(self.bincenter),scale,mu,sigma)-np.array(self.data))/np.array(self.errors))**2)
            # Creating nice smooth curve for plotting fit function
            if(plot):
                xnew = np.linspace(self.xlow,self.xhigh,300)
                l = plt.plot(xnew, gauss(xnew,scale,mu,sigma), 'r--', linewidth=2)
                maxfn=max(gauss(xnew,scale,mu,sigma))

        if(type=="norm+p0"):
            kguess=(self.data[0]+self.data[self.nbins-1])/2.
            index, value = max(enumerate(self.data), key=operator.itemgetter(1))
            initial=[max(self.data)-kguess,self.bincenter[index],1,kguess]
 #           print initial
            fit=optimization.curve_fit(gaussp0, np.array(self.bincenter), np.array(self.data), initial, np.array(self.errors))
            scale=fit[0][0]
            mu=fit[0][1]
            sigma=fit[0][2]
            const=fit[0][3]
            """
            scale=initial[0]
            mu=initial[1]
            sigma=initial[2]
            const=initial[3]
            """
            chi2=sum(((gaussp0(np.array(self.bincenter),scale,mu,sigma,const)-np.array(self.data))/np.array(self.errors))**2)
            # Creating nice smooth curve for plotting fit function
            if(plot):
                xnew = np.linspace(self.xlow,self.xhigh,300)
                l = plt.plot(xnew, gaussp0(xnew,scale,mu,sigma,const), 'r--', linewidth=2)
                maxfn=max(gaussp0(xnew,scale,mu,sigma,const))

        if(type=="norm+p1"):
            kguess=(self.data[0]+self.data[self.nbins-1])/2.
            index, value = max(enumerate(self.data), key=operator.itemgetter(1))
            initial=[max(self.data)-kguess,self.bincenter[index],1,0,kguess]
#            print initial
            fit=optimization.curve_fit(gaussp1, np.array(self.bincenter), np.array(self.data), initial, np.array(self.errors))
            scale=fit[0][0]
            mu=fit[0][1]
            sigma=fit[0][2]
            slope=fit[0][3]
            intercept=fit[0][4]
            """
            scale=initial[0]
            mu=initial[1]
            sigma=initial[2]
            const=initial[3]
            """
            chi2=sum(((gaussp1(np.array(self.bincenter),scale,mu,sigma,slope,intercept)-np.array(self.data))/np.array(self.errors))**2)
            # Creating nice smooth curve for plotting fit function
            if(plot):
                xnew = np.linspace(self.xlow,self.xhigh,300)
                l = plt.plot(xnew, gaussp1(xnew,scale,mu,sigma,slope,intercept), 'r--', linewidth=2)
                maxfn=max(gaussp1(xnew,scale,mu,sigma,slope,intercept))

        if (type=="exp"):
            # Using end bins to estimate scale and constant
            y1=self.data[0]
            y2=self.data[self.nbins-1]
            cguess=(y1/y2)**(self.xhigh/(self.xhigh-self.xlow))
            kguess=np.log(y1/y2)/(self.xlow-self.xhigh)
            initial=[cguess,kguess]
            fit=optimization.curve_fit(exp, np.array(self.bincenter), np.array(self.data), initial, np.array(self.errors))
            c=fit[0][0]
            k=fit[0][1]
            chi2=sum(((exp(np.array(self.bincenter),c,k)-np.array(self.data))/np.array(self.errors))**2)
            # Creating nice smooth curve for plotting fit function
            if(plot):
                xnew = np.linspace(self.xlow,self.xhigh,300)
                l = plt.plot(xnew, exp(xnew,c,k), 'r--', linewidth=2)
                maxfn=max(exp(xnew,c,k))

        if (type=="pow"):
            # Using end bins to estimate scale and constant
            y1=self.data[0]
            y2=self.data[self.nbins-1]
            aguess=np.log(y1/y2)/np.log(self.xlow/self.xhigh)
            bguess=y2**(-aguess)/self.xhigh
            initial=[bguess,aguess]
            fit=optimization.curve_fit(pow, np.array(self.bincenter), np.array(self.data), initial, np.array(self.errors))
            b=fit[0][0]
            a=fit[0][1]
            chi2=sum(((pow(np.array(self.bincenter),b,a)-np.array(self.data))/np.array(self.errors))**2)
            # Creating nice smooth curve for plotting fit function
            if(plot):
                xnew = np.linspace(self.xlow,self.xhigh,300)
                l = plt.plot(xnew, pow(xnew,b,a), 'r--', linewidth=2)
                maxfn=max(pow(xnew,b,a))



        if (type=="resonans"):
            index, value = max(enumerate(self.data), key=operator.itemgetter(1))
#            print index, value, self.bincenter[index]
            initial=[value,self.bincenter[index],10]
#            print initial
            fit=optimization.curve_fit(bw, np.array(self.bincenter), np.array(self.data), initial, np.array(self.errors))
            k=fit[0][0]
            m=fit[0][1]
            gamma=fit[0][2]
            chi2=sum(((bw(np.array(self.bincenter),k,m,gamma)-np.array(self.data))/np.array(self.errors))**2)
            # Creating nice smooth curve for plotting fit function
            if(plot):
                xnew = np.linspace(self.xlow,self.xhigh,300)
                l = plt.plot(xnew, bw(xnew,k,m,gamma), 'r--', linewidth=2)
                maxfn=max(bw(xnew,k,m,gamma))

        if (type=="resonans+bg"):
            print "Determining initial values:"
            seed=self.fit("resonans",False)["parameters"]
            initial=[
                seed[0],
                seed[1],
                seed[2],
                (self.data[self.nbins-1]-self.data[0])/(self.xhigh-self.xlow),
                (self.data[0]+self.data[self.nbins-1])/2.
                ]
            print initial
#            fit=optimization.curve_fit(bwlin, np.array(self.bincenter), np.array(self.data), initial, np.array(self.errors))
#            k=fit[0][0]
#            m=fit[0][1]
#            gamma=fit[0][2]
#            a=fit[0][3]
#            b=fit[0][4]
            k=initial[0]
            m=initial[1]
            gamma=initial[2]
            a=initial[3]
            b=initial[4]
            chi2=sum(((bwlin(np.array(self.bincenter),k,m,gamma,a,b)-np.array(self.data))/np.array(self.errors))**2)
            # Creating nice smooth curve for plotting fit function
            if(plot):
                xnew = np.linspace(self.xlow,self.xhigh,300)
                l = plt.plot(xnew, bwlin(xnew,k,m,gamma,a,b), 'r--', linewidth=2)
                maxfn=max(bwlin(xnew,k,m,gamma,a,b))


        if (type=="p0"):
            initial=[sum(self.data)/len(self.data)]
            fit=optimization.curve_fit(p0, np.array(self.bincenter), np.array(self.data), initial, np.array(self.errors))
            value=fit[0][0]
            chi2=sum(((p0(np.array(self.bincenter),value)-np.array(self.data))/np.array(self.errors))**2)
            # Creating nice smooth curve for plotting fit function
            if(plot):
                xnew = np.linspace(self.xlow,self.xhigh,300)
                l = plt.plot(xnew, p0(xnew,value), 'r--', linewidth=2)
                maxfn=value

        if (type=="p1"):
            # Assuming no slope
            initial=[0.,sum(self.data)/len(self.data)]
            fit=optimization.curve_fit(p1, np.array(self.bincenter), np.array(self.data), initial, np.array(self.errors))
            slope=fit[0][0]
            const=fit[0][1]
            chi2=sum(((p1(np.array(self.bincenter),slope,const)-np.array(self.data))/np.array(self.errors))**2)
            # Creating nice smooth curve for plotting fit function
            if(plot):
                xnew = np.linspace(self.xlow,self.xhigh,300)
                l = plt.plot(xnew, p1(xnew,slope,const), 'r--', linewidth=2)
                maxfn=max(p1(xnew,slope,const))

        if (type=="p2"):
            # Assuming p0
            initial=[0.,0.,sum(self.data)/len(self.data)]
            fit=optimization.curve_fit(p2, np.array(self.bincenter), np.array(self.data), initial, np.array(self.errors))
            a=fit[0][0]
            b=fit[0][1]
            c=fit[0][2]
            chi2=sum(((p2(np.array(self.bincenter),a,b,c)-np.array(self.data))/np.array(self.errors))**2)
            # Creating nice smooth curve for plotting fit function
            if(plot):
                xnew = np.linspace(self.xlow,self.xhigh,300)
                l = plt.plot(xnew, p2(xnew,a,b,c), 'r--', linewidth=2)
                maxfn=max(p2(xnew,a,b,c))


            
        
        if(plot):
            plt.axis([
                    self.xlow,self.xhigh,
                    min(0,min(np.array(self.data)-np.array(self.errors))),
                    max(max(np.array(self.data)+np.array(self.errors)),maxfn*1.02),
                    ])
            plt.show()
            plt.close()
#        else: plt.clf()
        # Calculating degrees of freedom
        ndof=self.nbins-len(fit[0])

        # Calculating fit probability
        prob=1-chi2dist.cdf(chi2,ndof)

        # Returning a dictionary of results
        result={}
        result["parameters"]=fit[0]
        result["chi2"]=chi2
        result["ndof"]=ndof
        result["prob"]=prob
        return result
