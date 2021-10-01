from scipy.special import gamma
from scipy.stats import norm
from numpy import *
from merge import *
from scipy.integrate import odeint



import ray

def dc(g,s):
    return int(ceil( 1. / ss / (g-1.) ))
def sc(g,d):
    return 1. / (g-1.) * 1. / float(d)


#divisional load at time t [in units of time required to produce a terminally differentiated cell] cf. Eqs. 12 & 13
def D(N,N0,tlife, n, g): return tlife/N0 * power(g,-(n - 1.)) + (n - 1)*(g - 1) + 1 
#probability of cancer by time t [in units of years]
def P(N,N0,mu, n, d, g, tlife): return N0*(2*mu**d/gamma(d + 1))*g**(n - 1)*((D(N,N0,tlife, n, g) - 1.5)**d - ( (n - 1)*(g - 1) + 1 - 1.5)**d)
#probability of cancer in years in ages [in units of years]
def dP(N,N0,mu, n, d, g,ages,tlife,tdev):
    tmp=[P(N,N0,mu,n,d,g,age*(N/(tlife-tdev))) for age in ages]
    return array(tmp[1:]+[P(N,N0,mu,n,d,g,(ages[-1]+1)*(N/(tlife-tdev)))])-array(tmp)

def dP1(N,N0,mu, n, d, g,age,tlife,tdev):
    x = P(N,N0,mu,n,d,g,(age+1.)*(N/(tlife-tdev))) - P(N,N0,mu,n,d,g,age*(N/(tlife-tdev)))
    if x<=0:
        x=1
    return x

@ray.remote
def Pir(mu,n,d,gamma,times,s_scheme):
    s = 1. / (gamma-1.) * 1. / float(d)
    if s_scheme == "neutral":
        s=0
    return array([iMQ_ivp_mc(gamma, t, s, d, mu, int(n), s_scheme) for t in times])     

@ray.remote
def Pisr(mu,n,s,gamma,times,s_scheme):

    return array([iMQ_ivp(gamma, t, s, mu, int(n), s_scheme) for t in times])     

def Pi(mu,n,d,gamma,times,s_scheme):
    s = 1. / (gamma-1.) * 1. / float(d)
    if s_scheme == "neutral":
        s=0
    return array([iMQ_ivp_mc(gamma, t, s, d, mu, int(n), s_scheme) for t in times])     

def Pis(mu,n,s,gamma,times,s_scheme="scd"):
    return array([iMQ_ivp(gamma, t, s, mu, int(n), s_scheme) for t in times])     

def Pisd(mu,n,s,d,gamma,times,s_scheme="scd"):
    return array([iMQ_ivp_mc(gamma, t, s, d, mu, n, s_scheme) for t in times])     

def dPi(N,N0,mu,n,d,gamma,ages,tlife,tdev,s_scheme="scd"):
    times=[age*(N/N0/(tlife-tdev)) for age in list(ages)+[list(ages)[-1]+1] ]
    P = Pi(mu,n,d,gamma,times,s_scheme)
    dPs = (P[1:] - P[:-1])*N0          
    return dPs





@ray.remote
def run_chain(N,N0,mu,n,d,g,steps,lldata,logprior):
    ll=lldata(N,N0,mu,n,g,d)+logprior(N,N0,mu,n,g,d)

    #accepted moves counter
    a=0     
    #mean step sizes
    sigma=[0.01,0.1,1.0]
    chain={}
    for var in ['log_mu','n','d','g','ll']:
        chain[var]=[]
    for i in range(steps):
        reject=False

        #pick a move
        r=random.random()
        if r<1.5/9.:
            mumu=mu
            nn=n+random.normal(0,sigma[random.randint(0,2)])
            #nn=n+(random.randint(0,1)-0.5)*2
            gg=g
            dd=d
            NN=N
            N0N0=N0

        elif r<3./9.:
            mumu=mu
            nn=n
            gg=g
            dd=d+random.normal(0,sigma[random.randint(0,2)])
            NN=N
            N0N0=N0

        elif r<4.5/9.:
            mumu=mu
            nn=n
            gg=g+random.normal(0,sigma[random.randint(0,2)])
            dd=d
            NN=N
            N0N0=N0

        elif r<6/9.:
            mumu=exp(log(mu)+random.normal(0,sigma[random.randint(0,2)]))
            nn=n
            gg=g
            dd=d
            NN=N
            N0N0=N0

        elif r<7/9.:
            mumu=mu            
            lg=log(g)
            lgg=lg+random.normal(0,sigma[random.randint(0,2)])
            gg=exp(lgg)
            nn=1+lg/lgg*(n-1)
            dd=d
            NN=N
            N0N0=N0

        elif r<8/9.:
            dd=d+random.normal(0,sigma[random.randint(0,2)])

            mumu=mu**(d/dd)
            nn=n+(d - dd)
            gg=g
            NN=N
            N0N0=N0

        else:
            lmu=log(mu)
            lmumu=lmu+random.normal(0,sigma[random.randint(0,2)])
            mumu=exp(lmumu)
            C=g**(1-n)
            CC=C*exp( (d*(lmu - lmumu))/(d-1) )        
            gg=CC**(1/(1 - n))
            nn=n         
            dd=d
            NN=N
            N0N0=N0

        #boundary conditions / uninformative prior
        if nn<2 or gg<2 or gg>10 or log(mumu)/log(10)<-9 or dd<1 or n>30:
            reject=True

        #Metropolis-Hastings step
        if not reject:
            llprior=logprior(NN,N0N0,mumu,nn,gg,dd)
            llll=lldata(NN,N0N0,mumu,nn,gg,dd)+llprior
            if (llll>ll or exp(llll-ll)>random.random()):
                ll=llll
                mu=mumu
                n=nn
                d=dd
                g=gg
                NN=N
                N0=N0N0

                a+=1.
                lprior=llprior

        #output chain state every 1000 iterations    
        if i%100==0:
            chain["log_mu"]=append(chain["log_mu"],log(mu)/log(10))
            chain["n"]=append(chain["n"],n)
            chain["g"]=append(chain["g"],g)
            chain["d"]=append(chain["d"],d)        
            chain["ll"]=append(chain["ll"],ll)        

    return {"chain":chain,"q":[mu,n,d,g]}

@ray.remote
def run_mug_chain(N,N0,mu,g,steps,lldata,logprior):
    ll=lldata(N,N0,mu,g)+logprior(N,N0,mu,g)

    #accepted moves counter
    a=0    
    #mean step sizes
    sigma=[0.01,0.1,1.0]
    chain={}
    for var in ['log_mu','g','ll']:
        chain[var]=[]
    for i in range(steps):
        reject=False

        #pick a move
        r=random.random()
        if r<1/2.:
            mumu=mu
            gg=2.+exp(log(g-2.)+random.normal(0,sigma[random.randint(0,2)]))
            NN=N
            N0N0=N0

        else:
            mumu=exp(log(mu)+random.normal(0,sigma[random.randint(0,2)]))
            gg=g
            NN=N
            N0N0=N0


        #boundary conditions / uninformative prior
        if  gg<2 or gg>10 or log(mumu)/log(10)<-9 or log(mumu)/log(10)>-2:
            reject=True

        #Metropolis-Hastings step
        if not reject:
            llll=lldata(NN,N0N0,mumu,gg)+logprior(NN,N0N0,mumu,gg)
            if (llll>ll or exp(llll-ll)>random.random()):
                ll=llll
                mu=mumu
                g=gg
                NN=N
                N0=N0N0

                a+=1.

        #output chain state every 1000 iterations    
        if i%1==0:
            chain["log_mu"]=append(chain["log_mu"],log(mu)/log(10))
            chain["g"]=append(chain["g"],g)
            chain["ll"]=append(chain["ll"],ll)        

    return {"chain":chain,"q":[mu,g]}

