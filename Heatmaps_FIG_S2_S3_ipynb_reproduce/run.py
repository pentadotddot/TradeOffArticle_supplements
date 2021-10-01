import ray
from numpy import *
from scipy.stats import poisson
from scipy.integrate import odeint, solve_ivp
######################Pmum#######################
from scipy.special import gamma as Gamma
def Di(tlife, n, g): return tlife * power(g,-(n - 1)) + (n - 1)*(g - 1) + 1 

def P(mu, n, d, g, tlife): return (2*mu**d/Gamma(d + 1))*g**(n - 1)*((Di(tlife, n, g) - 1.5)**d - ( (n - 1)*(g - 1) + 1 - 1.5)**d)

def cPmum(gamma, s, mu, tlife, n): return P(mu, n,  1. / s / (gamma-1.), gamma, tlife)
    
def Pmum(gamma, s, mu, tlife, n): return P(mu, n,  ceil(1. / s / (gamma-1.)), gamma, tlife)
#################################################

######################Q##########################
scalar_type = float64
int_type = int64

zero = scalar_type(0)
one = scalar_type(1)
two = scalar_type(2)

def one_m_exp_m(x): return x if x<1e-10 else one-exp(-x)

def gamma_sx(gamma,s,mc,x):
    return gamma / max(( 1 - s*(mc-x) * (gamma-1) ),finfo(scalar_type).tiny*1e6)
    #return (gamma * (gamma-1) / x)
def Qs(gamma, s, mc, tlife, mu, n, s_scheme, Nk):
    if s_scheme == "scd":
        gamma_s = array([nan]+[ gamma_sx(gamma,s,mc,x) for x in range(1,mc+1) ], dtype=scalar_type) #m = mc-x
    elif s_scheme == "neutral":
        gamma_s = array([nan]+[ gamma for x in range(1,mc+1) ], dtype=scalar_type) #m = mc-x
    else:
        gamma_s = array([nan]+[ (gamma - 2*s*(mc-x) * (gamma-1) ) / ( 1 - s*(mc-x) * (gamma-1) ) for x in range(1,mc+1) ], dtype=scalar_type) #m = mc-x
    #gamma_s = [None] + gamma_s

    tcorr=tlife
    #def W(x,k): return None if x==0 else power(gamma,k-n)*( one - s*(mc-x) * (gamma-one) )*tcorr/Nk[k]
    def W(x,k): return None if x==0 else power(gamma,k-n)*gamma_sx(gamma,s,mc,x)*tcorr/Nk[k]
        
    Q = array( [ [scalar_type(0.0) for k in range(n+1)] for x in range(mc+1) ], dtype=scalar_type )
    Q[:,n] = 0
    Q[0,:] = 1
    for x in range(1,mc+1):
        for k in range(0,n)[::-1]:
            Q[x, k] = one_m_exp_m( ( 
                mu * (gamma_s[x]-two) * Q[x-1,k  ] +
                mu *  gamma_s[x]      * Q[x-1,k+1] +
                      gamma_s[x]      * Q[x  ,k+1]           
                ) * one_m_exp_m(W(x,k))  )
    return Q

def PQ(gamma, tlife, s, mu, n, s_scheme="neutral", Nk = ones(0), rng = random.default_rng()):
    if Nk.size==0:
        Nk=ones(n)
        
    mc = int(ceil( 1. / s / (gamma-1.) ))

    Q = Qs(gamma, s, mc, tlife, mu, n, s_scheme, Nk)
                
    lambda0 = power( gamma, one-n ) * mu * tlife
    
    Pc = poisson.sf(mc-1,double(lambda0))

    Cm = array([poisson.sf(m,double(lambda0)) for m in range(mc)]) / mu

    Mm = array([ Cm[0] + poisson.pmf(0,double(lambda0)) - one if m==0 else Cm[m] + poisson.pmf(m,double(lambda0)) for m in range(mc)])

    Qm = array([ Q[mc-m,1] for m in range(mc) ])
              
    return two*Pc + sum( Mm * Qm ) + sum( Nk[1:n]*Q[mc,1:n] )

def iMQ_ivp_times(gamma, tlife, s, mu, n, times, s_scheme="scd", Nk = ones(0)):
    if Nk.size==0:
        Nk=ones(n)

    mc = int(ceil( 1. / s / (gamma-1.) ))
    Q = Qs(gamma, s, mc, tlife, mu, n, s_scheme, Nk)
    
    p0 = zeros(mc+1)
    p0[0] = one
    
    t = zero
    pi = zero

    def dpdt(t,p,gamma,s,mc,n):
        Q = Qs(gamma, s, mc, t, mu, n, s_scheme, Nk)
        r = power( gamma, one-n )
        #    d_m = r*( mu + mu * Q[mc-m-1,1] + (1-mu) * Q[mc-m,1] )
        d = zeros(mc, dtype=scalar_type)
        di = zeros(mc, dtype=scalar_type)
        for m in range(mc):
            d[m] = r * ( mu + mu * Q[mc-m-1,1] + (one-mu) * Q[mc-m,1] )
            di[m] = r * ( mu * Q[mc-m-1,1] + (one-mu) * Q[mc-m,1] )
        di[mc-1] += r * mu
        dp = zeros(p.size)
        dp[:mc] = -d * p[:mc]
        dp[1:mc] += d[:mc-1] * p[:mc-1]
        dp[mc] = sum(p[:mc] * di)
        return dp
        
    p = solve_ivp(dpdt,[0,tlife],p0,args=(gamma,s,mc,n),t_eval=times,method="Radau")
    return p.y[-1] + sum(Nk * Q[mc,1:])


def iMQ_ivp(gamma, tlife, s, mu, n, s_scheme="scd", Nk = ones(0)):
    mc = int(ceil( 1. / s / (gamma-1.) ))
    return iMQ_ivp_mc(gamma, tlife, s, mc, mu, n, s_scheme)

def iMQ_ivp_mc(gamma, tlife, s, mc, mu, n, s_scheme="scd", Nk = ones(0)):
    if Nk.size==0:
        Nk=ones(n)

    
    Q = Qs(gamma, s, mc, tlife, mu, n, s_scheme, Nk)
    di = zeros(mc, dtype=scalar_type)
    
    r = power( gamma, one-n )
#    d_m = r*( mu + mu * Q[mc-m-1,1] + (1-mu) * Q[mc-m,1] )
    d = zeros(mc, dtype=scalar_type)
    for m in range(mc):
        d[m] = r * ( mu + mu * Q[mc-m-1,1] + (one-mu) * Q[mc-m,1] )
        di[m] = r * ( mu * Q[mc-m-1,1] + (one-mu) * Q[mc-m,1] )
    di[mc-1] += r * mu
    
    p0 = zeros(mc+1)
    p0[0] = one
    
    t = zero
    pi = zero

    def dpdt(t,p,d,di):
        dp = zeros(p.size)
        dp[:mc] = -d * p[:mc]
        dp[1:mc] += d[:mc-1] * p[:mc-1]
        dp[mc] = sum(p[:mc] * di)
        return dp
        
    p = solve_ivp(dpdt,[0,tlife],p0,args=(d,di),method="Radau")
    return p.y[-1][-1] + sum(Nk * Q[mc,1:])


def iMQ_odeint(gamma, tlife, s, mu, n, s_scheme="scd", Nk = ones(0)):
    if Nk.size==0:
        Nk=ones(n)

    mc = int(ceil( 1. / s / (gamma-1.) ))
    
    Q = Qs(gamma, s, mc, tlife, mu, n, s_scheme, Nk)
    di = zeros(mc, dtype=scalar_type)
    
    r = power( gamma, one-n )
#    d_m = r*( mu + mu * Q[mc-m-1,1] + (1-mu) * Q[mc-m,1] )
    d = zeros(mc, dtype=scalar_type)
    for m in range(mc):
        d[m] = r * ( mu + mu * Q[mc-m-1,1] + (one-mu) * Q[mc-m,1] )
        di[m] = r * ( mu * Q[mc-m-1,1] + (one-mu) * Q[mc-m,1] )
    for m in range(mc):
        d[m] = r * ( mu + mu * Q[mc-m-1,1] + (one-mu) * Q[mc-m,1] )
        di[m] = r * ( mu * Q[mc-m-1,1] + (one-mu) * Q[mc-m,1] )
    di[mc-1] += r * mu
    
    p0 = zeros(mc+1)
    p0[0] = one
    
    t = zero
    pi = zero

    def dpdt(p,t,d,di):
        dp = zeros(p.size)
        dp[:mc] = -d * p[:mc]
        dp[1:mc] += d[:mc-1] * p[:mc-1]
        dp[mc] = sum(p[:mc] * di)
        return dp
        
    p = odeint(dpdt,p0,[0,tlife],args=(d,di))
    return p[1][mc] + sum(Nk * Q[mc,1:])


def iMQ_euler(gamma, tlife, s, mu, n, s_scheme="scd", delta=1., Nk = ones(0)):
    if Nk.size==0:
        Nk=ones(n)

    mc = int(ceil( 1. / s / (gamma-1.) ))
    
    Q = Qs(gamma, s, mc, tlife, mu, n, s_scheme, Nk)
    di = zeros(mc, dtype=scalar_type)
    
    r = power( gamma, one-n )
    d = zeros(mc, dtype=scalar_type)
    for m in range(mc):
        d[m] = r * ( mu + mu * Q[mc-m-1,1] + (one-mu) * Q[mc-m,1] )
        di[m] = r * ( mu * Q[mc-m-1,1] + (one-mu) * Q[mc-m,1] )
    for m in range(mc):
        d[m] = r * ( mu + mu * Q[mc-m-1,1] + (one-mu) * Q[mc-m,1] )
        di[m] = r * ( mu * Q[mc-m-1,1] + (one-mu) * Q[mc-m,1] )
    di[mc-1] += r * mu
    
    p=zeros(mc)
    dp=zeros(mc)
    p[0]=one
    delta_t=delta*one/r
    
    t=zero
    pi=zero
    while t<tlife:
        dp = -d * p
        dp[1:] += d[:mc-1] * p[:mc-1]
        pi += sum(p * di * delta_t)
        p = p + dp * delta_t        
        t += delta_t

    return pi + sum(Nk * Q[mc,1:])

def PMQ(gamma, tlife, s, mu, n, s_scheme="scd", epsilon=1e-1, Nk = ones(0), rng = random.default_rng()):
    if Nk.size==0:
        Nk=ones(n)

    mc = int(ceil( 1. / s / (gamma-1.) ))
    
    Q = Qs(gamma, s, mc, tlife, mu, n, s_scheme, Nk)
    
    r = power( gamma, one-n )
#    d_m = r*( mu + mu * Q[mc-m-1,1] + (1-mu) * Q[mc-m,1] )
    d = zeros(mc, dtype=scalar_type)
    for m in range(mc):
        d[m] = r * ( mu + mu * Q[mc-m-1,1] + (one-mu) * Q[mc-m,1] )
    for m in range(mc):
        d[m] = r * ( mu + mu * Q[mc-m-1,1] + (one-mu) * Q[mc-m,1] ) #+ (1.-rng.random())*epsilon*mean(d)
    for i in range(m+1):
        for j in range(m+1):
            if i!=j and min(d[i],d[j])*epsilon > abs(d[i]-d[j]):
                #print(i,j,d[i],d[j],d[i]-d[j])
                d[i]+= (1.-rng.random())*epsilon*min(d[i],d[j])
                d[j]+= (1.-rng.random())*epsilon*min(d[i],d[j])
                #print(i,j,d[i],d[j],d[i]-d[j])

#   C_m = r*(r*mu)^m * sum_{i=0}^{m} (1-exp(-d_i*t))/d_i / prod_{j=0,j!=i}^{m} (d_j-d_i)
    C = zeros(mc+1, dtype=scalar_type)
    for m in range(mc):
        positive_Cs=[]
        negative_Cs=[]
        prods=[]
        for i in range(m+1):
            prod = one
            for j in range(m+1):
                if i!=j:
                    #print(r * mu * Q[mc-j-1,1] , r * mu * Q[mc-i-1,1] , r * (one-mu) * Q[mc-j,1] , r * (one-mu) * Q[mc-i,1],r * mu * Q[mc-j-1,1] - r * mu * Q[mc-i-1,1] + r * (one-mu) * Q[mc-j,1] - r * (one-mu) * Q[mc-i,1])
                    prod *= d[j]- d[i] #- one #r * mu * Q[mc-j-1,1] - r * mu * Q[mc-i-1,1] + r * (one-mu) * Q[mc-j,1] - r * (one-mu) * Q[mc-i,1]
            #if prod!=1:
            #    prod *= power(d[i],m-1)
            C_term = r * power( r * mu , m ) * one_m_exp_m(d[i]*tlife) / d[i] / prod 
            prods+=[prod]
            if C_term>0:
                positive_Cs+=[C_term]
            else:
                negative_Cs+=[C_term]
                
            #print(one_m_exp_m(d[i]*tlife),d[i],prod, r * power( r * mu , m ) * one_m_exp_m(d[i]*tlife) / d[i] / prod)
        Csum = 0
        negative_Cs = sorted(negative_Cs)
        positive_Cs = sorted(positive_Cs)[::-1]
#        print(negative_Cs)
#        print(positive_Cs)
#        print(r * power( r * mu , m )*tlife/array(prods))
        Cs=[]
        for i in range(max(len(positive_Cs),len(negative_Cs))):
            Cs+=[0]
            if i<len(positive_Cs):
                Cs[i] += positive_Cs[i]
            if i<len(negative_Cs):
                Cs[i] += negative_Cs[i]
 #       print(Cs,"<<<",sum(Cs))
        C[m] = sum(Cs)

    M = zeros(mc+1, dtype=scalar_type)
#   M_0 = (1-mu)*C_0    
    M[0] = (one - mu) * C[0]
#   M_m = (1-mu)*C_m + mu*C_{m-1}
    for m in range(1,mc+1):
        M[m] = (one - mu) * C[m] + mu * C[m-1]
#   M_{mc} = mu*C_{mc-1}
    M[mc] =  mu * C[mc-1]
#   2*M_{mc} + sum_{m=0}^{mc-1} Q_{mc-m} * M_m
    PMQ = 2 * M[mc]
    for m in range(mc):
        PMQ += Q[mc-m,1] * M[m]
    return PMQ + sum(Nk * Q[mc,1:])
#################################################

#################merge###########################
scalar_type = float 

def merge(PA, PB, poi, mc):
    PAp = zeros(mc+1)
    PBp = zeros(mc+1)
    for k in range(mc+1): #eqs. 1-3
        PAp[k] = inner(PA[:k+1][::-1] , poi[:k+1])
        PBp[k] = inner(PB[:k+1][::-1] , poi[:k+1])
        
    PAB = zeros(mc+1)
    PAB[1:] = PAp[1:] * cumsum(PBp)[:mc] + PBp[1:] * cumsum(PAp)[:mc] + PAp[1:] * PBp[1:]  #eq. 4
    PAB[0] =  PAp[0] * PBp[0]
    return PAB

def PT0(mc):#inital condition at leaves of the cell lineage tree
    PT_0=zeros(mc+1, dtype=scalar_type)
    PT_0[0]=1
    return PT_0    
    
def PT(t,k,n,rates,tlife,mc,poi, rng):
    if k == n:
        return PT0(mc)

    t += 1./rates["R"][k]#rng.exponential(1./rates["R"][k])
    if t > tlife:
        return PT0(mc)
    
    if rates["scd"][k] > rng.random() * rates["R"][k]:
        PA = PT(t, k, n, rates, tlife, mc, poi, rng) 
        PB = PT(t, k, n, rates, tlife, mc, poi, rng)
    else:
        PA = PT(t, k+1, n, rates, tlife, mc, poi, rng) 
        PB = PT(t, k+1, n, rates, tlife, mc, poi, rng)
    PAp = zeros(mc+1)
    PBp = zeros(mc+1)
    return merge(PA, PB, poi, mc)

def PT0_flat(t,rates, n, tlife, mc, poi,rng):
    PAs=[]
    while t < tlife:
        PA=PT(t, 0, n, rates, tlife, mc, poi, rng) 
        PAs+=[PA]
        t += rng.exponential(1./rates["scdd"][0])
    PAB = PT0(mc)
    for PA in PAs:
        PAB = merge(PA, PAB, poi, mc)
    return PAB[mc] #eq. 5

#truncated PT0_flat
def PT0_flat_K(t,rates,tlife, mc, poi, K=100):
    if K>int(tlife*rates["scdd"][0]):
        K=int(ceil(tlife*rates["scdd"][0]))
    Krange=range(K)
    PAs=[PT(t, 0, rates, tlife, mc, poi) for i in Krange]
    PAB = PT0(mc)
    while t < tlife:
        PAB = merge(PAs[random.choice(Krange)], PAB, poi, mc)
        t += rng.exponential(1./rates["scdd"][0])
    return PAB[mc]#eq. 5

def T(gamma, n, tlife, Nk, s, mu, F=PT0_flat, rng=random.default_rng()):#rates and precalc.

    mc = int(ceil( 1. / s / (gamma-1.) ))
    
    poi = array([ poisson.pmf(m, mu) for m in range(mc+2)])
        
    delta = array([power( gamma, 1 - (n-k)  ) for k in range(n)])
    
    p = 1 # for now ..
    q = 2./gamma
    rates={}
    rates["scd"] = array([0 if k==0 or k==n else 0.5 * delta[k] / Nk[k] * (1-q) for k in range(n) ])
    rates["scdd"]= array([ delta[k] / Nk[k] if k==0 else 0.5 * delta[k] / Nk[k] for k in range(n) ])
    rates["R"] = rates["scd"] + rates["scdd"]
    
    return F(0,rates,n,tlife,mc,poi,rng)
#################################################

#################ray#############################
#utility functions for ray
@ray.remote
def PTray(t,k,rates,tlife,mc,poi):
    return PT(t,k,rates,tlife,mc,poi)
@ray.remote
def Tray(gamma, n, tlife, Nk, s, mu):
    return T(gamma, n, tlife, Nk, s, mu, PT0_flat_K)

#per progenitor subtree parallelization of PT0_flat
def PT0_flat_ray(t,rates,tlife, mc, poi):
    PAs=[]
    while t < tlife:
        PA=PTray.remote(t, 0, rates, tlife, mc, poi) 
        PAs+=[PA]
        t += rng.exponential(1./rates["scdd"][0])
    PAs = ray.get(PAs)
    PAB = PT0(mc)
    for PA in PAs:
        PAB = merge(PA, PAB, poi, mc)
    return PAB[mc]#eq. 5

#truncated per progenitor subtree parallelization of PT0_flat
def PT0_flat_ray_K(t,rates,tlife, mc, poi, K=100):
    Krange=range(K)
    PAs=ray.get([PTray.remote(t, 0, rates, tlife, mc, poi) for i in Krange])
    PAB = PT0(mc)
    while t < tlife:
        PAB = merge(PAs[random.choice(Krange)], PAB, poi, mc)
        t += rng.exponential(1./rates["scdd"][0])
    return PAB[mc]#eq. 5

#################################################

