
TITLE Stochastic Hodgkin and Huxley model incorporating channel noise (microscopic version).

COMMENT

This mod-file implementes a stochastic version of the HH model incorporating channel noise.
This version is the ``microscopic'' version, i.e. it employs a Markov model for the simulation
of the open-close kinetics of ion channels.

Author: Daniele Linaro - daniele.linaro@unige.it
Date: September 2010

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S)  = (siemens)
    (pS) = (picosiemens)
    (um) = (micron)
} : end UNITS


NEURON {
    SUFFIX HHmicro
    USEION na READ ena WRITE ina
    USEION k READ ek WRITE ik
    NONSPECIFIC_CURRENT il
    RANGE el
    RANGE gnabar, gkbar, gna, gk, gl
    RANGE m0h0,m1h0,m2h0,m3h0,m0h1,m1h1,m2h1,m3h1
    RANGE n0, n1, n2, n3, n4
    RANGE gamma_na, gamma_k
    RANGE Nna, Nk
    GLOBAL m_inf, h_inf, n_inf, tau_m, tau_h, tau_n
    RANGE seed
    THREADSAFE : assigned GLOBALs will be per thread
} : end NEURON


PARAMETER {
    gnabar = 0.12   (S/cm2)     : maximum sodium conductance
    gkbar  = 0.036  (S/cm2)     : maximum potassium conductance
    gl     = 0.0003 (S/cm2)     : maximum leakage conductance
    el       = -54.3 (mV)       : leakage reversal potential
    
    gamma_na = 10  (pS)		: single channel sodium conductance
    gamma_k  = 10  (pS)		: single channel potassium conductance
    seed = 5061983              : always use the same seed
} : end PARAMETER


STATE {
    m h n
} : end STATE


ASSIGNED {
    ina   (mA/cm2)
    ik    (mA/cm2)
    il    (mA/cm2)
    gna   (S/cm2)
    gk    (S/cm2)
    ena  (mV)
    ek   (mV)
    
    dt    (ms)
    area  (um2)
    celsius  (degC)
    v  (mV)
    
    Nna			 	: total number of sodium channels
    Nk			 	: total number of potassium channels

    m0h0			: inactivated state (sodium channels)
    m1h0			: inactivated state (sodium channels)
    m2h0			: inactivated state (sodium channels)
    m3h0			: inactivated state (sodium channels)
    m0h1			: closed state (sodium channels)
    m1h1			: closed state (sodium channels)
    m2h1			: closed state (sodium channels)
    m3h1			: open state (sodium channels)

    n1				: closed state (potassium channels)
    n2				: closed state (potassium channels)
    n3				: closed state (potassium channels)
    n4				: closed state (potassium channels)
    n5				: open state (potassium channels)
    
    m_inf h_inf n_inf
    tau_m (ms) tau_h (ms) tau_n (ms)
    
} : end ASSIGNED



INITIAL {
    m = 0
    h = 0
    n = 0
    
    Nna = ceil(((1e-8)*area)*(gnabar)/((1e-12)*gamma_na))   : area in um2 -> 1e-8*area in cm2; gnabar in S/cm2; gamma_na in pS -> 1e-12*gamma_na in S
    Nk = ceil(((1e-8)*area)*(gkbar)/((1e-12)*gamma_k))   : area in um2 -> 1e-8*area in cm2; gkbar in S/cm2; gamma_k in pS -> 1e-12*gamma_k in S
    
    VERBATIM
		/*
		fprintf(stdout, "HHmicro>> ");
    fprintf(stdout, "gamma_na = %.0f gamma_k = %.0f ", gamma_na, gamma_k);
    fprintf(stdout, "Nna = %.0f Nk = %.0f\n", Nna, Nk);
    fflush(stdout);
		*/
    ENDVERBATIM

    m0h0 = 0
    m1h0 = 0
    m2h0 = 0
    m3h0 = 0
    m0h1 = Nna		: therefore you should wait at the beginning of the simulation until the Na channels have reached steady state.
    m1h1 = 0
    m2h1 = 0
    m3h1 = 0

    n1   = Nk		: therefore you should wait at the beginning of the simulation until the K channels have reached steady state.
    n2   = 0
    n3   = 0
    n4   = 0
    n5   = 0
    
    set_seed(seed)
    
} : end INITIAL


BREAKPOINT {
    SOLVE states METHOD cnexp
    gna = m3h1*((1e-12)*gamma_na)/((1e-8)*area)
    gk = n5*((1e-12)*gamma_k)/((1e-8)*area)
    ina = gna * (v-ena)
    ik = gk * (v-ek)
    il = gl * (v-el)
} : end BREAKPOINT


DERIVATIVE states {   
    UNITSOFF
    m' = m
    h' = h
    n' = n
    UNITSON
    noise()
} : end DERIVATIVE


FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 

FUNCTION alpham(Vm (mV)) (/ms) {
    UNITSOFF
    alpham = .1 * vtrap(-(Vm+40),10)
    UNITSON
}


FUNCTION betam(Vm (mV)) (/ms) {
    UNITSOFF
    betam =  4 * exp(-(Vm+65)/18)
    UNITSON
}


FUNCTION alphah(Vm (mV)) (/ms) {
    UNITSOFF
    alphah = .07 * exp(-(Vm+65)/20)
    UNITSON
}


FUNCTION betah(Vm (mV)) (/ms) {
    UNITSOFF
    betah = 1 / (exp(-(Vm+35)/10) + 1)
    UNITSON
}


FUNCTION alphan(Vm (mV)) (/ms) {
    UNITSOFF
    alphan = .01*vtrap(-(Vm+55),10) 
    UNITSON
}


FUNCTION betan(Vm (mV)) (/ms) {
    UNITSOFF
    betan = .125*exp(-(Vm+65)/80)
    UNITSON
}


PROCEDURE noise() {
    LOCAL p,am,bm,ah,bh,an,bn,m0h0_merk,m1h0_merk,m2h0_merk,m3h0_merk,m0h1_merk,m1h1_merk,m2h1_merk,m3h1_merk,n1_merk,n2_merk,n3_merk,n4_merk,n5_merk,rnd,den
    
    : alpha_m and beta_m
    am = alpham(v)
    bm = betam(v)
    
    : alpha_h and beta_h
    ah = alphah(v)
    bh = betah(v)
    
    : alpha_n and beta_n
    an = alphan(v)
    bn = betan(v)

    m0h0_merk = m0h0
    m1h0_merk = m1h0
    m2h0_merk = m2h0
    m3h0_merk = m3h0
    m0h1_merk = m0h1
    m1h1_merk = m1h1
    m2h1_merk = m2h1
    m3h1_merk = m3h1

    n1_merk = n1
    n2_merk = n2
    n3_merk = n3
    n4_merk = n4
    n5_merk = n5
    
    : ------------- h0 -------------
    
    :m0h0
    den = ah+3*am
    p = 1-exp(-dt*den)
    FROM ii=1 TO m0h0_merk {
	:scop_random gives random number uniform between 0 and 1
	if (scop_random() <= p)	{					:probability that a channel in the state m0h0 goes to state m0h1 or m1h0
	    if (scop_random() <= ah/den) {       		        :probability that this channel goes to state m0h1 (via rate ah)
		m0h0=m0h0-1
		m0h1=m0h1+1
	    }
	    else {							:otherwise this channel goes to m1h0 (via rate 3*am)
		m0h0=m0h0-1
		m1h0=m1h0+1
	    }
	}
    }
    
    :m1h0
    den = ah+2*am+bm
    p = 1-exp(-dt*den)
    FROM ii=1 TO m1h0_merk {
	if (scop_random() <= p) {
	    rnd = scop_random()
	    if (rnd <= ah/den) {	
		m1h0=m1h0-1
		m1h1=m1h1+1
	    }
	    else if (rnd <= (ah+2*am)/den) {
		m1h0=m1h0-1
		m2h0=m2h0+1
	    }
	    else {
		m1h0=m1h0-1
		m0h0=m0h0+1
	    }
	}
    }
    
    :m2h0
    den = ah+am+2*bm
    p = 1-exp(-dt*den) 
    FROM ii=1 TO m2h0_merk {
	if (scop_random() <= p){
	    rnd = scop_random()
	    if (rnd <= ah/den)	{
		m2h0=m2h0-1
		m2h1=m2h1+1
	    }
	    else if (rnd <= (ah+am)/den) {
		m2h0=m2h0-1
		m3h0=m3h0+1
	    }
	    else {
		m2h0=m2h0-1
		m1h0=m1h0+1
	    }
	}
    }
    
    :m3h0
    den = ah+3*bm
    p = 1-exp(-dt*den)
    FROM ii=1 TO m3h0_merk {
	if (scop_random() <= p) {
	    rnd = scop_random()
	    if (rnd <= ah/den) {	
		m3h0=m3h0-1
		m3h1=m3h1+1
	    }
	    else {	
		m3h0=m3h0-1
		m2h0=m2h0+1
	    }
	}
    }
    
    
    : ------------- h1 -------------
    
    :m0h1
    den = bh+3*am
    p=1-exp(-dt*den) 
    FROM ii=1 TO m0h1_merk {
	if (scop_random() <= p) {
	    rnd = scop_random()
	    if (rnd <= bh/den)	{	
		m0h1=m0h1-1
		m0h0=m0h0+1
	    }
	    else {
		m0h1=m0h1-1
		m1h1=m1h1+1
	    }
	}
    }
    
    :m1h1
    den = bh+2*am+bm
    p = 1-exp(-dt*den) 
    FROM ii=1 TO m1h1_merk {
	if (scop_random() <= p) {
	    rnd = scop_random()
	    if (rnd <= bh/den)	{	
		m1h1=m1h1-1
		m1h0=m1h0+1
	    }
	    else if (rnd <= (bh+2*am)/den) {	
		m1h1=m1h1-1
		m2h1=m2h1+1
	    }
	    else {
		m1h1=m1h1-1
		m0h1=m0h1+1
	    }
	}
    }
    
    :m2h1
    den = bh+am+2*bm
    p = 1-exp(-dt*den) 
    FROM ii=1 TO m2h1_merk {
	if (scop_random()<= p){
	    rnd = scop_random()
	    if (rnd <= bh/den)	{	
		m2h1=m2h1-1
		m2h0=m2h0+1
	    }
	    else if (rnd <= (bh+am)/den) {	
		m2h1=m2h1-1
		m3h1=m3h1+1
	    }
	    else {
		m2h1=m2h1-1
		m1h1=m1h1+1
	    }
	}
    }
    
    :m3h1
    den = bh+3*bm
    p = 1-exp(-dt*den) 
    FROM ii=1 TO m3h1_merk {
	if (scop_random()<= p) {
	    rnd = scop_random()
	    if (rnd <= bh/den)	{	
		m3h1=m3h1-1
		m3h0=m3h0+1
	    }	
	    else {
		m3h1=m3h1-1
		m2h1=m2h1+1
	    }
	}
    }
    
    
    : ------------- n -------------
    
    :n1
    den = 4*an
    p = 1-exp(-dt*den) 
    FROM ii=1 TO n1_merk {
	if (scop_random() <= p) {
	    n1=n1-1
	    n2=n2+1	
	}
    }
    
    
    :n2
    den = 3*an+bn
    p = 1-exp(-dt*den) 
    FROM ii=1 TO n2_merk {
	if (scop_random() <= p) {
	    if (scop_random() <= 3*an/den) {	
		n2=n2-1
		n3=n3+1
	    }
	    else {
		n2=n2-1
		n1=n1+1
	    }
	}
    }
    
    
    :n3
    den = 2*an+2*bn
    p = 1-exp(-dt*den)
    FROM ii=1 TO n3_merk {
	if (scop_random() <= p) {
	    if (scop_random() <= 2*an/den) {
		n3=n3-1
		n4=n4+1
	    }
	    else {
		n3=n3-1
		n2=n2+1
	    }
	}
    }
    
    
    :n4
    den = an+3*bn
    p = 1-exp(-dt*den)
    FROM ii=1 TO n4_merk {
	if (scop_random() <= p) {
	    if (scop_random() <= an/den) {
		n4=n4-1
		n5=n5+1
	    }
	    else {
		n4=n4-1
		n3=n3+1
	    }
	}
    }
    
    
    :n5
    den = 4*bn
    p = 1-exp(-dt*den)
    FROM ii=1 TO n5_merk {
	if (scop_random() <= p) {
	    n5=n5-1
	    n4=n4+1	
	}
    }
}
