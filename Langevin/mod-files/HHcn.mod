
TITLE Stochastic Hodgkin and Huxley model incorporating channel noise (effective version).

COMMENT

This mod-file implementes a stochastic version of the HH model incorporating channel noise.
This version is the ``effective'' version, i.e. it employs a Ornstein-Uhlenbeck process with
the same stochastic properties (mean and covariance) as the Markov model.

Author: Daniele Linaro - daniele.linaro@unige.it
Date: September 2010

ENDCOMMENT

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S)  = (siemens)
    (pS) = (picosiemens)
    (um) = (micrometer)
} : end UNITS


NEURON {
    SUFFIX HHcn
    USEION na READ ena WRITE ina
    USEION k READ ek WRITE ik
    NONSPECIFIC_CURRENT il
    RANGE el, gl
    RANGE gnabar, gkbar
    RANGE gna, gk
    RANGE gamma_na, gamma_k
    RANGE m_inf, h_inf, n_inf
    RANGE tau_m, tau_h, tau_n
    RANGE tau_y1, tau_y2, tau_y3, tau_y4
    RANGE tau_z1, tau_z2, tau_z3, tau_z4, tau_z5, tau_z6, tau_z7
    RANGE var_y1, var_y2, var_y3, var_y4
    RANGE var_z1, var_z2, var_z3, var_z4, var_z5, var_z6, var_z7
    RANGE noise_y1, noise_y2, noise_y3, noise_y4
    RANGE noise_z1, noise_z2, noise_z3, noise_z4, noise_z5, noise_z6, noise_z7
    RANGE Nna, Nk
    RANGE seed    
    : these auxiliary variables are only needed if the exact method of solution is used    
    RANGE mu_y1, mu_y2, mu_y3, mu_y4
    RANGE mu_z1, mu_z2, mu_z3, mu_z4, mu_z5, mu_z6, mu_z7
    THREADSAFE
} : end NEURON


PARAMETER {
    gnabar  = 0.12   (S/cm2)
    gkbar   = 0.036  (S/cm2)
    gl      = 0.0003 (S/cm2)  
    el       = -54.3 (mV)       : leakage reversal potential
    
    gamma_na = 10  (pS)		: single channel sodium conductance
    gamma_k  = 10  (pS)		: single channel potassium conductance
    seed = 5061983              : always use the same seed
} : end PARAMETER


STATE {
    m h n
    y1 y2 y3 y4
    z1 z2 z3 z4 z5 z6 z7
} : end STATE


ASSIGNED {
    ina        (mA/cm2)
    ik         (mA/cm2)
    il         (mA/cm2)
    gna        (S/cm2)
    gk         (S/cm2)
    ena        (mV)
    ek         (mV)
    
    dt         (ms)
    area       (um2)
    celsius    (degC)
    v          (mV)
        
    Nna  (1) : number of potassium channels
    Nk   (1) : number of sodium channels
    
    m_inf h_inf n_inf
    noise_y1 noise_y2 noise_y3 noise_y4
    noise_z1 noise_z2 noise_z3 noise_z4 noise_z5 noise_z6 noise_z7
    var_y1 (ms2) var_y2 (ms2) var_y3 (ms2) var_y4 (ms2) 
    var_z1 (ms2) var_z2 (ms2) var_z3 (ms2) var_z4 (ms2) var_z5 (ms2) var_z6 (ms2) var_z7 (ms2)
    tau_m (ms) tau_h (ms) tau_n (ms)
    tau_y1 (ms) tau_y2 (ms) tau_y3 (ms) tau_y4 (ms) 
    tau_z1 (ms) tau_z2 (ms) tau_z3 (ms) tau_z4 (ms) tau_z5 (ms) tau_z6 (ms) tau_z7 (ms)
    : these auxiliary variables are only needed if the exact method of solution is used
    mu_y1 mu_y2 mu_y3 mu_y4
    mu_z1 mu_z2 mu_z3 mu_z4 mu_z5 mu_z6 mu_z7
} : end ASSIGNED

INITIAL {
    Nna = ceil(((1e-8)*area)*(gnabar)/((1e-12)*gamma_na))   : area in um2 -> 1e-8*area in cm2; gnabar in S/cm2; gamma_na in pS -> 1e-12*gamma_na in S
    Nk = ceil(((1e-8)*area)*(gkbar)/((1e-12)*gamma_k))   : area in um2 -> 1e-8*area in cm2; gkbar in S/cm2; gamma_k in pS -> 1e-12*gamma_k in S
    
    rates(v)
    m = m_inf
    h = h_inf
    n = n_inf
    y1 = 0.
    y2 = 0.
    y3 = 0.
    y4 = 0.
    z1 = 0.
    z2 = 0.
    z3 = 0.
    z4 = 0.
    z5 = 0.
    z6 = 0.
    z7 = 0.
    set_seed(seed)
} : end INITIAL


BREAKPOINT {
    SOLVE states
    gna = gnabar * (m*m*m*h + z1+z2+z3+z4+z5+z6+z7)
    if (gna < 0) {
	gna = 0
    }
    gk = gkbar * (n*n*n*n + y1+y2+y3+y4)
    if (gk < 0) {
	gk = 0
    }
    ina = gna * (v - ena)
    ik  = gk * (v - ek)
    il  = gl * (v - el)
} : end BREAKPOINT


PROCEDURE states() {
    rates(v)
    m = m + dt * (m_inf-m)/tau_m
    h = h + dt * (h_inf-h)/tau_h
    n = n + dt * (n_inf-n)/tau_n
    : Exact
    y1 = y1*mu_y1 + noise_y1
    y2 = y2*mu_y2 + noise_y2
    y3 = y3*mu_y3 + noise_y3
    y4 = y4*mu_y4 + noise_y4
    z1 = z1*mu_z1 + noise_z1
    z2 = z2*mu_z2 + noise_z2
    z3 = z3*mu_z3 + noise_z3
    z4 = z4*mu_z4 + noise_z4
    z5 = z5*mu_z5 + noise_z5
    z6 = z6*mu_z6 + noise_z6
    z7 = z7*mu_z7 + noise_z7
    : Euler-Maruyama
    :y1 = y1 - dt * y1 / tau_y1 + noise_y1
    :y2 = y2 - dt * y2 / tau_y2 + noise_y2
    :y3 = y3 - dt * y3 / tau_y3 + noise_y3
    :y4 = y4 - dt * y4 / tau_y4 + noise_y4
    :z1 = z1 - dt * z1 / tau_z1 + noise_z1
    :z2 = z2 - dt * z2 / tau_z2 + noise_z2
    :z3 = z3 - dt * z3 / tau_z3 + noise_z3
    :z4 = z4 - dt * z4 / tau_z4 + noise_z4
    :z5 = z5 - dt * z5 / tau_z5 + noise_z5
    :z6 = z6 - dt * z6 / tau_z6 + noise_z6
    :z7 = z7 - dt * z7 / tau_z7 + noise_z7
    
    VERBATIM
    return 0;
    ENDVERBATIM
} : end PROCEDURE states()


PROCEDURE rates(Vm (mV)) { 
    LOCAL a,b,m3_inf,n4_inf,q10,sum,one_minus_m,one_minus_h,one_minus_n
    
    UNITSOFF
    
    q10 = 3.^((celsius-6.3)/10.)
    
    ::: SODIUM :::
    
    : alpha_m and beta_m
    a = alpham(Vm)
    b = betam(Vm)
    sum = a+b
    tau_m = 1. / (q10*sum)
    m_inf = a / sum
    one_minus_m = 1. - m_inf
    m3_inf = m_inf*m_inf*m_inf

    : alpha_h and beta_h
    a = alphah(Vm)
    b = betah(Vm)
    sum = a+b
    tau_h = 1. / (q10*sum)
    h_inf = a / sum
    one_minus_h = 1. - h_inf
    
    tau_z1 = tau_h
    tau_z2 = tau_m
    tau_z3 = tau_m/2
    tau_z4 = tau_m/3
    tau_z5 = tau_m*tau_h/(tau_m+tau_h)
    tau_z6 = tau_m*tau_h/(tau_m+2*tau_h)
    tau_z7 = tau_m*tau_h/(tau_m+3*tau_h)
    var_z1 = 1.0 / Nna * m3_inf*m3_inf*h_inf * one_minus_h
    var_z2 = 3.0 / Nna * m3_inf*m_inf*m_inf*h_inf*h_inf * one_minus_m
    var_z3 = 3.0 / Nna * m3_inf*m_inf*h_inf*h_inf * one_minus_m*one_minus_m
    var_z4 = 1.0 / Nna * m3_inf*h_inf*h_inf * one_minus_m*one_minus_m*one_minus_m
    var_z5 = 3.0 / Nna * m3_inf*m_inf*m_inf*h_inf * one_minus_m*one_minus_h
    var_z6 = 3.0 / Nna * m3_inf*m_inf*h_inf * one_minus_m*one_minus_m*one_minus_h
    var_z7 = 1.0 / Nna * m3_inf*h_inf * one_minus_m*one_minus_m*one_minus_m*one_minus_h
    
    : Exact
    mu_z1 = exp(-dt/tau_z1)
    mu_z2 = exp(-dt/tau_z2)
    mu_z3 = exp(-dt/tau_z3)
    mu_z4 = exp(-dt/tau_z4)
    mu_z5 = exp(-dt/tau_z5)
    mu_z6 = exp(-dt/tau_z6)
    mu_z7 = exp(-dt/tau_z7)
    noise_z1 = sqrt(var_z1 * (1-mu_z1*mu_z1)) * normrand(0,1)
    noise_z2 = sqrt(var_z2 * (1-mu_z2*mu_z2)) * normrand(0,1)
    noise_z3 = sqrt(var_z3 * (1-mu_z3*mu_z3)) * normrand(0,1)
    noise_z4 = sqrt(var_z4 * (1-mu_z4*mu_z4)) * normrand(0,1)
    noise_z5 = sqrt(var_z5 * (1-mu_z5*mu_z5)) * normrand(0,1)
    noise_z6 = sqrt(var_z6 * (1-mu_z6*mu_z6)) * normrand(0,1)
    noise_z7 = sqrt(var_z7 * (1-mu_z7*mu_z7)) * normrand(0,1)
    : Euler-Maruyama
    :noise_z1 = sqrt(2 * dt * var_z1 / tau_z1) * normrand(0,1)
    :noise_z2 = sqrt(2 * dt * var_z2 / tau_z2) * normrand(0,1)
    :noise_z3 = sqrt(2 * dt * var_z3 / tau_z3) * normrand(0,1)
    :noise_z4 = sqrt(2 * dt * var_z4 / tau_z4) * normrand(0,1)
    :noise_z5 = sqrt(2 * dt * var_z5 / tau_z5) * normrand(0,1)
    :noise_z6 = sqrt(2 * dt * var_z6 / tau_z6) * normrand(0,1)
    :noise_z7 = sqrt(2 * dt * var_z7 / tau_z7) * normrand(0,1)
    
    ::: POTASSIUM :::
    
    : alpha_n and beta_n
    a = alphan(Vm)
    b = betan(Vm)
    sum = a+b
    tau_n = 1. / (q10*sum)
    n_inf = a / sum
    one_minus_n = 1. - n_inf
    n4_inf = n_inf * n_inf * n_inf * n_inf
    tau_y1 = tau_n
    tau_y2 = tau_n/2
    tau_y3 = tau_n/3
    tau_y4 = tau_n/4
    var_y1 = 4.0/Nk * n4_inf*n_inf*n_inf*n_inf * one_minus_n
    var_y2 = 6.0/Nk * n4_inf*n_inf*n_inf * one_minus_n*one_minus_n
    var_y3 = 4.0/Nk * n4_inf*n_inf * one_minus_n*one_minus_n*one_minus_n
    var_y4 = 1.0/Nk * n4_inf * one_minus_n*one_minus_n*one_minus_n*one_minus_n
    : Exact
    mu_y1 = exp(-dt/tau_y1)
    mu_y2 = exp(-dt/tau_y2)
    mu_y3 = exp(-dt/tau_y3)
    mu_y4 = exp(-dt/tau_y4)
    noise_y1 = sqrt(var_y1 * (1-mu_y1*mu_y1)) * normrand(0,1)
    noise_y2 = sqrt(var_y2 * (1-mu_y2*mu_y2)) * normrand(0,1)
    noise_y3 = sqrt(var_y3 * (1-mu_y3*mu_y3)) * normrand(0,1)
    noise_y3 = sqrt(var_y3 * (1-mu_y3*mu_y3)) * normrand(0,1)
    : Euler-Maruyama
    :noise_y1 = sqrt(2 * dt * var_y1 / tau_y1) * normrand(0,1)
    :noise_y2 = sqrt(2 * dt * var_y2 / tau_y2) * normrand(0,1)
    :noise_y3 = sqrt(2 * dt * var_y3 / tau_y3) * normrand(0,1)
    :noise_y4 = sqrt(2 * dt * var_y4 / tau_y4) * normrand(0,1)
    
    UNITSON
}

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

