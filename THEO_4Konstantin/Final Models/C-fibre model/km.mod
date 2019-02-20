: km.mod is the m type K+ current from
: Values taken fro Tigerholm et al (2015)

NEURON {
	SUFFIX km
	RANGE ek
	USEION k READ ek WRITE ik VALENCE 1
	RANGE gk, gkbar, ik, tau_n, ninf
}

UNITS {
	(S) = (siemens)
	(mV) = (millivolts)
	(mA) = (milliamp)
}

PARAMETER {
	gkbar = 6.9733e-3 (S/cm2)
	ek = -85 (mV)
	q10k = 3.3
}

ASSIGNED {
	v	(mV) : NEURON provides this
	ik	(mA/cm2)
	gk	(S/cm2)
	tau_ns	(ms)
	tau_nf  (ms)
	ninf
}

STATE { ns nf }

BREAKPOINT {
	SOLVE states METHOD cnexp
	gk = gkbar*(ns/4 + 3*nf/4)
	ik = gk*(v-ek)
}

INITIAL {
	ns = ninf
	nf = ninf
}

DERIVATIVE states {
	rates(v)
  ns' = (ninf - ns)/tau_ns
  nf' = (ninf - nf)/tau_nf
}

UNITSOFF

FUNCTION rates(v (mV)) (/ms) {
  ninf = 1/(1+exp(-(v+30)/6))

  tau_nf = 1/((0.00395*exp((v + 30)/40)) + (0.00395*exp(-(v + 30)/20)*q10k))

  if (v < -60){
		tau_ns = 219*q10k
	}
	else {
    tau_ns = 13*v + 1000*q10k
	}
}

UNITSON
