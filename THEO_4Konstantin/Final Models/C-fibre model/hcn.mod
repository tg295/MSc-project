: hcn.mod is the h current from
: Values taken from Tigerholm et al (2015)

NEURON {
	SUFFIX hcn
	USEION k READ ek WRITE ik VALENCE 1
  USEION na READ ena WRITE ina VALENCE 1
	RANGE gh, ghbar, ik, ina, tau_ns, tau_nf, ninf
}

UNITS {
	(S) = (siemens)
	(mV) = (millivolts)
	(mA) = (milliamp)
}

PARAMETER {
	ghbar = 2.5377e-3 (S/cm2)
	ek = -85 (mV)
  ena = 69 (mV)
	q10h = 3
}

ASSIGNED {
	v	(mV)
	ik	(mA/cm2)
  ina (mA/cm2)
	gh	(S/cm2)
	tau_ns (ms)
	tau_nf (ms)
	ninf
}

STATE { ns nf }

BREAKPOINT {
	SOLVE states METHOD cnexp
	gh = ghbar*(0.5*ns + 0.5*nf)
	ina = 0.5 * gh * (v+ena)
  ik = 0.5 * gh * (v+ek)
}

INITIAL {
	rates(v)
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

	ninf = 1/(1+exp((v+87.2)/9.7))

  if (v > -70){
		tau_ns = 300+542*exp((v + 25)/20)*q10h
  }
	if (v < -70){
    tau_ns = 2500+100*exp((v+240)/50)*q10h
	}
  if (v > -70){
    tau_nf = 140 + 50*exp((v+25)/(-20))*q10h
  }
  if (v < -70){
    tau_nf = 250 + 12*exp((v + 240)/50)*q10h
  }
}

UNITSON
