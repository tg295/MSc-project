: nav1p8.mod is the NaV1.8 Na+ current from
: Values taken from Tigerholm et al (2015)

NEURON {
	SUFFIX nav1p8
	RANGE ena
	USEION na READ ena WRITE ina VALENCE 1
	RANGE gnabar, gna, ina, tau_h, tau_m, tau_s, tau_u, minf, hinf, sinf, uinf
}

UNITS {
	(S) = (siemens)
	(mV) = (millivolts)
	(mA) = (milliamp)
}

PARAMETER {
	gnabar = 242.7124e-3 (S/cm2)
	ena = 69 (mV)
	q10na = 2.5
}

ASSIGNED {
	v	(mV) : NEURON provides this
	ina	(mA/cm2)
	gna	(S/cm2)
	tau_h	(ms)
	tau_m	(ms)
	tau_s (ms)
	tau_u (ms)
	minf
	hinf
	sinf
	uinf
}

STATE { m h s u }

BREAKPOINT {
	SOLVE states METHOD cnexp
	gna = gnabar*m*m*m*h*s*u
	ina = gna * (v-ena)
}

UNITSOFF

INITIAL {
	rates(v)
	m = alpham(v)/(alpham(v)+betam(v))
	h = 1/(1+exp((v + 32.2)/4))
	s = 1/(1+exp((v + 45)/8))
	u = 1/(1+exp((v + 51)/8))
}

DERIVATIVE states {
	rates(v)
  m' = (minf - m)/tau_m
  h' = (hinf - h)/tau_h
	s' = (sinf - s)/tau_s
	u' = (uinf - u)/tau_u
}


FUNCTION alpham(v (mV)) (/ms) {
	alpham = 2.85 - 2.839/(1 + exp((v - 1.159)/13.95))
}

FUNCTION alphas(v (mV)) (/ms) {
	alphas = 0.001*5.4203/(1 + exp((v + 79.816)/16.269))
}

FUNCTION alphau(v (mV)) (/ms) {
	alphau = 0.0002*2.0434/(1 + exp((v + 67.499)/19.51))
}

FUNCTION betam(v (mV)) (/ms) {
	betam= 7.6205/(1 + exp((v + 46.463)/8.8289))
}

FUNCTION betas(v (mV)) (/ms) {
	betas = 0.001*5.0757/(1+exp(-(v + 15.968)/11.542))
}

FUNCTION betau(v (mV)) (/ms) {
	betau = 0.0002*1.9952/(1+exp(-(v + 30.963)/14.792))
}

FUNCTION rates(v (mV)) (/ms) {
	tau_m = 1.0 / (alpham(v) + betam(v)) * q10na
	minf = alpham(v)/(alpham(v)+betam(v))

	tau_h = 1.218 + 42.043*exp(-((v + 38.1)^2)/(2*15.19^2))*q10na
	hinf = 1/(1+exp((v + 32.2)/4))

  tau_s = 1/(alphas(v) + betas(v)) * q10na
	sinf = 1/(1+exp((v + 45)/8))

	tau_u = 1/(alphau(v) + betau(v)) * q10na
	uinf = 1/(1+exp((v + 51)/8))
}

UNITSON 
