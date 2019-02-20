: nav1.7mod is the NaV1.7 Na+ current from
: Values taken from Tigerholm et al (2015)

NEURON {
	SUFFIX nav1p7
	RANGE ena
	USEION na READ ena WRITE ina VALENCE 1
	RANGE gna, gnabar, ina, tau_h, tau_m, tau_s, minf, hinf, sinf
}

UNITS {
	(S) = (siemens)
	(mV) = (millivolts)
	(mA) = (milliamp)
}

PARAMETER {
	gnabar = 106.6e-3 (S/cm2)
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
	minf
	hinf
	sinf
}

STATE { m h s }

BREAKPOINT {
	SOLVE states METHOD cnexp
	gna = gnabar * m * m * m * h * s
	ina = gna * (v-ena)
}


INITIAL {
rates (v)
	m = minf
	h = hinf
	s = sinf
}

DERIVATIVE states {
	rates(v)
	m' = (minf - m)/tau_m
	h' = (hinf - h)/tau_h
	s' = (sinf - s)/tau_s
}

UNITSOFF

FUNCTION alpham(v (mV)) (/ms) {
	alpham = 15.5/(1 + exp((v - 5)/(-12.08)))
}

FUNCTION alphas(v (mV)) (/ms) {
	alphas = 0.00003+0.00092/(1+exp((v+93.9)/16.6))
}

FUNCTION alphah(v (mV)) (/ms) {
	alphah = 0.38685/(1 + exp((v + 122.35)/15.29))
}

FUNCTION betam(v (mV)) (/ms) {
	betam= 35.2/(1 + exp((v + 72.7)/16.7))
}

FUNCTION betas(v (mV)) (/ms) {
	betas = 132.05 - 132.05/(1+exp((v - 384.9)/28.5))
}

FUNCTION betah(v (mV)) (/ms) {
	betah = -0.00283+2.00283/(1+exp((v + 5.5266)/(-12.70195)))
}

UNITSON

FUNCTION rates(v (mV)) (/ms) {
	tau_m = 1.0 / (alpham(v) + betam(v)) * q10na
	minf = alpham(v)/(alpham(v)+betam(v))

	tau_h = 1.0 / (alphah(v) + betah(v)) * q10na
	hinf = alphah(v)/(alphah(v)+betah(v))

  tau_s = 1/(alphas(v) + betas(v)) * q10na
	sinf = alphas(v)/(alphas(v)+betas(v))

}
