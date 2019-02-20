: ka.mod is the a type K+ current from
: Values taken fro Tigerholm et al (2015)

NEURON {
	SUFFIX ka
	RANGE ek
	USEION k READ ek WRITE ik VALENCE 1
	RANGE gk, gkbar, ik, tau_n, ninf, tau_h, hinf
}

UNITS {
	(S) = (siemens)
	(mV) = (millivolts)
	(mA) = (milliamp)
}

PARAMETER {
	gkbar = 12.7555e-3 (S/cm2)
	ek = -85 (mV)
	q10k = 3.3
	}

ASSIGNED {
	v	(mV) : NEURON provides this
	ik	(mA/cm2)
	gk	(S/cm2)
	tau_n	(ms)
  tau_h (ms)
  ninf
  hinf
}

STATE { n h }

BREAKPOINT {
	SOLVE states METHOD cnexp
	gk = gkbar*n*h
	ik = gk * (v-ek)
}

INITIAL {
	rates(v)
	n = ninf
  h = hinf
}

DERIVATIVE states {
	rates(v)
  n' = (ninf - n)/tau_n
  h' = (hinf - h)/tau_h
}

PROCEDURE rates(v) {
  ninf = (1/(1+exp(-(v+5.4+15)/16.4)))^4

  hinf = 1/(1+exp(-(v+49.9+15)/4.6))

  tau_n = 0.25 + 10.04*exp(-((v+24.67)^2)/(2*34.8^2))*q10k

  tau_h = 20 + 50*exp(-((v+40)^2)/(2*40^2))*q10k

  if (tau_h<5) {
  tau_h = 5
  }
}
