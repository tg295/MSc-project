: kdr.mod is the delayed-rectifier K+ current from
: Values taken fro Tigerholm et al (2015)

NEURON {
	SUFFIX kdr
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
	gkbar = 18.0017e-3 (S/cm2)
	eK = -85 (mV)
	q10k = 3.3
}

ASSIGNED {
	v	(mV) : NEURON provides this
	ik	(mA/cm2)
	gk	(S/cm2)
	ek (mV)
}

STATE { n }

BREAKPOINT {
	SOLVE states METHOD cnexp
	gk = gkbar*n*n*n*n
	ik = gk * (v-eK)
}

INITIAL {
	n = ninf(v)
}

DERIVATIVE states {
  n' = (ninf(v) - n)/tau_n(v)
}

FUNCTION ninf (Vm (mV)) () {
UNITSOFF
	ninf = 1/(1+exp(-(Vm+45)/15.4))
UNITSON
}

FUNCTION tau_n (Vm (mV)) (/ms) {

UNITSOFF
  if (v <= 31){
		tau_n = 1000*(0.000688 + 1/(exp((v +75.2)/6.5)+exp((v - 131.5)/-34.8)))*q10k
	}
	else{
    tau_n = 0.16 + 0.8*exp(-0.0267*(v + 11))*q10k
	}
UNITSON
}
