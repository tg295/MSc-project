TITLE sodium potassium pump
:  from Lindblad et al Am J Physiol 1996 275:H1666

: Original model has been modified to assume constant nai

NEURON {
	SUFFIX NaKpump
	USEION k READ ko WRITE ik
	USEION na READ nai WRITE ina
	RANGE  ik, ina, ko, nai
}

UNITS {
	(molar)  =  (1/liter)
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mM) = (millimolar)
	(S) = (siemens)


}

PARAMETER {
	gp = 0.0048e-3 (S/cm2)
}

ASSIGNED {
	v (mV)
	ko (mM)
	nai (mM)
	ik (mA/cm2)
	ina (mA/cm2)
}

UNITSOFF

BREAKPOINT {
	ik = gp/((1+1/ko)^2)*(1.62/(1+(6.7/(nai+8))^3)+1.0/(1+(67.6/(nai+8))^3))
	ina = -3/2*ik
}

UNITSON
