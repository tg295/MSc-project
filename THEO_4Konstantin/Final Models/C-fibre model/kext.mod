: Extracellular potassium ion accumulation

NEURON {
	SUFFIX kacc
	USEION k READ ik WRITE ko
  GLOBAL kbath0
	RANGE fhspace, txfer
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	FARADAY = (faraday) (coulombs)
	(molar) = (1/liter)
	(mM) = (millimolar)
  (S) = (siemens)
}

PARAMETER {
  kbath0 = 5.6 (mM)
  fhspace = 290 (angstrom)
	ko0	= 5.6  (mM)		:	Initial K conc in Extracellular space
  txfer = 50 (ms)
}

ASSIGNED {
	ik (mA/cm2)
}

STATE {
  ko (mM)
}


BREAKPOINT {
	SOLVE state METHOD derivimplicit
}

DERIVATIVE state {
  ko' = (1e8)*ik/(fhspace*FARADAY) + (kbath0 - ko)/txfer
}
