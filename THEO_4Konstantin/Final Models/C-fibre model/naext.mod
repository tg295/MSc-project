: Extracellular sodium ion accumulation

NEURON {
	SUFFIX naext
	USEION na READ ina WRITE nao
  GLOBAL nabath0
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
  nabath0 = 154 (mM)
  fhspace = 290 (angstrom)
  nao0	= 154  (mM)		:	Initial na conc in Extracellular space
  txfer = 50 (ms)
}

ASSIGNED {
	ina (mA/cm2)
}

STATE { nao (mM) }


BREAKPOINT {
	SOLVE state METHOD derivimplicit
}

DERIVATIVE state {
  nao' = (1e8)*ina/(fhspace*FARADAY) + (nabath0 - nao)/txfer
}
