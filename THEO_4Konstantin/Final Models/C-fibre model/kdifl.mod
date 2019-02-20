COMMENT
Longitudinal diffusion of potassium (no buffering)
(equivalent modified euler with standard method and
equivalent to diagonalized linear solver with CVODE )
ENDCOMMENT

NEURON {
	SUFFIX kdifl
	USEION k READ ik WRITE ki
	RANGE D
}

UNITS {
	(mM) = (milli/liter)
	(um) = (micron)
	FARADAY = (faraday) (coulomb)
	PI = (pi) (1)
}

PARAMETER {
	D = 1.85 (um2/ms)
}

ASSIGNED {
	ik (milliamp/cm2)
	diam (um)
}

STATE {
	ki (mM)
}

BREAKPOINT {
	SOLVE conc METHOD sparse
}

KINETIC conc {
	COMPARTMENT PI*diam*diam/4 {ki}
	LONGITUDINAL_DIFFUSION D*PI*diam*diam/4 {ki}
	~ ki << (-ik/(FARADAY)*PI*diam*(1e4))
}
