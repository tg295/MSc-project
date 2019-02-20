: il.mod is the leak current from
: Values taken fro Tigerholm et al (2015)

NEURON {
	SUFFIX leak
	USEION na READ ena WRITE ina
  USEION k READ ek  WRITE ik
	RANGE gk, gna, ina, ik
}

UNITS {
	(S) = (siemens)
	(mV) = (millivolts)
	(mA) = (milliamp)
}

PARAMETER {
	ek = -85 (mV)
  ena = 69 (mV)
	gnav1p8bar = 242.7124e-3 (S/cm2)
  gnav1p9bar = 0.0948e-3 (S/cm2)
	gnav1p7bar = 106.6e-3 (S/cm2)
  ghbar = 2.5377e-3 (S/cm2)
  gp =  0.0048e-3 (S/cm2)
	gkmbar = 6.9733e-3 (S/cm2)
	gkabar = 12.7555e-3 (S/cm2)
	gkdrbar = 18.0017e-3 (S/cm2)
	gknabar = 0.0012e-3 (S/cm2)
	v_init = -80 (mV)
}

ASSIGNED {
	v	(mV) : NEURON provides this
	ina_nav1p7 (mA/cm2)
	ina_nav1p8 (mA/cm2)
	ina_nav1p9 (mA/cm2)
	ina_NaKpump (mA/cm2)
	ina_h (mA/cm2)
	ik_km (mA/cm2)
	ik_ka (mA/cm2)
	ik_h (mA/cm2)
	ik_kdr (mA/cm2)
	ik_NaKpump (mA/cm2)
	ik_kna (mA/cm2)
	ik	(mA/cm2)
  ina (mA/cm2)
	gk	(S/cm2)
  gna (S/cm2)
}

BREAKPOINT {
	ina_nav1p7 = 	gnav1p7bar*(v-ena)
	ina_nav1p8 = gnav1p8bar*(v-ena)
	ina_nav1p9 = gnav1p9bar*(v-ena)
	ina_h = ghbar*(v-ena)
	ina_NaKpump = gp*(v-ena)
	ik_km = gkmbar*(v-ek)
	ik_ka = gkabar*(v-ek)
	ik_h = ghbar*(v-ek)
	ik_kdr = gkdrbar*(v-ek)
	ik_NaKpump = gp*(v-ek)
	ik_kna = gknabar*(v-ek)

  gna = -(ina_nav1p7 + ina_nav1p9 + ina_nav1p8 + ina_h + ina_NaKpump)/(v_init-ena)
  gk = -(ik_km + ik_ka + ik_h + ik_kdr + ik_NaKpump + ik_kna)/(v_init-ek)

	ik = gk*(v-ek)
	ina = gna*(v-ena)
}
