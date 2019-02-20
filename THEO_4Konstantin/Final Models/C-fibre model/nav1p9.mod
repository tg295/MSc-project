:nav1.7mod is the NaV1.9 Na+ current from
: Values taken from Tigerholm et al (2015)

NEURON {
 SUFFIX nav1p9
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
 gnabar = 0.0948e-3 (S/cm2)
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
 gna = gnabar * m * h * s
 ina = gna*(v-ena)
}


INITIAL {
rates (v)
 m = alpham(v)/(alpham(v)+betam(v))
 h = alphah(v)/(alphah(v)+betah(v))
 s = alphas(v)/(alphas(v)+betas(v))
}

DERIVATIVE states {
 rates(v)
 m' = (minf - m)/tau_m
 h' = (hinf - h)/tau_h
 s' = (sinf - s)/tau_s
}

UNITSOFF

FUNCTION alpham(v (mV)) (/ms) {
 alpham = 1.032/(1 + exp((v + 6.99)/(-14.87115)))
}

FUNCTION alphas(v (mV)) (/ms) {
 alphas = 0.00000016*exp(-v/12)
}

FUNCTION alphah(v (mV)) (/ms) {
 alphah = 0.06435/(1 + exp((v + 73.26415)/3.71928))
}

FUNCTION betam(v (mV)) (/ms) {
 betam= 5.79/(1 + exp((v + 130.4)/22.9))
}

FUNCTION betas(v (mV)) (/ms) {
 betas = 0.0005/(1+exp(-(v + 32)/23))
}

FUNCTION betah(v (mV)) (/ms) {
 betah = 0.13496/(1+exp((v + 10.27853)/(-9.09334)))
}


FUNCTION rates(v (mV)) (/ms) {
 tau_m = 1.0 / (alpham(v) + betam(v)) * q10na
 minf = alpham(v)/(alpham(v)+betam(v))

 tau_h = 1.0 / (alphah(v) + betah(v)) * q10na
 hinf = alphah(v)/(alphah(v)+betah(v))

 tau_s = 1/(alphas(v) + betas(v)) * q10na
 sinf = alphas(v)/(alphas(v)+betas(v))

}

FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = (y/x)/(-1-(x/2*y))
        }else{
                vtrap = 1/(1 - exp(x/y))
        }
}

UNITSON
