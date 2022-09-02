# ChIPanalyser

## ChIPanalyser
ChIPanalyser: Predicting Transcription Factor Binding Sites

## Authors
Patrick C.N. Martin <pm16057@essex.ac.uk>

and

Dr. Nicolae Radu Zabet <nzabet@essex.ac.uk>

## Description
Based on a statistical thermodynamic framework, ChIPanalyser
tries to produce ChIP-seq like profile.
The model relies on four consideration:
TF binding sites can be scored using a Position weight Matrix,
DNA accessibility plays a role in Transcription Factor binding,
binding profiles are dependant on	the number of transcription
factors bound to DNA and finally binding energy
(another way of describing PWM's) or binding specificity should be
modulated (hence the introduction of a binding specificity modulator).
The end result of ChIPanalyser is to produce profiles simulating
real ChIP-seq profile and provide accuracy measurements
of these predicted profiles after being compared to real ChIP-seq data.
The ultimate goal is to produce ChIP-seq like profiles predicting ChIP-seq
like profile to circumvent the need to produce costly ChIP-seq experiments.



## References

Zabet NR, Adryan B (2015) Estimating binding properties of transcription
factors from genome-wide binding profiles. Nucleic Acids Res., 43, 84–94.


Patrick C.N. Martin and Nicolae Radu Zabe (2020) Dissecting the binding 
mechanisms of transcription factors to DNA using a statistical 
thermodynamics framework. CSBJ, 18, 3590-3605.
