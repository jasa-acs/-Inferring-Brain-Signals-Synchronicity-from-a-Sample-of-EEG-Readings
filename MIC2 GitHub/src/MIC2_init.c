#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _MIC2_clustalign(SEXP, SEXP);
extern SEXP _MIC2_EigLap(SEXP, SEXP, SEXP);
extern SEXP _MIC2_MIC_mcmc(SEXP, SEXP, SEXP, SEXP);
extern SEXP _MIC2_SpecOnly(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _MIC2_specParzen(SEXP, SEXP, SEXP, SEXP);
extern SEXP _MIC2_SpecSim(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_MIC2_clustalign", (DL_FUNC) &_MIC2_clustalign, 2},
  {"_MIC2_EigLap",     (DL_FUNC) &_MIC2_EigLap,     3},
  {"_MIC2_MIC_mcmc",   (DL_FUNC) &_MIC2_MIC_mcmc,   4},
  {"_MIC2_SpecOnly",   (DL_FUNC) &_MIC2_SpecOnly,   6},
  {"_MIC2_specParzen", (DL_FUNC) &_MIC2_specParzen, 4},
  {"_MIC2_SpecSim",    (DL_FUNC) &_MIC2_SpecSim,    6},
  {NULL, NULL, 0}
};

void R_init_MIC2(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
