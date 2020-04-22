#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP _mcmcse_inseq(SEXP, SEXP);
extern SEXP _mcmcse_mbmC(SEXP, SEXP);
extern SEXP _mcmcse_mobmC(SEXP, SEXP);
extern SEXP _mcmcse_msveC(SEXP, SEXP, SEXP);
extern SEXP _mcmcse_lag(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_mcmcse_inseq", (DL_FUNC) &_mcmcse_inseq, 2},
    {"_mcmcse_mbmC",  (DL_FUNC) &_mcmcse_mbmC,  2},
    {"_mcmcse_mobmC",  (DL_FUNC) &_mcmcse_mobmC,  2},
    {"_mcmcse_msveC", (DL_FUNC) &_mcmcse_msveC, 3},
    {"_mcmcse_lag", (DL_FUNC) &_mcmcse_lag, 4},
    {NULL, NULL, 0}
};

void R_init_mcmcse(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
