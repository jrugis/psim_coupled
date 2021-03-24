/* -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2019, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------
 * Based on CVODE example problem.
 * -----------------------------------------------------------------*/

#include <fstream>
#include <string>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */

#include "cCVode.hpp"
#include "cCell_flow.hpp"
#include "global_defs.hpp"
#include "utils.hpp"

/*
 * f routine. Compute function f(t,y).
 */

static int f(realtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  // pointer to cell_flow object
  cCell_flow* pt_flow = static_cast<cCell_flow*>(user_data);

  // create input and output arrays for calling cell flow function
  Array1IC ymat, ydotmat;
  // copy input from sundials data structure to array for calling secretion function
  for (int i = 0; i < IONCOUNT; i++) { ymat(i) = NV_Ith_S(y, i); }

  // call secretion function
  pt_flow->secretion(t, ymat, ydotmat);

  // copy result back into sundials data structure
  for (int i = 0; i < IONCOUNT; i++) { NV_Ith_S(ydot, i) = ydotmat(i); }

  return (0);
}

/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

cCVode::cCVode(cCell_flow* flow_, std::ofstream& out_, realtype abstol_, realtype reltol_)
  : flow(flow_), out(out_), initialised(false), nvars(0), abstol(abstol_), reltol(reltol_)
{
  out << std::scientific;
  out << "<CVode> creating CVode solver" << std::endl;
  out << "<CVode> tolerances are " << abstol << " (absolute) and ";
  out << reltol << " (relative)" << std::endl;
  out << std::fixed;
  y = NULL;
  A = NULL;
  LS = NULL;
  cvode_mem = NULL;
}

cCVode::~cCVode()
{
  if (initialised) {
    /* Free y vector */
    N_VDestroy(y);

    /* Free integrator memory */
    CVodeFree(&cvode_mem);

    /* Free the linear solver memory */
    SUNLinSolFree(LS);

    /* Free the matrix memory */
    SUNMatDestroy(A);
  }
}

void cCVode::init(Array1IC& yini)
{
  int retval;
  nvars = IONCOUNT;
  realtype t0 = 0.0;

  /* Create serial vector of length NEQ for I.C. */
  y = N_VNew_Serial(nvars);
  check_retval(static_cast<void*>(y), "N_VNew_Serial", 0);

  /* Initialize y */
  for (sunindextype i = 0; i < nvars; i++) { NV_Ith_S(y, i) = yini(i); }

  /* Call CVodeCreate to create the solver memory and specify the
   * Backward Differentiation Formula */
  cvode_mem = CVodeCreate(CV_BDF);
  check_retval(static_cast<void*>(cvode_mem), "CVodeCreate", 0);

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the inital time T0, and
   * the initial dependent variable vector y. */
  retval = CVodeInit(cvode_mem, f, t0, y);
  check_retval(&retval, "CVodeInit", 1);

  /* Call CVodeSStolerances to specify the scalar relative tolerance
   * and scalar absolute tolerances */
  retval = CVodeSStolerances(cvode_mem, reltol, abstol);
  check_retval(&retval, "CVodeSStolerances", 1);

  /* Call CVodeSetUserData */
  retval = CVodeSetUserData(cvode_mem, static_cast<void*>(flow));
  check_retval(&retval, "CVodeSetUserData", 1);

  /* Create dense SUNMatrix for use in linear solves */
  A = SUNDenseMatrix(nvars, nvars);
  check_retval(static_cast<void*>(A), "SUNDenseMatrix", 0);

  /* Create dense SUNLinearSolver object for use by CVode */
  LS = SUNLinSol_Dense(y, A);
  check_retval(static_cast<void*>(LS), "SUNLinSol_Dense", 0);

  /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
  retval = CVodeSetLinearSolver(cvode_mem, LS, A);
  check_retval(&retval, "CVodeSetLinearSolver", 1);

  initialised = true;
}

void cCVode::run(realtype t, realtype tout, Array1IC& yout)
{
  int retval;

  // reinit CVode
  for (sunindextype i = 0; i < nvars; i++) { NV_Ith_S(y, i) = yout(i); }
  retval = CVodeReInit(cvode_mem, t, y);
  check_retval(&retval, "CVodeReinit", 1);

  // call CVode
  retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);

  // check for errors
  check_retval(&retval, "CVode", 1);
  if (retval != CV_SUCCESS) { utils::fatal_error("CVode did not succeed", out); }

  /* Print some final statistics */
  //  PrintFinalStatsDetailed();
  //  PrintFinalStatsBrief();

  /* store result for passing back */
  for (sunindextype i = 0; i < nvars; i++) { yout(i) = NV_Ith_S(y, i); }
}

/*
 * Get and print some final statistics
 */

void cCVode::PrintFinalStatsBrief()
{
  long int nst, nfe;
  int retval;

  retval = CVodeGetNumSteps(cvode_mem, &nst);
  check_retval(&retval, "CVodeGetNumSteps", 1);
  retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_retval(&retval, "CVodeGetNumRhsEvals", 1);

  out << "<cCVode> solver statistics: ";
  out << "num_steps = " << nfe << " num_func_calls = " << nfe << std::endl;
}

void cCVode::PrintFinalStatsDetailed()
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  int retval;

  retval = CVodeGetNumSteps(cvode_mem, &nst);
  check_retval(&retval, "CVodeGetNumSteps", 1);
  retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_retval(&retval, "CVodeGetNumRhsEvals", 1);
  retval = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_retval(&retval, "CVodeGetNumLinSolvSetups", 1);
  retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_retval(&retval, "CVodeGetNumErrTestFails", 1);
  retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_retval(&retval, "CVodeGetNumNonlinSolvIters", 1);
  retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_retval(&retval, "CVodeGetNumNonlinSolvConvFails", 1);

  retval = CVodeGetNumJacEvals(cvode_mem, &nje);
  check_retval(&retval, "CVodeGetNumJacEvals", 1);
  retval = CVodeGetNumLinRhsEvals(cvode_mem, &nfeLS);
  check_retval(&retval, "CVodeGetNumLinRhsEvals", 1);

  retval = CVodeGetNumGEvals(cvode_mem, &nge);
  check_retval(&retval, "CVodeGetNumGEvals", 1);

  out << "Final Statistics:" << std::endl;
  out << "nst = " << nfe << " nfe = " << nfe << " nsetups = " << nsetups << std::endl;
  out << "nfeLS = " << nfeLS << " nje = " << nje << " nni = " << nni << std::endl;
  out << "ncfn = " << ncfn << " netf = " << netf << " nge = " << nge << std::endl;
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns an integer value so check if
 *            retval < 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer
 */

void cCVode::check_retval(void* returnvalue, std::string funcname, int opt)
{
  int* retval;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && returnvalue == NULL) {
    utils::fatal_error("SUNDIALS_ERROR: " + funcname + "() failed - returned NULL pointer", out);
  }

  /* Check if retval < 0 */
  else if (opt == 1) {
    retval = (int*)returnvalue;
    if (*retval < 0) {
      utils::fatal_error("SUNDIALS_ERROR: " + funcname + "() failed with retval = " + std::to_string(*retval), out);
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && returnvalue == NULL) {
    utils::fatal_error("MEMORY_ERROR: " + funcname + "() failed - returned NULL pointer", out);
  }
}
