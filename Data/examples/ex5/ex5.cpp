/*
 * Create a very simple conic problem.
 * min x_2
 * s.t.
 *    3x_1+x_2=3
 *    (x_1,x_2) \in L^2
 *
 * Problem is unbounded without conic constraints
 *
 */
#include <iostream>
#include <mosek.h>

static void MSKAPI printstr(void *handle,
                            MSKCONST char str[]) {
  printf("%s",str);
} /* printstr */

struct Problem {
  MSKenv_t env_;
  MSKtask_t task_;
};

Problem * create_problem() {
  MSKrescodee  r;
  const MSKint32t numvar = 2;
  const MSKint32t numcon = 1;
  MSKboundkeye bkc[] = {MSK_BK_FX};
  double       blc[] = {3.0};
  double       buc[] = {3.0};
  MSKboundkeye bkx[] = {MSK_BK_FR,
                        MSK_BK_FR};
  double       blx[] = {-MSK_INFINITY,
                        -MSK_INFINITY};
  double       bux[] = {+MSK_INFINITY,
                        +MSK_INFINITY};
  double       c[]   = {0.0,
                        1.0};
  MSKint32t     aptrb[] = {0, 1};
  MSKint32t     aptre[] = {1, 2};
  MSKint32t     asub[] = {0,
                          0};
  double      aval[] = {3.0,
                        1.0};
  /* will be used for cones  */
  MSKint32t i, j;
  MSKint32t csub[2];
  MSKenv_t env  = NULL;
  MSKtask_t task = NULL;
  r = MSK_makeenv(&env,NULL);
  r = MSK_maketask(env,numcon,numvar,&task);
  r = MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,printstr);
  r = MSK_appendcons(task,numcon);
  r = MSK_appendvars(task,numvar);
  for(j=0; j<numvar && r == MSK_RES_OK; ++j) {
    r = MSK_putcj(task,j,c[j]);
    r = MSK_putvarbound(task,
			j,           /* Index of variable.*/
			bkx[j],      /* Bound key.*/
			blx[j],      /* Numerical value of lower bound.*/
			bux[j]);     /* Numerical value of upper bound.*/
    r = MSK_putacol(task,
		    j,                 /* Variable (column) index.*/
		    aptre[j]-aptrb[j], /* Number of non-zeros in column j.*/
		    asub+aptrb[j],     /* Pointer to row indexes of column j.*/
		    aval+aptrb[j]);    /* Pointer to Values of column j.*/
  }
  for(i=0; i<numcon && r==MSK_RES_OK; ++i)
    r = MSK_putconbound(task,
			i,           /* Index of constraint.*/
			bkc[i],      /* Bound key.*/
			blc[i],      /* Numerical value of lower bound.*/
			buc[i]);     /* Numerical value of upper bound.*/
  csub[0] = 0;
  csub[1] = 1;
  r = MSK_appendcone(task,
		     MSK_CT_QUAD,
		     0.0, /* For future use only, can be set to 0.0 */
		     2,
		     csub);
  r = MSK_putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE);
  Problem * p = new Problem;
  p->env_ = env;
  p->task_= task;
  return p;
}

int main() {
  Problem * p = create_problem();
  MSK_writedata(p->task_, "ex5.mps");
  MSK_deletetask(&(p->task_));
  MSK_deleteenv(&(p->env_));
  delete p;
  return 0;
}
