#include <Rcpp.h>

// solver
extern "C" {
  void solver_();
}


// ode system
extern "C" {
  void Rcppfct_(int* inp) {
    Rcpp::Rcout << *inp << std::endl;
  }
}


//' @export
// [[Rcpp::export]]
void solver_wrapper() {
  solver_();
}


typedef void (*rhs)(int* n, double* t, double* y, double* f);

// FATODE explicit Runge Kutta
extern "C" {
  void rk_(double* tin, double* tout, int* nvar, double* var,
                  double* rtol, double* atol, rhs fun);
}


void lv(int* n, double* t, double* y, double* f) {
  double a = 0.9;
  double b = 1.0;
  double c = 0.1;
  double d = 1.0;

  f[0] = y[0]*y[1]*c - y[0]*d;
  f[1] = y[1]*a - y[0]*y[1]*b;


  if(*t > 3 && *t < 5) {
    Rcpp::Rcout << f[0] << std::endl;
  }

}

//' @export
// [[Rcpp::export]]
void explicit_RK() {
  double tin;
  double tout;
  int nvar;
  double* var;
  double rtol;
  double* atol;
  rhs fun;
  fun = lv;

  tin = 0.0;
  tout = 10.0;
  nvar = 2;

  var = new double[2];
  var[0] = 10.0;
  var[1] = 10.0;
  rtol = 1e-2;

  atol = new double[2];
  atol[0] = 1e-2;
  atol[1] = 1e-2;

  rk_(&tin, &tout, &nvar, var,
             &rtol, atol, fun);

  delete[] var;
  delete[] atol;

}
