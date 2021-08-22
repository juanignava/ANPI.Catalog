#include <iostream>
#include <vector>
#include <ginac/ginac.h>
#include <matplotlibcpp.h>

using namespace std;
using namespace GiNaC;
namespace plt = matplotlibcpp;


int main()
  {

  symbol x("x");
  symtab table;
  table["x"] = x;
  parser reader(table);
  
  string ftext;
  cout << "Escriba la función: ";
  cin >> ftext;

  ex f=reader(ftext);
  ex fd=diff(f,x);

  string x_in;
  cout << "Escriba el valor inicial: ";
  cin >> x_in;
  ex x0=reader(x_in);

  ex xk=x0;
  ex tol=0.000001;
  ex error=tol+1;
  double k=0;
  std::vector< double > errors;
  std::vector< double > iter;

  while (error>=tol){
    k=k+1;
    ex num=evalf(subs(f,x==xk));
    ex den=evalf(subs(fd,x==xk));
    xk = xk -num/den;
    error=abs(evalf(subs(f,x==xk)));
    ex aux=evalf(error);
    double m = ex_to<numeric>(aux).to_double();
    errors.push_back(m);
    iter.push_back(k);
  }
  cout << "Iterations: " << k << endl;
  cout << "Approximation of Solution: "<< xk << endl;
  cout << "Error: " << error << endl;
  //cout << "Errors: " << errors << endl;

  plt::named_plot("Error |f(x_k)|", iter,errors);
  plt::title("Error Newton-Raphson");
  plt::legend();
  plt::show();

  return 0;
}