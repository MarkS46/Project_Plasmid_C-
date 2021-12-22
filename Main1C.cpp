#include <iostream>
#include <fstream>
#include <vector>
#include <exception>
#include <cstdlib>
#include <cmath> 

//*** parameters of the integration algorithm *********

    const int nvar = 3; // number of variables
    const double dt0 = 0.0005; // initial time step size
    const double dtsav = 0.05; // save data after time steps
    const double tEnd = 100000.0; // end time
    const double tolerance = 1.0e-6; // acceptable local error during numerical integration

//*** model parameters *********

    const double K0 = 4.0; // half saturation constant of plasmid bearing 
    const double K1 = 4.0; // half saturation constant of plasmid bearing
    const double e1 = 6.25 * 1e-7; // resource needed to divide once for plasmid bearing
    const double e0 = 6.25 * 1e-7; // resource needed to divide once for plasmid free
    const double l = 1e-3; // loss of plasmid 
    const double c = pow(10, -9); // conjugation factor
    const double D = 0.25; // flow rate
    const double r0 = 0.738; // growth rate of plasmid free
    const double r1 = 0.6642; // growth rate of plasmid bearing
    double Sin; //inflow concentration of resource   



//*** ODE description *********
void rhs(const double &t, const std::vector<double> &x, std::vector<double> &dxdt)
{

    double R = x[0]; // Resource
    double N0 = x[1]; // Plasmid free
    double N1 = x[2]; // Plasmid bearing 

    double psi0 = ((r0 * R) / (K0 + R)); // growth rate of plasmid free bacteria
    double psi1 = ((r1 * R) / (K1 + R)); // growth rate of plasmid bearing bacteria

    dxdt[0] = D * (Sin - R) - e0 * psi0 * N0 - e1 * psi1 * N1; // differential equation of the resource
    dxdt[1] = psi0 * N0  - D * N0 + l * N1 - c * N0 * N1; // differential equation of the plasmid free cells
    dxdt[2] = psi1 * N1 - D * N1 - l * N1 + c * N0 * N1; // differential equation of the plasmid bearing cells
 }
//*** ODE integration routine *********

    const double kdShrinkMax = 0.1; // decrease step size by no more than this factor
    const double kdGrowMax = 1.3; // increase step size by no more than this factor
    const double kdSafety = 0.9; // safety factor in adaptive stepsize control
    const double kdMinH = 1.0e-6; // minimum step size

//*** The Bogacki-Shampine stepper *********

bool BogackiShampineStepper(double &t, std::vector<double> &x,  std::vector<double>&dxdt, double &h)
{

    // step 2
    std::vector<double> xtmp(nvar);
    for(int i = 0; i < nvar; ++i)
    {
        xtmp[i] = x[i] + 0.5 * h * dxdt[i];
    }
    std::vector<double> dxdt2(nvar);
    rhs(t + 0.5 * h, xtmp, dxdt2);

    // step 3
    for(int i = 0; i < nvar; ++i)
    {
        xtmp[i] = x[i] + 0.75 * h * dxdt2[i];
    }
    std::vector<double> dxdt3(nvar);
    rhs(t + 0.75 * h, xtmp, dxdt3);

    // step 4
    for(int i = 0; i < nvar; ++i)
    {
        xtmp[i] = x[i] +(1.0/9.0) * h *(2 * dxdt[i] + 3 * dxdt2[i] + 4 * dxdt3[i]) ;
    }
    std::vector<double> dxdt4(nvar);
    rhs(t + h, xtmp, dxdt4);

    double errMax = 0.0;
    for(int i = 0; i < nvar; ++i) 
    {
        // compute error
        double erri = fabs(h * (5.0 * dxdt[i] / 72.0 - dxdt2[i] / 12.0 - dxdt3[i] / 9.0 + 0.125 * dxdt4[i])) / tolerance;
        if(erri > errMax)
        {
            errMax = erri;
        }
    }

    // adjust step size
    const double fct = errMax > 0.0 ? kdSafety / pow(errMax, 1.0/3.0) : kdGrowMax;
    if(errMax > 1.0) 
    {
        // reduce step size and reject step
        if(fct < kdShrinkMax) // is the factor to little?
        {
            h *= kdShrinkMax; // decrease step size by this factor
        }
        else
        {
            h *= fct; // else decrease by the fct factor
        }
        if(h < kdMinH)
        {
            throw std::runtime_error("step size underflow in eulerHeunAdaptiveStepper().");
        return false;
        }
    }
    else 
    {
        // update solution and increase step size
        x = xtmp;
        dxdt = dxdt4;
        t += h;
        if(fct > kdGrowMax) // is the factor to large? 
        {
            h *= kdGrowMax; // increase step size by this factor
        }
        else
        {
            h *= fct; // else increase by the fct factor 
        }
        return true;
    }
}

void do_analysis(std::string output_filename, const std::vector<double>& pars) 
{
  double initialN0 = 1.0;
  double initialN1 = 1.0e-5;

   // give initial values
  std::vector<double> x(nvar);
  x[0] = Sin;
  x[1] = initialN0;
  x[2] = initialN1;
  std::vector<double> dxdt(nvar);
  rhs(0.0, x, dxdt);

  // start numerical integration
  int nOK = 0, nStep = 0;
//  double dtMin = dt0, dtMax = kdMinH;
  double t;
  double dt;
  for(t = 0.0, dt = dt0; t < tEnd; ++nStep)
  {
      if(BogackiShampineStepper(t, x, dxdt, dt))
      {
          ++nOK;
      }
      if (fabs(dxdt[1]) < 1.0e-6 && fabs(dxdt[2]) < 1.0e-6)
      {
          break;
      }
  }

  if (fabs(dxdt[2]) < 1.0e-6 && x[2] < initialN1 * 1000)
  {
    x[2] = 0;
  }
  if (fabs(dxdt[1]) < 1.0e-6 && x[1] < initialN0 * 1000)
  {
    x[1] = 0;
  }

  std::ofstream ofs(output_filename.c_str(), std::ios::app);
   if(!ofs.is_open())
       throw std::runtime_error("unable to open file.\n");

  for(size_t i = 0; i < pars.size(); ++i ){
    ofs << pars[i] << ',';
  } ofs << x[1] << ',' << "N0" << ',' << dxdt[1] << '\n';

  for(size_t i = 0; i < pars.size(); ++i ){
    ofs << pars[i] << ',';
  } ofs << x[2] << ',' << "N1"  << ',' << dxdt[2] << '\n';


  ofs.close();
}


 //*** function main() *********

int main()
{
    try {
    std::string file_name = "results1comp0.25.csv";
    std::ofstream ofs(file_name.c_str());
    // give first row with variable names
    ofs  << "pars" << ',' << "popsize" << ',' << "population" << ',' << "DXDT" << "\n";
    ofs.close();

      for (double local_Sin = 5 ; local_Sin < 100; ++local_Sin ) 
      {
        Sin = local_Sin;
        std::vector<double> pars = {local_Sin};
        do_analysis(file_name, pars);
        std::cout << Sin << "\n";
      }
    }
    catch(std::exception &error)
    {
    std::cerr << "error: " << error.what();
    exit(EXIT_FAILURE);
    }
   
    return 0;
}