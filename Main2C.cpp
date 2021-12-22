#include <iostream>
#include <fstream>
#include <vector>
#include <exception>
#include <cstdlib>
#include <cmath>
#include <string>


//*** parameters of the integration algorithm *********

    const int nvar = 6; // number of variables
    const double dt0 = 0.0005; // initial time step size
    const double dtsav = 0.05; // save data after time steps
    const double tEnd = 10000.0; // end time
    const double tolerance = 1.0e-6; // acceptable local error during numerical integration
//*** model parameters *********

    // general
    const double K0 = 4.0; // half saturation constant of plasmid free
    const double K1 = 4.0; // half saturation constant of plasmid bearing
    const double e1 = 6.25 * 1e-7; // resource needed to divide once for plasmid bearing
    const double e0 = 6.25 * 1e-7; // resource needed to divide once for plasmid free
    const double l = 1e-3; // loss of plasmid
    const double r0 = 0.738; // max growth rate of plasmid free
    const double r1 = 0.6642; // max growth rate of plasmid bearing

    // in lumen
    const double VL = 10.0; // Volume of lumen chemostat
    const double WL = 4.5; // Rate at which the nutrient solution enters (and leaves) the chemostat (lumen)
    const double cL = pow(10, -9); // conjugation factor in lumen
    const double DL = WL/VL; // flow rate in lumen
    const double SLin = 33; // resource coming in to the lumen
    double KLW0 = 1e-9; // attaching of plasmid free cells to wall
    const double KLW1 = 1e-9; // attaching of plasmid bearing cells to wall
    const double d0L = DL; // death rate of plasmid free cells in  lumen
    const double d1L = DL; // death rate of plasmid bearing cells in lumen

    // at wall
    const double VW = 2.0; // Volume of wall chemostat 
    const double WW = 0.70; // Rate at which the nutrient solution enters (and leaves) the chemostat (wall)
    const double cW = pow(10, -9); // conjugation factor at wall
    const double DW = WW/VW; // flow rate at wall
    const double d0W = DW; // death rate of plasmide free at wall
    const double d1W = DW; // deat rate of plasmid bearing at wall
    const double SWin = 9; // resource coming to the wall
    const double KWL0 = 1e-9; // dettaching of plasmid free cells of wall
    const double KWL1 = 1e-9; // dettaching of plasmid bearing cells of wall

    const std::vector< double > D = {DW, DL};
    std::vector< double > Sin = {SWin, SLin};
    const std::vector< double > c = {cW, cL};

//*** ODE description *********

enum location {Wall, Lumen};

void update_state(const std::vector<double> &x, std::vector<double> &dxdt, location loc) {
  size_t init_index = 0;
  if (loc == Lumen) init_index = 3;

  double R =  x[init_index];
  double N0 = x[init_index + 1];
  double N1 = x[init_index + 2];

  const double psi0 = ((r0 * R) / (K0 + R));
  const double psi1 = ((r1 * R) / (K1 + R));


  dxdt[init_index] = D[loc] * (Sin[loc] - R) - e0 * psi0 * N0 - e1 * psi1 * N1;

  dxdt[init_index + 1] = psi0 * N0 - D[loc] * N0 + l * N1 - c[loc] * N0 * N1;

  dxdt[init_index + 2] = psi1 * N1 - D[loc] * N1 - l * N1 + c[loc] * N0 * N1;

  return;
}

void exchange_cells(const std::vector<double> &x, std::vector<double> &dxdt) 
{

  double N0W = x[1]; // plasmid free at wall
  double N1W = x[2]; // plasmid bearing at wall

  double N0L = x[4]; // plasmid free in lumen
  double N1L = x[5]; // plasmid bearing in lumen

  dxdt[1] += (KLW0 * N0L * VL)/VW - KWL0 * N0W; // Idk really but I thought this way I can compensate for the 
  dxdt[2] += (KLW1 * N1L * VL)/VW - KWL0 * N1W; // transferring of cells into the different volumes. 
  dxdt[4] += (KWL0 * N0W * VW)/VL - KLW0 * N0L; // and this is how i could portray the volumes 
  dxdt[5] += (KWL1 * N1W * VW)/VL - KLW1 * N1L;

  return;
}

void rhs(const double &t, const std::vector<double> &x, std::vector<double> &dxdt)
{
  update_state(x, dxdt, Wall);
  update_state(x, dxdt, Lumen);

  exchange_cells(x, dxdt);
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
        }
        return false;
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
  x[0] = SWin;
  x[1] = initialN0;
  x[2] = initialN1;
  x[3] = SLin;
  x[4] = initialN0;
  x[5] = initialN1;
  std::vector<double> dxdt(nvar);
  rhs(0.0, x, dxdt);

  // start numerical integration
  int nOK = 0, nStep = 0;
//  double dtMin = dt0, dtMax = kdMinH;
  double t;
  double tsav;
  double dt;
  for(t = 0.0, tsav = 0.0, dt = dt0; t < tEnd; ++nStep)
  {
      if(BogackiShampineStepper(t, x, dxdt, dt))
      {
          ++nOK;
      }
      if (fabs(dxdt[1]) < 1.0e-4 && fabs(dxdt[2]) < 1.0e-4 && fabs(dxdt[4]) < 1.0e-4 && fabs(dxdt[5]) < 1.0e-4)
      {
          break;
      }
  }

  if (dxdt[2] > -1.0e-4 && dxdt[2] < 0 && x[2] < initialN1)
  {
    x[2] = 0;
  }
    if (dxdt[5] > -1.0e-4 && dxdt[5] < 0 && x[5] < initialN1) 
  {
    x[5] = 0;
  }
  if (dxdt[1] > -1.0e-4 && dxdt[1] < 0 && x[1] < initialN0)
  {
    x[1] = 0;
  }
  if (dxdt[4] > -1.0e-4 && dxdt[4] < 0 && x[4] < initialN0) 
  {
    x[4] = 0;
  }

  std::ofstream ofs(output_filename.c_str(), std::ios::app);
   if(!ofs.is_open())
       throw std::runtime_error("unable to open file.\n");

  for(size_t i = 0; i < pars.size(); ++i ){
    ofs << pars[i] << ',';
  } ofs << x[1] << ',' << "N0" << ',' << "Wall" << '\n';

  for(size_t i = 0; i < pars.size(); ++i ){
    ofs << pars[i] << ',';
  } ofs << x[2] << ',' << "N1" << ',' << "Wall" << '\n';

  for(size_t i = 0; i < pars.size(); ++i ){
    ofs << pars[i] << ',';
  } ofs << x[4] << ',' << "N0" << ',' << "Lumen" << '\n';

  for(size_t i = 0; i < pars.size(); ++i ){
    ofs << pars[i] << ',';
  } ofs << x[5] << ',' << "N1" << ',' << "Lumen" << '\n';

  ofs.close();
}


//*** function main() *********

int main()
{
    try {

      std::string file_name = "resultsKLW0fy.csv";
      std::ofstream ofs(file_name.c_str());
      // give first row with variable names
      ofs  << "pars" << ',' << "popsize" << ',' << "population" << ',' << "location" << "\n";
      ofs.close();

      for (double local_KLW0 = 1e-10; local_KLW0 < pow(10, -0.5); local_KLW0 *= 2 ) 
      {
        KLW0 = local_KLW0;
        std::vector<double> pars = {local_KLW0};
        do_analysis(file_name, pars);
        std::cout << local_KLW0 << "\n";
      }
    }
    catch(std::exception &error)
    {
    std::cerr << "error: " << error.what();
    exit(EXIT_FAILURE);
    }
    return 0;
}
