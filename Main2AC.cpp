#include <iostream>
#include <fstream>
#include <vector>
#include <exception>
#include <cstdlib>
#include <cmath>
#include <string>

//*** THE EFFECT OF WALL ATTACHMENT ON PLASMID PERSISTENCE SCRIPT 4 *********

// This is the final script of the report, with which figure 3 and 4 were produced
// equillibrium values are recorded for different parameter values
// currently parameters are set to prodcue the figure 4 A2 
// that is stronger attachment of the plasmid free cells to the wall with also more conjugation at the wall
// Adjustments can be made to recreate the other figures


//*** model parameters *********

  // general
  const double K0 = 4.0; // half saturation constant of plasmid free
  const double K1 = 4.0; // half saturation constant of plasmid bearing
  const double e1 = 6.25 * 1e-7; // resource needed to divide once for plasmid bearing cells
  const double e0 = 6.25 * 1e-7; // resource needed to divide once for plasmid free cells
  const double l = 1e-3; // loss of plasmid by divison
  const double r0 = 0.738; // max growth rate of plasmid free 
  const double r1 = 0.6642; // max growth rate of plasmid bearing
  const double H0 = 4.0; // half saturation constant of plasmid free for antibiotics
  const double H1 = 4.0; // half saturation constant of plasmid bearing for antibiotics
  const double Adisp0 = 1e-5; // antibiotic needed for the death of a plasmid free cell
  const double Adisp1 = 1e-4; // antibiotic needed for the death of a plasmid bearing cell 
  const double d1 = 1e-9; // death of plasmid bearing cell by antibiotic
  const double d0 = 1e-2; // death of plasmid free cell by antibiotic


  // in lumen
  const double cL = pow(10, -9); // conjugation factor in lumen
  const double DL = 0.25; // turnover rate in lumen
  const double SLin = 10; // input resource concentration of the lumen
  const double KLW0 = 1e-9; // attaching rate of plasmid free cells to wall
  double KLW1; // attaching rate of plasmid bearing cells to wall
  const double d0L = DL; // death rate of plasmid free cells in lumen
  const double d1L = DL; // death rate of plasmid bearing cells in lumen
  const double AinL = 30.0; // input antibiotic concentraton in lumen

  // at wall
  const double cW = pow(10, -8.8); // conjugation factor at wall
  const double DW = 0.25; // turnover at wall
  const double SWin = 10; // input resource concentration at the wall
  const double KWL0 = 1e-9; // dettaching rate of plasmid free cells of wall
  const double KWL1 = 1e-9; // dettaching rate of plasmid bearing cells of wall
  const double d0W = DW; // death rate of plasmide free at wall
  const double d1W = DW; // deat rate of plasmid bearing at wall
  const double AinW = 0; // input antibiotic concecntration at wall 

  // vectors containg the wall and lumen specific parameters
  const std::vector< double > D = {DW, DL};
  const std::vector< double > Sin = {SWin, SLin};
  const std::vector< double > c = {cW, cL};
  const std::vector< double > Ain = {AinW, AinL};

//*** ODE description *********

  // enumerating wall and lumen in location
  enum location {Wall, Lumen};


// function to update the sate of the concentration, with the exception of migration which will be added seperately
void update_state(const std::vector<double> &x, std::vector<double> &dxdt, location loc) 
{
  // define both indexes
  size_t init_index = 0;
  if (loc == Lumen) init_index = 4;

  // define resource, plasmid free and plamsid bearing concentrations
  double R =  x[init_index];
  double N0 = x[init_index + 1];
  double N1 = x[init_index + 2];
  double A = x[init_index + 3];  

  // growth rates of plasmid free and plasmid bearing cells
  const double psi0 = ((r0 * R) / (K0 + R));
  const double psi1 = ((r1 * R) / (K1 + R));
  const double die0 = ((d0 * A) / (H0 + A));
  const double die1 = ((d1 * A) / (H1 + A));

  // acces the above vectors to calculate location specific difference and apply this to the correct population
  dxdt[init_index] = D[loc] * (Sin[loc] - R) - e0 * psi0 * N0 - e1 * psi1 * N1;

  dxdt[init_index + 1] = psi0 * N0 - D[loc] * N0 + l * N1 - c[loc] * N0 * N1 - die0 * N0;

  dxdt[init_index + 2] = psi1 * N1 - D[loc] * N1 - l * N1 + c[loc] * N0 * N1 - die1 * N1;

  dxdt[init_index + 3] = D[loc] * (Ain[loc] - A) - Adisp0 * N0 * die0 - Adisp1 * N1 * die1;

  return;
}

// function to add the change due to migration
void exchange_cells(const std::vector<double> &x, std::vector<double> &dxdt) 
{
  double N0W = x[1]; // plasmid free at wall
  double N1W = x[2]; // plasmid bearing at wall

  double N0L = x[5]; // plasmid free in lumen
  double N1L = x[6]; // plasmid bearing in lumen

  // compensate aswell for the difference in volume
  // by first measuring how many cells will transfer
  // and dividing this by their destination volume
  dxdt[1] += KLW0 * N0L  - KWL0 * N0W; 
  dxdt[2] += KLW1 * N1L - KWL0 * N1W; 
  dxdt[5] += KWL0 * N0W - KLW0 * N0L; 
  dxdt[6] += KWL1 * N1W - KLW1 * N1L;

  return;
}

// right hand site definition by accesing previously defined functions
void rhs(const double &t, const std::vector<double> &x, std::vector<double> &dxdt)
{
  update_state(x, dxdt, Wall);
  update_state(x, dxdt, Lumen);

  exchange_cells(x, dxdt);
}

//*** parameters of the integration algorithm *********

    const double kdShrinkMax = 0.1; // decrease step size by no more than this factor
    const double kdGrowMax = 1.3; // increase step size by no more than this factor
    const double kdSafety = 0.9; // safety factor in adaptive stepsize control
    const double kdMinH = 1.0e-6; // minimum step size
    const int nvar = 8; // number of variables
    const double dt0 = 0.0005; // initial time step size
    const double dtsav = 0.05; // save data after time steps
    const double tEnd = 100000.0; // end time
    const double tolerance = 1.0e-6; // acceptable local error during numerical integration

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

// analysis function to allow for easy looping over variables
void do_analysis(std::string output_filename, const std::vector<double>& pars) 
{
  // give initial values
  double initialN0 = 1.0;
  double initialN1 = 1.0e-5;

  std::vector<double> x(nvar);
  x[0] = SWin;
  x[1] = initialN0;
  x[2] = initialN1;
  x[3] = AinW;
  x[4] = SLin;
  x[5] = initialN0;
  x[6] = initialN1;
  x[7] = AinL;
  std::vector<double> dxdt(nvar);
  rhs(0.0, x, dxdt);

  // start numerical integration
  int nOK = 0, nStep = 0;
  for(double t = 0.0, tsav = 0.0, dt = dt0; t < tEnd; ++nStep)
  {
    if(BogackiShampineStepper(t, x, dxdt, dt))
    {
      ++nOK;
    }
    if (fabs(dxdt[1]) < 1.0e-6 && fabs(dxdt[2]) < 1.0e-6 && fabs(dxdt[5]) < 1.0e-6 && fabs(dxdt[6]) < 1.0e-6)
    {
       break; // if change very little, stop 
    }
  }

  std::cout << "At wall:  \n" << "Antibiotic: " << x[3] << "\nResource: " << x[0] << "\nPlasmid Free: " << x[1] << "\nPlasmid bearing: " << x[2] << "\n" << std::endl;
  std::cout << "At Lumen:  \n"<< "Antibiotic: " << x[7] << "\nResource: " << x[4] << "\nPlasmid Free: " << x[5] << "\nPlasmid bearing: " << x[6] << "\n" << std::endl;
  
  if (fabs(dxdt[2]) < 1.0e-6 && x[2] < initialN1 * 1000)  // if the almost equillibrium value is reached and the value is not 1000
  {                                                       // times larger than intial value, set it to 0 (cause it will go extinct) 
    x[2] = 0;
  }
  if (fabs(dxdt[6]) < 1.0e-6 && x[5] < initialN1 * 1000)
  {
    x[5] = 0;
  }
  if (fabs(dxdt[1]) < 1.0e-6 && x[1] < initialN0 * 1000)
  {
    x[1] = 0;
  }
  if (fabs(dxdt[5]) < 1.0e-6 && x[4] < initialN0 * 1000)
  {
    x[4] = 0;
  }

  std::ofstream ofs(output_filename.c_str(), std::ios::app);
  if(!ofs.is_open())
  {
    throw std::runtime_error("unable to open file.\n");
  }

  // write the paramater(s) followed by the concentration of the cell types
  for(size_t i = 0; i < pars.size(); ++i )
  {
    ofs << pars[i] << ',';
  } 
  ofs << x[1] << ',' << "N0" << ',' << "Wall" << '\n';

  for(size_t i = 0; i < pars.size(); ++i )
  {
    ofs << pars[i] << ',';
  } 
  ofs << x[2] << ',' << "N1" << ',' << "Wall" << '\n';

  for(size_t i = 0; i < pars.size(); ++i )
  {
    ofs << pars[i] << ',';
  } 
  ofs << x[5] << ',' << "N0" << ',' << "Lumen" << '\n';

  for(size_t i = 0; i < pars.size(); ++i )
  {
    ofs << pars[i] << ',';
  } 
  ofs << x[6] << ',' << "N1" << ',' << "Lumen" << '\n';
  
  ofs.close();
}


//*** function main() *********

int main()
{
  try 
  {
    // provide file name
    std::string file_name = "Main2AC4.csv";
    std::ofstream ofs(file_name.c_str());

    // give first row with variable names
    ofs  << "pars1"  << ',' << "popsize" << ',' << "population" << ',' << "location" << "\n";
    ofs.close();

    // do analysis for different values of the attachement rates, and increased conjugation  
    for (double local_KLW0 = 1e-10; local_KLW0 < pow(10, -3); local_KLW0 *= 3) 
    {
      KLW1 = local_KLW0; 
      std::vector<double> pars = {KLW1};
      std::cout << "migration of L to W of plasmid bearing at rate: " << KLW1 << "\n";
      do_analysis(file_name, pars);
      /*std::cout << KLW1 << ',' << KLW0 << "\n";*/
    }
  }
  catch(std::exception &error)
  {
  std::cerr << "error: " << error.what();
  exit(EXIT_FAILURE);
  }
  return 0;
}
