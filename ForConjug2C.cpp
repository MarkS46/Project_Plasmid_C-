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
    std::string name = "ForConjug2C3.csv";
//*** model parameters *********

    // general 
    const double K0 = 4.0; // half saturation constant of plasmid free
    const double K1 = 4.0; // half saturation constant of plasmid bearing
    const double e1 = 6.25 * 1e-7; // resource needed to divide once for plasmid bearing
    const double e0 = 6.25 * 1e-7; // resource needed to divide once for plasmid free
    const double l = 1e-3; // loss of plasmid 
    const double r0 = 0.738; // max growth rate of plasmid free
    const double r1 = 0.6642; // max growth rate of plasmid bearing
    const double alfa = 1 - (r1/r0); // selective advantage of plasmid free cells

    // in lumen
    const double cL = pow(10, -9); // conjugation factor in lumen
    const double DL = 0.45; // flow rate in lumen
    const double SLin = 33; // resource coming in to the lumen
    const double KLW0 = 1e-9; // attaching of plasmid free cells to wall
    const double KLW1 = 1e-9; // attaching of plasmid bearing cells to wall
    const double d0L = DL; // death rate of plasmid free cells in  lumen
    const double d1L = DL; // death rate of plasmid bearing cells in lumen 

    // at wall
    double cW; // conjugation factor at wall
    const double DW = 0.35; // flow rate at wall     
    const double d0W = DW; // death rate of plasmide free at wall 
    const double d1W = DW; // deat rate of plasmid bearing at wall
    const double SWin = 9; // resource coming to the wall
    const double KWL0 = 1e-9; // dettaching of plasmid free cells of wall
    const double KWL1 = 1e-9; // dettaching of plasmid bearing cells of wall

//*** ODE description *********

void rhs(const double &t, const std::vector<double> &x, std::vector<double> &dxdt, double cW)
{
    // main variables

    double RW = x[0]; // resource at wall
    double N0W = x[1]; // plasmid free at wall
    double N1W = x[2]; // plasmid bearing at wall 
    double RL = x[3]; // resource at wall
    double N0L = x[4]; // plasmid free in lumen
    double N1L = x[5]; // plasmid bearing in lumen
    
    // at wall

    double psi0W = ((r0 * RW) / (K0 + RW)); 
    double psi1W = ((r1 * RW) / (K1 + RW)); 

    dxdt[0] = DW * (SWin - RW) - e0 * psi0W * N0W - e1 * psi1W * N1W; // differential equation of the resource at the wall
    dxdt[1] = psi0W * N0W  - d0W * N0W + l * N1W - cW * N0W * N1W + KLW0 * N0L - KWL0 * N0W; // differential equation of the plasmid free cells at the wall
    dxdt[2] = psi1W * N1W - d1W * N1W - l * N1W + cW * N0W * N1W + KLW1 * N1L - KWL1 * N1W; // differential equation of the plasmid bearing cells at the wall
 
    // in lumen

    double psi0L = ((r0 * RL) / (K0 + RL)); 
    double psi1L = ((r1 * RL) / (K1 + RL)); 

    dxdt[3] = DL * (SLin - RL) - e0 * psi0L * N0L - e1 * psi1L * N1L; // differential equation of the resource at the wall
    dxdt[4] = psi0L * N0L  - d0L * N0L + l * N1L - cL * N0L * N1L + KWL0 * N0W - KLW0 * N0L; // differential equation of the plasmid free cells at the wall
    dxdt[5] = psi1L * N1L - d1L * N1L - l * N1L + cL * N0L * N1L + KWL1 * N1W - KLW1 * N1L; // differential equation of the plasmid bearing cells at the wall
}

//*** ODE integration routine *********

    const double kdShrinkMax = 0.1; // decrease step size by no more than this factor
    const double kdGrowMax = 1.3; // increase step size by no more than this factor
    const double kdSafety = 0.9; // safety factor in adaptive stepsize control
    const double kdMinH = 1.0e-6; // minimum step size

//*** The Bogacki-Shampine stepper *********

bool BogackiShampineStepper(double &t, std::vector<double> &x,  std::vector<double>&dxdt, double &h, double cW)
{

    // step 2
    std::vector<double> xtmp(nvar);
    for(int i = 0; i < nvar; ++i)
    {
        xtmp[i] = x[i] + 0.5 * h * dxdt[i];
    }
    std::vector<double> dxdt2(nvar);
    rhs(t + 0.5 * h, xtmp, dxdt2, cW);

    // step 3
    for(int i = 0; i < nvar; ++i)
    {
        xtmp[i] = x[i] + 0.75 * h * dxdt2[i];
    }
    std::vector<double> dxdt3(nvar);
    rhs(t + 0.75 * h, xtmp, dxdt3, cW);

    // step 4
    for(int i = 0; i < nvar; ++i)
    {
        xtmp[i] = x[i] +(1.0/9.0) * h *(2 * dxdt[i] + 3 * dxdt2[i] + 4 * dxdt3[i]) ;
    }
    std::vector<double> dxdt4(nvar);
    rhs(t + h, xtmp, dxdt4, cW);

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

//*** function main() *********

int main()
{
    try {
        // open data file
        std::ofstream ofs(name);
        if(!ofs.is_open())
            throw std::runtime_error("unable to open file.\n");
           

        // give first row with variable names
        ofs << "popsize" << ',' << "population" << ',' << "location" << ',' << "conjugationW" << "\n";


        double initialN0 = 1.0;
        double initialN1 = 1e-5;

        for (double i = 4.0; i < 12.0; i += 0.1)
        {
            std::cout << "x";
            double cWval = pow(10, -i);
            // give initial values
            std::vector<double> x(nvar);
            x[0] = SWin;
            x[1] = initialN0;
            x[2] = initialN1;
            x[3] = SLin;
            x[4] = initialN0;
            x[5] = initialN1;
            std::vector<double> dxdt(nvar);
            rhs(0.0, x, dxdt, cWval);

            // start numerical integration
            int nOK = 0, nStep = 0;
            double dtMin = dt0, dtMax = kdMinH;
            double t = 0.0;
            double dt;
            for(t, dt = dt0; t < tEnd; ++nStep) 
            {
                if(BogackiShampineStepper(t, x, dxdt, dt, cWval))
                {
                    ++nOK;
                }
                if (fabs(dxdt[1]) < 1.0e-4 && fabs(dxdt[2]) < 1.0e-4 && fabs(dxdt[4]) < 1.0e-4 && fabs(dxdt[5]) < 1.0e-4)
                {
                    break;
                }
            }

            if (dxdt[2] > -1.0e-4 && x[2] < 1)
            {
                x[2] = 0;
            }
            if (dxdt[5] > -1.0e-4 && x[5] < 1) 
            {
                x[5] = 0;
            }
            if (dxdt[1] > -1.0e-4 && x[1] < 1)
            {
                x[1] = 0;
            }
            if (dxdt[4] > -1.0e-4 && x[4] < 1) 
            {
                x[4] = 0;
            }



            ofs << x[1] << ',' << "N0" << ',' << "Wall" << ',' << cWval << '\n' 
                << x[2] << ',' << "N1" << ',' << "Wall" << ',' << cWval << '\n'
                << x[4] << ',' << "N0" << ',' << "Lumen" << ',' << cWval << '\n'
                << x[5] << ',' << "N1" << ',' << "Lumen" << ',' << cWval << '\n';  

        }
        ofs.close();
    }
    catch(std::exception &error) 
    {
    std::cerr << "error: " << error.what();
    exit(EXIT_FAILURE);
    }
    return 0;
}