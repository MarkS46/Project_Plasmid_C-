#include <iostream>
#include <fstream>
#include <vector>
#include <exception>
#include <cstdlib>
#include <cmath> 
#include <string>

//*** THE EFFECT OF WALL ATTACHMENT ON PLASMID PERSISTENCE SCRIPT 2 *********

// this script provides a time simulation of the two compartment model
// right hand side definition is here very bulky but will be improved in a more advanced script later
// this script again was not used for visualisation in the report but was examined to check the dynamics
// Currently the parameters are set to a scenario with a net flux of 0 
// in which in both compartments the plasmid will not persist 
// increasing the resource input concentration in the lumen to 40 and in the wall to 25 will show persistence
// In this script it is assumed that both compartments are of the same size 
// in Main2C.cpp this will be accounted for more properly 

//*** model parameters *********

    // general 
    const double K0 = 4.0; // half saturation constant of plasmid free
    const double K1 = 4.0; // half saturation constant of plasmid bearing
    const double e1 = 1; // resource needed to divide once for plasmid bearing cells
    const double e0 = 1; // resource needed to divide once for plasmid free cells
    const double l = 0; // loss of plasmid by divsion
    const double r0 = 0.4; // max growth rate of plasmid free 
    const double r1 = 0.3; // max growth rate of plasmid bearing
    const double alfa = 1 - (r1/r0); // selective advantage of plasmid free cells

    // in lumen
    const double cL = 0.0; // conjugation rate in lumen
    const double DL = 0.45; // turnover rate in lumen
    const double SLin = 25000; // resource coming in to the lumen
    const double KLW0 = 0.004; // attaching of plasmid free cells to wall
    const double KLW1 = 0.004; // attaching of plasmid bearing cells to wall
    const double d0L = 0.99; // death rate of plasmid free cells in  lumen
    const double d1L = 0.99; // death rate of plasmid bearing cells in lumen 

    // at wall
    const double cW = 0.0; // conjugation rate at wall
    const double DW = 0.45; // turnover rate at wall     
    const double SWin = 25000; // resource coming to the wall
    const double KWL0 = 0.004; // dettaching of plasmid free cells of wall
    const double KWL1 = 0.004; // dettaching of plasmid bearing cells of wall
    const double d0W = 0.99; // death rate of plasmide free at wall 
    const double d1W = 0.99; // death rate of plasmid bearing at wall

//*** ODE description *********

void rhs(const double &t, const std::vector<double> &x, std::vector<double> &dxdt)
{
    // main variables

    double RW = x[0]; // resource concentration at wall
    double N0W = x[1]; // plasmid free concentration at wall
    double N1W = x[2]; // plasmid bearing concentraiton at wall 
    double RL = x[3]; // resource concentration at wall
    double N0L = x[4]; // plasmid free concentration in lumen
    double N1L = x[5]; // plasmid bearing concentration in lumen
    
    // at wall

    double psi0W = ((r0 * RW) / (K0 + RW)); // growth rate of plasmid free cells at the wall
    double psi1W = ((r1 * RW) / (K1 + RW)); // growth rate of plasmid bearing cells at the wall

    dxdt[0] = DW * (SWin - RW) - e0 * psi0W * N0W - e1 * psi1W * N1W; // differential equation of the resource at the wall
    dxdt[1] = psi0W * N0W  - d0W * N0W + l * N1W - cW * N0W * N1W + KLW0 * (N0L/(N0L + N1L)) * N0L - KWL0 * (N0W/(N0W + N1W)) * N0W; // differential equation of the plasmid free cell concentration at the wall
    dxdt[2] = psi1W * N1W - d1W * N1W - l * N1W + cW * N0W * N1W + KLW1 * (N1L/(N0L + N1L)) * N1L - KWL1 * (N1W/(N0W + N1W)) * N1W; // differential equation of the plasmid bearing cell concentration at the wall
 
    // in lumen

    double psi0L = ((r0 * RL) / (K0 + RL)); // growth rate of plasmid free cells at the wall
    double psi1L = ((r1 * RL) / (K1 + RL)); // growth rate of plasmid bearing cells at the wall 

    dxdt[3] = DL * (SLin - RL) - e0 * psi0L * N0L - e1 * psi1L * N1L; // differential equation of the resource concentration at the wall
    dxdt[4] = psi0L * N0L  - d0L * N0L + l * N1L - cL * N0L * N1L + KWL0 * (N0W/(N0W + N1W)) * N0W - KLW0 *(N0L/(N0L + N1L)) * N0L; // differential equation of the plasmid free cell concenctration at the wall
    dxdt[5] = psi1L * N1L - d1L * N1L - l * N1L + cL * N0L * N1L + KWL1 * (N1W/(N0W + N1W)) * N1W - KLW1 * (N1L/(N0L + N1L)) * N1L; // differential equation of the plasmid bearing cell concentration at the wall
}
//*** parameters of the integration algorithm *********

    const double kdShrinkMax = 0.1; // decrease step size by no more than this factor
    const double kdGrowMax = 1.3; // increase step size by no more than this factor
    const double kdSafety = 0.9; // safety factor in adaptive stepsize control
    const double kdMinH = 1.0e-6; // minimum step size
    const int nvar = 6; // number of variables
    const double dt0 = 0.0005; // initial time step size
    const double dtsav = 0.05; // save data after time steps
    const double tEnd = 5000.0; // end time
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
        return true;
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
    try 
    {
        // file name
        std::string name = "Basic2C2.csv";

        // open data file
        std::ofstream ofs(name);
        if(!ofs.is_open())
            throw std::runtime_error("unable to open file.\n");
           
        // give first row with variable names
        ofs << "t" << ',' << "popsize" << ',' << "population" << ',' << "location" << "\n";

        // give initial values
        double initialN0 = 250;
        double initialN1 = 250;

        std::vector<double> x(nvar);
        x[0] = 62500;
        x[1] = initialN0;
        x[2] = 0;
        x[3] = 62500;
        x[4] = initialN0;
        x[5] = 0;
        std::vector<double> dxdt(nvar);
        rhs(0.0, x, dxdt);

        // initiate parameters for the integration
        int nOK = 0, nStep = 0;
        double dtMin = dt0, dtMax = kdMinH; 

        bool input = false; 
        // start numerical integration
        for(double t = 0.0, tsav = 0.0, dt = dt0; t < tEnd; ++nStep) 
        {
            
            std::cout << x[3] << "  " << x[0] << std::endl; 

            
            if(BogackiShampineStepper(t, x, dxdt, dt))
            {
                ++nOK;
            }
            if (t > 1000 && input == false)
            {
                std::cout << (x[1]) <<std::endl;
                std::cout << (x[4]) <<std::endl;
                x[2] = 0.01* x[1];
                x[5] = 0.01* x[4];
                std::cout << (x[2]) <<std::endl;
                std::cout << (x[5]) <<std::endl;
                input = true;
            }
            /*if (fabs(dxdt[1]) < 1.0e-6 && fabs(dxdt[2]) < 1.0e-6 && fabs(dxdt[4]) < 1.0e-6 && fabs(dxdt[5]) < 1.0e-6 )
            {
                break;
            }*/
            if(dt < dtMin) // keep track of the smallest step size
            {
                dtMin = dt;
            }
            else if(dt > dtMax) // keep track of the largest step size
            {
                dtMax = dt;
            }


            double input1;
            double input2;
            double input4;
            double input5;

            if (x[1] < 1e-10)
            {
                input1 = 1e-10;
            } 
            else
            {
               input1 = x[1];
            }

            if (x[2] < 1e-10)
            {
                input2 = 1e-10;
            } 
            else
            {
                input2 = x[2];
            }

            if (x[4] < 1e-10)
            {
                input4 = 1e-10;
            } 
            else
            {
                input4 = x[4];
            }

            if (x[5] < 1e-10)
            {
                input5 = 1e-10;
            } 
            else
            {
                input5 = x[5];
            }

            if(t > tsav) 
            {
                ofs << t << ',' << input1 << ',' << "N0" << ',' << "Wall" << '\n' 
                    << t << ',' << input2 << ',' << "N1" << ',' << "Wall" << '\n'
                    << t << ',' << input4 << ',' << "N0" << ',' << "Lumen" << '\n'
                    << t << ',' << input5 << ',' << "N1" << ',' << "Lumen" << '\n';  
                tsav += dtsav;
            }
        }

        // return final concentration
        std::cout << " \n at the wall: \n" 
                  << "plasmid free = " << x[1]  << "   plasmid bearing = " << x[2] << "   resource = " << x[0] 
                  << "\n in the lumen: \n" 
                  << "plasmid free = " << x[4]  << "   plasmid bearing = " << x[5] << "   resource = " << x[3] << "\n"; 
        
        // report integration data
        std::cout << "\nintegration complete.\n"
        << "number of steps : " << nStep << '\n'
        << "proportion bad steps : " << 1.0 - nOK * 1.0 / nStep << '\n'
        << "average step size : " << tEnd / nStep << '\n'
        << "min step size : " << dtMin << '\n'
        << "max step size : " << dtMax << '\n';

        // return alfa
        std::cout << "alpha = "  << alfa << "\n"; 
        
        ofs.close();
    }
    catch(std::exception &error) 
    {
    std::cerr << "error: " << error.what();
    exit(EXIT_FAILURE);
    }
    return 0;
}