#include <iostream>
#include <fstream>
#include <vector>
#include <exception>
#include <cstdlib>
#include <cmath> 

//*** parameters of the integration algorithm *********

    const int nvar = 3; // number of variables
    const double dt0 = 0.001; // initial time step size
    const double dtsav = 0.05; // save data after time steps
    const double tEnd = 1000.0; // end time
    const double tolerance = 1.0e-6; // acceptable local error during numerical integration

//*** model parameters *********

    const double K0 = 4.0; // half saturation constant of plasmid bearing 
    const double K1 = 4.0; // half saturation constant of plasmid bearing
    const double e1 = 6.25 * 1e-7; // resource needed to divide once for plasmid bearing
    const double e0 = 6.25 * 1e-7; // resource needed to divide once for plasmid free
    const double l = 1e-3; // loss of plasmid 
    double c; // conjugation factor
    const double D = 0.55; // flow rate
    const double r0 = 0.738; // growth rate of plasmid free
    const double r1 = 0.6642; // growth rate of plasmid bearing
    double Sin; //inflow concentration of resource   
    const double alfa = 1 - (r1/r0); // selective advantage of plasmid free cells


//*** ODE description *********
void rhs(const double &t, const std::vector<double> &x, std::vector<double> &dxdt, double S, double c)
{

    double R = x[0]; // Resource
    double N0 = x[1]; // Plasmid free
    double N1 = x[2]; // Plasmid bearing 

    double psi0 = ((r0 * R) / (K0 + R)); // growth rate of plasmid free bacteria
    double psi1 = ((r1 * R) / (K1 + R)); // growth rate of plasmid bearing bacteria

    dxdt[0] = D * (S - R) - e0 * psi0 * N0 - e1 * psi1 * N1; // differential equation of the resource
    dxdt[1] = psi0 * N0  - D * N0 + l * N1 - c * N0 * N1; // differential equation of the plasmid free cells
    dxdt[2] = psi1 * N1 - D * N1 - l * N1 + c * N0 * N1; // differential equation of the plasmid bearing cells
 }
//*** ODE integration routine *********

    const double kdShrinkMax = 0.1; // decrease step size by no more than this factor
    const double kdGrowMax = 1.3; // increase step size by no more than this factor
    const double kdSafety = 0.9; // safety factor in adaptive stepsize control
    const double kdMinH = 1.0e-6; // minimum step size

//*** The Bogacki-Shampine stepper *********

bool BogackiShampineStepper(double &t, std::vector<double> &x,  std::vector<double>&dxdt, double &h, double S, double c)
{

    // step 2
    std::vector<double> xtmp(nvar);
    for(int i = 0; i < nvar; ++i)
    {
        xtmp[i] = x[i] + 0.5 * h * dxdt[i];
    }
    std::vector<double> dxdt2(nvar);
    rhs(t + 0.5 * h, xtmp, dxdt2, S, c);

    // step 3
    for(int i = 0; i < nvar; ++i)
    {
        xtmp[i] = x[i] + 0.75 * h * dxdt2[i];
    }
    std::vector<double> dxdt3(nvar);
    rhs(t + 0.75 * h, xtmp, dxdt3, S, c);

    // step 4
    for(int i = 0; i < nvar; ++i)
    {
        xtmp[i] = x[i] +(1.0/9.0) * h *(2 * dxdt[i] + 3 * dxdt2[i] + 4 * dxdt3[i]) ;
    }
    std::vector<double> dxdt4(nvar);
    rhs(t + h, xtmp, dxdt4, S, c);

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

 //*** function main() *********

int main()
{
    try {
        // open data file
        std::ofstream ofs("ForConjug&Sin1C.csv");
        if(!ofs.is_open())
        {
            throw std::runtime_error("unable to open file.\n");
        }   

        // give first row with variable names
        ofs << "popsize" << ',' << "popid" << ',' << "Sin" << ',' << "c" << "\n";

        //loop over the c values to be tested
        for (double i = 5; i < 12;  i += 0.5 )
        {
            std::cout << "x";
            c = pow(10, -i);

            // loop over the Sin values to be tested 
            for (Sin = 10.0; Sin < 151.0; ++Sin)
            {  
                // give initial values
                std::vector<double> x(nvar);
                x[0] = Sin;
                x[1] = 1.0;
                x[2] = 1.0;
                std::vector<double> dxdt(nvar);
                rhs(0.0, x, dxdt, Sin, c);

                // start numerical integration
                int nOK = 0, nStep = 0;
                double dtMin = dt0, dtMax = kdMinH;
                double t;   
                for(double t = 0.0, tsav = 0.0, dt = dt0; t < tEnd; ++nStep) 
                {
                    if(BogackiShampineStepper(t, x, dxdt, dt, Sin, c))
                    {
                        ++nOK;
                    }

                    if (fabs(dxdt[1]) < 1.0e-6 && fabs(dxdt[2]) < 1.0e-6)
                    {
                        break; 
                    } 
           
                    if(dt < dtMin)
                    {
                        dtMin = dt;
                    }
                    else if(dt > dtMax)
                    {
                        dtMax = dt;
                    }
                }

                /* 
                // report integration data
                std::cout << "integration complete.\n"
                << "number of steps : " << nStep << '\n'
                << "proportion bad steps : " << 1.0 - nOK * 1.0 / nStep << '\n'
                << "average step size : " << tEnd / nStep << '\n'
                << "min step size : " << dtMin << '\n'
                << "max step size : " << dtMax << "\n\n";

                // return final population sizes
                std::cout << "plasmid free = " << x[1] << " plasmid bearing = " << x[2] <<"\n";

                // return alpha
                std::cout << "alpha = "  << alfa << "\n";

                // return concentration of plasmid bearing bacteria
                std::cout << "F+ = " << (x[2] / (x[1] + x[2]))*100 << "\n"; 

                // report wether equilibrium (as indicated above) was found
                if (fabs(dxdt[1]) < 1.0e-6 && fabs(dxdt[2]) < 1.0e-6)
                {
                    std::cout << "equilibrium was reached for all populations" <<  "\n\n";
                }
                else
                {
                    std::cout << "for one or more populations the equillibrium was not reached" << "\n\n";
                } */

                 // write data to file
                ofs << x[1] << ',' << "N0" << ',' << Sin << ',' << c << '\n'
                    << x[2] << ',' << "N1" << ',' << Sin << ',' << c << '\n';
            }
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