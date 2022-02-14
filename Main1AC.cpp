#include <iostream>
#include <fstream>
#include <vector>
#include <exception>
#include <cstdlib>
#include <cmath> 

//*** THE EFFECT OF WALL ATTACHMENT ON PLASMID PERSISTENCE SCRIPT 3 *********

// This script provides a way to record equillibrium values while varying a parameter value
// the below found script produces the left graph in figure 1 in the report
// But it was modified to create figure 1 and 2 
// This script has an extinction criteria so that when the plasmid won't persist it's value 
// will be set to 0 

//*** model parameters *********

    const double K0 = 4.0; // half saturation constant of plasmid bearing 
    const double K1 = 4.0; // half saturation constant of plasmid bearing
    const double e1 = 6.25 * 1e-7; // resource needed to divide once for plasmid bearing cells
    const double e0 = 6.25 * 1e-7; // resource needed to divide once for plasmid free cells
    const double l = 1e-3; // loss of plasmid by divsion
    const double c = 1e-9; // conjugation rate
    const double D = 0.25; // turnover rate
    const double r0 = 0.738; // max growth rate of plasmid free
    const double r1 = 0.6642; // max growth rate of plasmid bearing
    double Sin;  // input resource concentration  
    const double d1 = 1e-9; // death of plasmid bearing cell by antibiotic
    const double d0 = 1e-2; // death of plasmid free cell by antibiotic
    const double Ain = 0; // input antibiotic concentration 
    const double H0 = 4.0; // half saturation constant of plasmid free for antibiotics
    const double H1 = 4.0; // half saturation constant of plasmid bearing for antibiotics
    const double Adisp = 1e-5; // antibiotic needed for the death of a cell

//*** ODE description *********

void rhs(const double &t, const std::vector<double> &x, std::vector<double> &dxdt)
{
    double R = x[0]; // Resource concentration
    double N0 = x[1]; // Plasmid free concentration
    double N1 = x[2]; // Plasmid bearing concentration
    double A = x[3]; // Antibiotic concentration 

    double psi0 = ((r0 * R) / (K0 + R)); // growth rate of plasmid free bacteria
    double psi1 = ((r1 * R) / (K1 + R)); // growth rate of plasmid bearing bacteria
    double die0 = ((d0 * A) / (H0 + A)); // die rate due to antibiotics for plasmid free
    double die1 = ((d1 * A) / (H1 + A)); // die rate due to antibiotics for plasmid bearing
    
    dxdt[0] = D * (Sin - R) - e0 * psi0 * N0 - e1 * psi1 * N1; // differential equation of the resource
    dxdt[1] = psi0 * N0  - D * N0 + l * N1 - c * N0 * N1 - die0 * N0; // differential equation of the plasmid free cells
    dxdt[2] = psi1 * N1 - D * N1 - l * N1 + c * N0 * N1 - die1 * N1; // differential equation of the plasmid bearing cells
    dxdt[3] = D * (Ain - A) - Adisp * die0 * N0 - Adisp * die1 * N1; // differential equation of the antibiotic 
 }

//*** parameters of the integration algorithm *********

    const double kdShrinkMax = 0.1; // decrease step size by no more than this factor
    const double kdGrowMax = 1.3; // increase step size by no more than this factor
    const double kdSafety = 0.9; // safety factor in adaptive stepsize control
    const double kdMinH = 1.0e-6; // minimum step size
    const int nvar = 4; // number of variables
    const double dt0 = 0.0005; // initial time step size
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

// analysis function to allow for easy looping over variables
void do_analysis(std::string output_filename, const std::vector<double>& pars) 
{
    // give initial values 
    double initialN0 = 1.0;
    double initialN1 = 1.0e-5;

    std::vector<double> x(nvar);
    x[0] = Sin;
    x[1] = initialN0;
    x[2] = initialN1;
    x[3] = Ain;
    std::vector<double> dxdt(nvar);
    rhs(0.0, x, dxdt);

    // start numerical integration
    int nOK = 0, nStep = 0;
    for(double t = 0.0, dt = dt0; t < tEnd; ++nStep)
    {
        if(BogackiShampineStepper(t, x, dxdt, dt))
        {
            ++nOK;
        }
        if (fabs(dxdt[1]) < 1.0e-6 && fabs(dxdt[2]) < 1.0e-6) // if change very little, stop (equilibrium most likely reached)
        {
            break;
        }
    }


    if (fabs(dxdt[2]) < 1.0e-6 && x[2] < initialN1 * 1000) // if the equillibrium value is reached and the value is not 1000
    {                                                      // times larger than intial value, set it to 0 (cause it will go extinct) 
        x[2] = 0;
    }
    if (fabs(dxdt[1]) < 1.0e-6 && x[1] < initialN0 * 1000)
    {
        x[1] = 0;
    }

    // open file
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
    ofs << x[1] << ',' << "N0" << ',' << dxdt[1] << '\n';

    for(size_t i = 0; i < pars.size(); ++i )
    {
        ofs << pars[i] << ',';
    } 
    ofs << x[2] << ',' << "N1"  << ',' << dxdt[2] << '\n';


    ofs.close();
}

//*** function main() *********

int main()
{
    try 
    {   
        // provide file name
        std::string file_name = "Main1AC9.csv";
        std::ofstream ofs(file_name.c_str());

        // give first row with variable names
        ofs  << "pars" << ',' << "popsize" << ',' << "population" << ',' << "DXDT" << "\n";
        ofs.close();

        // do analysis for different values of the resource input concentration 
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