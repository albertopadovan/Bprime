#include <Eigen/Dense>
#include <iostream>
#include <math.h>
 
#include "mutation++.h"
 
#ifdef _GNU_SOURCE
#include <fenv.h>
#endif
 
using namespace std;
using namespace Mutation;
using namespace Mutation::Thermodynamics;
using namespace Mutation::Utilities;
 
// Simply stores the command line options
typedef struct {
    double T1;
    double T2;
    double dT;
 
    double P1;
    double P2; 
    double dP;

    double Bg1;
    double Bg2;
    double dBg;
 
    std::string mixture;
    std::string boundary_layer_comp;
    std::string pyrolysis_composition;
 
    bool pyrolysis_exist = false;
} Options;
 
// Checks if an option is present
bool optionExists(int argc, char** argv, const std::string& option) {
    return (std::find(argv, argv + argc, option) != argv + argc);
}
 
// Returns the value associated with a particular option
std::string getOption(int argc, char** argv, const std::string& option) {
    std::string value;
    char** ptr = std::find(argv, argv + argc, option);
 
    if (ptr == argv + argc || ptr + 1 == argv + argc)
        value = "";
    else
        value = *(ptr + 1);
 
    return value;
}
 
// Prints the program's usage information and exits.
void printHelpMessage(const char* const name) {
    std::string tab("    ");
 
    cout.setf(std::ios::left, std::ios::adjustfield);
 
    cout << endl;
    cout << "Usage: " << name << " [OPTIONS] mixture" << endl;
    cout << "Compute the non-dimensional mass blowing rate for mixture "
         << "over a set of temperatures and pressure using the Mutation++ "
            "library."
         << endl;
    cout << endl;
    cout << tab << "-h, --help          prints this help message" << endl;
    cout << tab
         << "-T                  temperature range in K \"T1:dT:T2\" or simply "
            "T (default = 300:100:5000 K)"
         << endl;
    cout << tab << "-P                  pressure in Pa P (default = 1 atm)"
         << endl;
    cout << tab
         << "-b                  pyrolysis non-dimensional mass blowing rate "
            "(default = 0)"
         << endl;
    cout << tab << "-m                  mixture name" << endl;
    cout << tab << "-bl                 boundary layer edge composition name"
         << endl;
    cout << tab
         << "-py                 pyrolysis composition name (default = null)"
         << endl;
 
    cout << endl;
    cout << "Example:" << endl;
    cout
        << tab << name
        << " -T 300:100:5000 -P 101325 -b 10 -m carbonPhenol -bl BLedge -py Gas"
        << endl;
    cout << endl;
    cout << "Mixture file:" << endl;
    cout << tab << "carbonPhenol - corresponds to the name of the mixture"
         << endl;
    cout << tab
         << "BLedge - corresponds to the boundary layer edge elemental "
            "composition"
         << endl;
    cout << tab
         << "Gas - corresponds to the pyrolysis elemental gas composition"
         << endl;
    cout << endl;
 
    exit(0);
}
 
// Parses a temperature or pressure range
bool parseRange(const std::string& range, double& x1, double& x2, double& dx) {
    std::vector<std::string> tokens;
    String::tokenize(range, tokens, ":");
 
    if (!String::isNumeric(tokens)) return false;
 
    switch (tokens.size()) {
        case 1:
            x1 = atof(tokens[0].c_str());
            x2 = x1;
            dx = 1.0;
            break;
        case 3:
            x1 = atof(tokens[0].c_str());
            x2 = atof(tokens[2].c_str());
            dx = atof(tokens[1].c_str());
            break;
        default:
            return false;
    }
 
    if (dx == 0.0) {
        x2 = x1;
        dx = 1.0;
    }
 
    return true;
}
 
bool parseRange(const std::string& range, double& x1) {
    std::vector<std::string> tokens;
    String::tokenize(range, tokens, ":");
 
    if (!String::isNumeric(tokens)) return false;
 
    x1 = atof(tokens[0].c_str());
 
    return true;
}
 
// Parse the command line options to determine what the user wants to do
Options parseOptions(int argc, char** argv) {
    Options opts;
 
    // Print the help message and exit if desired
    if (argc < 2) printHelpMessage(argv[0]);
 
    if (optionExists(argc, argv, "-h") || optionExists(argc, argv, "--help"))
        printHelpMessage(argv[0]);
 
    // Get the temperature range
    if (optionExists(argc, argv, "-T")) {
        if (!parseRange(getOption(argc, argv, "-T"), opts.T1, opts.T2,
                        opts.dT)) {
            cout << "Bad format for temperature range!" << endl;
            printHelpMessage(argv[0]);
        }
    } else {
        cout << "Bad format for temperature range!" << endl;
        printHelpMessage(argv[0]);
    }
 
    // Get the pressure range
    if (optionExists(argc, argv, "-P")) {
        if (!parseRange(getOption(argc, argv, "-P"), opts.P1, opts.P2, opts.dP)) {
            cout << "Bad format for pressure range!" << endl;
            printHelpMessage(argv[0]);
        }
    } else {
        cout << "Bad format for pressure range!" << endl;
        printHelpMessage(argv[0]);
    }
 
    // Get the blowing rate range
    if (optionExists(argc, argv, "-b")) {
        if (!parseRange(getOption(argc, argv, "-b"), opts.Bg1, opts.Bg2, opts.dBg)) {
            cout << "Bad format for blowing rate range!" << endl;
            printHelpMessage(argv[0]);
        }
    } else {
        cout << "Bad format for blowing rate range!" << endl;
        printHelpMessage(argv[0]);
    }
 
    if (optionExists(argc, argv, "-m")) {
        opts.mixture = getOption(argc, argv, "-m");
    } else {
        printHelpMessage(argv[0]);
    }
 
    if (optionExists(argc, argv, "-bl")) {
        opts.boundary_layer_comp = getOption(argc, argv, "-bl");
    } else {
        printHelpMessage(argv[0]);
    }
 
    if (optionExists(argc, argv, "-py")) {
        opts.pyrolysis_composition = getOption(argc, argv, "-py");
        opts.pyrolysis_exist = true;
    }
 
    return opts;
}

 
int main(int argc, char* argv[]) {
#ifdef _GNU_SOURCE
    feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif
 
    Options opts = parseOptions(argc, argv);
 
    Mixture mix(opts.mixture);
 
    const int ne = mix.nElements();
    const int ns = mix.nSpecies();
 
    std::vector<double> Yke(ne, 0);
    std::vector<double> Ykg(ne, 0);
    std::vector<double> Xw(ns, 0);
 
    // Run conditions
    double T1 = opts.T1;
    double T2 = opts.T2;
    double dT = opts.dT;
    double P1 = opts.P1;
    double P2 = opts.P2;
    double dP = opts.dP;
    double Bg1 = opts.Bg1;
    double Bg2 = opts.Bg2;
    double dBg = opts.dBg;
    double Bc, hw;
 
    mix.getComposition(opts.boundary_layer_comp, Yke.data(), Composition::MASS);
 
    if (opts.pyrolysis_exist)
        mix.getComposition(opts.pyrolysis_composition, Ykg.data(),
                           Composition::MASS);
 
    // cout << setw(10) << T1 << setw(10) << T2 << setw(10) << dT << endl;
    // cout << setw(10) << P1 << setw(10) << P2 << setw(10) << dP << endl;
    // cout << setw(10) << Bg1 << setw(10) << Bg2 << setw(10) << dBg << endl;
    cout << setw(10) <<"\"P[bar]\"" << setw(15) << "\"B'g\"" << setw(15) << "\"Tw[K]\"" << setw(15) << "\"B'c\"" << setw(15)
         << "\"hw[MJ/kg]\"";
    for (int i = 0; i < ns; ++i)
        cout << setw(25) << "\"" + mix.speciesName(i) + "\"";
    cout << endl;
    
    for (double P = P1; P < P2 + 1.0e-6; P *= dP) {
        for (double Bg = Bg1; Bg < Bg2 + 1.0e-6; Bg += dBg) {
            for (double T = T1; T < T2 + 1.0e-6; T += dT) {
                // replace mix.surfaceMassBalance()
                mix.surfaceMassBalance(Yke.data(), Ykg.data(), T, P, abs(Bg), Bc, hw, Xw.data());
                cout << setw(10) << P*1.0e-5 << setw(15) << Bg << setw(15) << T << setw(15) << Bc << setw(15) << hw / 1.0e6;
                for (int i = 0; i < ns; ++i) cout << setw(25) << Xw[i];
                cout << endl;
            } // end temperature loop
        } // end blowing rate loop
    } // end pressure loop
}
