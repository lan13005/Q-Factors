#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooClassFactory.h"
#include "TROOT.h"
 
using namespace RooFit;

void makeBivariateGaus()
{
    ////////////////////
    ////////////////////
    // For example if we want to fit the 2D distribution M(pi0) vs M(eta)
    // x = M(pi0)
    // y = M(eta)
    // sx = width of gaussian in M(pi0)
    // sy = width of gaussian in M(eta)
    // px = peak of gaussian in M(pi0)
    // py = peak of gaussian in M(eta)
    // r = correlation
    ////////////////////
    ////////////////////
    string coeff="1/(2*3.14157*sx*sy*sqrt(1-pow(r,2)))";
    string argCoeff="-1/(2*(1-pow(r,2)))";
    string xterm="((x-px)/sx)";
    string yterm="((y-py)/py)";
    string crossterm="-2*r*"+xterm+"*"+yterm;
    string expression=coeff+"*exp("+argCoeff+"*(pow("+xterm+",2)"+crossterm+"+pow("+yterm+",2)))";
    cout << "Signal expression: " << expression << endl;
    RooClassFactory::makePdf("bivariateGaus", "x,y,px,sx,py,sy,r", "", expression.c_str());
}
