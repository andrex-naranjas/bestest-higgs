//ChromoDipole class includes
#ifndef CHROMODIPOLETOPQUARK_CXX
#define CHROMODIPOLETOPQUARK_CXX

#include "ChromoDipoleTopQuark.h"

#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>
//#include <complex>
ChromoDipoleTopQuark::ChromoDipoleTopQuark()
{
}

ChromoDipoleTopQuark::~ChromoDipoleTopQuark(){}

double ChromoDipoleTopQuark::execute(){
  std::cout<<"testing.."<<std::endl;

  double test = MagneticDipole();
  std::cout<<test<<std::endl;


  return 0.;
}



double ChromoDipoleTopQuark::MagneticDipole(){

  double P1 = -0.238605;
  double P2 = 0.238605;
  double S1 = 0.238605;
  double S2 = 0.238605;

  double qq = 91*91;
  double mfi = 173;
  double ms = 125;

  double pi_val = 3.14192;

  double C0 = 1;

  double mu = 1.;

  double mfj = 283.35;

  double value1 = (mfj*(P1*S2+P2*S1))/(16*pi_val*pi_val*mfi*qq*(4*mfi*mfi-qq));
  double value2 = 2.*std::pow(mfi, 2)*qq*(std::pow(mfi, 2) - std::pow(mfj, 2) + std::pow(ms, 2))*C0;
  double value3 = qq*(std::pow(mfi, 2) + std::pow(mfj, 2) - std::pow(ms, 2) )*std::log(std::pow(mfj, 2)/std::pow(ms, 2));
  double value4 = 2.*std::pow(mfi, 2)*qq*std::log(std::pow(mu, 2)/std::pow(mfj, 2));
  double value5 = 2.*std::pow(mfi, 2)*std::pow(qq*(qq - 4*std::pow(mfj, 2)), 0.5) * std::log((std::pow(qq*(qq - 4*std::pow(mfj, 2)), 0.5) + 2*std::pow(mfj, 2) - qq) / (2.* std::pow(mfj, 2)));
  double value6 = 2.*std::pow(mfi, 2)*qq*std::log(std::pow(mu, 2)/std::pow(ms, 2));

  double value7p = 2*qq* std::pow( std::pow(mfi, 4) - 2*std::pow(mfi, 2)*( std::pow(mfj, 2) + std::pow(ms, 2)) + std::pow( std::pow(mfj, 2) - std::pow(ms, 2), 2),   0.5);

  double value7pp = std::log(((-std::pow(mfi, 2) + std::pow(std::pow(mfi, 4) - 2*std::pow(mfi, 2)* (std::pow(mfj, 2) + std::pow(ms, 2)) + std::pow(std::pow(mfj, 2) - std::pow(ms, 2), 2), 0.5) + std::pow(mfj, 2) + std::pow(ms, 2))/(2*mfj*ms)));
  double value7 = value7p * value7pp;


  return 1.;
}

#endif
