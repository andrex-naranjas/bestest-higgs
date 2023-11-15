//ChromoDipoleTopQuark includes
#ifndef CHROMODIPOLETOPQUARK_H
#define CHROMODIPOLETOPQUARK_H

#include <string>
#include <vector>

class ChromoDipoleTopQuark{

public:
  ChromoDipoleTopQuark();
  virtual ~ChromoDipoleTopQuark();
  virtual double execute(double MZ);

private:
  double pi_val = 3.1415926536;

  // Integral functions
  virtual double MagneticDipole();

};


//to talk to python
extern "C"{
  double loopint_execute(double MZ){
    ChromoDipoleTopQuark* m_loopint = new ChromoDipoleTopQuark();
    return m_loopint->execute(MZ);
  }
}

#endif //> !CHROMODIPOLETOPQUARK_H
