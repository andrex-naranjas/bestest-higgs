//ChromoDipoleTopQuark includes
#ifndef CHROMODIPOLETOPQUARK_H
#define CHROMODIPOLETOPQUARK_H

#include <string>
#include <vector>

class ChromoDipoleTopQuark{

public:
  ChromoDipoleTopQuark();
  virtual ~ChromoDipoleTopQuark();
  virtual double execute();

private:
  double pi_val = 3.1415926536;

  // Integral functions
  virtual double MagneticDipole();

};

//to talk to python
extern "C"{
  double loopint_execute(){
    ChromoDipoleTopQuark* m_loopint = new ChromoDipoleTopQuark();
    return m_loopint->execute();
  }
}


#endif //> !CHROMODIPOLETOPQUARK_H
