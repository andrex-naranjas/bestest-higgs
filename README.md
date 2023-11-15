# Bestest little higgs

Code to compute loop integrals with the results of FeynCalc

## Framework installation

To install the framework you need anaconda and git on a linux machine. In a terminal type:
1. Clone the repository:
  ```
  git clone git@github.com:andrex-naranjas/bestest-higgs.git
  ```
2. Access the code:
  ```
  cd bestest-higgs
  ```
3. Install the conda enviroment:
  ```
  conda env create -f config.yml
  conda activate higgs
  conda develop .
  ```
3.1 Update the conda enviroment:
   ```
   conda env update --file config.yml --prune
   ```
4. Compile the LoopTools package (here we use gfortran):
  ```
  cd ./feyncalc/LoopTools
  python3 compile_loop_tools.py
  cd ../..
  ```