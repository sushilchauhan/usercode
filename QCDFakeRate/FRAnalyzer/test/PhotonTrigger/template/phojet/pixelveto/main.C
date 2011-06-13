#include "myPlot.C"
#include "TROOT.h"


int main(){
  gROOT->ProcessLine("#include <vector>");
  gROOT->ProcessLine("#include <map>");
  myPlot m;
  m.Loop();
}
