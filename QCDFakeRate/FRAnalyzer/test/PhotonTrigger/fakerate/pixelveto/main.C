#include "myPlot.C"
#include "TROOT.h"

int main(){

  gROOT->ProcessLine("#include <vector>");
  gROOT->ProcessLine("#include <map>");
  gSystem->Load("libDCache.so");
  myPlot m;
  m.Loop();

}
