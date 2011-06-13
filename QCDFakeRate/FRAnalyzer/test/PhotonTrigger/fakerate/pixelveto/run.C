void run(){
gROOT->ProcessLine("#include <vector>");
gROOT->ProcessLine("#include <map>");
gSystem->Load("libDCache.so");         
gSystem->Load("libmyPlot.so");
myPlot t;
t.Loop();
}
