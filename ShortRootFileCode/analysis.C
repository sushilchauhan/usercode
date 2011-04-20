#include "myCuts.C" 
#include <ctime>
int main(){ 
  time_t start, end ; 
  double timeConsumed ; 
  time(&start) ; 
  myCuts *a = new myCuts();
  a->Loop();
  a->myCuts::~myCuts() ; 
  time(&end) ; 
  timeConsumed = difftime (end,start) ; 
  cout << " time for processing job : " << timeConsumed << " seconds " << endl ; 
  return 0; 

}

