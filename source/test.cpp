#include <fstream>
#include <iostream>
#include <cstdlib>
#include <ctime>
using namespace std;

int main(){
  //cout << "i = " << endl;
  //fflush(stdout);
  srand(time(nullptr));
  int i;
  int j;
  int lines=7;
  int columns=7;
  //ifstream ifs;
  //ifs.open ("waterfall_grey_1920_2520.raw", std::ifstream::in);
  int matrix[lines][columns];
  char c = ifs.get();

  //while (ifs.good()) {
    //cout << (unsigned int)c << " ";;
    for (i = 1; i < lines-1; i++)
    {
      for (j = 1; j < columns-1; j++)
      {
        matrix[i][j]= rand() % 255;
       // c = ifs.get();
      }
      cout << "i = " << i << endl;
    }
    //break;
    //c = ifs.get();
   // i++;
  //}

  for (i = 0; i < lines; i++)
  {
    for (j = 0; j < columns; j++)
    {
      cout<< matrix[i][j] << " " ;
      
    }
    cout << endl;
  }
  cout << endl;

  
//arxh loop
//stelnoyme grammes sthles klp...




//lamvabnoyme



//eswterikh epe3ergasia 2-> (lines-2) kai (2-> columns-2)




//wait() gia receives
//epe3ergasia




//wait() gia sents

//loop kai antistrofh pinakwn



  return 0;
}

/*
gia na gracw se raw arxeio
{
   MyData myData = GetData();
   ofstream binaryFile ("file.raw", ios::out | ios::binary);
   binaryFile.write ((char*)&myData, sizeof (MyData));
   binaryFile.close();
}
*/