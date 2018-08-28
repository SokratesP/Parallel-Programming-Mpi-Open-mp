#include <fstream>
#include <iostream>

using namespace std;

int main(){
  cout << "i = " << endl;
  int i;
  int j;
  ifstream ifs;
  ifs.open ("waterfall_grey_1920_2520.raw", std::ifstream::in);
  unsigned int matrix[1920][2550];
  char c = ifs.get();

  while (ifs.good()) {
    //cout << (unsigned int)c << " ";;
    for (i = 0; i < 1920; ++i)
    {
      for (j = 0; i < 2550; ++j)
      {
        cout << "i = " << i << endl;
        matrix[i][j]= (unsigned int)c;
        c = ifs.get();
      }
    }
    //c = ifs.get();
   // i++;
  }
    cout << endl;
    cout << "i = " << i << endl;
    ifs.close();
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