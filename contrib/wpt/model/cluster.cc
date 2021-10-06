#include <iostream>
#include <vector>
#include <math.h>
#include <unordered_map>
#include "cluster.h"
#include <sstream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <limits>
#include <cmath>
#include <stack>
// #include<eigen3/Eigen/Dense>
// using Eigen::MatrixXd;
// using namespace Eigen;
using namespace std;

//MCL
void
Expansion (vector<vector<double>> &matrix, int e)
{
  vector<vector<double>> temp = matrix;
  for (int i = 0; i < matrix.size (); i++)
    {
      for (int j = 0; j < matrix[0].size (); j++)
        {
          double sum = 0;
          for (int k = 0; k < matrix[0].size (); k++)
            {
              sum += (matrix[i][k] * matrix[k][j]);
            }
          temp[i][j] = sum;
        }
    }
  matrix = temp;
  // MatrixXd re(matrix.size(),matrix.size());
  // for(int i=0; i<matrix.size(); i++)
  //     re. row(i)=  VectorXd::Map(&matrix[i][0],matrix[i].size());
  // m<<
}
void
print (vector<vector<double>> matrix)
{
  for (int i = 0; i < matrix.size (); i++)
    {
      for (int j = 0; j < matrix[0].size (); j++)
        {
          cout << matrix[i][j] << " ";
        }
      cout << endl;
    }
}
void
Inflation (vector<vector<double>> &matrix, double r)
{

  for (int j = 0; j < matrix[0].size (); j++)
    {
      double sum = 0;
      for (int i = 0; i < matrix.size (); i++)
        {
          matrix[i][j] = pow (matrix[i][j], r);
          sum += matrix[i][j];
        }
      for (int i = 0; i < matrix.size (); i++)
        {
          matrix[i][j] /= sum;
        }
    }
}
void
MCL (vector<vector<double>> &matrix, int e, int r)
{
  // Expansion(matrix,e);
  // Expansion(matrix,e);
  Inflation (matrix, r);
}
void
Init (vector<vector<double>> &matrix)
{
  for (int j = 0; j < matrix[0].size (); j++)
    {
      double sum = 0;
      for (int i = 0; i < matrix.size (); i++)
        {
          sum += matrix[i][j];
        }
      for (int i = 0; i < matrix.size (); i++)
        {
          matrix[i][j] /= sum;
        }
    }
}
// bool judge(vector<vector<double>>& matrix)
// {
//     unordered_map<int,int> h;
//     bool f=true;
//     for(int i=0;i<matrix.size();i++)
//     {
//         for(int j=0;j<matrix[0].size();j++)
//         {
//             if(matrix[i][j]>0)
//             {
//                 double father=matrix[i][j];
//                 for(int k=j+1;k<matrix[0].size();k++)
//                 {
//                     if(matrix[i][k]>0&&abs(matrix[i][k]-father)>0.0001)
//                     {
//                         f=false;
//                     }
//                 }
//                 break;
//             }
//         }
//     }
//     return f;
// }
void
divide ()
{
  for (int i = 0; i < m.size (); i++)
    {
      for (int j = 0; j < m[0].size (); j++)
        {
          if (m[i][j] > 0)
            {
              for (int k = j + 1; k < m[0].size (); k++)
                {
                  if (m[i][k] > 0)
                    {
                      clu[i].push_back (k);
                    }
                }
              num++;
              break;
            }
        }
    }
}

//
// int main()
// {
//     // vector<vector<double>> m={{1,1,1,1},{1,1,0,1},{0.01,0,1,0},{1,1,0,1}};
//     // vector<vector<double>> m(10,vector<double>(10,rand() % (1000)));
//     for(int i=0;i<10;i++)
//     {
//         for(int j=0;j<10;j++)
//         {
//             m[i][j]=rand()%1000;
//         }
//     }
//     Init(m);
//     // print(m);
//     // for(int i=0;i<10;i++)
//     unordered_map<int,int> h;
//     for(int i=0;i<m.size();i++)
//     {
//         h[i]=i;
//     }
//     while(1)
//     {
//         MCL(m,2,10);
//         if(judge(m))
//         {
//             break;
//         }
//         // print(m);
//         // cout<<endl;
//     }
//     print(m);
// }

//DBSCAN

class point
{
public:
  float x;
  float y;
  int cluster = 0;
  int pointType = 1; //1 noise 2 border 3 core
  int pts = 0; //points in MinPts
  vector<int> corepts;
  // vector<int> borderpts;
  int visited = 0;
  point ()
  {
  }
  point (float a, float b, int c)
  {
    x = a;
    y = b;
    cluster = c;
  }
};
vector<vector<point *>> clus;
int range = 2;
float
stringToFloat (string i)
{
  stringstream sf;
  float score = 0;
  sf << i;
  sf >> score;
  return score;
}
vector<point>
openFile (const char *dataset)
{
  fstream file;
  file.open (dataset, ios::in);
  if (!file)
    {
      cout << "Open File Failed!" << endl;
      vector<point> a;
      return a;
    }
  vector<point> data;
  int i = 1;
  while (!file.eof ())
    {
      string temp;
      file >> temp;
      int split = temp.find (',', 0);
      point p (stringToFloat (temp.substr (0, split)),
               stringToFloat (temp.substr (split + 1, temp.length () - 1)), i++);
      vector<point *> t;
      clus.push_back (t);
      clus[i - 2].push_back (&p);
      data.push_back (p);
    }
  file.close ();
  cout << "successful!" << endl;
  return data;
}
float
squareDistance (point a, point b)
{
  return sqrt ((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

bool
judge (int cluster)
{
  for (int i = 0; i < clus[cluster].size (); i++)
    {
      for (int j = i + 1; j < clus[cluster].size (); j++)
        {
          point *a = clus[cluster][i];
          point *b = clus[cluster][j];
          double temp = (a->x - b->x) * (a->x - b->x) + (a->y - b->y) * (a->y - b->y);
          if (temp > range * range)
            {
              return false;
            }
        }
    }
  return true;
}
void
DBSCAN (vector<point> dataset, float Eps, int MinPts)
{
  int len = dataset.size ();
  //calculate pts
  cout << "calculate pts" << endl;
  for (int i = 0; i < len; i++)
    {
      for (int j = i + 1; j < len; j++)
        {
          if (squareDistance (dataset[i], dataset[j]) < Eps)
            {
              dataset[i].pts++;
              dataset[j].pts++;
            }
        }
    }
  //core point
  cout << "core point " << endl;
  vector<point> corePoint;
  for (int i = 0; i < len; i++)
    {
      if (dataset[i].pts >= MinPts)
        {
          dataset[i].pointType = 3;
          corePoint.push_back (dataset[i]);
        }
    }
  cout << "joint core point" << endl;
  //joint core point
  for (int i = 0; i < corePoint.size (); i++)
    {
      for (int j = i + 1; j < corePoint.size (); j++)
        {
          if (squareDistance (corePoint[i], corePoint[j]) < Eps)
            {
              corePoint[i].corepts.push_back (j);
              corePoint[j].corepts.push_back (i);
            }
        }
    }
  for (int i = 0; i < corePoint.size (); i++)
    {
      stack<point *> ps;
      if (corePoint[i].visited == 1)
        continue;
      ps.push (&corePoint[i]);
      point *v;
      while (!ps.empty ())
        {
          v = ps.top ();
          v->visited = 1;
          ps.pop ();
          for (int j = 0; j < v->corepts.size (); j++)
            {
              if (corePoint[v->corepts[j]].visited == 1)
                continue;
              int temp = corePoint[v->corepts[j]].cluster;
              corePoint[v->corepts[j]].cluster = corePoint[i].cluster;
              clus[corePoint[v->corepts[j]].cluster].push_back (&corePoint[i]);
              if (!judge (corePoint[i].cluster))
                {
                  corePoint[v->corepts[j]].cluster = temp;
                  clus[corePoint[v->corepts[j]].cluster].pop_back ();
                }
              corePoint[v->corepts[j]].visited = 1;
              ps.push (&corePoint[v->corepts[j]]);
            }
        }
    }
  cout << "border point,joint border point to core point" << endl;
  //border point,joint border point to core point
  for (int i = 0; i < len; i++)
    {
      if (dataset[i].pointType == 3)
        continue;
      for (int j = 0; j < corePoint.size (); j++)
        {
          if (squareDistance (dataset[i], corePoint[j]) < Eps)
            {
              dataset[i].pointType = 2;
              dataset[i].cluster = corePoint[j].cluster;
              break;
            }
        }
    }

  cout << "noise point,joint noise point to core point" << endl;
  //noise point,joint noise point to core point
  for (int i = 0; i < len; i++)
    {
      if (dataset[i].pointType == 3 || dataset[i].pointType == 2)
        continue;
      double minn = 10000000;
      for (int j = 0; j < corePoint.size (); j++)
        {
          if (squareDistance (dataset[i], corePoint[j]) < minn)
            {
              minn = squareDistance (dataset[i], corePoint[j]);
              dataset[i].pointType = 1;
              dataset[i].cluster = corePoint[j].cluster;
            }
        }
    }
  cout << "output" << endl;
  //output
  fstream clustering;
  clustering.open ("clustering.txt", ios::out);
  for (int i = 0; i < len; i++)
    {
      if (dataset[i].pointType == 2)
        clustering << dataset[i].x << "," << dataset[i].y << "," << dataset[i].cluster << "\n";
    }
  clustering << endl;
  for (int i = 0; i < corePoint.size (); i++)
    {
      clustering << corePoint[i].x << "," << corePoint[i].y << "," << corePoint[i].cluster << "\n";
    }

  clustering << endl;
  for (int i = 0; i < len; i++)
    {
      if (dataset[i].pointType == 1)
        clustering << dataset[i].x << "," << dataset[i].y << "," << dataset[i].cluster << "\n";
    }
  clustering.close ();
}
int
main (int argc, char **argv)
{
  vector<point> dataset = openFile ("dataset3.txt");
  DBSCAN (dataset, 1.5, 2);
  return 0;
}