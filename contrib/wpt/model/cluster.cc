#include<iostream>
#include<vector>
#include <math.h>
#include<unordered_map>
#include "cluster.h"
// #include<eigen3/Eigen/Dense>
// using Eigen::MatrixXd;
// using namespace Eigen;
using namespace std;
void Expansion(vector<vector<double>>& matrix,int e)
{
    vector<vector<double>> temp=matrix;
    for(int i=0;i<matrix.size();i++)
    {
        for(int j=0;j<matrix[0].size();j++)
        {
            double sum=0;
            for(int k=0;k<matrix[0].size();k++)
            {
                sum+=(matrix[i][k]*matrix[k][j]);
            }
            temp[i][j]=sum;
        }
    }
    matrix=temp;
    // MatrixXd re(matrix.size(),matrix.size());
    // for(int i=0; i<matrix.size(); i++)
    //     re. row(i)=  VectorXd::Map(&matrix[i][0],matrix[i].size());
    // m<<
}
void print(vector<vector<double>> matrix)
{
    for(int i=0;i<matrix.size();i++)
    {
        for(int j=0;j<matrix[0].size();j++)
        {
            cout<<matrix[i][j]<<" ";
        }
        cout<<endl;
    }
}
void Inflation(vector<vector<double>>& matrix,double r)
{
    
    for(int j=0;j<matrix[0].size();j++)
    {
        double sum=0;
        for(int i=0;i<matrix.size();i++)
        {
            matrix[i][j]=pow(matrix[i][j],r);
            sum+=matrix[i][j];
        }
        for(int i=0;i<matrix.size();i++)
        {
            matrix[i][j]/=sum;
        }
    }
}
void MCL(vector<vector<double>>& matrix,int e,int r)
{
    Expansion(matrix,e);
    Inflation(matrix,r);
}
void Init(vector<vector<double>>& matrix)
{
    for(int j=0;j<matrix[0].size();j++)
    {
        double sum=0;
        for(int i=0;i<matrix.size();i++)
        {
            sum+=matrix[i][j];
        }
        for(int i=0;i<matrix.size();i++)
        {
            matrix[i][j]/=sum;
        }
    }
}
bool judge(vector<vector<double>>& matrix)
{
    unordered_map<int,int> h;
    bool f=true;
    for(int i=0;i<matrix.size();i++)
    {
        for(int j=0;j<matrix[0].size();j++)
        {
            if(matrix[i][j]>0)
            {
                double father=matrix[i][j];
                for(int k=j+1;k<matrix[0].size();k++)
                {
                    if(matrix[i][k]>0&&abs(matrix[i][k]-father)>0.0001)
                    {
                        f=false;
                    }
                }
                break;
            }
        }
    }
    return f;
}
void divide()
{
    for(int i=0;i<m.size();i++)
    {
        for(int j=0;j<m[0].size();j++)
        {
            if(m[i][j]>0)
            {
                for(int k=j+1;k<m[0].size();k++)
                {
                    if(m[i][k]>0)
                    {
                        clu[i].push_back(k);
                    }
                }
                num++;
                break;
            }
        }
    }
}
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