#ifndef CLUSTER_H
#define CLUSTER_H


#include<iostream>
#include<vector>
using namespace std;
struct anchor
{
    double x,y,prio,tmin;
    anchor(double x,double y)
    {
        this->x=x;
        this->y=y;
    }
};
int num=0;
vector<vector<int>> clu;
vector<vector<double>> m;
void Expansion(vector<vector<double>>& matrix,int e);
void print(vector<vector<double>> matrix);
void Inflation(vector<vector<double>>& matrix,double r);
void MCL(vector<vector<double>>& matrix,int e,int r);
void Init(vector<vector<double>>& matrix);
bool judge(vector<vector<double>>& matrix);
void divide();



#endif /* CLUSTER_H */


// #ifdef __cplusplus
// extern "C"
// {
// #endif
 
// #include<iostream>
// #include<vector>
// using namespace std;
// struct anchor
// {
//     double x,y,prio;
//     anchor(double x,double y)
//     {
//         this->x=x;
//         this->y=y;
//     }
// };
// int num=0;
// vector<vector<int>> clu;
// vector<vector<double>> m;
// void Expansion(vector<vector<double>>& matrix,int e);
// void print(vector<vector<double>> matrix);
// void Inflation(vector<vector<double>>& matrix,double r);
// void MCL(vector<vector<double>>& matrix,int e,int r);
// void Init(vector<vector<double>>& matrix);
// bool judge(vector<vector<double>>& matrix);
// void divide();
 
// #ifdef __cplusplus
// }
#endif