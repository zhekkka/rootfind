//
// Created by zhest on 19.10.2022.
//
#ifndef SLAE_FUNCTIONS_H
#define SLAE_FUNCTIONS_H
#include<vector>
std::vector<double> operator*(double a, std::vector<double>& vector);
std::vector<double> operator+(std::vector<double>&a,std::vector<double>&b);
std::vector<double> operator-(std::vector<double>&a,std::vector<double>&b);
double Norm(std::vector<double>&a);
double func(double &x);
void contEnum(double &leftBorder,double &rightBorder,double &eps);
double hybridNewtonMethod(double &leftBorder, double &rightBorder);
std::vector<double> newtoneSystemMethod(std::vector<double> &vector,int &k);
std::vector<double> phi();
std::vector<double> systemFunc(std::vector<double>& vector);
#endif //SLAE_FUNCTIONS_H
