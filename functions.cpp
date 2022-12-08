#include "iostream"
#include "cmath"
#include <vector>

std::vector<double> operator*(double a, std::vector<double>& vector){
    std::vector<double> c;
    c.resize(vector.size());
    for (int i = 0; i < vector.size(); ++i) {
        c[i]=vector[i]*a;
    }
    return c;
}
std::vector<double> operator+(std::vector<double>&a,std::vector<double>&b){
    std::vector<double> c;
    c.resize(a.size());
    for (int i = 0; i < a.size(); ++i) {
        c[i]=a[i]+b[i];
    }
    return c;
}
std::vector<double> operator-(std::vector<double>&a,std::vector<double>&b){
    std::vector<double> c;
    c.resize(a.size());
    for (int i = 0; i < a.size(); ++i) {
        c[i]=a[i]-b[i];
    }
    return c;
}
double Norm(std::vector<double>&a){
    double s=0;
    for (int i = 0; i < a.size(); ++i) {
        s+=(a[i]*a[i]);
    }
    return std::sqrt(s);
}
double func(double &x){
    return (log10(x)-(7/(2*x+6)));
}
double derivedFunc(double &x){
    return ((7/(2*((x+3)*(x+3))))+(1/x));
}
void newtoneIteration(double &xk,double&xK){
    xK=xk-(func(xk)/ (derivedFunc(xk)));
}
void halfDivision(double &xk,double& xK,double &leftBorder, double& rightBorder){
    if (func(xk)* func(rightBorder)>0){
        rightBorder=xk;
    }
    if (func(xk)* func(leftBorder)>0){
        leftBorder=xk;
    }
    xK=(rightBorder+leftBorder)/2;
}
double hybridNewtonMethod(double &leftBorder, double &rightBorder){
    double xk,xK;
    xK=(leftBorder+rightBorder)/2;
    do{
        xk=xK;
        newtoneIteration(xk,xK);
        if (xK>rightBorder or xK<leftBorder){
            halfDivision(xk,xK,leftBorder,rightBorder);
        }
    }
    while(std::abs(xk-xK)>=(0.0001));
    return xK;
}
void contEnum(double &leftBorder,double &rightBorder,double &eps){
    double a=leftBorder;
    double b=leftBorder+eps;
    while(func(a)*func(b)>=0 and (b<=rightBorder)){
        a+=0.05;
        b+=0.05;
    }
    leftBorder=a;
    rightBorder=b;
}
std::vector<double> systemFunc(std::vector<double>& vector){
    std::vector<double> x;
    x.resize(2);
    x[0]=sin(vector[0])+2*vector[1]-2;
    x[1]=vector[0]+cos(vector[1]-1)-0.7;
    return x;
}
std::vector<double> iteration(std::vector<double> &vector){
    std::vector<double> *matrix=new std::vector<double>[2];
    matrix[0].resize(2);
    matrix[1].resize(2);
    std::vector<double> delta;
    delta.resize(2);
    matrix[0][0]= cos(vector[0]);
    matrix[0][1]=2;
    matrix[1][0]=0;
    matrix[1][1]=(-sin(vector[1]-1)-(2/(cos(vector[0]))));
    delta[1]=((-systemFunc(vector)[1])/matrix[1][1]);
    delta[0]=(((-systemFunc(vector)[0])-(matrix[0][1]*delta[1]))/(matrix[0][0]));
    delta[0]+=vector[0];
    delta[1]+=vector[1];
    delete[] matrix;
    return delta;
}
std::vector<double> newtoneSystemMethod(std::vector<double> &vector,int &k){
    std::vector<double> xk,xK,difference,delta;
    xK=vector;
    do{
        ++k;
        xk=xK;
        xK= iteration(xk);
        difference=xK-xk;
    } while(Norm(difference)>0.0001);
    return xK;
}
std::vector<double> systemPhiFunc(double &lamda,std::vector<double>& vector){
    std::vector<double> x;
    x.resize(2);
    x[0]=(lamda*sin(vector[0]))+2*vector[1]-2;
    x[1]=vector[0]+(lamda*cos(vector[1]-1))-0.7;
    return x;
}
std::vector<double>phiFunc(double &lamda,std::vector<double>&vector){
    std::vector<double> *matrix=new std::vector<double>[2];
    matrix[0].resize(2);
    matrix[1].resize(2);
    std::vector<double> delta;
    delta.resize(2);
    matrix[0][0]= (lamda*cos(vector[0]));
    matrix[0][1]=2;
    matrix[1][0]=0;
    matrix[1][1]=(((-lamda)*sin(vector[1]-1))-(2/((lamda)*cos(vector[0]))));
    delta[1]=((-systemPhiFunc(lamda,vector)[1])/matrix[1][1]);
    delta[0]=(((-systemPhiFunc(lamda,vector)[0])-(matrix[0][1]*delta[1]))/(matrix[0][0]));
    delta[0]+=vector[0];
    delta[1]+=vector[1];
    delete[] matrix;
    return delta;
}
std::vector<double> phi(){
    std::vector<double> x;
    x.resize(2);
    int n=10;
    x[1]=1;
    x[0]=0.7;
    for (int i = 1; i < n+1; ++i) {
        double coef=(((double)i)/((double)n));
        x= phiFunc(coef,x);
    }
    return x;
}