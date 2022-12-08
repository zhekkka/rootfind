#include <iostream>
#include "functions.h"
#include <iomanip>
int main() {
    int cond;
    std::cout<<"0-equation 1-system"<<std::endl;
    std::cin>>cond;
    switch (cond) {
        case 0:{
            double leftBorder=1;
            double rightBorder=50;
            double eps=1;
            //Для того, чтобы точка вылетела за промежуток нужно закомментить следующую строчку(14-ю)
            //contEnum(leftBorder,rightBorder,eps);
            std::cout<<leftBorder<<" "<<rightBorder<<std::endl;
            double x=hybridNewtonMethod(leftBorder,rightBorder);
            std::cout<<func(x)<<std::endl;
            std::cout<<std::setprecision(10)<<x;
            break;
        }
        case 1:{
            std::vector<double> x;
            x.resize(2);
            x[0]=-0.2;
            x[1]=1.2;
            int k=0;
            std::cout<<"---Graphic method---"<<std::endl;
            std::cout<<x[0]<<' '<<x[1]<<std::endl;
            x=newtoneSystemMethod(x,k);
            std::cout<<"Solve: ";
            std::cout<<x[0]<<" "<<x[1]<<" "<<k<<std::endl;
            std::vector<double> y= systemFunc(x);
            std::cout<<y[0]<<" "<<y[1]<<std::endl;
            std::cout<<"---Next method---"<<std::endl;
            x=phi();
            std::cout<<x[0]<<" "<<x[1]<<std::endl;
            std::cout<<"Solve: "<<std::endl;
            k=0;
            x=newtoneSystemMethod(x,k);
            std::cout<<x[0]<<" "<<x[1]<<" "<<k;
            break;
        }
        default:{
            std::cout<<"Wrong option";
        }
    }

    return 0;
}
