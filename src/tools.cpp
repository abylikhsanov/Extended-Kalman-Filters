#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd toReturn(4);
    toReturn << 0,0,0,0;

  if(estimations.size() != ground_truth.size()){
      return toReturn;
  }

  for(unsigned int i = 0; i<estimations.size();++i){
      VectorXd temp = estimations[i] - ground_truth[i];
      temp = temp.array()*temp.array();
      toReturn += temp;
  }
    toReturn = toReturn/estimations.size();
    toReturn = toReturn.array().sqrt();
    std::cout<<toReturn<<std::endl;
    return toReturn;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

    double a = pow(x_state[0],2)+pow(x_state[1],2);
    double b = (x_state[1]*(x_state[2]*x_state[1]-x_state[3]*x_state[0]))/(pow(a,3/2));
    double c = (x_state[0]*(x_state[3]*x_state[0]-x_state[2]*x_state[1]))/(pow(a,3/2));

    MatrixXd toReturn(3,4);

    if(fabs(a) < 0.0001){
        cout<<"Illegal: Divison by zero"<<endl;
        a = 0.0001;
    }

    toReturn << x_state[0]/sqrt(a), x_state[1]/sqrt(a), 0, 0,
            -(x_state[1]/(a)), x_state[0]/(a), 0, 0,
            b, c, x_state[0]/sqrt(a), x_state[1]/sqrt(a);

    return toReturn;

}
