
#include "kalman_filter.h"
#include<iostream>
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;

}

void KalmanFilter::Predict() {
    x_ = F_*x_;
    P_ = F_*P_*(F_.transpose())+Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    std::cout<<H_.size()<<std::endl;
    VectorXd y = z - H_*x_;
    MatrixXd S = H_*P_*(H_.transpose()) + R_;
    MatrixXd K = P_*(H_.transpose())*(S.inverse());
    x_ = x_ + K*y;
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size,x_size);
    P_ = (I-K*H_)*P_;

}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    VectorXd x_temp(3);
    x_temp << sqrt(x_(0)*x_(0)+x_(1)*x_(1)),
              atan2(x_(1),x_(0)),
            (x_(0)*x_(2)+x_(1)*x_(3))/(sqrt(pow(x_(0),2)+pow(x_(1),2)));
    if(x_temp(0) == 0){
        x_temp(0) = 0.0001;
    }
    VectorXd y = z - x_temp;
    y(1) = atan2(sin(y(1)), cos(y(1)));
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_*P_*Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_*Ht;
    MatrixXd K = PHt*Si;

    x_ = x_ + K*y;
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size,x_size);
    std::cout<<H_.size()<<std::endl;
    P_ = (I-K*H_)*P_;
    std::cout<<"P is done"<<std::endl;



  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
}
