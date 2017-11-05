#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;


  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);
  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  H_laser_ << 1,0,0,0,
             0,1,0,0;


    ekf_.R_ = MatrixXd(2, 2);
    ekf_.R_ << 0.0225, 0,
            0, 0.0225;

    //measurement matrix
    ekf_.H_ = MatrixXd(2, 4);
    ekf_.H_ << 1, 0, 0, 0,
            0, 1, 0, 0;

    //the initial transition matrix F_
    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ << 1, 0, 1, 0,
            0, 1, 0, 1,
            0, 0, 1, 0,
            0, 0, 0, 1;

    //set the acceleration noise components
    noise_ax = 9.0;
    noise_ay = 9.0;

    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;//dummy entry

    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1000, 0,
            0, 0, 0, 1000;

    ekf_.Q_ = MatrixXd(4,4);
    ekf_.Q_ << 0,0,0,0,
               0,0,0,0,
               0,0,0,0,
               0,0,0,0;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */


}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {

    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {


        float px = measurement_pack.raw_measurements_[0]*cos(measurement_pack.raw_measurements_[1]);
        float py = measurement_pack.raw_measurements_[0]*sin(measurement_pack.raw_measurements_[1]);
        ekf_.x_ << px,py,0,0;

      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
        cout<<"Problem is here "<<endl;
        ekf_.x_ << measurement_pack.raw_measurements_(0),measurement_pack.raw_measurements_(1),0,0;

    }
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  double delta_t = (measurement_pack.timestamp_-previous_timestamp_ )/1000000.0;

  ekf_.F_ << 1,0,delta_t,0,
             0,1,0,delta_t,
             0,0,1,0,
             0,0,0,1;

  ekf_.Q_ << (pow(delta_t,4)/4)*pow(noise_ax,2),0,(pow(delta_t,3)/2)*pow(noise_ax,2),0,
             0,(pow(delta_t,4)/4)*pow(noise_ay,2),0,(pow(delta_t,3)/2)*pow(noise_ay,2),
             (pow(delta_t,3)/2)*pow(noise_ax,2),0,pow(delta_t,2)*pow(noise_ax,2),0,
             0,(pow(delta_t,3)/2)*pow(noise_ay,2),0,pow(delta_t,2)*pow(noise_ay,2);

   previous_timestamp_ = measurement_pack.timestamp_;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      cout<<"recieved Radar"<<endl;
      Hj_ = tools.CalculateJacobian(ekf_.x_);
      ekf_.H_ = Hj_;
      ekf_.R_ = R_radar_;
      cout<<"Calling UpdateEKF"<<endl;
      ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
      ekf_.H_ = H_laser_;
      ekf_.R_ = R_laser_;
      ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
