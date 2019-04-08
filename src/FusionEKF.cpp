#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
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
  ekf_.Q_ = Eigen::MatrixXd(4,4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */

   H_laser_ << 1.0, 0.0, 0.0, 0.0,
               0.0, 1.0, 0.0, 0.0;

   ekf_.P_ = MatrixXd(4, 4);
   ekf_.P_ << 1.0, 0.0, 0.0, 0.0,
              0.0, 1.0, 0.0, 0.0,
              0.0, 0.0, 1000.0, 0.0,
              0.0, 0.0, 0.0, 1000.0;

   ekf_.F_ = MatrixXd(4, 4);
   ekf_.F_ << 1.0, 0.0, 0.0, 0.0,
              0.0, 1.0, 0.0, 0.0,
              0.0, 0.0, 1.0, 0.0,
              0.0, 0.0, 0.0, 1.0; 

   ekf_.H_ = MatrixXd(4, 4);
   ekf_.H_ << 1.0, 0.0, 0.0, 0.0,
              0.0, 1.0, 0.0, 0.0,
              0.0, 0.0, 1.0, 0.0,
              0.0, 0.0, 0.0, 1.0;

}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    previous_timestamp_ = measurement_pack.timestamp_;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates 
      //         and initialize state.
      cout << "Init Radar Start: " << endl;
      const double rho =  measurement_pack.raw_measurements_[0];
      const double phi =  measurement_pack.raw_measurements_[1];

      double px = rho * cos(phi);
      double py = rho * sin(phi);
      
      if(fabs(px) < 0.0001){
         px = 0.1;
      }

      if(fabs(py) < 0.001){
        py = 0.01;
      }

      ekf_.x_ << px, py, 0, 0;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state.
      cout << "Init Laser Start: " << endl;
      double px = measurement_pack.raw_measurements_[0];
      double py = measurement_pack.raw_measurements_[1];

      if(fabs(px) < 0.0001){
         px = 0.1;
      }

      if(fabs(py) < 0.001){
        py = 0.01;
      }

      ekf_.x_ << px, py, 0, 0;
      cout << "Init Laser End " << endl;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    cout << "ekf_.x_: " << ekf_.x_ << endl;
    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
   cout << "Start Prediction" << endl;
   double delta = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
   previous_timestamp_ = measurement_pack.timestamp_;

   double delta_2 = delta * delta;
   double delta_3 = delta_2 * delta;
   double delta_4 = delta_3 * delta;

   ekf_.F_(0, 2) = delta;
   ekf_.F_(1, 3) = delta;

   cout << ekf_.F_ << endl;

   const double noise_ax = 9;
   const double noise_ay = 9;

   ekf_.Q_ <<   delta_4/4.0 * noise_ax, 0, delta_3/2.0 * noise_ax, 0,
                0, delta_4/4.0 * noise_ay, 0, delta_3/2.0 * noise_ay,
                delta_3/2.0 * noise_ax, 0, delta_2 * noise_ax, 0,
                0, delta_3/2.0 * noise_ay, 0, delta_2 * noise_ay;

   cout << "ekf Prediction run" << endl;
   ekf_.Predict();
   cout << "ekf Prediction end" << endl; 
  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Radar updates
      
      cout << "Radar Measurement update start" << endl;
      ekf_.hx_ = VectorXd(3);
      double px = ekf_.x_[0];
      double py = ekf_.x_[1];
      double vx = ekf_.x_[2];
      double vy = ekf_.x_[3];

      double rho;
      double phi;
      double drho;

      if(fabs(px) < 0.0001 or fabs(py) < 0.0001){
        
        if(fabs(px) < 0.0001){
          px = 0.0001;
        }

        if(fabs(py) < 0.0001){
          py = 0.0001;
        }
        
        rho = sqrt(px*px + py*py);
        phi = 0.0;
        drho = 0.0;
  
      }
      else{
        rho = sqrt(px*px + py*py);
        phi = atan2(py, px);
        drho = (px*vx + py*vy) / rho;
      }

      ekf_.hx_ << rho, phi, drho;
      Hj_ = tools.CalculateJacobian(ekf_.x_);

      if (Hj_.isZero(0)){
        return;
      }

      ekf_.H_ = Hj_;
      ekf_.R_ = R_radar_;

      cout << "Radar Measurement updateEKF start" << endl;
      ekf_.UpdateEKF(measurement_pack.raw_measurements_);
      cout << "Radar Measurement updateEKF end" << endl;
  } else {
    // TODO: Laser updates
    cout << "Laser Measurement update start" << endl;
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
    cout << "Laser Measurement update end" << endl;    
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
