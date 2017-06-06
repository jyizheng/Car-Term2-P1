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

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  H_laser << 1,0,0,0,
		0,1,0,0; 
  noise_ax = 9;
  noise_ay = 9;
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
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double rho = measurement_pack.raw_measurements[0];
      double phi = measurement_pack.raw_measurements[1];
      double dot_rho = measurement_pack.raw_measurements[2];

      ekf_.x_ << rho*cos(phi), rho*sin(phi), dot_rho*cos(phi), dot_rho*sin(phi);  
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_ << measurement_pack.raw_measurements[0], measurement_pack.raw_measurements[1], 0, 0; 
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /* Prepare the space for matrix in ekf_ */
  VectorXd x(4);   
  MatrixXd P(4,4);   
  MatrixXd F(4,4);
  MatrixXd Q(4,4);
  MatrixXd H;
  MatrixXd R;

  P << 1,0,0,0,
	0,1,0,0,
	0,0,1000,0,
	0,0,0,1000; 

  F << 1,0,1,0,
   	0,1,0,1,
    	0,0,1,0,
     	0,0,0,1;
  

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    //Calculate Jacobian Matrix
    Hj_ = CalculateJacobian(ekf_.x_);  
    H = Hj_;
    R = R_radar;
  } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    H = H_laser;
    R = R_laser;
  }

  efk_.Init(x, P, F, H, R, Q);  

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
  float dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
  previous_timestamp = measurement_timestamp;

  float dt_2 = dt*dt;
  float dt_3 = dt_2*dt;
  float dt_4 = dt_3*dt;  

  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;
  ekf_.Q_ = Matrix(4,4);
  ekf_.Q_ << dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
             	0, dt_4/4*noise_ay, 0, dt_3/2*noisy_ay,
             	dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
             	0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

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
    // Radar updates
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
