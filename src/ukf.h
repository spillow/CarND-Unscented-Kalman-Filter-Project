#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_ = false;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_ = true;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_ = true;

  ///* time when the state is true, in us
  long long time_us_ = 0;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  const double std_a_ = 3.0;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  const double std_yawdd_ = 0.5;

  ///* Laser measurement noise standard deviation position1 in m
  const double std_laspx_ = 0.15;

  ///* Laser measurement noise standard deviation position2 in m
  const double std_laspy_ = 0.15;

  ///* Radar measurement noise standard deviation radius in m
  const double std_radr_ = 0.3;

  ///* Radar measurement noise standard deviation angle in rad
  const double std_radphi_ = 0.03;

  ///* Radar measurement noise standard deviation radius change in m/s
  const double std_radrd_ = 0.3;

  ///* State dimension
  const int n_x_ = 5;

  ///* Augmented state dimension
  const int n_aug_ = 7;

  ///* Sigma point spreading parameter
  const double lambda_ = (double)(3 - n_aug_);

  ///* Weights of sigma points
  VectorXd weights_ = VectorXd(2 * n_aug_ + 1);

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_ = VectorXd(5);

  ///* state covariance matrix
  MatrixXd P_ = MatrixXd(5, 5);

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);

  void InitFirstMeasurement(const MeasurementPackage &meas_package);
};

#endif /* UKF_H */
