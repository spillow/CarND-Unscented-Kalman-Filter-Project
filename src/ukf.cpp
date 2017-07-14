#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
    /**
    TODO:

    Hint: one or more values initialized above might be wildly off...
    */

    // set weights
    double weight_0 = lambda_ / (lambda_ + n_aug_);
    weights_(0) = weight_0;
    for (int i = 1; i < 2 * n_aug_ + 1; i++)
    {
        double weight = 0.5 / (n_aug_ + lambda_);
        weights_(i) = weight;
    }
}

UKF::~UKF() {}

void UKF::InitFirstMeasurement(const MeasurementPackage &meas_package)
{
    VectorXd measures = meas_package.raw_measurements_;

    P_ <<
        1, 0, 0, 0, 0,       // px
        0, 1, 0, 0, 0,       // py
        0, 0, 15 * 15, 0, 0, // v   (assume bikes mostly stay under 15 m/s)
        0, 0, 0, 6 * 6, 0,   // yaw (normalized angle |yaw| < 2pi ~ 6)
        0, 0, 0, 0, 2 * 2;   // yaw rate (2pi rad/3 s ~ 2 max)?
          
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
        float rho = measures(0);
        float theta = measures(1);
        x_ << rho*cos(theta), rho*sin(theta), 0, 0, 0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
        x_ << measures(0), measures(1), 0, 0, 0;
    }
    else
    {
        assert(0 && "unknown measurement type!");
    }
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{
    long long currTime = meas_package.timestamp_;
    long long prevTime = time_us_;

    time_us_ = currTime;

    if (!is_initialized_)
    {
        InitFirstMeasurement(meas_package);

        is_initialized_ = true;
        return;
    }

    const double dt = (double)(currTime - prevTime) / 1000000.0;

    Prediction(dt);

    switch (meas_package.sensor_type_)
    {
    case MeasurementPackage::RADAR:
        UpdateRadar(meas_package);
        break;
    case MeasurementPackage::LASER:
        UpdateLidar(meas_package);
        break;
    default:
        assert(0 && "unknown measurement type!");
    }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}
