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

static MatrixXd AugmentedSigmaPoints(
    int n_x, int n_aug, double lambda,
    double std_a, double std_yawdd, VectorXd x, MatrixXd P)
{
    //create augmented mean vector
    VectorXd x_aug = VectorXd(n_aug);

    //create augmented state covariance
    MatrixXd P_aug = MatrixXd(n_aug, n_aug);

    //create sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);

    //create augmented mean state
    x_aug.head(n_x) = x;
    x_aug(n_x+0) = 0.0;
    x_aug(n_x+1) = 0.0;

    //create augmented covariance matrix
    P_aug.fill(0.0);
    P_aug.topLeftCorner(n_x, n_x) = P;
    P_aug(n_x, n_x) = std_a*std_a;
    P_aug(n_x+1, n_x+1) = std_yawdd*std_yawdd;

    //create square root matrix
    MatrixXd L = P_aug.llt().matrixL();

    //create augmented sigma points
    Xsig_aug.col(0) = x_aug;
    for (int i = 0; i < n_aug; i++)
    {
        Xsig_aug.col(i + 1) = x_aug + sqrt(lambda + n_aug) * L.col(i);
        Xsig_aug.col(i + 1 + n_aug) = x_aug - sqrt(lambda + n_aug) * L.col(i);
    }

    return Xsig_aug;
}

static MatrixXd SigmaPointPrediction(
    int n_x, int n_aug, const MatrixXd &Xsig_aug, double delta_t)
{
    //create matrix with predicted sigma points as columns
    MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);

    //predict sigma points
    for (int i = 0; i < 2 * n_aug + 1; i++)
    {
        double p_x      = Xsig_aug(0, i);
        double p_y      = Xsig_aug(1, i);
        double v        = Xsig_aug(2, i);
        double yaw      = Xsig_aug(3, i);
        double yawd     = Xsig_aug(4, i);
        double nu_a     = Xsig_aug(5, i);
        double nu_yawdd = Xsig_aug(6, i);

        //predicted state values
        double px_p, py_p;

        //avoid division by zero
        if (fabs(yawd) > 0.001) {
            px_p = p_x + v / yawd * (sin(yaw + yawd*delta_t) - sin(yaw));
            py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd*delta_t));
        }
        else {
            px_p = p_x + v*delta_t*cos(yaw);
            py_p = p_y + v*delta_t*sin(yaw);
        }

        double v_p = v;
        double yaw_p = yaw + yawd*delta_t;
        double yawd_p = yawd;

        //add noise
        px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
        py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
        v_p = v_p + nu_a*delta_t;

        yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
        yawd_p = yawd_p + nu_yawdd*delta_t;

        //write predicted sigma point into right column
        Xsig_pred(0, i) = px_p;
        Xsig_pred(1, i) = py_p;
        Xsig_pred(2, i) = v_p;
        Xsig_pred(3, i) = yaw_p;
        Xsig_pred(4, i) = yawd_p;
    }

    return Xsig_pred;
}

static std::pair<VectorXd, MatrixXd> PredictMeanAndCovariance(
    int n_x, int n_aug, double lambda,
    const MatrixXd &Xsig_pred, const MatrixXd &weights)
{
    VectorXd x = VectorXd(n_x);

    MatrixXd P = MatrixXd(n_x, n_x);

    //predicted state mean
    x.fill(0.0);
    //iterate over sigma points
    for (int i = 0; i < 2 * n_aug + 1; i++)
    {
        x += weights(i) * Xsig_pred.col(i);
    }

    //predicted state covariance matrix
    P.fill(0.0);
    //iterate over sigma points
    for (int i = 0; i < 2 * n_aug + 1; i++)
    {
        // state difference
        VectorXd x_diff = Xsig_pred.col(i) - x;
        //angle normalization
        while (x_diff(3) > M_PI) x_diff(3) -= 2.*M_PI;
        while (x_diff(3) < -M_PI) x_diff(3) += 2.*M_PI;

        P = P + weights(i) * x_diff * x_diff.transpose();
    }

    return std::make_pair(x, P);
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

    MatrixXd sigmas      = AugmentedSigmaPoints(n_x_, n_aug_, lambda_, std_a_, std_yawdd_, x_, P_);
    MatrixXd predictions = SigmaPointPrediction(n_x_, n_aug_, sigmas, delta_t);
    std::tie(x_, P_)     = PredictMeanAndCovariance(n_x_, n_aug_, lambda_, predictions, weights_);

    //std::cout << "x = " << x_ << ", P = " << P_ << std::endl;
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
