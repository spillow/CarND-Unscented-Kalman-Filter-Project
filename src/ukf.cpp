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

    double v = 8;
    double yaw = 4;
    double yawd = 2;

    P_ <<
        .8, 0, 0, 0, 0,         // px
        0, .8, 0, 0, 0,         // py
        0, 0, v*v, 0, 0,       // v   (assume bikes mostly stay under 15 m/s)
        0, 0, 0, yaw*yaw, 0,   // yaw (normalized angle |yaw| < 2pi ~ 6)
        0, 0, 0, 0, yawd*yawd; // yaw rate (2pi rad/3 s ~ 2 max)?
          
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
        if (use_radar_)
            UpdateRadar(meas_package);
        break;
    case MeasurementPackage::LASER:
        if (use_laser_)
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

static double Normalize(double angle)
{
    int k = (int)(angle / M_2_PI);
    double shift = angle - (double)k * M_2_PI;

    if (shift > M_PI)
        shift -= M_2_PI;
    else if (shift < -M_PI)
        shift += M_2_PI;

    return shift;
};

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
        x_diff(3) = Normalize(x_diff(3));

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
    // Estimate the object's location. Modify the state vector, x_.
    // Predict sigma points, the state, and the state covariance matrix.

    MatrixXd sigmas  = AugmentedSigmaPoints(n_x_, n_aug_, lambda_, std_a_, std_yawdd_, x_, P_);
    Xsig_pred_       = SigmaPointPrediction(n_x_, n_aug_, sigmas, delta_t);
    std::tie(x_, P_) = PredictMeanAndCovariance(n_x_, n_aug_, lambda_, Xsig_pred_, weights_);
}

static void UpdateState(
    int n_x, int n_aug, double lambda, const VectorXd &weights,
    MatrixXd Xsig_pred, MatrixXd Zsig, VectorXd z_pred, MatrixXd S, VectorXd z,
    VectorXd &x, MatrixXd &P)
{
    int n_z = z_pred.size();

    //create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x, n_z);

    //calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        //angle normalization
        z_diff(1) = Normalize(z_diff(1));

        // state difference
        VectorXd x_diff = Xsig_pred.col(i) - x;
        //angle normalization
        x_diff(3) = Normalize(x_diff(3));

        Tc += weights(i) * x_diff * z_diff.transpose();
    }

    //Kalman gain K;
    MatrixXd K = Tc * S.inverse();

    //residual
    VectorXd z_diff = z - z_pred;

    //angle normalization
    z_diff(1) = Normalize(z_diff(1));

    //update state mean and covariance matrix
    x += K * z_diff;
    P += -K*S*K.transpose();
}

static std::pair<VectorXd, MatrixXd> PredictMeasurement(
    MatrixXd &Zsig,
    int n_x, int n_aug, int n_z, double lambda,
    const VectorXd &weights, MatrixXd Xsig_pred,
    std::function<VectorXd(VectorXd)> h)
{
    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug + 1; i++)
    {  //2n+1 simga points
        Zsig.col(i) = h(Xsig_pred.col(i));
    }

    //mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    for (int i = 0; i < 2 * n_aug + 1; i++)
    {
        z_pred += weights(i) * Zsig.col(i);
    }

    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z, n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug + 1; i++)
    {   //2n+1 simga points
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        //angle normalization
        z_diff(1) = Normalize(z_diff(1));

        S += weights(i) * z_diff * z_diff.transpose();
    }

    return std::make_pair(z_pred, S);
}

double calcNIS(VectorXd z, VectorXd z_pred, MatrixXd S)
{
    VectorXd diff = z - z_pred;
    return diff.transpose() * S.inverse() * diff;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  // Use lidar data to update the belief about the object's
  // position. Modify the state vector, x_, and covariance, P_.

  // You'll also need to calculate the lidar NIS.

    auto h = [](VectorXd Xpred)
    {
        double p_x = Xpred(0);
        double p_y = Xpred(1);

        VectorXd zpred(2);
        zpred << p_x, p_y;

        return zpred;
    };

    VectorXd z_pred;
    MatrixXd S;
    MatrixXd Zsig = MatrixXd(2, 2 * n_aug_ + 1);

    std::tie(z_pred, S) = PredictMeasurement(Zsig, n_x_, n_aug_, 2, lambda_, weights_, Xsig_pred_, h);

    //add measurement noise covariance matrix
    MatrixXd R = MatrixXd(2, 2);
    R << std_laspx_*std_laspx_, 0,
         0, std_laspy_*std_laspy_;
    S += R;

    VectorXd z = meas_package.raw_measurements_;

    UpdateState(
        n_x_, n_aug_, lambda_, weights_,
        Xsig_pred_, Zsig, z_pred, S, z,
        x_, P_);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    // Use radar data to update the belief about the object's
    // position. Modify the state vector, x_, and covariance, P_.

    // You'll also need to calculate the radar NIS.

    auto h = [](VectorXd Xpred)
    {
        double p_x = Xpred(0);
        double p_y = Xpred(1);
        double v   = Xpred(2);
        double yaw = Xpred(3);

        double denom = sqrt(p_x*p_x + p_y*p_y);

        if (fabs(denom) < 0.0001)
        {
            std::cout << "h() - Division by Zero" << std::endl;
            return VectorXd(3, 1);
        }

        double v1 = cos(yaw)*v;
        double v2 = sin(yaw)*v;

        // measurement model
        VectorXd zpred(3);
        zpred << denom,                     // r
                 atan2(p_y, p_x),           // phi
                 (p_x*v1 + p_y*v2) / denom; // r_dot

        return zpred;
    };

    VectorXd z_pred;
    MatrixXd S;
    MatrixXd Zsig = MatrixXd(3, 2 * n_aug_ + 1);

    std::tie(z_pred, S) = PredictMeasurement(Zsig, n_x_, n_aug_, 3, lambda_, weights_, Xsig_pred_, h);

    //add measurement noise covariance matrix
    MatrixXd R = MatrixXd(3, 3);
    R << std_radr_*std_radr_, 0, 0,
         0, std_radphi_*std_radphi_, 0,
         0, 0, std_radrd_*std_radrd_;
    S += R;

    VectorXd z = meas_package.raw_measurements_;

    UpdateState(
        n_x_, n_aug_, lambda_, weights_,
        Xsig_pred_, Zsig, z_pred, S, z,
        x_, P_);
}
