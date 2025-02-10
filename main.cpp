#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <chrono>

#include "mpc_utils.hpp"
#include "mpc_controller.hpp"
#include "OsqpEigen/OsqpEigen.h"
#include "matplotlibcpp.h"

using namespace std::chrono;
namespace plt = matplotlibcpp;

int main() {
    std::cout << "Starting MPC Simulation..." << std::endl;

    // allocate the initial and the reference state space
    Eigen::Matrix<double, 4, 1> x0;
    Eigen::Matrix<double, 2, 1> ctr;
    x0 << 0, 0, 0, 0;
    ctr << 0, 0;

    double dl = 0.1;  // Course tick
    std::vector<double> cx, cy, cyaw, ck;

    //getStraightCourse(dl, cx, cy, cyaw, ck);
    //getStraightCourse2(dl, cx, cy, cyaw, ck);
    //getStraightCourse3(dl, cx, cy, cyaw, ck);
    //getForwardCourse(dl, cx, cy, cyaw, ck);

    getSwitchBackCourse(dl, cx, cy, cyaw, ck);
    x0 << 2, 0, 0, 0; // start and start is the same

    std::vector<double> sp = calcSpeedProfile(cx, cy, cyaw, TARGET_SPEED);
    // set the preview window
    int mpcWindow = 5;

    // allocate the dynamics matrices
    Eigen::Matrix<double, 4, 4> a;
    Eigen::Matrix<double, 4, 2> b;
    Eigen::Matrix<double, 4, 1> c;

    // allocate the constraints vector
    Eigen::Matrix<double, 4, 1> xMax;
    Eigen::Matrix<double, 4, 1> xMin;
    Eigen::Matrix<double, 2, 1> uMax;
    Eigen::Matrix<double, 2, 1> uMin;

    // allocate the weight matrices
    Eigen::DiagonalMatrix<double, 4> Q;
    Eigen::DiagonalMatrix<double, 2> R;


    // allocate QP problem matrices and vectors
    Eigen::SparseMatrix<double> hessian;
    Eigen::VectorXd gradient;
    Eigen::SparseMatrix<double> linearMatrix;
    Eigen::VectorXd lowerBound;
    Eigen::VectorXd upperBound;

    // set the initial and the desired states
    auto [target_ind, mind] = calcNearestIndex(x0, cx, cy, cyaw, 0, N_IND_SEARCH);

    if (x0(3) - cyaw[0] >= M_PI) {
        x0(3) -= 2 * M_PI;
    } else if (x0(3) - cyaw[0] <= -M_PI) {
        x0(3) += 2 * M_PI;
    }

    double goal_x = cx.back(), goal_y = cy.back();

    smoothYaw(cyaw);

    Eigen::MatrixXd xRef(4, mpcWindow + 1);
    for (int t = 0; t <= mpcWindow; t++) {
        xRef(0, t) = cx[t];
        xRef(1, t) = cy[t];
        xRef(2, t) = sp[t];
        xRef(3, t) = cyaw[t];
    }

    // set MPC problem quantities for direct linearized model
    //double v = x0(2);
    //double phi = x0(3);
    //double delta = 0;
    //setDynamicsMatrices(a, b, c, v, phi, delta);

    setDynamicsMatrices(a, b, c, x0, ctr, DT);
    setInequalityConstraints(xMax, xMin, uMax, uMin);
    setWeightMatrices(Q, R);

    // cast the MPC problem as QP problem
    castMPCToQPHessian(Q, R, mpcWindow, hessian);
    castMPCToQPGradient(Q, xRef, mpcWindow, gradient);
    castMPCToQPConstraintMatrix(a, b, mpcWindow, linearMatrix);
    castMPCToQPConstraintVectors(xMax, xMin, uMax, uMin, x0, mpcWindow, lowerBound, upperBound);

    // instantiate the solver
    OsqpEigen::Solver solver;

    // settings
    // solver.settings()->setVerbosity(false);
    solver.settings()->setWarmStart(true);

    // set the initial data of the QP solver
    solver.data()->setNumberOfVariables(4 * (mpcWindow + 1) + 2 * mpcWindow);
    solver.data()->setNumberOfConstraints(2 * 4 * (mpcWindow + 1) + 2 * mpcWindow);
    solver.data()->setHessianMatrix(hessian);
    solver.data()->setGradient(gradient);
    solver.data()->setLinearConstraintsMatrix(linearMatrix);
    solver.data()->setLowerBound(lowerBound);
    solver.data()->setUpperBound(upperBound);

    // instantiate the solver
    solver.initSolver();
    Eigen::VectorXd QPSolution;

    // number of iteration steps
    int numberOfSteps = 20000;

    std::vector<double> x_data, y_data, delta_data; // Store trajectory data

    for (int i = 0; i < numberOfSteps; i++)
    {
        auto loop_start = high_resolution_clock::now();
        // direct linearization
        //setDynamicsMatrices(a, b, c, v, phi, delta);

        //RK4
        setDynamicsMatrices(a, b, c, x0, ctr, DT);

        auto [xRef, target_ind_new, dref] = calcRefTrajectory(x0, cx, cy, cyaw, ck, sp, dl, target_ind, mpcWindow, NX, N_IND_SEARCH, DT);
        target_ind = target_ind_new;

        castMPCToQPGradient(Q, xRef, mpcWindow, gradient);
        solver.updateGradient(gradient);

        castMPCToQPConstraintMatrix(a, b, mpcWindow, linearMatrix);
        solver.updateLinearConstraintsMatrix(linearMatrix);

        castMPCToQPConstraintVectors(xMax, xMin, uMax, uMin, x0, mpcWindow, lowerBound, upperBound);
        updateConstraintVectors(x0, lowerBound, upperBound);
        solver.updateBounds(lowerBound, upperBound);

        // solve the QP problem
        if (solver.solveProblem() != OsqpEigen::ErrorExitFlag::NoError)
            return 1;

        // get the controller input
        QPSolution = solver.getSolution();
        ctr = QPSolution.block(4 * (mpcWindow + 1), 0, 2, 1);

        // propagate the model
        double v_new = x0(2) + ctr(0) * DT;

        x0(3) += (v_new / WB) * tan(ctr(1)) * DT;
        x0(0) += v_new * cos(x0(3)) * DT;
        x0(1) += v_new * sin(x0(3)) * DT;
        x0(2) = v_new;

        //x0 = a * x0 + b * ctr + c; linearized model

        // for using in direct linearized model without nonlinearity
        //v = x0(2);
        //phi = x0(3);
        //delta = ctr(1);

        // Goal condition
        if (std::hypot(x0(0) - goal_x, x0(1) - goal_y) <= GOAL_DIS && std::abs(x0(2)) < STOP_SPEED) {
            std::cout << "Goal reached!" << std::endl;
            break;
        }

        // Store x and y positions
        x_data.push_back(x0(0));
        y_data.push_back(x0(1));
        delta_data.push_back(ctr(1));


        auto loop_end = high_resolution_clock::now();
        duration<double> loop_duration = loop_end - loop_start;
        double loop_freq = 1.0 / loop_duration.count();
        std::cout << "Loop Frequency: " << loop_freq << " Hz" << std::endl;

    }
    plt::figure();

    // First subplot: Vehicle trajectory vs. reference trajectory
    plt::subplot(2, 1, 1);  // (rows, cols, index)
    plt::named_plot("MPC Trajectory", x_data, y_data, "r-");
    plt::named_plot("Reference", cx, cy, "b--");

    // Labels and formatting
    plt::xlabel("X Position");
    plt::ylabel("Y Position");
    plt::title("MPC Trajectory vs. Reference");
    plt::legend();
    plt::grid(true);

    // Second subplot: Steering angle (delta) over time
    plt::subplot(2, 1, 2);  // Second subplot for delta
    plt::plot(delta_data, "g-"); // Green line for delta

    // Labels and formatting
    plt::xlabel("Time Steps");
    plt::ylabel("Steering Angle (delta) [rad]");
    plt::title("Steering Control Over Time");
    plt::grid(true);

    // Show all subplots
    plt::show();
    return 0;


}
