/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "calibration.h"
#include "matrix_algo.h"


using namespace easy3d;



/**
 * TODO: Finish this function for calibrating a camera from the corresponding 3D-2D point pairs.
 *       You may define a few functions for some sub-tasks.
 * @return True on success, otherwise false. On success, the camera parameters are returned by fx, fy, cx, cy, skew, R, and t).
 */
bool Calibration::calibration(
        const std::vector<Vector3D>& points_3d, /// input: An array of 3D points.
        const std::vector<Vector2D>& points_2d, /// input: An array of 2D image points.
        double& fx,  /// output: focal length (i.e., K[0][0]).
        double& fy,  /// output: focal length (i.e., K[1][1]).
        double& cx,  /// output: x component of the principal point (i.e., K[0][2]).
        double& cy,  /// output: y component of the principal point (i.e., K[1][2]).
        double& s,   /// output: skew factor (i.e., K[0][1]), which is s = -alpha * cot(theta).
        Matrix33& R, /// output: the 3x3 rotation matrix encoding camera rotation.
        Vector3D& t) /// output：a 3D vector encoding camera translation.
{
    std::cout << "\nTODO: implement the 'calibration()' function in the file 'Calibration/calibration_method.cpp'\n\n";

    std::cout << "[Liangliang]:\n"
                 "\tIn this assignment, two essential data structures, 'Matrix' and 'Vector', are provided for the\n"
                 "\tmanipulation and storage of matrices and vectors. These data structures are defined in:\n"
                 "\t    - Calibration/matrix.h: handles matrices of arbitrary dimensions and related functions.\n"
                 "\t    - Calibration/vector.h: manages vectors of arbitrary sizes and related functions.\n"
                 "\tCamera calibration requires computing the SVD and inverse of matrices. These functions, along\n"
                 "\twith several other relevant ones, are provided in:\n"
                 "\t    - Calibration/matrix_algo.h: contains functions for determinant, inverse, SVD, linear least-squares...\n"
                 "\tIn the 'Calibration::calibration(...)' function, code snippets are provided for your reference.\n"
                 "\tFor more details about these data structures and a complete list of related functions, please\n"
                 "\trefer to the header files mentioned above.\n\n"
                 "\tFor your final submission, adhere to the following guidelines:\n"
                 "\t    - submit ONLY the 'Calibration/calibration_method.cpp' file.\n"
                 "\t    - remove ALL unrelated test code, debugging code, and comments.\n"
                 "\t    - ensure that your code compiles and can reproduce your results WITHOUT ANY modification.\n\n" << std::flush;

    // TODO: the above code just demonstrates some useful data structures and APIs. Please remove all above code in your
    //       final submission.

    //--------------------------------------------------------------------------------------------------------------
    // implementation starts ...

    std::cout << "\n[Liangliang]:\n"
                 "\tThis function takes two arrays as input parameters:\n"
                 "\t\t- points_3d: An array of 3D points representing the scene\n"
                 "\t\t- points_2d: An array of 2D image points corresponding to the 3D points\n"
                 "\tThe function should return either 'true' upon successful calibration or 'false' otherwise.\n"
                 "\tUpon success, the following parameters must be stored in the specified variables:\n"
                 "\t\t- fx and fy: focal lengths along the x and y axes, respectively\n"
                 "\t\t- cx and cy: coordinates of the principal point\n"
                 "\t\t- s: the skew factor, i.e., s = -alpha * cot(theta)\n"
                 "\t\t- R: the 3x3 rotation matrix encoding camera orientation\n"
                 "\t\t- t: a 3D vector encoding camera location.\n"
                 "\tIMPORTANT: don't forget to write your recovered parameters to the above variables." << std::endl;

    // TODO: check if input is valid (e.g., number of correspondences >= 6, sizes of 2D/3D points must match)

    if (points_3d.size()<6 || points_2d.size()<6 || points_3d.size() != points_2d.size()) {
        std::cout << "Wrong inputs!" << "\n";
        std::cout << "3d: " << points_3d.size() << "\n";
        std::cout << "2d: " << points_2d.size() << std::endl;
    } else {
        std::cout << "Points_3d size: " << points_3d.size() << "\n";
        std::cout << "Points_2d size: " << points_2d.size() << std::endl;
    }

    // TODO: construct the P matrix (so P * m = 0).

    const int n = points_3d.size();
    const int matrix_length = n * 2;

    Matrix P(matrix_length, 12, 0.0); // 12x12

    int j = 0;

    for (int i = 0; i < matrix_length; ++i)
    {
        double u = points_2d[j][0];
        double v = points_2d[j][1];

        const double X = points_3d[j][0];
        const double Y = points_3d[j][1];
        const double Z = points_3d[j][2];

        j = j + 1;

        u = u * -1.0;
        v = v * -1.0;

        P[i][0] = X;
        P[i][1] = Y;
        P[i][2] = Z;
        P[i][3] = 1.0;

        P[i][8] = u * X;
        P[i][9] = u * Y;
        P[i][10] = u * Z;
        P[i][11] = u;

        i = i + 1;

        P[i][4] = X;
        P[i][5] = Y;
        P[i][6] = Z;
        P[i][7] = 1.0;

        P[i][8] = v * X;
        P[i][9] = v * Y;
        P[i][10] = v * Z;
        P[i][11] = v;
    }

    std::cout << "matrix looks as follows: " << std::endl;
    for (int i=0; i<matrix_length; ++i) {
        std::cout << "\t" << i << ": " << P.get_row(i) << std::endl;
    }


    // TODO: solve for M (the whole projection matrix, i.e., M = K * [R, t]) using SVD decomposition.
    //   Optional: you can check if your M is correct by applying M on the 3D points. If correct, the projected point
    //             should be very close to your input images points.

    // SVD --> A = UDV^T --> P = UDV^T = (2n x 2n) (2n x 12) (12 x 12)
    // p = MP = K [R t]P

    Matrix U(matrix_length, matrix_length, 0.0);   // initialized with 0s
    Matrix D(matrix_length, 12, 0.0);   // initialized with 0s
    Matrix V(12, 12, 0.0);   // initialized with 0s

    svd_decompose(P,U,D,V);

    Vector m = V.get_column(V.cols()-1); // get m from last column of V

    Matrix34 M(m[0],m[1],m[2],m[3],
                m[4],m[5],m[6],m[7],
                m[8],m[9],m[10],m[11]);

    std::cout << "M: \n" << M << std::endl;

    Vector3D a1(M[0][0],M[0][1],M[0][2]);
    Vector3D a2(M[1][0],M[1][1],M[1][2]);
    Vector3D a3(M[2][0],M[2][1],M[2][2]);

    Vector3D b(M[0][3],M[1][3],M[2][3]);

    // TODO: extract intrinsic parameters from M.

    // DOES NOT HANGLE SIGN OF p YET !!!!
    double reau = 1 / a3.norm();
    double reau2 = reau * reau;

    std::cout << "reau: " << reau << "\n" << std::endl;

    cx = reau2 * dot(a1,a3);
    cy = reau2 * dot(a2,a3);

    std::cout << "cx: " << cx << "\n" << std::endl;
    std::cout << "cy: " << cy << "\n" << std::endl;

    // Formula says dot product, but program does not take dot product on 1x1 vector, so that is why * is used
    double cos_angle = -1 * (dot(cross(a1,a3),cross(a2,a3)) / (cross(a1,a3).norm() * cross(a2,a3).norm()));
    double sin_angle = sqrt(1-(cos_angle * cos_angle));

    std::cout << "cos_angle: " << cos_angle << "\n" << std::endl;
    std::cout << "sin_angle: " << sin_angle << "\n" << std::endl;

    double alpha = reau2 * cross(a1,a3).norm() * sin_angle;
    double beta = reau2 * cross(a2,a3).norm() * sin_angle;

    std::cout << "alpha: " << alpha << "\n" << std::endl;
    std::cout << "beta: " << beta << "\n" << std::endl;

    // TODO: extract extrinsic parameters from M.

    Vector3D r1 = cross(a2,a3) / cross(a2,a3).norm();
    Vector3D r3 = reau * a3;
    Vector3D r2 = cross(r3,r1);

    std::cout << "r1: " << r1 << "\n" << std::endl;
    std::cout << "r2: " << r2 << "\n" << std::endl;
    std::cout << "r3: " << r3 << "\n" << std::endl;

    // TODO: make sure the recovered parameters are passed to the corresponding variables (fx, fy, cx, cy, s, R, and t)

    fx = alpha;
    fy = beta/sin_angle;
    s = -1 * alpha * cos_angle/sin_angle;

    R.set_row(0,{r1[0],r1[1],r1[2]});
    R.set_row(1,{r2[0],r2[1],r2[2]});
    R.set_row(2,{r3[0],r3[1],r3[2]});

    std::cout << "R " << R << "\n" <<  std::endl;

    Matrix33 K(fx,s,cx,
                0.0,fy,cy,
                0.0,0.0,1.0);

    std::cout << "K " << K << std::endl;

    t= reau * inverse(K) * b;

    std::cout << "t: " << t << std::endl;


    std::cout << "\n\tTODO: After you implement this function, please return 'true' - this will trigger the viewer to\n"
                 "\t\tupdate the rendering using your recovered camera parameters. This can help you to visually check\n"
                 "\t\tif your calibration is successful or not.\n\n" << std::flush;
    return true;
}