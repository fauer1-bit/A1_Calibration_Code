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

    /// Below are a few examples showing some useful data structures and functions.

    // This is a 1D array of 'double' values. Alternatively, you can use 'double mat[25]' but you cannot change it
    // length. With 'std::vector', you can append/delete/insert elements, and much more. The 'std::vector' can store
    // not only 'double', but also any other types of objects. In case you may want to learn more about 'std::vector'
    // check here: https://en.cppreference.com/w/cpp/container/vector
    // std::vector<double> array = {1, 3, 3, 4, 7, 6, 2, 8, 2, 8, 3, 2, 4, 9, 1, 7, 3, 23, 2, 3, 5, 2, 1, 5, 8, 9, 22};
    // array.push_back(5); // append 5 to the array (so the size will increase by 1).
    // array.insert(array.end(), 10, 3);  // append ten 3 (so the size will grow by 10).
    //
    // /// To access the value of an element.
    // double a = array[2];
    //
    // /// define a 2D vector/point
    // Vector2D b(1.1, 2.2);
    //
    // /// define a 3D vector/point
    // Vector3D c(1.1, 2.2, 3.3);
    //
    // /// get the Cartesian coordinates of a (a is treated as Homogeneous coordinates)
    // Vector2D p = c.cartesian();
    //
    // /// get the Homogeneous coordinates of p
    // Vector3D q = p.homogeneous();
    //
    // /// the length of a vector
    // double len = p.length();
    // /// the squared length of a vector
    // double sqr_len = p.length2();
    //
    // /// the dot product of two vectors
    // double dot_prod = dot(p, q);
    //
    // /// the cross product of two vectors
    // Vector cross_prod = cross(c, q);
    //
    // /// normalize this vector
    // cross_prod.normalize();
    //
    // // Define an m-by-n double valued matrix.
    // // Here I use the above array to initialize it. You can also use A(i, j) to initialize/modify/access its elements.
    // const int m = 6, n = 5;
    // Matrix A(m, n, array.data());    // 'array.data()' returns a pointer to the array.
    // std::cout << "M: \n" << A << std::endl;
    //
    // /// define a 3 by 4 matrix (and all elements initialized to 0.0)
    // Matrix M(3, 4, 0.0);
    //
    // /// set first row by a vector
    // M.set_row(0, Vector4D(1.1, 2.2, 3.3, 4.4));
    //
    // /// set second column by a vector
    // M.set_column(1, Vector3D(5.5, 5.5, 5.5));
    //
    // /// define a 3 by 3 matrix (and all elements initialized to 0.0)
    // Matrix33 B;
    //
    // /// define and initialize a 3 by 3 matrix
    // Matrix33 T(1.1, 2.2, 3.3,
    //            0, 2.2, 3.3,
    //            0, 0, 1);
    //
    // /// define and initialize a 3 by 4 matrix
    // Matrix34 P(1.1, 2.2, 3.3, 0,
    //            0, 2.2, 3.3, 1,
    //            0, 0, 1, 1);
    //
    // /// define a 15 by 9 matrix (and all elements initialized to 0.0)
    // Matrix W(15, 9, 0.0);
    // /// set the first row by a 9-dimensional vector
    // W.set_row(0, {0, 1, 2, 3, 4, 5, 6, 7, 8}); // {....} is equivalent to a std::vector<double>
    //
    // /// get the number of rows.
    // int num_rows = W.rows();
    //
    // /// get the number of columns.
    // int num_cols = W.cols();
    //
    // /// get the the element at row 1 and column 2
    // double value = W(1, 2);
    //
    // /// get the last column of a matrix
    // Vector last_column = W.get_column(W.cols() - 1);
    //
    // /// define a 3 by 3 identity matrix
    // Matrix33 I = Matrix::identity(3, 3, 1.0);
    //
    // /// matrix-vector product
    // Vector3D v = M * Vector4D(1, 2, 3, 4); // M is 3 by 4
    //
    // Matrix U(m, m, 0.0);   // initialized with 0s
    // Matrix S(m, n, 0.0);   // initialized with 0s
    // Matrix V(n, n, 0.0);   // initialized with 0s
    //
    // // Compute the SVD decomposition of A
    // svd_decompose(A, U, S, V);
    //
    // // Now let's check if the SVD result is correct
    //
    // // Check 1: U is orthogonal, so U * U^T must be identity
    // std::cout << "U*U^T: \n" << U * transpose(U) << std::endl;
    //
    // // Check 2: V is orthogonal, so V * V^T must be identity
    // std::cout << "V*V^T: \n" << V * transpose(V) << std::endl;
    //
    // // Check 3: S must be a diagonal matrix
    // std::cout << "S: \n" << S << std::endl;
    //
    // // Check 4: according to the definition, A = U * S * V^T
    // std::cout << "M - U * S * V^T: \n" << A - U * S * transpose(V) << std::endl;
    //
    // // Compute the inverse of a matrix
    // Matrix invT;
    // inverse(T, invT);
    // // Let's check if the inverse is correct
    // std::cout << "T * invT: \n" << T * invT << std::endl;
    //
    // // TODO: the above code just demonstrates some useful data structures and APIs. Please remove all above code in your
    // //       final submission.

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

    if (points_3d.size() < 6 || points_2d.size() < 6 || points_3d.size() != points_2d.size()) {
        std::cerr << "Invalid input: need at least 6 correspondences and 3D and 2D size must match" << std::endl;
        return false;
    }

    // TODO: construct the P matrix (so P * m = 0).

    // The p matrix is a 12*12 matrix if we use 6 correspondences, otherwise it yields 2n homogeneous linear equations, which make it a 2n*12 where n = correspondences

    const int correspondences_p = points_3d.size(), mp_length = 2 * correspondences_p;

    Matrix P(mp_length, 12, 0.0); // matrix P has length of m = 2*points and width of 12. (12 *12)

    for (int i = 0; i < correspondences_p; i++) {

        // set point variables for 3d (X,Y,Z) and for 2d (image plane) (u,v)
        const double u = -points_2d[i][0]; //set to negative
        const double v = -points_2d[i][1]; //set to negative
        const double X = points_3d[i][0];
        const double Y = points_3d[i][1];
        const double Z = points_3d[i][2];

        //assign the values of the points to the variables in the matrix
        // set row and then loop over in pairs of 2

        //row 1 (i * 2) at 0 = 0, at i = 1, skips row1 and moves automatically to row 2...
        P.set_row(i*2, {X, Y, Z, 1.0, 0.0, 0.0, 0.0, 0.0, u * X, u * Y, u * Z, u});

        //row 2 (i * 2 + 1), always next after row 1..
        P.set_row(i*2+1, {X, Y, Z, 1.0, 0.0, 0.0, 0.0, 0.0, v * X, v * Y, v * Z, v});

    }

    // Check matrix
    std::cout << "Matrix P:" << std::endl;
    for (int i = 0; i < mp_length; ++i) {
        std::cout << "\tRow " << i << ": " << P.get_row(i) << std::endl;
    }


    // TODO: solve for M (the whole projection matrix, i.e., M = K * [R, t]) using SVD decomposition.
    //   Optional: you can check if your M is correct by applying M on the 3D points. If correct, the projected point
    //             should be very close to your input images points.
    // let P = UDV^T, set m equal to the lac column of V. (where

    // Matrix U, D and V initialise for SVD on matrix P
    Matrix U(mp_length, mp_length);
    Matrix D(mp_length, 12); // diagonal matrix
    Matrix V(12, 12);

    svd_decompose(P, U, D, V);

    //get last vector for m
    Vector m = V.get_column(V.cols() -1);

    //Construct Matrix M from m
    Matrix34 M(
        m[0], m[1], m[2], m[3],
        m[4], m[5], m[6], m[7],
        m[8], m[9], m[10], m[11]
    );

    //check matrix
    std::cout << "Projection Matrix M: \n" << M << std::endl;

    // extract the R (rotation) and translation vector t from M

    Vector3D a1(M[0][0], M[0][1], M[0][2]);
    Vector3D a2(M[1][0], M[1][1], M[1][2]);
    Vector3D a3(M[2][0], M[2][1], M[2][2]);

    Vector3D b(M[0][3], M[1][3], M[2][3]);


    // TODO: extract intrinsic parameters from M.



    // TODO: extract extrinsic parameters from M.

    // TODO: make sure the recovered parameters are passed to the corresponding variables (fx, fy, cx, cy, s, R, and t)

    std::cout << "\n\tTODO: After you implement this function, please return 'true' - this will trigger the viewer to\n"
                 "\t\tupdate the rendering using your recovered camera parameters. This can help you to visually check\n"
                 "\t\tif your calibration is successful or not.\n\n" << std::flush;
    return false;
}





