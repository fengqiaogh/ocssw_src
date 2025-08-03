/* --------------------------------------------------------------------------------------- */
/* T.S. Kostadinov, D.A. Siegel, and S. Maritorena Particle Size Distribution Algorithm    */
/*                                                                                         */
/* References:                                                                             */
/*                                                                                         */
/*    Kostadinoc, T.S. , D.A. Siegel, and S. Maritorena (2009), Retrieval of the particle  */
/*    size distribution from satellite ocean color observations, J. Geophys. Res., 114,    */
/*    C09015, doi:10.1029/2009JC005303H.                                                   */
/*                                                                                         */
/* Implementation:                                                                         */
/* Program written by Robert Lossing September 2014 for NASA GSFC OBPG             */
/*                                                             */
/* --------------------------------------------------------------------------------------- */
/**********************************************************************************************************************************************************/
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "l12_proto.h"

// The five (5) look up tables "LUTs" provided by Dr. Tihomir Kostadinov are below:

float bbp_slope[449] = {-1.5, -1.49, -1.48, -1.47, -1.46, -1.45, -1.44, -1.43,
    -1.42, -1.41, -1.4, -1.39, -1.38, -1.37, -1.36, -1.35, -1.34, -1.33,
    -1.32, -1.31, -1.3, -1.29, -1.28, -1.27, -1.26, -1.25, -1.24, -1.23,
    -1.22, -1.21, -1.2, -1.19, -1.18, -1.17, -1.16, -1.15, -1.14, -1.13,
    -1.12, -1.11, -1.1, -1.09, -1.08, -1.07, -1.06, -1.05, -1.04, -1.03,
    -1.02, -1.01, -1, -0.99, -0.98, -0.97, -0.96, -0.95, -0.94, -0.93,
    -0.92, -0.91, -0.9, -0.89, -0.88, -0.87, -0.86, -0.85, -0.84, -0.83,
    -0.82, -0.81, -0.8, -0.79, -0.78, -0.77, -0.76, -0.75, -0.74, -0.73,
    -0.72, -0.71, -0.7, -0.69, -0.68, -0.67, -0.66, -0.65, -0.64, -0.63,
    -0.62, -0.61, -0.6, -0.59, -0.58, -0.57, -0.56, -0.55, -0.54, -0.53,
    -0.52, -0.51, -0.5, -0.49, -0.48, -0.47, -0.46, -0.45, -0.44, -0.43,
    -0.42, -0.41, -0.4, -0.39, -0.38, -0.37, -0.36, -0.35, -0.34, -0.33,
    -0.32, -0.31, -0.3, -0.29, -0.28, -0.27, -0.26, -0.25, -0.24, -0.23,
    -0.22, -0.21, -0.2, -0.19, -0.18, -0.17, -0.16, -0.15, -0.14, -0.13,
    -0.12, -0.11, -0.1, -0.09, -0.08, -0.07, -0.06, -0.05, -0.04, -0.03,
    -0.02, -0.01, 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
    0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21,
    0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33,
    0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45,
    0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57,
    0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69,
    0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81,
    0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93,
    0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1, 1.01, 1.02, 1.03, 1.04, 1.05,
    1.06, 1.07, 1.08, 1.09, 1.1, 1.11, 1.12, 1.13, 1.14, 1.15, 1.16, 1.17,
    1.18, 1.19, 1.2, 1.21, 1.22, 1.23, 1.24, 1.25, 1.26, 1.27, 1.28, 1.29,
    1.3, 1.31, 1.32, 1.33, 1.34, 1.35, 1.36, 1.37, 1.38, 1.39, 1.4, 1.41,
    1.42, 1.43, 1.44, 1.45, 1.46, 1.47, 1.48, 1.49, 1.5, 1.51, 1.52, 1.53,
    1.54, 1.55, 1.56, 1.57, 1.58, 1.59, 1.6, 1.61, 1.62, 1.63, 1.64, 1.65,
    1.66, 1.67, 1.68, 1.69, 1.7, 1.71, 1.72, 1.73, 1.74, 1.75, 1.76, 1.77,
    1.78, 1.79, 1.8, 1.81, 1.82, 1.83, 1.84, 1.85, 1.86, 1.87, 1.88, 1.89,
    1.9, 1.91, 1.92, 1.93, 1.94, 1.95, 1.96, 1.97, 1.98, 1.99, 2, 2.01,
    2.02, 2.03, 2.04, 2.05, 2.06, 2.07, 2.08, 2.09, 2.1, 2.11, 2.12, 2.13,
    2.14, 2.15, 2.16, 2.17, 2.18, 2.19, 2.2, 2.21, 2.22, 2.23, 2.24, 2.25,
    2.26, 2.27, 2.28, 2.29, 2.3, 2.31, 2.32, 2.33, 2.34, 2.35, 2.36, 2.37,
    2.38, 2.39, 2.4, 2.41, 2.42, 2.43, 2.44, 2.45, 2.46, 2.47, 2.48, 2.49,
    2.5, 2.51, 2.52, 2.53, 2.54, 2.55, 2.56, 2.57, 2.58, 2.59, 2.6, 2.61,
    2.62, 2.63, 2.64, 2.65, 2.66, 2.67, 2.68, 2.69, 2.7, 2.71, 2.72, 2.73,
    2.74, 2.75, 2.76, 2.77, 2.78, 2.79, 2.8, 2.81, 2.82, 2.83, 2.84, 2.85,
    2.86, 2.87, 2.88, 2.89, 2.9, 2.91, 2.92, 2.93, 2.94, 2.95, 2.96, 2.97,
    2.98};

float PSD_slope[449] = {3.04, 3.042, 3.046, 3.05, 3.053, 3.056, 3.059, 3.062,
    3.064, 3.066, 3.07, 3.073, 3.076, 3.079, 3.081, 3.084, 3.087, 3.089,
    3.093, 3.095, 3.096, 3.1, 3.103, 3.106, 3.109, 3.112, 3.114, 3.117,
    3.12, 3.123, 3.126, 3.129, 3.132, 3.134, 3.137, 3.14, 3.143, 3.145,
    3.148, 3.15, 3.153, 3.156, 3.159, 3.162, 3.165, 3.168, 3.171, 3.174,
    3.177, 3.18, 3.182, 3.185, 3.189, 3.192, 3.195, 3.198, 3.2, 3.204,
    3.208, 3.211, 3.214, 3.216, 3.219, 3.223, 3.225, 3.229, 3.232, 3.235,
    3.237, 3.241, 3.244, 3.247, 3.249, 3.252, 3.256, 3.259, 3.263, 3.266,
    3.269, 3.272, 3.275, 3.278, 3.281, 3.284, 3.288, 3.291, 3.295, 3.298,
    3.302, 3.306, 3.309, 3.313, 3.315, 3.319, 3.321, 3.325, 3.328, 3.331,
    3.335, 3.338, 3.341, 3.344, 3.348, 3.352, 3.356, 3.36, 3.363, 3.366,
    3.369, 3.373, 3.377, 3.379, 3.383, 3.386, 3.39, 3.394, 3.398, 3.402,
    3.405, 3.409, 3.412, 3.414, 3.418, 3.422, 3.426, 3.43, 3.433, 3.437,
    3.441, 3.445, 3.449, 3.452, 3.457, 3.46, 3.463, 3.468, 3.471, 3.474,
    3.478, 3.482, 3.485, 3.49, 3.493, 3.497, 3.501, 3.506, 3.51, 3.513,
    3.518, 3.522, 3.526, 3.53, 3.535, 3.539, 3.544, 3.548, 3.553, 3.556,
    3.561, 3.566, 3.57, 3.575, 3.579, 3.584, 3.588, 3.592, 3.596, 3.601,
    3.606, 3.61, 3.614, 3.619, 3.624, 3.628, 3.633, 3.638, 3.642, 3.647,
    3.651, 3.656, 3.662, 3.666, 3.671, 3.676, 3.682, 3.687, 3.693, 3.698,
    3.703, 3.709, 3.714, 3.719, 3.725, 3.73, 3.736, 3.742, 3.747, 3.753,
    3.757, 3.763, 3.769, 3.775, 3.78, 3.786, 3.792, 3.798, 3.805, 3.811,
    3.817, 3.823, 3.829, 3.835, 3.841, 3.847, 3.853, 3.86, 3.866, 3.873,
    3.879, 3.885, 3.892, 3.898, 3.905, 3.911, 3.918, 3.924, 3.931, 3.938,
    3.944, 3.951, 3.958, 3.965, 3.972, 3.979, 3.986, 3.993, 4, 4.007, 4.014,
    4.021, 4.028, 4.036, 4.043, 4.05, 4.058, 4.065, 4.073, 4.08, 4.088,
    4.096, 4.103, 4.111, 4.119, 4.127, 4.134, 4.142, 4.15, 4.158, 4.166,
    4.174, 4.182, 4.19, 4.199, 4.207, 4.215, 4.223, 4.231, 4.24, 4.248,
    4.257, 4.265, 4.273, 4.282, 4.29, 4.299, 4.308, 4.316, 4.325, 4.334,
    4.342, 4.351, 4.36, 4.369, 4.378, 4.386, 4.395, 4.404, 4.413, 4.422,
    4.431, 4.44, 4.449, 4.458, 4.467, 4.476, 4.486, 4.495, 4.504, 4.513,
    4.522, 4.532, 4.541, 4.55, 4.559, 4.569, 4.578, 4.587, 4.597, 4.606,
    4.616, 4.625, 4.634, 4.644, 4.653, 4.663, 4.672, 4.682, 4.691, 4.701,
    4.71, 4.72, 4.729, 4.739, 4.749, 4.758, 4.768, 4.777, 4.787, 4.797,
    4.806, 4.816, 4.826, 4.835, 4.845, 4.855, 4.864, 4.874, 4.884, 4.893,
    4.903, 4.913, 4.922, 4.932, 4.942, 4.952, 4.961, 4.971, 4.981, 4.991,
    5.001, 5.01, 5.02, 5.03, 5.04, 5.05, 5.059, 5.069, 5.079, 5.089, 5.099,
    5.108, 5.118, 5.128, 5.138, 5.148, 5.158, 5.168, 5.177, 5.187, 5.197,
    5.207, 5.217, 5.227, 5.237, 5.247, 5.256, 5.266, 5.276, 5.286, 5.296,
    5.306, 5.316, 5.326, 5.336, 5.346, 5.356, 5.366, 5.375, 5.385, 5.395,
    5.405, 5.415, 5.425, 5.435, 5.445, 5.455, 5.465, 5.475, 5.485, 5.495,
    5.505, 5.515, 5.525, 5.535, 5.545, 5.555, 5.565, 5.575, 5.585, 5.595,
    5.605, 5.615, 5.625, 5.635, 5.645, 5.656, 5.666, 5.676, 5.686, 5.696,
    5.706, 5.716, 5.726, 5.736, 5.747, 5.757, 5.767, 5.777, 5.787, 5.797,
    5.808, 5.818, 5.828, 5.838, 5.849, 5.859, 5.869, 5.879, 5.89, 5.9, 5.91,
    5.921, 5.931, 5.941, 5.952, 5.962, 5.973, 5.983, 5.993};

float PSD_slope_sigma[449] = {0.209, 0.209, 0.207, 0.206, 0.205, 0.205, 0.204,
    0.203, 0.203, 0.204, 0.202, 0.202, 0.201, 0.2, 0.201, 0.2, 0.2, 0.2,
    0.199, 0.199, 0.2, 0.199, 0.198, 0.198, 0.196, 0.197, 0.196, 0.196,
    0.195, 0.195, 0.195, 0.194, 0.194, 0.194, 0.194, 0.193, 0.192, 0.193,
    0.193, 0.193, 0.193, 0.193, 0.192, 0.192, 0.191, 0.19, 0.19, 0.189,
    0.189, 0.188, 0.19, 0.189, 0.188, 0.187, 0.186, 0.186, 0.187, 0.186,
    0.185, 0.185, 0.184, 0.185, 0.185, 0.183, 0.183, 0.182, 0.182, 0.181,
    0.182, 0.182, 0.181, 0.181, 0.182, 0.183, 0.181, 0.18, 0.179, 0.18,
    0.178, 0.179, 0.179, 0.179, 0.179, 0.179, 0.178, 0.177, 0.177, 0.175,
    0.174, 0.173, 0.173, 0.173, 0.173, 0.173, 0.174, 0.174, 0.174, 0.173,
    0.172, 0.173, 0.172, 0.173, 0.172, 0.17, 0.169, 0.167, 0.167, 0.168,
    0.169, 0.168, 0.167, 0.169, 0.167, 0.169, 0.167, 0.166, 0.164, 0.163,
    0.164, 0.163, 0.163, 0.165, 0.164, 0.163, 0.162, 0.161, 0.161, 0.16,
    0.16, 0.159, 0.158, 0.158, 0.156, 0.157, 0.159, 0.156, 0.156, 0.157,
    0.156, 0.155, 0.156, 0.154, 0.154, 0.155, 0.153, 0.151, 0.15, 0.151,
    0.149, 0.149, 0.149, 0.147, 0.145, 0.145, 0.143, 0.143, 0.141, 0.141,
    0.139, 0.138, 0.136, 0.135, 0.134, 0.133, 0.134, 0.134, 0.135, 0.135,
    0.133, 0.134, 0.134, 0.132, 0.13, 0.131, 0.129, 0.13, 0.131, 0.129,
    0.131, 0.13, 0.127, 0.128, 0.126, 0.126, 0.123, 0.122, 0.12, 0.12,
    0.118, 0.116, 0.114, 0.113, 0.111, 0.109, 0.107, 0.105, 0.104, 0.102,
    0.108, 0.106, 0.104, 0.104, 0.102, 0.1, 0.098, 0.095, 0.093, 0.091,
    0.09, 0.088, 0.086, 0.084, 0.086, 0.083, 0.081, 0.079, 0.078, 0.076,
    0.075, 0.074, 0.072, 0.071, 0.07, 0.07, 0.068, 0.067, 0.066, 0.066,
    0.064, 0.063, 0.063, 0.062, 0.061, 0.059, 0.058, 0.057, 0.056, 0.056,
    0.055, 0.054, 0.053, 0.052, 0.051, 0.05, 0.049, 0.048, 0.047, 0.046,
    0.045, 0.045, 0.044, 0.043, 0.042, 0.041, 0.041, 0.04, 0.039, 0.038,
    0.038, 0.037, 0.036, 0.035, 0.035, 0.034, 0.034, 0.033, 0.032, 0.032,
    0.031, 0.03, 0.03, 0.029, 0.029, 0.028, 0.028, 0.027, 0.027, 0.026,
    0.026, 0.025, 0.025, 0.024, 0.024, 0.023, 0.023, 0.022, 0.022, 0.022,
    0.021, 0.021, 0.02, 0.02, 0.02, 0.019, 0.019, 0.019, 0.018, 0.018,
    0.018, 0.017, 0.017, 0.017, 0.016, 0.016, 0.016, 0.016, 0.015, 0.015,
    0.015, 0.015, 0.014, 0.014, 0.014, 0.014, 0.014, 0.013, 0.013, 0.013,
    0.013, 0.013, 0.012, 0.012, 0.012, 0.012, 0.012, 0.012, 0.011, 0.011,
    0.011, 0.011, 0.011, 0.011, 0.011, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
    0.01, 0.01, 0.01, 0.01, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009,
    0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.009, 0.008, 0.008,
    0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008,
    0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008,
    0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008, 0.008,
    0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007,
    0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007,
    0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007,
    0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007,
    0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007,
    0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007, 0.007};

float log10_of_bbp440_over_No[449] = {-18.899, -18.9, -18.902, -18.903,
    -18.904, -18.905, -18.906, -18.907, -18.908, -18.909, -18.91, -18.911,
    -18.912, -18.913, -18.913, -18.914, -18.914, -18.915, -18.915, -18.916,
    -18.916, -18.917, -18.918, -18.919, -18.921, -18.92, -18.921, -18.922,
    -18.924, -18.924, -18.925, -18.925, -18.926, -18.926, -18.927, -18.928,
    -18.928, -18.93, -18.93, -18.932, -18.932, -18.932, -18.933, -18.934,
    -18.935, -18.935, -18.936, -18.936, -18.937, -18.937, -18.939, -18.939,
    -18.94, -18.941, -18.941, -18.941, -18.941, -18.942, -18.943, -18.944,
    -18.945, -18.944, -18.944, -18.944, -18.944, -18.943, -18.944, -18.944,
    -18.945, -18.947, -18.947, -18.947, -18.948, -18.949, -18.949, -18.949,
    -18.95, -18.95, -18.951, -18.95, -18.95, -18.951, -18.952, -18.952,
    -18.952, -18.952, -18.953, -18.952, -18.952, -18.951, -18.952, -18.952,
    -18.953, -18.952, -18.952, -18.952, -18.952, -18.951, -18.951, -18.951,
    -18.951, -18.951, -18.95, -18.949, -18.949, -18.948, -18.947, -18.947,
    -18.948, -18.947, -18.947, -18.947, -18.946, -18.946, -18.945, -18.944,
    -18.944, -18.944, -18.942, -18.942, -18.942, -18.942, -18.942, -18.941,
    -18.94, -18.939, -18.938, -18.937, -18.936, -18.935, -18.934, -18.933,
    -18.932, -18.93, -18.929, -18.928, -18.927, -18.928, -18.926, -18.926,
    -18.926, -18.924, -18.923, -18.921, -18.919, -18.918, -18.917, -18.916,
    -18.914, -18.913, -18.911, -18.909, -18.907, -18.905, -18.903, -18.901,
    -18.899, -18.898, -18.896, -18.895, -18.892, -18.89, -18.888, -18.886,
    -18.884, -18.88, -18.878, -18.875, -18.873, -18.87, -18.867, -18.865,
    -18.862, -18.86, -18.857, -18.855, -18.852, -18.848, -18.844, -18.842,
    -18.838, -18.835, -18.832, -18.828, -18.824, -18.821, -18.818, -18.814,
    -18.81, -18.806, -18.802, -18.8, -18.796, -18.792, -18.787, -18.783,
    -18.779, -18.775, -18.771, -18.767, -18.762, -18.757, -18.753, -18.748,
    -18.744, -18.739, -18.734, -18.729, -18.724, -18.72, -18.714, -18.709,
    -18.704, -18.699, -18.693, -18.688, -18.683, -18.677, -18.671, -18.666,
    -18.66, -18.654, -18.648, -18.643, -18.637, -18.631, -18.625, -18.619,
    -18.612, -18.606, -18.6, -18.594, -18.587, -18.581, -18.574, -18.567,
    -18.56, -18.553, -18.546, -18.539, -18.532, -18.525, -18.518, -18.51,
    -18.503, -18.495, -18.488, -18.48, -18.472, -18.465, -18.457, -18.449,
    -18.441, -18.433, -18.424, -18.416, -18.408, -18.399, -18.391, -18.382,
    -18.374, -18.365, -18.356, -18.348, -18.339, -18.33, -18.321, -18.312,
    -18.303, -18.294, -18.284, -18.275, -18.266, -18.256, -18.247, -18.237,
    -18.227, -18.218, -18.208, -18.198, -18.188, -18.178, -18.168, -18.158,
    -18.148, -18.138, -18.128, -18.117, -18.107, -18.097, -18.086, -18.076,
    -18.065, -18.055, -18.044, -18.033, -18.023, -18.012, -18.001, -17.99,
    -17.979, -17.968, -17.957, -17.946, -17.935, -17.924, -17.912, -17.901,
    -17.89, -17.879, -17.867, -17.856, -17.844, -17.833, -17.821, -17.81,
    -17.798, -17.786, -17.775, -17.763, -17.751, -17.739, -17.728, -17.716,
    -17.704, -17.692, -17.68, -17.668, -17.656, -17.644, -17.632, -17.619,
    -17.607, -17.595, -17.583, -17.571, -17.558, -17.546, -17.534, -17.521,
    -17.509, -17.496, -17.484, -17.471, -17.459, -17.446, -17.433, -17.421,
    -17.408, -17.395, -17.383, -17.37, -17.357, -17.344, -17.331, -17.319,
    -17.306, -17.293, -17.28, -17.267, -17.254, -17.241, -17.228, -17.215,
    -17.201, -17.188, -17.175, -17.162, -17.149, -17.135, -17.122, -17.109,
    -17.095, -17.082, -17.069, -17.055, -17.042, -17.028, -17.015, -17.001,
    -16.987, -16.974, -16.96, -16.946, -16.933, -16.919, -16.905, -16.891,
    -16.878, -16.864, -16.85, -16.836, -16.822, -16.808, -16.794, -16.78,
    -16.766, -16.752, -16.738, -16.723, -16.709, -16.695, -16.681, -16.666,
    -16.652, -16.638, -16.623, -16.609, -16.594, -16.58, -16.565, -16.55,
    -16.536, -16.521, -16.506, -16.492, -16.477, -16.462, -16.447, -16.432,
    -16.417, -16.402, -16.387, -16.372, -16.357, -16.341, -16.326, -16.311,
    -16.296, -16.28, -16.265, -16.249, -16.234, -16.218, -16.202, -16.187,
    -16.171, -16.155, -16.139, -16.123, -16.107, -16.091, -16.075, -16.059,
    -16.043, -16.026, -16.01, -15.993, -15.984};

float log10_of_bbp440_over_No_sigma_std[449] = {0.422, 0.422, 0.422, 0.421,
    0.421, 0.421, 0.42, 0.42, 0.42, 0.421, 0.421, 0.421, 0.421, 0.421, 0.42,
    0.421, 0.421, 0.421, 0.421, 0.421, 0.422, 0.421, 0.421, 0.42, 0.42,
    0.42, 0.42, 0.42, 0.42, 0.42, 0.42, 0.421, 0.421, 0.421, 0.421, 0.421,
    0.421, 0.422, 0.422, 0.422, 0.423, 0.423, 0.424, 0.423, 0.423, 0.423,
    0.423, 0.423, 0.423, 0.423, 0.423, 0.424, 0.424, 0.424, 0.424, 0.425,
    0.426, 0.426, 0.427, 0.427, 0.428, 0.428, 0.428, 0.428, 0.428, 0.429,
    0.429, 0.429, 0.429, 0.43, 0.43, 0.431, 0.431, 0.432, 0.432, 0.432,
    0.433, 0.433, 0.434, 0.435, 0.435, 0.435, 0.436, 0.437, 0.438, 0.438,
    0.439, 0.439, 0.439, 0.439, 0.44, 0.441, 0.442, 0.442, 0.443, 0.443,
    0.444, 0.444, 0.444, 0.445, 0.445, 0.446, 0.447, 0.447, 0.447, 0.447,
    0.448, 0.448, 0.45, 0.45, 0.451, 0.451, 0.452, 0.453, 0.453, 0.453,
    0.453, 0.454, 0.455, 0.455, 0.456, 0.457, 0.458, 0.458, 0.459, 0.46,
    0.461, 0.461, 0.461, 0.462, 0.462, 0.462, 0.463, 0.463, 0.464, 0.464,
    0.465, 0.466, 0.466, 0.467, 0.468, 0.468, 0.469, 0.471, 0.471, 0.471,
    0.471, 0.472, 0.472, 0.472, 0.473, 0.473, 0.473, 0.473, 0.473, 0.473,
    0.473, 0.474, 0.474, 0.474, 0.474, 0.475, 0.475, 0.475, 0.475, 0.475,
    0.475, 0.475, 0.475, 0.475, 0.475, 0.475, 0.475, 0.475, 0.475, 0.475,
    0.475, 0.475, 0.475, 0.475, 0.474, 0.474, 0.474, 0.474, 0.473, 0.473,
    0.472, 0.472, 0.471, 0.471, 0.47, 0.471, 0.47, 0.47, 0.469, 0.468,
    0.468, 0.468, 0.467, 0.467, 0.466, 0.465, 0.465, 0.464, 0.463, 0.463,
    0.462, 0.461, 0.461, 0.46, 0.46, 0.459, 0.458, 0.458, 0.457, 0.456,
    0.455, 0.455, 0.454, 0.453, 0.453, 0.452, 0.451, 0.451, 0.451, 0.45,
    0.449, 0.448, 0.447, 0.447, 0.446, 0.446, 0.445, 0.444, 0.444, 0.443,
    0.442, 0.441, 0.441, 0.44, 0.439, 0.438, 0.438, 0.437, 0.436, 0.435,
    0.435, 0.434, 0.433, 0.432, 0.432, 0.431, 0.43, 0.429, 0.429, 0.428,
    0.427, 0.426, 0.426, 0.425, 0.424, 0.424, 0.423, 0.422, 0.421, 0.421,
    0.42, 0.419, 0.419, 0.418, 0.417, 0.417, 0.416, 0.415, 0.415, 0.414,
    0.414, 0.413, 0.412, 0.412, 0.411, 0.411, 0.41, 0.409, 0.409, 0.408,
    0.408, 0.407, 0.407, 0.406, 0.406, 0.405, 0.405, 0.404, 0.404, 0.403,
    0.403, 0.402, 0.402, 0.402, 0.401, 0.401, 0.4, 0.4, 0.4, 0.399, 0.399,
    0.398, 0.398, 0.398, 0.397, 0.397, 0.397, 0.396, 0.396, 0.396, 0.395,
    0.395, 0.395, 0.394, 0.394, 0.394, 0.394, 0.393, 0.393, 0.393, 0.393,
    0.392, 0.392, 0.392, 0.392, 0.391, 0.391, 0.391, 0.391, 0.391, 0.39,
    0.39, 0.39, 0.39, 0.39, 0.39, 0.389, 0.389, 0.389, 0.389, 0.389, 0.389,
    0.388, 0.388, 0.388, 0.388, 0.388, 0.388, 0.388, 0.388, 0.387, 0.387,
    0.387, 0.387, 0.387, 0.387, 0.387, 0.387, 0.387, 0.386, 0.386, 0.386,
    0.386, 0.386, 0.386, 0.386, 0.386, 0.386, 0.386, 0.386, 0.386, 0.386,
    0.385, 0.385, 0.385, 0.385, 0.385, 0.385, 0.385, 0.385, 0.385, 0.385,
    0.385, 0.385, 0.385, 0.385, 0.385, 0.385, 0.385, 0.385, 0.385, 0.384,
    0.384, 0.384, 0.384, 0.384, 0.384, 0.384, 0.384, 0.384, 0.384, 0.384,
    0.384, 0.384, 0.384, 0.384, 0.384, 0.384, 0.384, 0.384, 0.384, 0.384,
    0.384, 0.384, 0.384, 0.384, 0.384, 0.384, 0.384, 0.384, 0.384, 0.384,
    0.384, 0.384, 0.384, 0.384, 0.384, 0.384, 0.384, 0.384, 0.384, 0.384,
    0.384, 0.384, 0.384, 0.384, 0.384, 0.384, 0.384, 0.384, 0.384, 0.384,
    0.377};

/**********************************************************************************************************************************************************/
void calc_psd_ksm(float eta, float bbp_443, int ip, float *abundance_micro_ksm,
        float *abundance_nano_ksm, float *abundance_pico_ksm,
        float *volume_micro_ksm, float *volume_nano_ksm, float *volume_pico_ksm,
        float *ratio_micro_ksm, float *ratio_nano_ksm, float *ratio_pico_ksm) {

    /* This bio-optical algorithm computes the particle size distribution "PSD" via the particle backscattering coefficient spectrum "bbp(lambda)",
     * which provides important biogeochemical and ecological information about the pelagic ocean ecosystem structure and function.

     The PSD allows us to determine the abundance and volume of differing size classes of phytoplankton in the oceans.
     N_micro_KSM= the total particle abundance "N" in the size distribution of microplankton (20-50 um "microns" in diameter).
     N_nano_KSM= the particle size distribution of nanoplankton (2-20 um "microns" in diameter).
     N_pico_KSM= the particle size distribution of picoplankton (0.5-2 um "microns" in diameter).
     V_micro_KSM= the particle size distribution of microplankton (20-50 um "microns" in diameter).
     V_nano_KSM= the particle size distribution of nanoplankton (2-20 um "microns" in diameter).
     V_pico_KSM= the particle size distribution of picoplankton (0.5-2 um "microns" in diameter).
     */

    // Use the retrieved "eta" value from las_iop.c for determining the index of the look up tables "LUTs"
    int LUT_index = 0;

    while (eta > bbp_slope[LUT_index]) {
        LUT_index += 1;
    }

    /**********************************************************************************************************************************************************/
    //Using the LUT_index to determine the location of the lower and upper values bounding "eta".
    float bbp_slope_low, bbp_slope_high;

    bbp_slope_low = bbp_slope[LUT_index];
    bbp_slope_high = bbp_slope[LUT_index + 1];

    //Linearly interpolate the fractional location of "eta" between the bbp_slope LUT lower and upper bounding values.
    float bbp_over_no, f;
    float bbp_over_no_low, bbp_over_no_high;

    bbp_over_no_low = log10_of_bbp440_over_No[LUT_index];
    bbp_over_no_high = log10_of_bbp440_over_No[LUT_index + 1];

    f = (eta - bbp_slope_low) / (bbp_slope_high - bbp_slope_low);
    bbp_over_no = ((1 - f) * bbp_over_no_low) + (f * bbp_over_no_high);

    //set the lower and upper bound of the PSD_slope using the LUT_index and the PSD_slope LUT. Interpolate using the previously calculated "f" to determine "xi".
    float PSD_slope_low, PSD_slope_high, xi;

    PSD_slope_low = PSD_slope[LUT_index];
    PSD_slope_high = PSD_slope[LUT_index + 1];

    xi = ((1 - f) * PSD_slope_low) + (f * PSD_slope_high);

    /**********************************************************************************************************************************************************/
    // No is the reference particle number concentration at (2um) , No [units of M^-4].
    float No_log10;


    No_log10 = log10(bbp_443) - (bbp_over_no);

    float No;

    No = pow(10, No_log10);
    /**********************************************************************************************************************************************************/

    // CALUCULATING ABUNDANCE "N" AND VOLUME "V" using the PSD slope "xi".
    //With the above variables computed it is now possible to compute a variety of relevant biogeochemical and ecological parameters using the derived PSD slope, "xi", and reference abundance, N0.
    //Here, particle abundance and volume in the size ranges corresponding to picoplankton,nanoplankton, and microplankton can be calculated. See the size classes immediately below:
    // micro= the particle size distribution of microplankton which are 20-50 um "microns" in diameter.
    // nano= the particle size distribution of nanoplankton which are 2-20 um "microns" in diameter.
    // pico= the particle size distribution of picoplankton which are 0.5-2 um "microns" in diameter.
    //The total particle abundance "N", [particles/m^3], in a given size range Dmin to Dmax [m] is calculated using the equation below:
    float Do, micro_Dmin, micro_Dmax, nano_Dmin, nano_Dmax, pico_Dmin,
            pico_Dmax;

    //Reference diameter (2 um) in meters.
    Do = 2 * pow(10, -6);

    //Taken from Table 1. The diameter size range of the phytoplankton types.
    micro_Dmin = 2 * pow(10, -5);
    micro_Dmax = 5 * pow(10, -5);
    nano_Dmin = 2 * pow(10, -6);
    nano_Dmax = 2 * pow(10, -5);
    pico_Dmin = 5 * pow(10, -7);
    pico_Dmax = 2 * pow(10, -6);

    //N= integration(Dmin,Dmax) No (D/Do)^-xi dD=1/(1-xi)NoDo^xi(Dmax^(1-xi)-Dmin^(1-xi));
    //float abundance_micro_ksm,abundance_nano_ksm,abundance_pico_ksm;

    //N=1/(1-xi)NoDo^xi(Dmax^(1-xi)-Dmin^(1-xi));
    *abundance_micro_ksm = (1 / (1 - xi)) * No * pow(Do, xi)
            * (pow(micro_Dmax, 1 - xi) - pow(micro_Dmin, 1 - xi));
    *abundance_nano_ksm = (1 / (1 - xi)) * No * pow(Do, xi)
            * (pow(nano_Dmax, 1 - xi) - pow(nano_Dmin, 1 - xi));
    *abundance_pico_ksm = (1 / (1 - xi)) * No * pow(Do, xi)
            * (pow(pico_Dmax, 1 - xi) - pow(pico_Dmin, 1 - xi));

    /**********************************************************************************************************************************************************/

    //Below the particle volume "V" is computed.
    //However when xi = 4, the integral used is not defined and V is instead calculated as :
    //float volume_micro_ksm,volume_nano_ksm,volume_pico_ksm;
    if (xi == 4) {
        //V=integration(Dmin^Dmax) (PIE/6)D^3No(D/Do)^-xi dD =(PIE/6)NoDo^xi(Dmax/Dmin)

        //V= (OEL_PI/6)No*Do^xi(Dmax/Dmin);
        *volume_micro_ksm = (OEL_PI / 6) * No * pow(Do, 4)
                * log(micro_Dmax / micro_Dmin);
        *volume_nano_ksm = (OEL_PI / 6) * No * pow(Do, 4)
                * log(nano_Dmax / nano_Dmin);
        *volume_pico_ksm = (OEL_PI / 6) * No * pow(Do, 4)
                * log(pico_Dmax / pico_Dmin);
    }// Otherwise it is calculated as:

    else {
        //V=integration(Dmin^Dmax) (PIE/6)D^3No(D/Do)^-xi dD =(1/(4-xi))(PIE/6)NoDo^xi(Dmax^(4-xi)-Dmin^(4-xi))
        double piOverSix = OEL_PI / 6;
        //V= (1/(4-xi))(OEL_PI/6)No*Do^xi(Dmax^(4-xi)-Dmin^(4-xi));
        *volume_micro_ksm = (1 / (4 - xi)) * piOverSix * No * pow(Do, xi)
                * (pow(micro_Dmax, 4 - xi) - pow(micro_Dmin, 4 - xi));
        *volume_nano_ksm = (1 / (4 - xi)) * piOverSix * No * pow(Do, xi)
                * (pow(nano_Dmax, 4 - xi) - pow(nano_Dmin, 4 - xi));
        *volume_pico_ksm = (1 / (4 - xi)) * piOverSix * No * pow(Do, xi)
                * (pow(pico_Dmax, 4 - xi) - pow(pico_Dmin, 4 - xi));
    }
    /**********************************************************************************************************************************************************/

    //Calculating the ratio of micro,nano,pico vs. volume. Will need to relate to sea water volume as well...!!!!!!!!
    float V_ratio_micro_KSM, V_ratio_nano_KSM, V_ratio_pico_KSM;
    //float *ratio_micro_ksm,*ratio_nano_ksm,*ratio_pico_ksm;
    float total_volume_KSM;

    total_volume_KSM = *volume_micro_ksm + *volume_nano_ksm + *volume_pico_ksm;

    V_ratio_micro_KSM = *volume_micro_ksm / total_volume_KSM;
    V_ratio_nano_KSM = *volume_nano_ksm / total_volume_KSM;
    V_ratio_pico_KSM = *volume_pico_ksm / total_volume_KSM;

    // If the ratio is computed at greater than 1. Set upper boundary to 1.
    if (V_ratio_micro_KSM > 1) {
        *ratio_micro_ksm = 1;
    } else {
        *ratio_micro_ksm = V_ratio_micro_KSM;
    }
    if (V_ratio_nano_KSM > 1) {
        *ratio_nano_ksm = 1;
    } else {
        *ratio_nano_ksm = V_ratio_nano_KSM;
    }
    if (V_ratio_pico_KSM > 1) {
        *ratio_pico_ksm = 1;
    } else {
        *ratio_pico_ksm = V_ratio_pico_KSM;
    }
    // If the ratio is computed at less than 1. Set lower boundary to 0.
    if (V_ratio_micro_KSM < 0) {
        *ratio_micro_ksm = 0;
    } else {
        *ratio_micro_ksm = V_ratio_micro_KSM;
    }
    if (V_ratio_nano_KSM < 0) {
        *ratio_nano_ksm = 0;
    } else {
        *ratio_nano_ksm = V_ratio_nano_KSM;
    }
    if (V_ratio_pico_KSM < 0) {
        *ratio_pico_ksm = 0;
    } else {
        *ratio_pico_ksm = V_ratio_pico_KSM;
    }
}

/**********************************************************************************************************************************************************/

void get_psd_ksm(l2str *l2rec, l2prodstr *p, float prod[]) {
    int ip;
    float abundance_micro_ksm, abundance_nano_ksm, abundance_pico_ksm,
            volume_micro_ksm, volume_nano_ksm, volume_pico_ksm, ratio_micro_ksm,
            ratio_nano_ksm, ratio_pico_ksm;

    // The slope of the particle back scatter "bbp", "eta" is calculated using a log-transform and a linear regression from the Loisel et al. 2005 algorithm (las_iop.c).
    // The computed eta value is then applied to the look up tables "LUT" included above.

    // Use the retrieved "eta" value from las_iop.c for determining the index of the look up tables "LUTs"

    // the particle backscatter at 443 nm "bbp_440" is returned from the Loisel et al. 2005 algorithm (las_iop.c).

    for (ip = 0; ip < l2rec->l1rec->npix; ip++) {


        // The slope of the particle back scatter "bbp", "eta" is calculated using a log-transform and a linear regression from the Loisel et al. 2005 algorithm (las_iop.c).
        // The computed eta value is then applied to the look up tables "LUT" included above.
        float eta;

        //old "full spectrum" calculated eta
        eta = get_bbp_las_eta(l2rec, ip);

        //adding different eta calculations for different wavelengths
        //        if (input->ksd_opt=0)
        //            eta = get_bbp_las_eta(l2rec, ip);
        //        else

        //interface with las_iop.c which calculates eta based on ONLY 490,510,550 nm wavelengths only for SeaWiFS
        //eta = get_bbp_las_eta_ksm(l2rec, ip);


        if (eta == -32767) {
            prod[ip] = BAD_FLT;
            continue;
        }
        // If the returned eta is outside of the range of the Kostadinov algorithm's scope then this algorithm program quits and returns a BAD_FLT for the outputted PSD parameters.

        if (eta < -1.5 || eta > 2.8) {
            prod[ip] = BAD_FLT;
            continue;

        }
        // Use the retrieved "eta" value from las_iop.c for determining the index of the look up tables "LUTs"

        // the particle backscatter at 443 nm "bbp_440" is returned from the Loisel et al. 2005 algorithm (las_iop.c).
        //bbp_440=get_bbp_las(l2str *l2rec, int ip, float tab_wave[], float tab_bbp[], int tab_nwave);

        float wave = 443;
        float bbp_443;

        if (!get_bbp_las(l2rec, ip, &wave, &bbp_443, 1)) {
            prod[ip] = BAD_FLT;
            continue;

        }

        //Values of bbp(443)  outside the range (.00001 to .1) are rejected according to algorithm methodology
        //if (bbp_443 < .00001 || bbp_443 > .1) {
        //prod[ip] = BAD_FLT;
        //continue;
        //}

        /**********************************************************************************************************************************************************/


        switch (p->cat_ix) {
        case CAT_microplankton_abundanceksm:
            calc_psd_ksm(eta, bbp_443, ip, &prod[ip], &abundance_nano_ksm,
                    &abundance_pico_ksm, &volume_micro_ksm, &volume_nano_ksm,
                    &volume_pico_ksm, &ratio_micro_ksm, &ratio_nano_ksm,
                    &ratio_pico_ksm);
            break;
        case CAT_nanoplankton_abundanceksm:
            calc_psd_ksm(eta, bbp_443, ip, &abundance_micro_ksm, &prod[ip],
                    &abundance_pico_ksm, &volume_micro_ksm, &volume_nano_ksm,
                    &volume_pico_ksm, &ratio_micro_ksm, &ratio_nano_ksm,
                    &ratio_pico_ksm);
            break;
        case CAT_picoplankton_abundanceksm:
            calc_psd_ksm(eta, bbp_443, ip, &abundance_micro_ksm,
                    &abundance_nano_ksm, &prod[ip], &volume_micro_ksm,
                    &volume_nano_ksm, &volume_pico_ksm, &ratio_micro_ksm,
                    &ratio_nano_ksm, &ratio_pico_ksm);
            break;
        case CAT_microplankton_volumeksm:
            calc_psd_ksm(eta, bbp_443, ip, &abundance_micro_ksm,
                    &abundance_nano_ksm, &abundance_pico_ksm, &prod[ip],
                    &volume_nano_ksm, &volume_pico_ksm, &ratio_micro_ksm,
                    &ratio_nano_ksm, &ratio_pico_ksm);
            break;
        case CAT_nanoplankton_volumeksm:
            calc_psd_ksm(eta, bbp_443, ip, &abundance_micro_ksm,
                    &abundance_nano_ksm, &abundance_pico_ksm, &volume_micro_ksm,
                    &prod[ip], &volume_pico_ksm, &ratio_micro_ksm,
                    &ratio_nano_ksm, &ratio_pico_ksm);
            break;
        case CAT_picoplankton_volumeksm:
            calc_psd_ksm(eta, bbp_443, ip, &abundance_micro_ksm,
                    &abundance_nano_ksm, &abundance_pico_ksm, &volume_micro_ksm,
                    &volume_nano_ksm, &prod[ip], &ratio_micro_ksm,
                    &ratio_nano_ksm, &ratio_pico_ksm);
            break;
        case CAT_microplankton_ratioksm:
            calc_psd_ksm(eta, bbp_443, ip, &abundance_micro_ksm,
                    &abundance_nano_ksm, &abundance_pico_ksm, &volume_micro_ksm,
                    &volume_nano_ksm, &volume_pico_ksm, &prod[ip],
                    &ratio_nano_ksm, &ratio_pico_ksm);
            break;
        case CAT_nanoplankton_ratioksm:
            calc_psd_ksm(eta, bbp_443, ip, &abundance_micro_ksm,
                    &abundance_nano_ksm, &abundance_pico_ksm, &volume_micro_ksm,
                    &volume_nano_ksm, &volume_pico_ksm, &ratio_micro_ksm,
                    &prod[ip], &ratio_pico_ksm);
            break;
        case CAT_picoplankton_ratioksm:
            calc_psd_ksm(eta, bbp_443, ip, &abundance_micro_ksm,
                    &abundance_nano_ksm, &abundance_pico_ksm, &volume_micro_ksm,
                    &volume_nano_ksm, &volume_pico_ksm, &ratio_micro_ksm,
                    &ratio_nano_ksm, &prod[ip]);
            break;

        default:
            printf("get_psd_ksm.c can not produce product %s\n",
                    p->algorithm_id);
            exit(1);
        }
        if (isnan(prod[ip])) {
            prod[ip] = BAD_FLT;
            l2rec->l1rec->flags[ip] |= PRODFAIL;
        }
    }
}

