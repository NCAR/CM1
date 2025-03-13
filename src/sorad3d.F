  MODULE sorad3d_module

  implicit none

  private
  public :: sorad

  integer :: i,j

!-----------------------------------------------------------------------
!
!  The following DATA statements originally came from file
!  "cai.dat", which define pre-computed tables cai
!
!-----------------------------------------------------------------------
!
!  integer   nm,nt,na
!  parameter (nm=11,nt=9,na=11)

  REAL :: caib(11,9,11),caif(9,11)

!!!  COMMON /radtab004/ caib,caif

  DATA ((caib(1,i,j),j=1,11),i=1,9)/                                    &
      0.,0.084,0.170,0.258,0.349,0.442,0.540,0.643,0.752,0.870,1.,      &
      0.,0.076,0.153,0.233,0.315,0.401,0.492,0.589,0.698,0.826,1.,      &
      0.,0.069,0.139,0.210,0.283,0.359,0.438,0.523,0.619,0.741,1.,      &
      0.,0.065,0.130,0.196,0.263,0.332,0.403,0.478,0.561,0.671,1.,      &
      0.,0.062,0.126,0.189,0.253,0.318,0.385,0.454,0.533,0.681,1.,      &
      0.,0.062,0.123,0.185,0.248,0.311,0.375,0.444,0.539,0.738,1.,      &
      0.,0.061,0.122,0.183,0.244,0.307,0.370,0.444,0.577,0.752,1.,      &
      0.,0.060,0.121,0.182,0.243,0.304,0.368,0.457,0.574,0.732,1.,      &
      0.,0.060,0.120,0.181,0.242,0.302,0.367,0.451,0.551,0.690,1./
  DATA ((caib(2,i,j),j=1,11),i=1,9)/                                    &
      0.,0.094,0.190,0.286,0.384,0.482,0.583,0.684,0.788,0.893,1.,      &
      0.,0.090,0.180,0.274,0.368,0.466,0.566,0.668,0.774,0.885,1.,      &
      0.,0.082,0.167,0.254,0.343,0.436,0.533,0.636,0.746,0.866,1.,      &
      0.,0.074,0.151,0.229,0.310,0.395,0.485,0.582,0.692,0.822,1.,      &
      0.,0.068,0.138,0.209,0.282,0.358,0.439,0.527,0.630,0.769,1.,      &
      0.,0.065,0.130,0.197,0.265,0.336,0.410,0.493,0.599,0.766,1.,      &
      0.,0.063,0.126,0.190,0.256,0.323,0.394,0.480,0.604,0.772,1.,      &
      0.,0.062,0.124,0.186,0.249,0.314,0.385,0.475,0.594,0.749,1.,      &
      0.,0.061,0.122,0.183,0.245,0.309,0.378,0.462,0.564,0.704,1./
  DATA ((caib(3,i,j),j=1,11),i=1,9)/                                    &
      0.,0.097,0.194,0.292,0.392,0.491,0.591,0.692,0.794,0.897,1.,      &
      0.,0.094,0.189,0.286,0.383,0.482,0.582,0.684,0.788,0.893,1.,      &
      0.,0.090,0.180,0.274,0.368,0.466,0.566,0.668,0.775,0.885,1.,      &
      0.,0.083,0.167,0.254,0.344,0.438,0.535,0.638,0.748,0.868,1.,      &
      0.,0.075,0.152,0.232,0.314,0.400,0.492,0.592,0.703,0.835,1.,      &
      0.,0.069,0.140,0.212,0.288,0.366,0.451,0.546,0.659,0.805,1.,      &
      0.,0.066,0.132,0.200,0.270,0.344,0.424,0.518,0.638,0.794,1.,      &
      0.,0.063,0.127,0.192,0.259,0.329,0.406,0.499,0.616,0.769,1.,      &
      0.,0.062,0.124,0.187,0.251,0.318,0.391,0.477,0.580,0.721,1./
  DATA ((caib(4,i,j),j=1,11),i=1,9)/                                    &
      0.,0.098,0.196,0.295,0.394,0.494,0.594,0.695,0.796,0.898,1.,      &
      0.,0.096,0.193,0.290,0.389,0.488,0.589,0.690,0.792,0.896,1.,      &
      0.,0.093,0.187,0.282,0.379,0.478,0.578,0.680,0.784,0.891,1.,      &
      0.,0.087,0.176,0.268,0.362,0.458,0.557,0.661,0.768,0.881,1.,      &
      0.,0.080,0.163,0.248,0.336,0.428,0.524,0.628,0.739,0.862,1.,      &
      0.,0.074,0.149,0.226,0.308,0.393,0.485,0.586,0.700,0.835,1.,      &
      0.,0.068,0.138,0.210,0.285,0.365,0.452,0.550,0.668,0.814,1.,      &
      0.,0.065,0.131,0.199,0.269,0.344,0.426,0.522,0.638,0.786,1.,      &
      0.,0.063,0.126,0.191,0.257,0.327,0.404,0.491,0.596,0.737,1./
  DATA ((caib(5,i,j),j=1,11),i=1,9)/                                    &
      0.,0.098,0.197,0.296,0.396,0.496,0.596,0.696,0.797,0.898,1.,      &
      0.,0.097,0.195,0.293,0.392,0.492,0.592,0.693,0.794,0.897,1.,      &
      0.,0.095,0.190,0.287,0.385,0.484,0.584,0.686,0.789,0.894,1.,      &
      0.,0.090,0.182,0.276,0.372,0.470,0.570,0.672,0.778,0.887,1.,      &
      0.,0.084,0.171,0.259,0.351,0.446,0.545,0.649,0.758,0.875,1.,      &
      0.,0.077,0.157,0.239,0.325,0.415,0.511,0.614,0.728,0.855,1.,      &
      0.,0.071,0.145,0.220,0.300,0.384,0.476,0.578,0.694,0.832,1.,      &
      0.,0.067,0.135,0.206,0.280,0.358,0.444,0.542,0.658,0.802,1.,      &
      0.,0.064,0.129,0.195,0.264,0.337,0.416,0.504,0.610,0.751,1./
  DATA ((caib(6,i,j),j=1,11),i=1,9)/                                    &
      0.,0.099,0.198,0.298,0.397,0.497,0.597,0.698,0.798,0.899,1.,      &
      0.,0.098,0.196,0.295,0.394,0.494,0.594,0.695,0.796,0.898,1.,      &
      0.,0.096,0.193,0.290,0.389,0.488,0.588,0.690,0.792,0.896,1.,      &
      0.,0.093,0.187,0.282,0.379,0.478,0.578,0.680,0.784,0.891,1.,      &
      0.,0.088,0.177,0.268,0.362,0.460,0.559,0.663,0.770,0.882,1.,      &
      0.,0.081,0.164,0.250,0.340,0.433,0.531,0.636,0.747,0.868,1.,      &
      0.,0.074,0.151,0.231,0.314,0.402,0.497,0.600,0.716,0.846,1.,      &
      0.,0.069,0.140,0.214,0.291,0.372,0.462,0.561,0.675,0.815,1.,      &
      0.,0.065,0.132,0.200,0.271,0.346,0.427,0.517,0.624,0.764,1./
  DATA ((caib(7,i,j),j=1,11),i=1,9)/                                    &
      0.,0.099,0.199,0.298,0.398,0.498,0.598,0.698,0.798,0.899,1.,      &
      0.,0.098,0.197,0.296,0.396,0.496,0.596,0.696,0.797,0.898,1.,      &
      0.,0.097,0.195,0.293,0.392,0.492,0.592,0.693,0.794,0.897,1.,      &
      0.,0.094,0.190,0.287,0.384,0.484,0.584,0.686,0.789,0.894,1.,      &
      0.,0.090,0.182,0.276,0.372,0.470,0.570,0.673,0.779,0.888,1.,      &
      0.,0.084,0.171,0.260,0.353,0.448,0.548,0.652,0.761,0.877,1.,      &
      0.,0.078,0.158,0.241,0.328,0.419,0.516,0.620,0.733,0.858,1.,      &
      0.,0.072,0.145,0.222,0.302,0.386,0.478,0.578,0.692,0.827,1.,      &
      0.,0.067,0.135,0.206,0.278,0.355,0.438,0.530,0.637,0.775,1./
  DATA ((caib(8,i,j),j=1,11),i=1,9)/                                    &
      0.,0.100,0.199,0.299,0.399,0.498,0.599,0.699,0.799,0.900,1.,      &
      0.,0.099,0.198,0.298,0.397,0.497,0.597,0.698,0.798,0.899,1.,      &
      0.,0.098,0.196,0.295,0.395,0.494,0.594,0.695,0.796,0.898,1.,      &
      0.,0.096,0.193,0.291,0.390,0.489,0.589,0.690,0.793,0.896,1.,      &
      0.,0.093,0.188,0.283,0.380,0.479,0.579,0.682,0.786,0.892,1.,      &
      0.,0.088,0.178,0.270,0.365,0.462,0.562,0.665,0.772,0.884,1.,      &
      0.,0.082,0.165,0.252,0.341,0.434,0.532,0.636,0.747,0.867,1.,      &
      0.,0.074,0.151,0.230,0.313,0.400,0.492,0.593,0.706,0.837,1.,      &
      0.,0.069,0.139,0.211,0.286,0.364,0.448,0.541,0.649,0.786,1./
  DATA ((caib(9,i,j),j=1,11),i=1,9)/                                    &
      0.,0.100,0.200,0.299,0.399,0.499,0.599,0.700,0.800,0.900,1.,      &
      0.,0.100,0.199,0.299,0.399,0.498,0.599,0.699,0.799,0.900,1.,      &
      0.,0.099,0.198,0.298,0.397,0.497,0.597,0.698,0.798,0.899,1.,      &
      0.,0.098,0.196,0.295,0.394,0.494,0.594,0.695,0.796,0.898,1.,      &
      0.,0.096,0.192,0.290,0.388,0.487,0.588,0.689,0.791,0.895,1.,      &
      0.,0.092,0.185,0.280,0.376,0.474,0.574,0.676,0.781,0.889,1.,      &
      0.,0.086,0.173,0.262,0.354,0.449,0.547,0.650,0.759,0.875,1.,      &
      0.,0.078,0.157,0.239,0.324,0.412,0.506,0.607,0.719,0.846,1.,      &
      0.,0.071,0.143,0.217,0.293,0.373,0.459,0.552,0.660,0.796,1./
  DATA ((caib(10,i,j),j=1,11),i=1,9)/                                   &
      0.,0.100,0.200,0.300,0.400,0.500,0.600,0.700,0.800,0.900,1.,      &
      0.,0.100,0.200,0.300,0.400,0.500,0.600,0.700,0.800,0.900,1.,      &
      0.,0.100,0.200,0.300,0.400,0.500,0.600,0.700,0.800,0.900,1.,      &
      0.,0.100,0.200,0.299,0.399,0.499,0.599,0.699,0.799,0.900,1.,      &
      0.,0.099,0.198,0.297,0.396,0.496,0.595,0.696,0.796,0.898,1.,      &
      0.,0.096,0.192,0.289,0.387,0.485,0.585,0.686,0.788,0.893,1.,      &
      0.,0.090,0.181,0.273,0.366,0.462,0.560,0.662,0.769,0.881,1.,      &
      0.,0.082,0.164,0.248,0.335,0.424,0.519,0.620,0.730,0.854,1.,      &
      0.,0.073,0.147,0.223,0.301,0.382,0.468,0.563,0.671,0.805,1./
  DATA ((caib(11,i,j),j=1,11),i=1,9)/                                   &
      0.,0.100,0.200,0.300,0.400,0.500,0.600,0.700,0.800,0.900,1.,      &
      0.,0.100,0.200,0.300,0.400,0.500,0.600,0.700,0.800,0.900,1.,      &
      0.,0.100,0.200,0.300,0.400,0.500,0.600,0.700,0.800,0.900,1.,      &
      0.,0.100,0.200,0.300,0.400,0.500,0.600,0.700,0.800,0.900,1.,      &
      0.,0.100,0.200,0.300,0.400,0.500,0.600,0.700,0.800,0.900,1.,      &
      0.,0.100,0.200,0.299,0.398,0.496,0.595,0.695,0.795,0.897,1.,      &
      0.,0.096,0.190,0.284,0.378,0.474,0.572,0.674,0.778,0.886,1.,      &
      0.,0.086,0.171,0.257,0.345,0.436,0.531,0.632,0.740,0.861,1.,      &
      0.,0.076,0.152,0.229,0.308,0.390,0.478,0.573,0.681,0.813,1./
!
  DATA ((caif(i,j),j=1,11),i=1,9)/                                      &
      0.,0.099,0.198,0.298,0.397,0.497,0.597,0.698,0.798,0.899,1.,      &
      0.,0.098,0.196,0.295,0.395,0.495,0.595,0.695,0.796,0.898,1.,      &
      0.,0.096,0.193,0.291,0.390,0.489,0.589,0.690,0.793,0.896,1.,      &
      0.,0.093,0.187,0.283,0.380,0.479,0.580,0.681,0.786,0.892,1.,      &
      0.,0.087,0.177,0.268,0.363,0.460,0.561,0.665,0.772,0.884,1.,      &
      0.,0.079,0.160,0.246,0.335,0.430,0.529,0.635,0.747,0.869,1.,      &
      0.,0.068,0.141,0.217,0.298,0.385,0.481,0.586,0.704,0.840,1.,      &
      0.,0.059,0.120,0.185,0.256,0.333,0.419,0.518,0.637,0.789,1.,      &
      0.,0.051,0.104,0.160,0.219,0.284,0.357,0.443,0.551,0.707,1./
!
!-----------------------------------------------------------------------
!
!  The following DATA statements originally came from file
!  "cah.dat", which define pre-computed tables used for co2.
!
!-----------------------------------------------------------------------
!
  REAL :: cah(22,19)

!!!  COMMON /radtab005/ cah

  DATA ((cah(i,j),i=1,22),j= 1, 5)/                                     &
      0.9923, 0.9922, 0.9921, 0.9920, 0.9916, 0.9910, 0.9899, 0.9882,   &
      0.9856, 0.9818, 0.9761, 0.9678, 0.9558, 0.9395, 0.9188, 0.8945,   &
      0.8675, 0.8376, 0.8029, 0.7621, 0.7154, 0.6647, 0.9876, 0.9876,   &
      0.9875, 0.9873, 0.9870, 0.9864, 0.9854, 0.9837, 0.9811, 0.9773,   &
      0.9718, 0.9636, 0.9518, 0.9358, 0.9153, 0.8913, 0.8647, 0.8350,   &
      0.8005, 0.7599, 0.7133, 0.6627, 0.9808, 0.9807, 0.9806, 0.9805,   &
      0.9802, 0.9796, 0.9786, 0.9769, 0.9744, 0.9707, 0.9653, 0.9573,   &
      0.9459, 0.9302, 0.9102, 0.8866, 0.8604, 0.8311, 0.7969, 0.7565,   &
      0.7101, 0.6596, 0.9708, 0.9708, 0.9707, 0.9705, 0.9702, 0.9697,   &
      0.9687, 0.9671, 0.9647, 0.9612, 0.9560, 0.9483, 0.9372, 0.9221,   &
      0.9027, 0.8798, 0.8542, 0.8253, 0.7916, 0.7515, 0.7054, 0.6551,   &
      0.9568, 0.9568, 0.9567, 0.9565, 0.9562, 0.9557, 0.9548, 0.9533,   &
      0.9510, 0.9477, 0.9428, 0.9355, 0.9250, 0.9106, 0.8921, 0.8700,   &
      0.8452, 0.8171, 0.7839, 0.7443, 0.6986, 0.6486/

  DATA ((cah(i,j),i=1,22),j= 6,10)/                                     &
      0.9377, 0.9377, 0.9376, 0.9375, 0.9372, 0.9367, 0.9359, 0.9345,   &
      0.9324, 0.9294, 0.9248, 0.9181, 0.9083, 0.8948, 0.8774, 0.8565,   &
      0.8328, 0.8055, 0.7731, 0.7342, 0.6890, 0.6395, 0.9126, 0.9126,   &
      0.9125, 0.9124, 0.9121, 0.9117, 0.9110, 0.9098, 0.9079, 0.9052,   &
      0.9012, 0.8951, 0.8862, 0.8739, 0.8579, 0.8385, 0.8161, 0.7900,   &
      0.7585, 0.7205, 0.6760, 0.6270, 0.8809, 0.8809, 0.8808, 0.8807,   &
      0.8805, 0.8802, 0.8796, 0.8786, 0.8770, 0.8747, 0.8712, 0.8659,   &
      0.8582, 0.8473, 0.8329, 0.8153, 0.7945, 0.7697, 0.7394, 0.7024,   &
      0.6588, 0.6105, 0.8427, 0.8427, 0.8427, 0.8426, 0.8424, 0.8422,   &
      0.8417, 0.8409, 0.8397, 0.8378, 0.8350, 0.8306, 0.8241, 0.8148,   &
      0.8023, 0.7866, 0.7676, 0.7444, 0.7154, 0.6796, 0.6370, 0.5897,   &
      0.7990, 0.7990, 0.7990, 0.7989, 0.7988, 0.7987, 0.7983, 0.7978,   &
      0.7969, 0.7955, 0.7933, 0.7899, 0.7846, 0.7769, 0.7664, 0.7528,   &
      0.7357, 0.7141, 0.6866, 0.6520, 0.6108, 0.5646/

  DATA ((cah(i,j),i=1,22),j=11,15)/                                     &
      0.7515, 0.7515, 0.7515, 0.7515, 0.7514, 0.7513, 0.7511, 0.7507,   &
      0.7501, 0.7491, 0.7476, 0.7450, 0.7409, 0.7347, 0.7261, 0.7144,   &
      0.6992, 0.6793, 0.6533, 0.6203, 0.5805, 0.5357, 0.7020, 0.7020,   &
      0.7020, 0.7019, 0.7019, 0.7018, 0.7017, 0.7015, 0.7011, 0.7005,   &
      0.6993, 0.6974, 0.6943, 0.6894, 0.6823, 0.6723, 0.6588, 0.6406,   &
      0.6161, 0.5847, 0.5466, 0.5034, 0.6518, 0.6518, 0.6518, 0.6518,   &
      0.6518, 0.6517, 0.6517, 0.6515, 0.6513, 0.6508, 0.6500, 0.6485,   &
      0.6459, 0.6419, 0.6359, 0.6273, 0.6151, 0.5983, 0.5755, 0.5458,   &
      0.5095, 0.4681, 0.6017, 0.6017, 0.6017, 0.6017, 0.6016, 0.6016,   &
      0.6016, 0.6015, 0.6013, 0.6009, 0.6002, 0.5989, 0.5967, 0.5932,   &
      0.5879, 0.5801, 0.5691, 0.5535, 0.5322, 0.5043, 0.4700, 0.4308,   &
      0.5518, 0.5518, 0.5518, 0.5518, 0.5518, 0.5518, 0.5517, 0.5516,   &
      0.5514, 0.5511, 0.5505, 0.5493, 0.5473, 0.5441, 0.5393, 0.5322,   &
      0.5220, 0.5076, 0.4878, 0.4617, 0.4297, 0.3929/

  DATA ((cah(i,j),i=1,22),j=16,19)/                                     &
      0.5031, 0.5031, 0.5031, 0.5031, 0.5031, 0.5030, 0.5030, 0.5029,   &
      0.5028, 0.5025, 0.5019, 0.5008, 0.4990, 0.4960, 0.4916, 0.4850,   &
      0.4757, 0.4624, 0.4441, 0.4201, 0.3904, 0.3564, 0.4565, 0.4565,   &
      0.4565, 0.4564, 0.4564, 0.4564, 0.4564, 0.4563, 0.4562, 0.4559,   &
      0.4553, 0.4544, 0.4527, 0.4500, 0.4460, 0.4400, 0.4315, 0.4194,   &
      0.4028, 0.3809, 0.3538, 0.3227, 0.4122, 0.4122, 0.4122, 0.4122,   &
      0.4122, 0.4122, 0.4122, 0.4121, 0.4120, 0.4117, 0.4112, 0.4104,   &
      0.4089, 0.4065, 0.4029, 0.3976, 0.3900, 0.3792, 0.3643, 0.3447,   &
      0.3203, 0.2923, 0.3696, 0.3696, 0.3696, 0.3696, 0.3696, 0.3696,   &
      0.3695, 0.3695, 0.3694, 0.3691, 0.3687, 0.3680, 0.3667, 0.3647,   &
      0.3615, 0.3570, 0.3504, 0.3409, 0.3279, 0.3106, 0.2892, 0.2642/
!
!-----------------------------------------------------------------------

  CONTAINS
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SORAD                      ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######               University of Oklahoma                 ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE sorad(idim,jdim,m,n,np,                                      &
           pl,ta,wa,oa,co2,                                             &
           taucldi,taucldl,reffi,reffl,fcld,ict,icb,                    &
           taual,rsirbm,rsirdf,rsuvbm,rsuvdf,cosz,                      &
           flx,flc,fdirir,fdifir,fdirpar,fdifpar,                       &
           fdiruv,fdifuv,                                               & ! MS add
           sdf,sclr,csm,cc,                                             &
           tauclb,tauclf,dp,wh,oh,scal,swh,so2,df,                      &
           tem2d1, tem2d2, tem2d3, tem2d4, tem2d5,                      &
           tem2d6, tem2d7, tem2d8, tem2d9, tem2d10,                     &
           tem2d11,tem2d12,tem2d13,tem2d14,tem2d15,                     &
           tem2d16,tem2d17,tem2d18,tem2d19,                             &
           tem3d1, tem3d2, tem3d3, tem3d4, tem3d5,                      &
           tem4d1, tem4d2, tem4d3, tem4d4, tem4d5,                      &
           tem5d1, tem5d2, tem5d3, tem5d4, tem5d5)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate solar fluxes due to the absoption by water vapor, ozone,
!  co2, o2, clouds, and aerosols and due to the scattering by clouds,
!  aerosols, and gases.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (a) Radiative Transfer Model: M.-D. Chou and M. Suarez
!          (b) Cloud Optics:Tao, Lang, Simpson, Sui, Ferrier and
!              Chou (1996)
!
!  MODIFICATION HISTORY:
!
!  03/15/1996 (Yuhe Liu)
!  Adopted the original code and formatted it in accordance with the
!  ARPS coding standard.
!
!-----------------------------------------------------------------------
!
!  ORIGINAL COMMENTS:
!
!****************  charney  /silo/z1mhy/new/solar.f  ***10/24/95******
!********************************************************************
!
! This routine computes solar fluxes due to the absoption by water
!  vapor, ozone, co2, o2, clouds, and aerosols and due to the
!  scattering by clouds, aerosols, and gases.
!
! This is a vectorized code. It computes the fluxes simultaneous for
!  (m x n) soundings, which is a subset of the (m x ndim) soundings.
!  In a global climate model, m and ndim correspond to the numbers of
!  grid boxes in the zonal and meridional directions, respectively.
!
! Ice and liquid cloud particles are allowed to co-exist in any of the
!  np layers. Two sets of cloud parameters are required as inputs, one
!  for ice paticles and the other for liquid particles.  These parameters
!  are optical thickness (taucld) and effective particle size (reff).
!
! If no information is available for reff, a default value of
!  10 micron for liquid water and 75 micron for ice can be used.
!
! Clouds are grouped into high, middle, and low clouds separated by the
!  level indices ict and icb.  For detail, see the subroutine cldscale.
!
!----- Input parameters:
!                                                units      size
!   number of soundings in zonal direction (m)       n/d        1
!   number of soundings in meridional direction (n)  n/d        1
!   maximum number of soundings in                   n/d        1
!        meridional direction (ndim)
!   number of atmospheric layers (np)                n/d        1
!   level pressure (pl)                              mb       m*ndim*(np+1)
!   layer temperature (ta)                           k        m*ndim*np
!   layer specific humidity (wa)                     gm/gm    m*ndim*np
!   layer ozone concentration (oa)                   gm/gm    m*ndim*np
!   co2 mixing ratio by volumn (co2)               parts/part   1
!   cloud optical thickness (taucld)                 n/d      m*ndim*np*2
!             index 1 for ice particles
!             index 2 for liquid drops
!   effective cloud-particle size (reff)           micrometer m*ndim*np*2
!             index 1 for ice particles
!             index 2 for liquid drops
!
!   ARPS notes: taucld were changed to taucldi and taucldl
!            reff   were changed to reffi   and reffl
!
!   cloud amount (fcld)                            fraction   m*ndim*np
!   level index separating high and middle           n/d        1
!             clouds (ict)
!   level index separating middle and low clouds     n/d        1
!             clouds (icb)
!   aerosol optical thickness (taual)                n/d      m*ndim*np
!   solar ir surface albedo for beam                fraction   m*ndim
!             radiation (rsirbm)
!   solar ir surface albedo for diffuse             fraction   m*ndim
!             radiation (rsirdf)
!   uv + par surface albedo for beam                     fraction   m*ndim
!             radiation (rsuvbm)
!   uv + par surface albedo for diffuse                  fraction   m*ndim
!             radiation (rsuvdf)
!   cosine of solar zenith angle (cosz)            n/d        m*ndim
!
!----- Output parameters
!
!   all-sky flux (downward minus upward) (flx)     fraction   m*ndim*(np+1)
!   clear-sky flux (downward minus upward) (flc)   fraction   m*ndim*(np+1)
!   all-sky direct downward ir (0.7-10 micron)
!             flux at the surface (fdirir)      fraction   m*ndim
!   all-sky diffuse downward ir flux at
!             the surface (fdifir)              fraction   m*ndim
!   all-sky direct downward par (0.4-0.7 micron)
!             flux at the surface (fdirpar)     fraction   m*ndim
!   all-sky diffuse downward par flux at
!             the surface (fdifpar)             fraction   m*ndim
!
! MS - add new variables to track the uv flux at the ground
!   all-sky direct downward uv (<0.4 micron)
!             flux at the surface (fdiruv)     fraction   m*ndim
!   all-sky diffuse downward par flux at
!             the surface (fdifuv)             fraction   m*ndim
!
!----- Notes:
!
!    (1) The unit of flux is fraction of the incoming solar radiation
!     at the top of the atmosphere.  Therefore, fluxes should
!     be equal to flux multiplied by the extra-terrestrial solar
!     flux and the cosine of solar zenith angle.
!    (2) Clouds and aerosols can be included in any layers by specifying
!     fcld(i,j,k), taucld(i,j,k,*) and taual(i,j,k), k=1,np.
!     For an atmosphere without clouds and aerosols,
!     set fcld(i,j,k)=taucld(i,j,k,*)=taual(i,j,k)=0.0.
!    (3) Aerosol single scattering albedos and asymmetry
!     factors are specified in the subroutines solir and soluv.
!    (4) pl(i,j,1) is the pressure at the top of the model, and
!     pl(i,j,np+1) is the surface pressure.
!
!    ARPS note: pl was replaced by pa at scalar points (layers)
!
!    (5) the pressure levels ict and icb correspond approximately
!     to 400 and 700 mb.
!
!**************************************************************************
!
!fpp$ expand (expmn)
!!dir$ inline always expmn
!*$*  inline routine (expmn)
!
!-----------------------------------------------------------------------
!
  use radlib3d_module, only : expmn
  IMPLICIT NONE

  INTEGER :: idim,jdim             ! ARPS nx, ny. NX is used in this
                                   ! subroutine for other purpose

  INTEGER :: m,n,np
  INTEGER :: ict,icb

  REAL :: pl(idim,jdim,np+1)       ! pressure at scalar points (layers)
  REAL :: ta(idim,jdim,np)
  REAL :: wa(idim,jdim,np)
  REAL :: oa(idim,jdim,np)
  REAL :: taucldi(idim,jdim,np)
  REAL :: taucldl(idim,jdim,np)
  REAL :: reffi(idim,jdim,np)
  REAL :: reffl(idim,jdim,np)
  REAL :: fcld(idim,jdim,np)
  REAL :: taual(idim,jdim,np)
  REAL :: rsirbm(idim,jdim)
  REAL :: rsirdf(idim,jdim)
  REAL :: rsuvbm(idim,jdim)
  REAL :: rsuvdf(idim,jdim)
  REAL :: cosz(idim,jdim)
  REAL :: co2

  REAL :: flx(idim,jdim,np+1)
  REAL :: flc(idim,jdim,np+1)
  REAL :: fdirir(idim,jdim)
  REAL :: fdifir(idim,jdim)
  REAL :: fdirpar(idim,jdim)
  REAL :: fdifpar(idim,jdim)
  REAL :: fdiruv(idim,jdim) ! MS add
  REAL :: fdifuv(idim,jdim) ! MS add
!
!-----------------------------------------------------------------------
!
!  Local temporary arrays
!
!-----------------------------------------------------------------------
!
  REAL :: sdf(m,n)
  REAL :: sclr(m,n)
  REAL :: csm(m,n)
  REAL :: cc(m,n,3)

  REAL :: tauclb(m,n,np)
  REAL :: tauclf(m,n,np)
  REAL :: dp(m,n,np)
  REAL :: wh(m,n,np)
  REAL :: oh(m,n,np)
  REAL :: scal(m,n,np)
  REAL :: swh(m,n,np+1)
  REAL :: so2(m,n,np+1)
  REAL :: df(m,n,np+1)

!-----temporary arrays used in subroutine cldflx

  REAL :: tem2d1 (m,n)
  REAL :: tem2d2 (m,n)
  REAL :: tem2d3 (m,n)
  REAL :: tem2d4 (m,n)
  REAL :: tem2d5 (m,n)
  REAL :: tem2d6 (m,n)
  REAL :: tem2d7 (m,n)
  REAL :: tem2d8 (m,n)
  REAL :: tem2d9 (m,n)
  REAL :: tem2d10(m,n)
  REAL :: tem2d11(m,n)
  REAL :: tem2d12(m,n)
  REAL :: tem2d13(m,n)
  REAL :: tem2d14(m,n)
  REAL :: tem2d15(m,n)
  REAL :: tem2d16(m,n)
  REAL :: tem2d17(m,n)
  REAL :: tem2d18(m,n)
  REAL :: tem2d19(m,n)

  REAL :: tem3d1(m,n,np)
  REAL :: tem3d2(m,n,np)
  REAL :: tem3d3(m,n,np+1)
  REAL :: tem3d4(m,n,np+1)
  REAL :: tem3d5(m,n,np+1)

  REAL :: tem4d1(m,n,np+1,2)
  REAL :: tem4d2(m,n,np+1,2)
  REAL :: tem4d3(m,n,np+1,2)
  REAL :: tem4d4(m,n,np+1,2)
  REAL :: tem4d5(m,n,np+1,2)

  REAL :: tem5d1(m,n,np+1,2,2)
  REAL :: tem5d2(m,n,np+1,2,2)
  REAL :: tem5d3(m,n,np+1,2,2)
  REAL :: tem5d4(m,n,np+1,2,2)
  REAL :: tem5d5(m,n,np+1,2,2)

!-----temporary variables

  INTEGER :: i,j,k
  REAL :: x
!
!-----------------------------------------------------------------------
!
!  Functions:
!
!-----------------------------------------------------------------------
!
!!!  REAL :: expmn
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO j=1,n
    DO i=1,m

      swh(i,j,1)=0.
      so2(i,j,1)=0.

!-----csm is the effective secant of the solar zenith angle
!  see equation (12) of Lacis and Hansen (1974, JAS)

      csm(i,j)=35./SQRT(1224.*cosz(i,j)*cosz(i,j)+1.)

    END DO
  END DO


  DO k=1,np
    DO j=1,n
      DO i=1,m

!-----compute layer thickness and pressure-scaling function.
!  indices for the surface level and surface layer
!  are np+1 and np, respectively.

        dp(i,j,k)=pl(i,j,k+1)-pl(i,j,k)
        scal(i,j,k)=dp(i,j,k)*(.5*(pl(i,j,k)+pl(i,j,k+1))/300.)**.8

!-----compute scaled water vapor amount, unit is g/cm**2

        wh(i,j,k)=1.02*wa(i,j,k)*scal(i,j,k)*                           &
                  (1.-0.00135*(ta(i,j,k)-240.))
        swh(i,j,k+1)=swh(i,j,k)+wh(i,j,k)

!-----compute ozone amount, unit is (cm-atm)stp.

        oh(i,j,k)=1.02*oa(i,j,k)*dp(i,j,k)*466.7

      END DO
    END DO
  END DO


!-----scale cloud optical thickness in each layer from taucld (with
!  cloud amount fcld) to tauclb and tauclf (with cloud amount cc).
!  tauclb is the scaled optical thickness for beam radiation and
!  tauclf is for diffuse radiation.

  CALL cldscale(idim,jdim,                                              &
                m,n,np,cosz,fcld,taucldi,taucldl,ict,icb,               &
                cc,tauclb,tauclf)

!-----initialize fluxes for all-sky (flx), clear-sky (flc), and
!  flux reduction (df)

  DO k=1,np+1
    DO j=1,n
      DO i=1,m
        flx(i,j,k)=0.
        flc(i,j,k)=0.
        df(i,j,k)=0.
      END DO
    END DO
  END DO

!-----compute solar ir fluxes

  CALL solir (idim,jdim,                                                &
              m,n,np,wh,taucldi,taucldl,tauclb,tauclf,                  &
              reffi,reffl,ict,icb,                                      &
              fcld,cc,taual,csm,rsirbm,rsirdf,                          &
              flx,flc,fdirir,fdifir,                                    &
              tem2d1, tem2d2, tem2d3, tem2d4,                           &
              tem2d5, tem2d6, tem2d7, tem2d8, tem2d9,                   &
              tem2d10,tem2d11,tem2d12,tem2d13,tem2d14,                  &
              tem3d1, tem3d2, tem3d3, tem3d4,                           &
              tem4d1, tem4d2, tem4d3, tem4d4, tem4d5,                   &
              tem2d15,tem2d16,tem2d17,tem2d18,tem2d19,                  &
              tem3d5,                                                   &
              tem5d1,tem5d2,tem5d3,tem5d4,tem5d5)

!  DO k=1,np
!    write (6,'(a,i2,a,i2,a,e20.10)')
!    :  'IR:  flx(2,2,',k+1,')-flx(2,2,',k,')=',
!    :  flx(2,2,k+1)-flx(2,2,k)
!  ENDDO

!-----compute and update uv and par fluxes

  CALL soluv (idim,jdim,                                                &
              m,n,np,oh,dp,taucldi,taucldl,tauclb,tauclf,               &
              reffi,reffl,ict,icb,                                      &
              fcld,cc,taual,csm,rsuvbm,rsuvdf,                          &
              flx,flc,fdirpar,fdifpar,                                  &
              fdiruv,fdifuv,                                            & ! MS add
              tem2d1, tem2d2, tem2d3,                                   &
              tem2d5, tem2d6, tem2d7, tem2d8, tem2d9,                   &
              tem2d10,tem2d11,tem2d12,tem2d13,tem2d14,                  &
              tem3d1, tem3d3, tem3d4,                                   &
              tem4d1, tem4d2, tem4d3, tem4d4, tem4d5,                   &
              tem2d15,tem2d16,tem2d17,tem2d18,tem2d19,                  &
              tem3d5,                                                   &
              tem5d1,tem5d2,tem5d3,tem5d4,tem5d5)

!  write (6,'(a/a3,2a15)') 'Total flux for PAR',
!    :  'k','dflxPAR','flxPAR'
!  DO k=1,np
!    write (6,'(i3,2e15.7)')
!    :  k,flx(2,2,k+1)-flx(2,2,k),flx(2,2,k)
!  ENDDO

!-----compute scaled amount of o2 (so2), unit is (cm-atm)stp.

  DO k=1,np
    DO j=1,n
      DO i=1,m
        so2(i,j,k+1)=so2(i,j,k)+165.22*scal(i,j,k)
      END DO
    END DO
  END DO

!-----compute flux reduction due to oxygen following
!   chou (J. climate, 1990). The fraction 0.0287 is the
!   extraterrestrial solar flux in the o2 bands.

  DO k=2,np+1
    DO j=1,n
      DO i=1,m
        x=so2(i,j,k)*csm(i,j)
        df(i,j,k)=df(i,j,k)+0.0287*(1.-expmn(-0.00027*SQRT(x)))
      END DO
    END DO
  END DO

!  write (6,'(a/a3,3a15)')
!    : 'Flux reduction due to oxygen, i=j=3',
!    : 'k','ddf','df','so2'
!  DO k=2,np+1
!    write (6,'(i3,3e15.7)')
!    :  k,df(3,3,k)-df(3,3,k-1),df(3,3,k),so2(3,3,k)
!  ENDDO

!-----compute scaled amounts for co2 (so2). unit is (cm-atm)stp.

  DO k=1,np
    DO j=1,n
      DO i=1,m
        so2(i,j,k+1)=so2(i,j,k)+co2*789.*scal(i,j,k)
      END DO
    END DO
  END DO

!-----compute and update flux reduction due to co2 following
!  chou (J. Climate, 1990)

  CALL flxco2(m,n,np,so2,swh,csm,df)

!-----adjust for the effect of o2 cnd co2 on clear-sky fluxes.

  DO k=2,np+1
    DO j=1,n
      DO i=1,m
        flc(i,j,k)=flc(i,j,k)-df(i,j,k)
      END DO
    END DO
  END DO

!-----adjust for the all-sky fluxes due to o2 and co2.  It is
!  assumed that o2 and co2 have no effects on solar radiation
!  below clouds. clouds are assumed randomly overlapped.

  DO j=1,n
    DO i=1,m
      sdf(i,j)=0.0
      sclr(i,j)=1.0
    END DO
  END DO

  DO k=1,np
    DO j=1,n
      DO i=1,m

!-----sclr is the fraction of clear sky.
!  sdf is the flux reduction below clouds.

        IF(fcld(i,j,k) > 0.01) THEN
          sdf(i,j)=sdf(i,j)+df(i,j,k)*sclr(i,j)*fcld(i,j,k)
          sclr(i,j)=sclr(i,j)*(1.-fcld(i,j,k))
        END IF
        flx(i,j,k+1)=flx(i,j,k+1)-sdf(i,j)-df(i,j,k+1)*sclr(i,j)

      END DO
    END DO
  END DO

!-----adjust for the direct downward ir flux.
  DO j=1,n
    DO i=1,m
!!!      fdirir(i,j)=fdirir(i,j)-sdf(i,j)-df(i,j,np+1)*sclr(i,j)
      ! MS - this line was altered to keep it consistent with the net flux calculation
      fdirir(i,j)=fdirir(i,j)-( sdf(i,j) + df(i,j,np+1)*sclr(i,j) )/(1. - rsirbm(i,j))
    END DO
  END DO

!  write (6,'(a3,2a15)') 'k','flx','flc'
!  DO k=1,np
!    write (6,'(i3,2e15.8)') k,flx(1,1,k),flc(1,1,k)
!  ENDDO
!
!  write (6,'(4a15)')
!    : 'fdirir','fdifir','fdirpar','fdifpar'
!
!  write (6,'(4e15.8)')
!    :  fdirir(1,1),fdifir(1,1),fdirpar(1,1),fdifpar(1,1)

  RETURN
END SUBROUTINE sorad
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SOLIR                      ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE solir(idim,jdim,                                             &
           m,n,np,wh,taucldi,taucldl,tauclb,tauclf,reffi,reffl,         &
           ict,icb,fcld,cc,taual,csm,rsirbm,rsirdf,                     &
           flx,flc,fdirir,fdifir,                                       &
           fsdir,fsdif, ssaclt,asyclt,                                  &
           rr1t,tt1t,td1t,rs1t,ts1t,                                    &
           rr2t,tt2t,td2t,rs2t,ts2t,                                    &
           ssacl,asycl,fall,fclr,                                       &
           rr,tt,td,rs,ts,                                              &
           tem2d1,tem2d2,tem2d3,tem2d4,tem2d5,                          &
           tem3d,                                                       &
           tem5d1,tem5d2,tem5d3,tem5d4,tem5d5)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate solar flux in the infrared region. The spectrum is
!  divided into three bands. (See original comments.)
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (a) Radiative Transfer Model: M.-D. Chou and M. Suarez
!          (b) Cloud Optics:Tao, Lang, Simpson, Sui, Ferrier and
!              Chou (1996)
!
!  MODIFICATION:
!
!  03/11/1996 (Yuhe Liu)
!  Modified the original code from 1-D to 3-D
!
!  12/14/1998 (Keith Brewster)
!  Modified values of ssaal and asymal in data statement.
!
!-----------------------------------------------------------------------
!
!  ORIGINAL COMMENTS:
!
!************************************************************************
!  compute solar flux in the infrared region. The spectrum is divided
!   into three bands:
!
!       band   wavenumber(/cm)  wavelength (micron)
!        1       1000-4400         2.27-10.0
!        2       4400-8200         1.22-2.27
!        3       8200-14300        0.70-1.22
!
!----- Input parameters:                            units      size
!
!   number of soundings in zonal direction (m)       n/d        1
!   number of soundings in meridional direction (n)  n/d        1
!   maximum number of soundings in                   n/d        1
!       meridional direction (ndim)
!   number of atmospheric layers (np)                n/d        1
!   layer water vapor content (wh)                 gm/cm^2    m*n*np
!   cloud optical thickness (taucld)                 n/d      m*ndim*np*2
!       index 1 for ice paticles
!       index 2 for liquid particles
!   scaled cloud optical thickness                   n/d      m*n*np
!       for beam radiation (tauclb)
!   scaled cloud optical thickness                   n/d      m*n*np
!       for diffuse radiation  (tauclf)
!   effective cloud-particle size (reff)           micrometer m*ndim*np*2
!       index 1 for ice paticles
!       index 2 for liquid particles
!   level index separating high and                  n/d      m*n
!       middle clouds (ict)
!   level index separating middle and                n/d      m*n
!       low clouds (icb)
!   cloud amount (fcld)                            fraction   m*ndim*np
!   cloud amount of high, middle, and                n/d      m*n*3
!       low clouds (cc)
!   aerosol optical thickness (taual)                n/d      m*ndim*np
!   cosecant of the solar zenith angle (csm)         n/d      m*n
!   near ir surface albedo for beam                fraction   m*ndim
!             radiation (rsirbm)
!   near ir surface albedo for diffuse             fraction   m*ndim
!             radiation (rsirdf)
!
!----- output (updated) parameters:
!
!   all-sky flux (downward-upward) (flx)           fraction   m*ndim*(np+1)
!   clear-sky flux (downward-upward) (flc)         fraction   m*ndim*(np+1)
!   all-sky direct downward ir flux at
!       the surface (fdirir)                    fraction   m*ndim
!   all-sky diffuse downward ir flux at
!       the surface (fdifir)                    fraction   m*ndim
!
!----- note: the following parameters must be specified by users:
!   aerosol single scattering albedo (ssaal)         n/d      nband
!   aerosol asymmetry factor (asyal)                 n/d      nband
!
!*************************************************************************
!
!fpp$ expand (expmn)
!fpp$ expand (deledd)
!fpp$ expand (sagpol)
!!dir$ inline always expmn,deledd,sagpol
!*$*  inline routine (expmn,deledd,sagpol)
!
!-----------------------------------------------------------------------
!
!  Variable declarations
!
!-----------------------------------------------------------------------
!
  use radlib3d_module, only : expmn,deledd,sagpol
  IMPLICIT NONE

!-----input parameters

  INTEGER :: idim,jdim
  INTEGER :: m,n,np,ict,icb

  REAL :: taucldi(idim,jdim,np)
  REAL :: taucldl(idim,jdim,np)
  REAL :: reffi(idim,jdim,np)
  REAL :: reffl(idim,jdim,np)
  REAL :: fcld(idim,jdim,np)
  REAL :: rsirbm(idim,jdim)
  REAL :: rsirdf(idim,jdim)
  REAL :: taual(idim,jdim,np)

  REAL :: tauclb(m,n,np)
  REAL :: tauclf(m,n,np)
  REAL :: cc(m,n,3)
  REAL :: wh(m,n,np)
  REAL :: csm(m,n)

!-----output (updated) parameters

  REAL :: flx(idim,jdim,np+1)
  REAL :: flc(idim,jdim,np+1)
  REAL :: fdirir(idim,jdim)
  REAL :: fdifir(idim,jdim)

!-----static parameters

  INTEGER :: nk,nband
  PARAMETER (nk=10,nband=3)
  REAL :: xk(nk),hk(nband,nk),ssaal(nband),asyal(nband)
  REAL :: aia(nband,3),awa(nband,3),aig(nband,3),awg(nband,3)

!-----temporary array

  REAL :: fsdir(m,n)
  REAL :: fsdif(m,n)

  REAL :: ssaclt(m,n)
  REAL :: asyclt(m,n)
  REAL :: rr1t(m,n)
  REAL :: tt1t(m,n)
  REAL :: td1t(m,n)
  REAL :: rs1t(m,n)
  REAL :: ts1t(m,n)
  REAL :: rr2t(m,n)
  REAL :: tt2t(m,n)
  REAL :: td2t(m,n)
  REAL :: rs2t(m,n)
  REAL :: ts2t(m,n)

  REAL :: ssacl(m,n,np)
  REAL :: asycl(m,n,np)
  REAL :: fall(m,n,np+1)
  REAL :: fclr(m,n,np+1)

  REAL :: rr(m,n,np+1,2)
  REAL :: tt(m,n,np+1,2)
  REAL :: td(m,n,np+1,2)
  REAL :: rs(m,n,np+1,2)
  REAL :: ts(m,n,np+1,2)

!-----temporary arrays used in subroutine cldflx

  REAL :: tem2d1(m,n)
  REAL :: tem2d2(m,n)
  REAL :: tem2d3(m,n)
  REAL :: tem2d4(m,n)
  REAL :: tem2d5(m,n)

  REAL :: tem3d (m,n,np+1)

  REAL :: tem5d1(m,n,np+1,2,2)
  REAL :: tem5d2(m,n,np+1,2,2)
  REAL :: tem5d3(m,n,np+1,2,2)
  REAL :: tem5d4(m,n,np+1,2,2)
  REAL :: tem5d5(m,n,np+1,2,2)

!-----temporary variables

  INTEGER :: ib,ik,i,j,k

  REAL :: tauwv,tausto,ssatau,asysto,tauto,ssato,asyto
  REAL :: taux,reff1,reff2,w1,w2,g1,g2

!-----function

!!!  REAL :: expmn

!-----water vapor absorption coefficient for 10 k-intervals.
!  unit: cm^2/gm

  DATA xk/                                                              &
      0.0010, 0.0133, 0.0422, 0.1334, 0.4217,                           &
      1.334,  5.623,  31.62,  177.8,  1000.0/

!-----water vapor k-distribution function,
!  the sum of hk is 0.52926. unit: fraction

  DATA hk/                                                              &
      .01074,.08236,.20673,  .00360,.01157,.03497,                      &
      .00411,.01133,.03011,  .00421,.01143,.02260,                      &
      .00389,.01240,.01336,  .00326,.01258,.00696,                      &
      .00499,.01381,.00441,  .00465,.00650,.00115,                      &
      .00245,.00244,.00026,  .00145,.00094,.00000/

!-----aerosol single-scattering albedo and asymmetry factor

!  data ssaal,asyal/nband*0.999,nband*0.850/
  DATA ssaal/0.75,0.55,0.90/
  DATA asyal/0.40,0.50,0.60/

!  data ssaal,asyal/nband*0.400,nband*0.500/

!-----coefficients for computing the single scattering albedo of
!  ice clouds from ssa=1-(aia(*,1)+aia(*,2)*reff+aia(*,3)*reff**2)

  DATA aia/                                                             &
      .08938331, .00215346,-.00000260,                                  &
      .00299387, .00073709, .00000746,                                  &
      -.00001038,-.00000134, .00000000/

!-----coefficients for computing the single scattering albedo of
!  liquid clouds from ssa=1-(awa(*,1)+awa(*,2)*reff+awa(*,3)*reff**2)

  DATA awa/                                                             &
      .01209318,-.00019934, .00000007,                                  &
      .01784739, .00088757, .00000845,                                  &
      -.00036910,-.00000650,-.00000004/

!-----coefficients for computing the asymmetry factor of ice clouds
!  from asycl=aig(*,1)+aig(*,2)*reff+aig(*,3)*reff**2

  DATA aig/                                                             &
      .84090400, .76098937, .74935228,                                  &
      .00126222, .00141864, .00119715,                                  &
      -.00000385,-.00000396,-.00000367/

!-----coefficients for computing the asymmetry factor of liquid clouds
!  from asycl=awg(*,1)+awg(*,2)*reff+awg(*,3)*reff**2

  DATA awg/                                                             &
      .83530748, .74513197, .79375035,                                  &
      .00257181, .01370071, .00832441,                                  &
      .00005519,-.00038203,-.00023263/
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!-----initialize surface fluxes, reflectances, and transmittances

  DO j= 1, n
    DO i= 1, m
      fdirir(i,j)=0.0
      fdifir(i,j)=0.0
      rr(i,j,np+1,1)=rsirbm(i,j)
      rr(i,j,np+1,2)=rsirdf(i,j)
      rs(i,j,np+1,1)=rsirbm(i,j)
      rs(i,j,np+1,2)=rsirdf(i,j)
      td(i,j,np+1,1)=0.0
      td(i,j,np+1,2)=0.0
      tt(i,j,np+1,1)=0.0
      tt(i,j,np+1,2)=0.0
      ts(i,j,np+1,1)=0.0
      ts(i,j,np+1,2)=0.0
    END DO
  END DO

!-----integration over spectral bands

  DO ib=1,nband

!-----compute cloud single scattering albedo and asymmetry factor
!  for a mixture of ice and liquid particles.

    DO k= 1, np

      DO j= 1, n
        DO i= 1, m

          ssaclt(i,j)=1.0
          asyclt(i,j)=1.0

          taux=taucldi(i,j,k)+taucldl(i,j,k)
          IF (taux > 0.05 .AND. fcld(i,j,k) > 0.01) THEN

            reff1=MIN(reffi(i,j,k),130.)
            reff2=MIN(reffl(i,j,k),20.0)

            w1=(1.-(aia(ib,1)+(aia(ib,2)+                               &
                aia(ib,3)*reff1)*reff1))*taucldi(i,j,k)
            w2=(1.-(awa(ib,1)+(awa(ib,2)+                               &
                awa(ib,3)*reff2)*reff2))*taucldl(i,j,k)
            ssaclt(i,j)=(w1+w2)/taux

            g1=(aig(ib,1)+(aig(ib,2)+aig(ib,3)*reff1)*reff1)*w1
            g2=(awg(ib,1)+(awg(ib,2)+awg(ib,3)*reff2)*reff2)*w2
            asyclt(i,j)=(g1+g2)/(w1+w2)

          END IF

        END DO
      END DO

      DO j=1,n
        DO i=1,m
          ssacl(i,j,k)=ssaclt(i,j)
        END DO
      END DO
      DO j=1,n
        DO i=1,m
          asycl(i,j,k)=asyclt(i,j)
        END DO
      END DO

    END DO

!-----integration over the k-distribution function

    DO ik=1,nk

      DO k= 1, np

        DO j= 1, n
          DO i= 1, m

            tauwv=xk(ik)*wh(i,j,k)

!-----compute total optical thickness, single scattering albedo,
!  and asymmetry factor.

            tausto=tauwv+taual(i,j,k)+1.0E-8
            ssatau=ssaal(ib)*taual(i,j,k)
            asysto=asyal(ib)*ssaal(ib)*taual(i,j,k)

!-----compute reflectance and transmittance for cloudless layers

            tauto=tausto
            ssato=ssatau/tauto+1.0E-8

            IF (ssato > 0.001) THEN

              ssato=MIN(ssato,0.999999)
              asyto=asysto/(ssato*tauto)

              CALL deledd (tauto,ssato,asyto,csm(i,j),                  &
                           rr1t(i,j),tt1t(i,j),td1t(i,j))

              CALL sagpol (tauto,ssato,asyto,rs1t(i,j),ts1t(i,j))

            ELSE

              td1t(i,j)=expmn(-tauto*csm(i,j))
              ts1t(i,j)=expmn(-1.66*tauto)
              tt1t(i,j)=0.0
              rr1t(i,j)=0.0
              rs1t(i,j)=0.0

            END IF

!-----compute reflectance and transmittance for cloud layers

            IF (tauclb(i,j,k) < 0.01) THEN

              rr2t(i,j)=rr1t(i,j)
              tt2t(i,j)=tt1t(i,j)
              td2t(i,j)=td1t(i,j)
              rs2t(i,j)=rs1t(i,j)
              ts2t(i,j)=ts1t(i,j)

            ELSE

              tauto=tausto+tauclb(i,j,k)
              ssato=(ssatau+ssacl(i,j,k)*tauclb(i,j,k))/tauto+1.0E-8
              ssato=MIN(ssato,0.999999)
              asyto=(asysto+asycl(i,j,k)*ssacl(i,j,k)*tauclb(i,j,k))/   &
                    (ssato*tauto)

              CALL deledd (tauto,ssato,asyto,csm(i,j),                  &
                           rr2t(i,j),tt2t(i,j),td2t(i,j))

              tauto=tausto+tauclf(i,j,k)
              ssato=(ssatau+ssacl(i,j,k)*tauclf(i,j,k))/tauto+1.0E-8
              ssato=MIN(ssato,0.999999)
              asyto=(asysto+asycl(i,j,k)*ssacl(i,j,k)*tauclf(i,j,k))/   &
                    (ssato*tauto)

              CALL sagpol (tauto,ssato,asyto,rs2t(i,j),ts2t(i,j))

            END IF

          END DO
        END DO

        DO j=1,n
          DO i=1,m
            rr(i,j,k,1)=rr1t(i,j)
          END DO
        END DO
        DO j=1,n
          DO i=1,m
            tt(i,j,k,1)=tt1t(i,j)
          END DO
        END DO
        DO j=1,n
          DO i=1,m
            td(i,j,k,1)=td1t(i,j)
          END DO
        END DO
        DO j=1,n
          DO i=1,m
            rs(i,j,k,1)=rs1t(i,j)
          END DO
        END DO
        DO j=1,n
          DO i=1,m
            ts(i,j,k,1)=ts1t(i,j)
          END DO
        END DO

        DO j=1,n
          DO i=1,m
            rr(i,j,k,2)=rr2t(i,j)
          END DO
        END DO
        DO j=1,n
          DO i=1,m
            tt(i,j,k,2)=tt2t(i,j)
          END DO
        END DO
        DO j=1,n
          DO i=1,m
            td(i,j,k,2)=td2t(i,j)
          END DO
        END DO
        DO j=1,n
          DO i=1,m
            rs(i,j,k,2)=rs2t(i,j)
          END DO
        END DO
        DO j=1,n
          DO i=1,m
            ts(i,j,k,2)=ts2t(i,j)
          END DO
        END DO

      END DO

!-----flux calculations

      CALL cldflx (m,n,np,ict,icb,cc,rr,tt,td,rs,ts,                    &
                   fclr,fall,fsdir,fsdif,                               &
                   tem2d1,tem2d2,tem2d3,tem2d4,tem2d5,                  &
                   tem3d, tem5d1,tem5d2,tem5d3,tem5d4,tem5d5)

      DO k= 1, np+1
        DO j= 1, n
          DO i= 1, m
            flx(i,j,k) = flx(i,j,k)+fall(i,j,k)*hk(ib,ik)
          END DO
        END DO
        DO j= 1, n
          DO i= 1, m
            flc(i,j,k) = flc(i,j,k)+fclr(i,j,k)*hk(ib,ik)
          END DO
        END DO
      END DO

!  write (6,'(a/2a3,2a15)') 'Flux for IR',
!    :  'ib','k','dflxIR','flxIR'
!  DO k=1,np
!    write (6,'(2i3,2e15.7)')
!    :  ib,k,fall(2,2,k+1)-fall(2,2,k),fall(2,2,k)
!  ENDDO

      DO j= 1, n
        DO i= 1, m
          fdirir(i,j) = fdirir(i,j)+fsdir(i,j)*hk(ib,ik)
          fdifir(i,j) = fdifir(i,j)+fsdif(i,j)*hk(ib,ik)
        END DO
      END DO

    END DO
  END DO

  RETURN
END SUBROUTINE solir
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SOLUV                      ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######               University of Oklahoma                 ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE soluv(idim,jdim,                                             &
           m,n,np,oh,dp,taucldi,taucldl,tauclb,tauclf,reffi,reffl,      &
           ict,icb,fcld,cc,taual,csm,rsuvbm,rsuvdf,                     &
           flx,flc,fdirpar,fdifpar,                                     &
           fdiruv,fdifuv,                                               & ! MS add
           fsdir,fsdif,asyclt,                                          &
           rr1t,tt1t,td1t,rs1t,ts1t,                                    &
           rr2t,tt2t,td2t,rs2t,ts2t,                                    &
           asycl,fall,fclr,                                             &
           td,rr,tt,rs,ts,                                              &
           tem2d1,tem2d2,tem2d3,tem2d4,tem2d5,                          &
           tem3d,                                                       &
           tem5d1,tem5d2,tem5d3,tem5d4,tem5d5)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate solar fluxes in the uv+visible region. the spectrum is
!  grouped into 8 bands. (See original comments.)
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (a) Radiative Transfer Model: M.-D. Chou and M. Suarez
!          (b) Cloud Optics:Tao, Lang, Simpson, Sui, Ferrier and
!              Chou (1996)
!
!  MODIFICATION HISTORY:
!
!  03/15/1996 (Yuhe Liu)
!  Adapted the original code and formatted it in accordance with the
!  ARPS coding standard.
!
!  12/07/1998 (Keith Brewster)
!  Modified values of ssaal and asymal in data statement.
!
!-----------------------------------------------------------------------
!
!  ORIGINAL COMMENTS:
!
!************************************************************************
!  compute solar fluxes in the uv+visible region. the spectrum is
!  grouped into 8 bands:
!
!           Band     Micrometer
!
!    UV-C    1.     .175 - .225
!            2.     .225 - .245
!                   .260 - .280
!            3.     .245 - .260
!
!    UV-B    4.     .280 - .295
!            5.     .295 - .310
!            6.     .310 - .320
!
!    UV-A    7.     .320 - .400
!
!    PAR     8.     .400 - .700
!
!----- Input parameters:                            units      size
!
!   number of soundings in zonal direction (m)       n/d        1
!   number of soundings in meridional direction (n)  n/d        1
!   maximum number of soundings in                   n/d        1
!        meridional direction (ndim)
!   number of atmospheric layers (np)                n/d        1
!   layer ozone content (oh)                      (cm-atm)stp m*n*np
!   layer pressure thickness (dp)                    mb       m*n*np
!   cloud optical thickness (taucld)                 n/d      m*ndim*np*2
!       index 1 for ice paticles
!       index 2 for liquid particles
!   scaled cloud optical thickness                   n/d      m*n*np
!       for beam radiation (tauclb)
!   scaled cloud optical thickness                   n/d      m*n*np
!       for diffuse radiation  (tauclf)
!   effective cloud-particle size (reff)           micrometer m*ndim*np*2
!       index 1 for ice paticles
!       index 2 for liquid particles
!   level indiex separating high and                 n/d      m*n
!       middle clouds (ict)
!   level indiex separating middle and               n/d      m*n
!       low clouds (icb)
!   cloud amount (fcld)                            fraction   m*ndim*np
!   cloud amounts of high, middle, and               n/d      m*n*3
!       low clouds (cc)
!   aerosol optical thickness (taual)                n/d      m*ndim*np
!   cosecant of the solar zenith angle (csm)         n/d      m*n
!   uv+par surface albedo for beam                 fraction   m*ndim
!        radiation (rsuvbm)
!   uv+par surface albedo for diffuse              fraction   m*ndim
!        radiation (rsuvdf)
!
!----- output (updated) parameters:
!
!   all-sky net downward flux (flx)                fraction   m*ndim*(np+1)
!   clear-sky net downward flux (flc)              fraction   m*ndim*(np+1)
!   all-sky direct downward par flux at
!       the surface (fdirpar)                   fraction   m*ndim
!   all-sky diffuse downward par flux at
!       the surface (fdifpar)                   fraction   m*ndim
!
! MS - add new vars to track uv fluxes at surface
!   all-sky direct downward par flux at
!       the surface (fdiruv)                   fraction   m*ndim
!   all-sky diffuse downward par flux at
!       the surface (fdifuv)                   fraction   m*ndim
!
!----- note: the following parameters must be specified by users:
!
!   aerosol single scattering albedo (ssaal)         n/d        1
!   aerosol asymmetry factor (asyal)                 n/d        1
!
!
!***********************************************************************
!
!fpp$ expand (deledd)
!fpp$ expand (sagpol)
!!dir$ inline always deledd,sagpol
!*$*  inline routine (deledd,sagpol)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  use radlib3d_module, only : deledd,sagpol
  IMPLICIT NONE

!-----input parameters

  INTEGER :: idim,jdim
  INTEGER :: m,n,np,ict,icb

  REAL :: taucldi(idim,jdim,np)
  REAL :: taucldl(idim,jdim,np)
  REAL :: reffi(idim,jdim,np)
  REAL :: reffl(idim,jdim,np)
  REAL :: fcld(idim,jdim,np)
  REAL :: taual(idim,jdim,np)
  REAL :: rsuvbm(idim,jdim)
  REAL :: rsuvdf(idim,jdim)

  REAL :: tauclb(m,n,np)
  REAL :: tauclf(m,n,np)
  REAL :: cc(m,n,3)
  REAL :: oh(m,n,np)
  REAL :: dp(m,n,np)
  REAL :: csm(m,n)

!-----output (updated) parameter

  REAL :: flx(idim,jdim,np+1)
  REAL :: flc(idim,jdim,np+1)
  REAL :: fdirpar(idim,jdim)
  REAL :: fdifpar(idim,jdim)
  REAL :: fdiruv(idim,jdim) ! MS add
  REAL :: fdifuv(idim,jdim) ! MS add

!-----static parameters

  INTEGER :: nband
  PARAMETER (nband=8)
  REAL :: hk(nband)
  REAL :: xk(nband)
  REAL :: ry(nband)
  REAL :: asyal(nband)
  REAL :: ssaal(nband)
  REAL :: aig(3)
  REAL :: awg(3)

!-----temporary array

  REAL :: fsdir(m,n)
  REAL :: fsdif(m,n)
  REAL :: asyclt(m,n)
  REAL :: rr1t(m,n)
  REAL :: tt1t(m,n)
  REAL :: td1t(m,n)
  REAL :: rs1t(m,n)
  REAL :: ts1t(m,n)
  REAL :: rr2t(m,n)
  REAL :: tt2t(m,n)
  REAL :: td2t(m,n)
  REAL :: rs2t(m,n)
  REAL :: ts2t(m,n)

  REAL :: asycl(m,n,np)
  REAL :: fall(m,n,np+1)
  REAL :: fclr(m,n,np+1)

  REAL :: td(m,n,np+1,2)
  REAL :: rr(m,n,np+1,2)
  REAL :: tt(m,n,np+1,2)
  REAL :: rs(m,n,np+1,2)
  REAL :: ts(m,n,np+1,2)

!-----temporary arrays used in subroutine cldflx

  REAL :: tem2d1(m,n)
  REAL :: tem2d2(m,n)
  REAL :: tem2d3(m,n)
  REAL :: tem2d4(m,n)
  REAL :: tem2d5(m,n)

  REAL :: tem3d (m,n,np+1)

  REAL :: tem5d1(m,n,np+1,2,2)
  REAL :: tem5d2(m,n,np+1,2,2)
  REAL :: tem5d3(m,n,np+1,2,2)
  REAL :: tem5d4(m,n,np+1,2,2)
  REAL :: tem5d5(m,n,np+1,2,2)

!-----temporary variables

  INTEGER :: i,j,k,ib
  REAL :: taurs,tauoz,tausto,ssatau,asysto,tauto,ssato,asyto
  REAL :: taux,reff1,reff2,g1,g2

!-----hk is the fractional extra-terrestrial solar flux.
!  the sum of hk is 0.47074.

  DATA hk/.00057, .00367, .00083, .00417,                               &
          .00600, .00556, .05913, .39081/

!-----xk is the ozone absorption coefficient. unit: /(cm-atm)stp

  DATA xk /30.47, 187.2,  301.9,   42.83,                               &
           7.09,  1.25,   0.0345,  0.0539/

!-----ry is the extinction coefficient for Rayleigh scattering.
!  unit: /mb.

  DATA ry /.00604, .00170, .00222, .00132,                              &
           .00107, .00091, .00055, .00012/

!-----aerosol single-scattering albedo and asymmetry factor

!  data ssaal,asyal/nband*0.999,nband*0.850/
  DATA ssaal,asyal/nband*0.960,nband*0.700/

!-----coefficients for computing the asymmetry factor of ice clouds
!  from asycl=aig(*,1)+aig(*,2)*reff+aig(*,3)*reff**2

  DATA aig/.74625000,.00105410,-.00000264/

!-----coefficients for computing the asymmetry factor of liquid
!  clouds from asycl=awg(*,1)+awg(*,2)*reff+awg(*,3)*reff**2

  DATA awg/.82562000,.00529000,-.00014866/

!-----initialize surface reflectances and transmittances

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO j= 1, n
    DO i= 1, m
      fdirpar(i,j)=0.0
      fdifpar(i,j)=0.0
      fdiruv(i,j)=0.0 ! MS add
      fdifuv(i,j)=0.0 ! MS add
      rr(i,j,np+1,1)=rsuvbm(i,j)
      rr(i,j,np+1,2)=rsuvdf(i,j)
      rs(i,j,np+1,1)=rsuvbm(i,j)
      rs(i,j,np+1,2)=rsuvdf(i,j)
      td(i,j,np+1,1)=0.0
      td(i,j,np+1,2)=0.0
      tt(i,j,np+1,1)=0.0
      tt(i,j,np+1,2)=0.0
      ts(i,j,np+1,1)=0.0
      ts(i,j,np+1,2)=0.0
    END DO
  END DO

!-----compute cloud asymmetry factor for a mixture of
!  liquid and ice particles.  unit of reff is micrometers.

  DO k= 1, np

    DO j= 1, n
      DO i= 1, m

        asyclt(i,j)=1.0

        taux=taucldi(i,j,k)+taucldl(i,j,k)
        IF (taux > 0.05 .AND. fcld(i,j,k) > 0.01) THEN

          reff1=MIN(reffi(i,j,k),130.)
          reff2=MIN(reffl(i,j,k),20.0)

          g1=(aig(1)+(aig(2)+aig(3)*reff1)*reff1)*taucldi(i,j,k)
          g2=(awg(1)+(awg(2)+awg(3)*reff2)*reff2)*taucldl(i,j,k)
          asyclt(i,j)=(g1+g2)/taux

        END IF

      END DO
    END DO

    DO j=1,n
      DO i=1,m
        asycl(i,j,k)=asyclt(i,j)
      END DO
    END DO

  END DO

!-----integration over spectral bands

  DO ib=1,nband

    DO k= 1, np

      DO j= 1, n
        DO i= 1, m

!-----compute ozone and rayleigh optical thicknesses

          taurs=ry(ib)*dp(i,j,k)
          tauoz=xk(ib)*oh(i,j,k)

!-----compute total optical thickness, single scattering albedo,
!  and asymmetry factor

          tausto=taurs+tauoz+taual(i,j,k)+1.0E-8
          ssatau=ssaal(ib)*taual(i,j,k)+taurs
          asysto=asyal(ib)*ssaal(ib)*taual(i,j,k)

!-----compute reflectance and transmittance for cloudless layers

          tauto=tausto
          ssato=ssatau/tauto+1.0E-8
          ssato=MIN(ssato,0.999999)
          asyto=asysto/(ssato*tauto)

          CALL deledd (tauto,ssato,asyto,csm(i,j),                      &
                       rr1t(i,j),tt1t(i,j),td1t(i,j))

          CALL sagpol (tauto,ssato,asyto,rs1t(i,j),ts1t(i,j))

!-----compute reflectance and transmittance for cloud layers

          IF (tauclb(i,j,k) < 0.01) THEN

            rr2t(i,j)=rr1t(i,j)
            tt2t(i,j)=tt1t(i,j)
            td2t(i,j)=td1t(i,j)
            rs2t(i,j)=rs1t(i,j)
            ts2t(i,j)=ts1t(i,j)

          ELSE

            tauto=tausto+tauclb(i,j,k)
            ssato=(ssatau+tauclb(i,j,k))/tauto+1.0E-8
            ssato=MIN(ssato,0.999999)
            asyto=(asysto+asycl(i,j,k)*tauclb(i,j,k))/(ssato*tauto)

            CALL deledd (tauto,ssato,asyto,csm(i,j),                    &
                         rr2t(i,j),tt2t(i,j),td2t(i,j))

            tauto=tausto+tauclf(i,j,k)
            ssato=(ssatau+tauclf(i,j,k))/tauto+1.0E-8
            ssato=MIN(ssato,0.999999)
            asyto=(asysto+asycl(i,j,k)*tauclf(i,j,k))/(ssato*tauto)

            CALL sagpol (tauto,ssato,asyto,rs2t(i,j),ts2t(i,j))

          END IF

        END DO
      END DO

      DO j=1,n
        DO i=1,m
          rr(i,j,k,1)=rr1t(i,j)
        END DO
      END DO
      DO j=1,n
        DO i=1,m
          tt(i,j,k,1)=tt1t(i,j)
        END DO
      END DO
      DO j=1,n
        DO i=1,m
          td(i,j,k,1)=td1t(i,j)
        END DO
      END DO
      DO j=1,n
        DO i=1,m
          rs(i,j,k,1)=rs1t(i,j)
        END DO
      END DO
      DO j=1,n
        DO i=1,m
          ts(i,j,k,1)=ts1t(i,j)
        END DO
      END DO

      DO j=1,n
        DO i=1,m
          rr(i,j,k,2)=rr2t(i,j)
        END DO
      END DO
      DO j=1,n
        DO i=1,m
          tt(i,j,k,2)=tt2t(i,j)
        END DO
      END DO
      DO j=1,n
        DO i=1,m
          td(i,j,k,2)=td2t(i,j)
        END DO
      END DO
      DO j=1,n
        DO i=1,m
          rs(i,j,k,2)=rs2t(i,j)
        END DO
      END DO
      DO j=1,n
        DO i=1,m
          ts(i,j,k,2)=ts2t(i,j)
        END DO
      END DO

    END DO

!-----flux calculations

    CALL cldflx (m,n,np,ict,icb,cc,rr,tt,td,rs,ts,                      &
                 fclr,fall,fsdir,fsdif,                                 &
                 tem2d1,tem2d2,tem2d3,tem2d4,tem2d5,                    &
                 tem3d, tem5d1,tem5d2,tem5d3,tem5d4,tem5d5)

    DO k= 1, np+1
      DO j= 1, n
        DO i= 1, m
          flx(i,j,k)=flx(i,j,k)+fall(i,j,k)*hk(ib)
        END DO
      END DO
      DO j= 1, n
        DO i= 1, m
          flc(i,j,k)=flc(i,j,k)+fclr(i,j,k)*hk(ib)
        END DO
      END DO
    END DO

!  write (6,'(a/2a3,2a15)') 'Flux for PAR',
!    :  'ib','k','dflxpa','flxpa'
!  DO k=1,np
!    write (6,'(2i3,2e15.7)')
!    :  ib,k,fall(2,2,k+1)-fall(2,2,k),fall(2,2,k)
!  ENDDO

   ! MS change for energy conservation  
    IF(ib == 8) THEN
      DO j=1,n
        DO i=1,m
!!!          fdirpar(i,j)=fsdir(i,j)*hk(ib)
!!!          fdifpar(i,j)=fsdif(i,j)*hk(ib)
          fdirpar(i,j)= fdirpar(i,j) + fsdir(i,j)*hk(ib)
          fdifpar(i,j)= fdifpar(i,j) + fsdif(i,j)*hk(ib)
        END DO
      END DO
    ELSE ! MS addition to count the UV radiation hitting the ground
      DO j=1,n
        DO i=1,m
          fdiruv(i,j)= fdiruv(i,j) + fsdir(i,j)*hk(ib)
          fdifuv(i,j)= fdifuv(i,j) + fsdif(i,j)*hk(ib)
        END DO
      END DO
    END IF

  END DO

  RETURN
END SUBROUTINE soluv
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CLDSCALE                   ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cldscale(idim,jdim,                                          &
           m,n,np,cosz,fcld,taucldi,taucldl,ict,icb,                    &
           cc,tauclb,tauclf)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the covers of high, middle, and low clouds and scales
!  the cloud optical thickness.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (a) Radiative Transfer Model: M.-D. Chou and M. Suarez
!          (b) Cloud Optics:Tao, Lang, Simpson, Sui, Ferrier and
!              Chou (1996)
!
!  MODIFICATION:
!
!  03/11/1996 (Yuhe Liu)
!  Formatted code to ARPS standard format
!
!-----------------------------------------------------------------------
!
!  ORIGINAL COMMENTS:
!
!********************************************************************
!
!   This subroutine computes the covers of high, middle, and
!    low clouds and scales the cloud optical thickness.
!
!   To simplify calculations in a cloudy atmosphere, clouds are
!    grouped into high, middle and low clouds separated by the levels
!    ict and icb (level 1 is the top of the atmosphere).
!
!   Within each of the three groups, clouds are assumed maximally
!    overlapped, and the cloud cover (cc) of a group is the maximum
!    cloud cover of all the layers in the group.  The optical thickness
!    (taucld) of a given layer is then scaled to new values (tauclb and
!    tauclf) so that the layer reflectance corresponding to the cloud
!    cover cc is the same as the original reflectance with optical
!    thickness taucld and cloud cover fcld.
!
!---input parameters
!
!    number of grid intervals in zonal direction (m)
!    number of grid intervals in meridional direction (n)
!    maximum number of grid intervals in meridional direction (ndim)
!    number of atmospheric layers (np)
!    cosine of the solar zenith angle (cosz)
!    fractional cloud cover (fcld)
!    cloud optical thickness (taucld)
!    index separating high and middle clouds (ict)
!    index separating middle and low clouds (icb)
!
!---output parameters
!
!    fractional cover of high, middle, and low clouds (cc)
!    scaled cloud optical thickness for beam radiation (tauclb)
!    scaled cloud optical thickness for diffuse radiation (tauclf)
!
!********************************************************************
!
!-----------------------------------------------------------------------
!
!  Variable declarations
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

!-----input parameters

  INTEGER :: idim,jdim
  INTEGER :: m,n,np,ict,icb

  REAL :: cosz(idim,jdim)
  REAL :: fcld(idim,jdim,np)
  REAL :: taucldi(idim,jdim,np)
  REAL :: taucldl(idim,jdim,np)

!-----output parameters

  REAL :: cc(m,n,3)
  REAL :: tauclb(m,n,np)
  REAL :: tauclf(m,n,np)

!-----temporary variables

  INTEGER :: i,j,k,im,it,ia,kk
  REAL :: fm,ft,fa,xai,taux

!-----pre-computed table

  INTEGER :: nm,nt,na
  PARAMETER (nm=11,nt=9,na=11)
  REAL :: dm,dt,da,t1
!!!  REAL :: dm,dt,da,t1,caib(nm,nt,na),caif(nt,na)
  PARAMETER (dm=0.1,dt=0.30103,da=0.1,t1=-0.9031)

!-----include the pre-computed table for cai

!  include "cai.dat"
!  save caib,caif

!!!  COMMON /radtab004/ caib,caif
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!-----clouds within each of the high, middle, and low clouds are
!  assumed maximally overlapped, and the cloud cover (cc)
!  for a group is the maximum cloud cover of all the layers
!  in the group

!!!  print *,'  ict,icb = ',ict,icb

    DO j=1,n
    DO i=1,m
      cc(i,j,1)=0.0
      DO k=1,ict-1
        cc(i,j,1)=MAX(cc(i,j,1),fcld(i,j,k))
      END DO
    END DO
    END DO

    DO j=1,n
    DO i=1,m
      cc(i,j,2)=0.0
      DO k=ict,icb-1
        cc(i,j,2)=MAX(cc(i,j,2),fcld(i,j,k))
      END DO
    END DO
    END DO

    DO j=1,n
    DO i=1,m
      cc(i,j,3)=0.0
      DO k=icb,np
        cc(i,j,3)=MAX(cc(i,j,3),fcld(i,j,k))
      END DO
    END DO
    END DO

!!!  print *,cc(1,1,1),cc(1,1,2),cc(1,1,3)

!-----scale the cloud optical thickness.
!  taucldi(i,j,k) is the optical thickness for ice particles, and
!  taucldl(i,j,k) is the optical thickness for liquid particles.

  DO k=1,np

    IF(k < ict) THEN
      kk=1
    ELSE IF(k >= ict .AND. k < icb) THEN
      kk=2
    ELSE
      kk=3
    END IF

    DO j=1,n
      DO i=1,m

        tauclb(i,j,k) = 0.0
        tauclf(i,j,k) = 0.0

        taux=taucldi(i,j,k)+taucldl(i,j,k)
        IF (taux > 0.05 .AND. fcld(i,j,k) > 0.01) THEN

!-----normalize cloud cover

          fa=fcld(i,j,k)/cc(i,j,kk)

!-----table look-up

          taux=MIN(taux,32.)

          fm=cosz(i,j)/dm
          ft=(LOG10(taux)-t1)/dt
          fa=fa/da

          im=INT(fm+1.5)
          it=INT(ft+1.5)
          ia=INT(fa+1.5)

          im=MAX(im,2)
          it=MAX(it,2)
          ia=MAX(ia,2)

          im=MIN(im,nm-1)
          it=MIN(it,nt-1)
          ia=MIN(ia,na-1)

          fm=fm-FLOAT(im-1)
          ft=ft-FLOAT(it-1)
          fa=fa-FLOAT(ia-1)

!-----scale cloud optical thickness for beam radiation.
!  the scaling factor, xai, is a function of the solar zenith
!  angle, optical thickness, and cloud cover.

          xai=    (-caib(im-1,it,ia)*(1.-fm)+                           &
              caib(im+1,it,ia)*(1.+fm))*fm*.5+caib(im,it,ia)*(1.-fm*fm)

          xai=xai+(-caib(im,it-1,ia)*(1.-ft)+                           &
              caib(im,it+1,ia)*(1.+ft))*ft*.5+caib(im,it,ia)*(1.-ft*ft)

          xai=xai+(-caib(im,it,ia-1)*(1.-fa)+                           &
              caib(im,it,ia+1)*(1.+fa))*fa*.5+caib(im,it,ia)*(1.-fa*fa)

          xai= xai-2.*caib(im,it,ia)
          xai=MAX(xai,0.0)

          tauclb(i,j,k) = taux*xai

!-----scale cloud optical thickness for diffuse radiation.
!  the scaling factor, xai, is a function of the cloud optical
!  thickness and cover but not the solar zenith angle.

          xai=    (-caif(it-1,ia)*(1.-ft)+                              &
              caif(it+1,ia)*(1.+ft))*ft*.5+caif(it,ia)*(1.-ft*ft)

          xai=xai+(-caif(it,ia-1)*(1.-fa)+                              &
              caif(it,ia+1)*(1.+fa))*fa*.5+caif(it,ia)*(1.-fa*fa)

          xai= xai-caif(it,ia)
          xai=MAX(xai,0.0)

          tauclf(i,j,k) = taux*xai

        END IF

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE cldscale
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE FLXCO2                     ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######               University of Oklahoma                 ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE flxco2(m,n,np,swc,swh,csm,df)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the reduction of clear-sky downward solar flux
!  due to co2 absorption.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (a) Radiative Transfer Model: M.-D. Chou and M. Suarez
!          (b) Cloud Optics:Tao, Lang, Simpson, Sui, Ferrier and
!              Chou (1996)
!
!  MODIFICATION HISTORY:
!
!  03/15/1996 (Yuhe Liu)
!  Adopted the original code and formatted it in accordance with the
!  ARPS coding standard.
!
!-----------------------------------------------------------------------
!
!  ORIGINAL COMMENTS:
!
!*****************************************************************
!
!-----compute the reduction of clear-sky downward solar flux
!  due to co2 absorption.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

!-----input parameters

  INTEGER :: m,n,np

  REAL :: csm(m,n)
  REAL :: swc(m,n,np+1)
  REAL :: swh(m,n,np+1)
!!!  REAL :: cah(22,19)

!-----output (undated) parameter

  REAL :: df(m,n,np+1)

!-----temporary variables

  INTEGER :: i,j,k,ic,iw
  REAL :: xx,CLOG,wlog,dc,dw,x1,x2,y2

!********************************************************************
!-----include co2 look-up table
!
!  include "cah.dat"
!  save cah

!!!  COMMON /radtab005/ cah
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!********************************************************************
!-----table look-up for the reduction of clear-sky solar
!  radiation due to co2. The fraction 0.0343 is the
!  extraterrestrial solar flux in the co2 bands.

!  write (6,'(a/3a3,7a10)') 'Flux reduction due to co2, i=j=3',
!    :'k','ic','iw','x1','y2','x1-y2','ddf','df','swc','swh'

  DO k= 2, np+1
    DO j= 1, n
      DO i= 1, m
        xx=1./.3
        CLOG=LOG10(swc(i,j,k)*csm(i,j))
        wlog=LOG10(swh(i,j,k)*csm(i,j))
        ic=INT( (CLOG+3.15)*xx+1.)
        iw=INT( (wlog+4.15)*xx+1.)
        IF(ic < 2)ic=2
        IF(iw < 2)iw=2
        IF(ic > 22)ic=22
        IF(iw > 19)iw=19
        dc=CLOG-FLOAT(ic-2)*.3+3.
        dw=wlog-FLOAT(iw-2)*.3+4.
        x1=cah(1,iw-1)+(cah(1,iw)-cah(1,iw-1))*xx*dw
        x2=cah(ic-1,iw-1)+(cah(ic-1,iw)-cah(ic-1,iw-1))*xx*dw
        y2=x2+(cah(ic,iw-1)-cah(ic-1,iw-1))*xx*dc
        IF (x1 < y2) x1=y2
        df(i,j,k)=df(i,j,k)+0.0343*(x1-y2)
      END DO
    END DO

!    write (6,'(3i3,7f10.5)')
!    :  k,ic,iw,x1,y2,x1-y2,
!    :  df(3,3,k)-df(3,3,k-1),df (3,3,k),swc(3,3,k),swh(3,3,k)

  END DO

  RETURN
END SUBROUTINE flxco2
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CLDFLX                     ######
!######                                                      ######
!######                     Developed by                     ######
!######                                                      ######
!######    Goddard Cumulus Ensemble Modeling Group, NASA     ######
!######                                                      ######
!######     Center for Analysis and Prediction of Storms     ######
!######               University of Oklahoma                 ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cldflx(m,n,np,ict,icb,cc,rr,tt,td,rs,ts,                     &
           fclr,fall,fsdir,fsdif,                                       &
           ch,cm,ct,fdndir,fdndif,                                      &
           flxdn,                                                       &
           rra,tta,tda,rsa,rxa)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate upward and downward fluxes using a two-stream adding
!  method.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: (a) Radiative Transfer Model: M.-D. Chou and M. Suarez
!          (b) Cloud Optics:Tao, Lang, Simpson, Sui, Ferrier and
!              Chou (1996)
!
!  MODIFICATION HISTORY:
!
!  03/15/1996 (Yuhe Liu)
!  Adopted the original code and formatted it in accordance with the
!  ARPS coding standard.
!
!-----------------------------------------------------------------------
!
!  ORIGINAL COMMENTS:
!
!*******************************************************************
!  compute upward and downward fluxes using a two-stream adding method
!  following equations (3)-(5) of Chou (1992, JAS).
!
!  clouds are grouped into high, middle, and low clouds which are
!  assumed randomly overlapped. It involves eight sets of calculations.
!  In each set of calculations, each atmospheric layer is homogeneous,
!  either with clouds or without clouds.

!  input parameters:
!  m:   number of soundings in zonal direction
!  n:   number of soundings in meridional direction
!  np:  number of atmospheric layers
!  ict: the level separating high and middle clouds
!  icb: the level separating middle and low clouds
!  cc:  effective cloud covers for high, middle and low clouds
!  tt:  diffuse transmission of a layer illuminated by beam radiation
!  td:  direct beam tranmssion
!  ts:  transmission of a layer illuminated by diffuse radiation
!  rr:  reflection of a layer illuminated by beam radiation
!  rs:  reflection of a layer illuminated by diffuse radiation
!
!  output parameters:
!  fclr:  clear-sky flux (downward minus upward)
!  fall:  all-sky flux (downward minus upward)
!  fdndir: surface direct downward flux
!  fdndif: surface diffuse downward flux
!*********************************************************************c
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

!-----input parameters

  INTEGER :: m,n,np,ict,icb

  REAL :: rr(m,n,np+1,2)
  REAL :: tt(m,n,np+1,2)
  REAL :: td(m,n,np+1,2)
  REAL :: rs(m,n,np+1,2)
  REAL :: ts(m,n,np+1,2)
  REAL :: cc(m,n,3)

!-----output parameters

  REAL :: fclr(m,n,np+1)
  REAL :: fall(m,n,np+1)
  REAL :: fsdir(m,n)
  REAL :: fsdif(m,n)

!-----temporary array

  REAL :: ch(m,n)
  REAL :: cm(m,n)
  REAL :: ct(m,n)
  REAL :: fdndir(m,n)
  REAL :: fdndif(m,n)

  REAL :: flxdn(m,n,np+1)

  REAL :: rra(m,n,np+1,2,2)
  REAL :: tta(m,n,np+1,2,2)
  REAL :: tda(m,n,np+1,2,2)
  REAL :: rsa(m,n,np+1,2,2)
  REAL :: rxa(m,n,np+1,2,2)

!-----temporary variables

  INTEGER :: i,j,k,ih,im,is
  REAL :: fupdif
  REAL :: denm,xx
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

    jloop:  DO j=1,n
    iloop:  DO i=1,m

!-----initialize all-sky flux (fall) and surface downward fluxes

  DO k=1,np+1
        fall(i,j,k)=0.0
  END DO

      fsdir(i,j)=0.0
      fsdif(i,j)=0.0

!-----compute transmittances and reflectances for a composite of
!  layers. layers are added one at a time, going down from the top.
!  tda is the composite transmittance illuminated by beam radiation
!  tta is the composite diffuse transmittance illuminated by
!      beam radiation
!  rsa is the composite reflectance illuminated from below
!      by diffuse radiation
!  tta and rsa are computed from eqs. (4b) and (3b) of Chou

!-----for high clouds. indices 1 and 2 denote clear and cloudy
!  situations, respectively.

  DO ih=1,2

        tda(i,j,1,ih,1)=td(i,j,1,ih)
        tta(i,j,1,ih,1)=tt(i,j,1,ih)
        rsa(i,j,1,ih,1)=rs(i,j,1,ih)
        tda(i,j,1,ih,2)=td(i,j,1,ih)
        tta(i,j,1,ih,2)=tt(i,j,1,ih)
        rsa(i,j,1,ih,2)=rs(i,j,1,ih)
        DO k= 2, ict-1
          denm = ts(i,j,k,ih)/( 1.-rsa(i,j,k-1,ih,1)*rs(i,j,k,ih))
          tda(i,j,k,ih,1)= tda(i,j,k-1,ih,1)*td(i,j,k,ih)
          tta(i,j,k,ih,1)= tda(i,j,k-1,ih,1)*tt(i,j,k,ih)               &
                        +(tda(i,j,k-1,ih,1)*rr(i,j,k,ih)                &
                        *rsa(i,j,k-1,ih,1)+tta(i,j,k-1,ih,1))*denm
          rsa(i,j,k,ih,1)= rs(i,j,k,ih)+ts(i,j,k,ih)                    &
                        *rsa(i,j,k-1,ih,1)*denm
          tda(i,j,k,ih,2)= tda(i,j,k,ih,1)
          tta(i,j,k,ih,2)= tta(i,j,k,ih,1)
          rsa(i,j,k,ih,2)= rsa(i,j,k,ih,1)
        END DO

!-----for middle clouds

    DO im=1,2

      DO k= ict, icb-1
            denm = ts(i,j,k,im)/( 1.-rsa(i,j,k-1,ih,im)*rs(i,j,k,im))
            tda(i,j,k,ih,im)= tda(i,j,k-1,ih,im)*td(i,j,k,im)
            tta(i,j,k,ih,im)= tda(i,j,k-1,ih,im)*tt(i,j,k,im)           &
                          +(tda(i,j,k-1,ih,im)*rr(i,j,k,im)             &
                          *rsa(i,j,k-1,ih,im)+tta(i,j,k-1,ih,im))*denm
            rsa(i,j,k,ih,im)= rs(i,j,k,im)+ts(i,j,k,im)                 &
                          *rsa(i,j,k-1,ih,im)*denm
      END DO

    END DO
  END DO

!-----layers are added one at a time, going up from the surface.
!  rra is the composite reflectance illuminated by beam radiation
!  rxa is the composite reflectance illuminated from above
!      by diffuse radiation
!  rra and rxa are computed from eqs. (4a) and (3a) of Chou

!-----for the low clouds

  DO is=1,2

        rra(i,j,np+1,1,is)=rr(i,j,np+1,is)
        rxa(i,j,np+1,1,is)=rs(i,j,np+1,is)
        rra(i,j,np+1,2,is)=rr(i,j,np+1,is)
        rxa(i,j,np+1,2,is)=rs(i,j,np+1,is)
        DO k=np,icb,-1
          denm=ts(i,j,k,is)/( 1.-rs(i,j,k,is)*rxa(i,j,k+1,1,is) )
          rra(i,j,k,1,is)=rr(i,j,k,is)+(td(i,j,k,is)                    &
              *rra(i,j,k+1,1,is)+tt(i,j,k,is)*rxa(i,j,k+1,1,is))*denm
          rxa(i,j,k,1,is)= rs(i,j,k,is)+ts(i,j,k,is)                    &
              *rxa(i,j,k+1,1,is)*denm
          rra(i,j,k,2,is)=rra(i,j,k,1,is)
          rxa(i,j,k,2,is)=rxa(i,j,k,1,is)
        END DO

!-----for middle clouds

    DO im=1,2

      DO k= icb-1,ict,-1
            denm=ts(i,j,k,im)/( 1.-rs(i,j,k,im)*rxa(i,j,k+1,im,is) )
            rra(i,j,k,im,is)= rr(i,j,k,im)+(td(i,j,k,im)                &
                *rra(i,j,k+1,im,is)+tt(i,j,k,im)*rxa(i,j,k+1,im,is))*denm
            rxa(i,j,k,im,is)= rs(i,j,k,im)+ts(i,j,k,im)                 &
                *rxa(i,j,k+1,im,is)*denm
      END DO

    END DO
  END DO

!-----integration over eight sky situations.
!  ih, im, is denotes high, middle and low cloud groups.

  DO ih=1,2

!-----clear portion

    IF(ih == 1) THEN
          ch(i,j)=1.0-cc(i,j,1)

    ELSE

!-----cloudy portion

          ch(i,j)=cc(i,j,1)

    END IF

    DO im=1,2

!-----clear portion

      IF(im == 1) THEN

            cm(i,j)=ch(i,j)*(1.0-cc(i,j,2))

      ELSE

!-----cloudy portion

            cm(i,j)=ch(i,j)*cc(i,j,2)

      END IF

      DO is=1,2

!-----clear portion

        IF(is == 1) THEN

              ct(i,j)=cm(i,j)*(1.0-cc(i,j,3))

        ELSE

!-----cloudy portion

              ct(i,j)=cm(i,j)*cc(i,j,3)

        END IF

!-----add one layer at a time, going down.

        DO k= icb, np
              denm = ts(i,j,k,is)/( 1.-rsa(i,j,k-1,ih,im)*rs(i,j,k,is) )
              tda(i,j,k,ih,im)= tda(i,j,k-1,ih,im)*td(i,j,k,is)
              tta(i,j,k,ih,im)=  tda(i,j,k-1,ih,im)*tt(i,j,k,is)        &
                   +(tda(i,j,k-1,ih,im)*rr(i,j,k,is)                    &
                   *rsa(i,j,k-1,ih,im)+tta(i,j,k-1,ih,im))*denm
              rsa(i,j,k,ih,im)= rs(i,j,k,is)+ts(i,j,k,is)               &
                   *rsa(i,j,k-1,ih,im)*denm
        END DO

!-----add one layer at a time, going up.

        DO k= ict-1,1,-1
              denm =ts(i,j,k,ih)/(1.-rs(i,j,k,ih)*rxa(i,j,k+1,im,is))
              rra(i,j,k,im,is)= rr(i,j,k,ih)+(td(i,j,k,ih)              &
                  *rra(i,j,k+1,im,is)+tt(i,j,k,ih)*rxa(i,j,k+1,im,is))*denm
              rxa(i,j,k,im,is)= rs(i,j,k,ih)+ts(i,j,k,ih)               &
                  *rxa(i,j,k+1,im,is)*denm
        END DO

!-----compute fluxes following eq (5) of Chou (1992)

!  fdndir is the direct  downward flux
!  fdndif is the diffuse downward flux
!  fupdif is the diffuse upward flux

        DO k=2,np+1
              denm= 1./(1.- rxa(i,j,k,im,is)*rsa(i,j,k-1,ih,im))
              fdndir(i,j)= tda(i,j,k-1,ih,im)
              xx = tda(i,j,k-1,ih,im)*rra(i,j,k,im,is)
              fdndif(i,j)= (xx*rsa(i,j,k-1,ih,im)+tta(i,j,k-1,ih,im))*denm
              fupdif= (xx+tta(i,j,k-1,ih,im)*rxa(i,j,k,im,is))*denm
              flxdn(i,j,k)=fdndir(i,j)+fdndif(i,j)-fupdif
        END DO

            flxdn(i,j,1)=1.0-rra(i,j,1,im,is)

!-----summation of fluxes over all (eight) sky situations.

        DO k=1,np+1
              IF(ih == 1 .AND. im == 1 .AND. is == 1) THEN
                fclr(i,j,k)=flxdn(i,j,k)
              END IF
              fall(i,j,k)=fall(i,j,k)+flxdn(i,j,k)*ct(i,j)
        END DO

            fsdir(i,j)=fsdir(i,j)+fdndir(i,j)*ct(i,j)
            fsdif(i,j)=fsdif(i,j)+fdndif(i,j)*ct(i,j)

      END DO
    END DO
  END DO

  END DO  iloop
  END DO  jloop

  RETURN
END SUBROUTINE cldflx

  END MODULE sorad3d_module
