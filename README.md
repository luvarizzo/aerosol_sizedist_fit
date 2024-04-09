# aerosol_sizedist_fit
Fit one to three lognormal curves to submicrometer aerosol size distribution samples, using Ordinary Least Squares minization.
It can be run in the fully automatic mode, so that the algorithm decides between one, two or three modes.
As output, it provides the lognormal parameters for each fitted mode: N (total number concentration), Dpg (geometric mean diameter), sg (geometric standard deviation).
As an option, the size distribution can be extrapolated for greater diameter based on the fitted lognormals.
Extrapolation is only recommended up to 600-700 nm.

Funcionalities: 
* fully automatic mode option, without the need of user intervention
* manual mode option, in which the user marks the start guess of mode diameters in the size distribution plot
* partially automatic mode option, in which the user is prompted for a decision only if the algorithm is not succesfull to find a good fitting
* option to interpolate to a common diameter sequence (shape-preserving interpolation)
* option to extrapolate to greater diameters (based on the lognormal fits)
* automatically correct fit order if Dp1>Dp2 or Dp2>Dp3
* Uses the last sample fitting parameters as a first guess to the current sample
* automatically classifies fitted modes into: nucleation(Dpg<30nm); Aitken(30<=Dpg<=90nm); Accum1(90<Dpg<160nm); Accum2(Dpg>=160nm)
* saves the sample indexes that eventually did not obtain a satisfactory fit in the automatic mode

Main script: smps_fit_lognormal_v8.m (matlab)
Required functions: ffitgaussian3_notoolbox.m; rsquared.m; rmse.m (built-in matlab function since version R2022b) 

INPUT: 
data (1st column=index or timeline; 1st row=diam; rows=dN/dlogDp)
extradiam (extra diameter array to extrapolate smps data - optional)   

OUTPUT: 
extrasmps (extrapolated dN/dlogDp values for the defined extradiam array)
smpsnew (if interpolation was made, this output holds the new interpolated dN/dlogDp values)
fits col 1-7: matlab time and datevec
	col 8-16: fit parameters (N1,N2,N3,Dpg1,Dpg2,Dpg3,sg1,sg2,sg3)
	col 17-21: goodness of fitting (sse,rsquare,dfe,adjrsquare,rmse)           
	col 22-23: Ndata | Nfit (total concentration comparison)
fits_into_5modes col 1-7: matlab time and datevec
	col 8-22: fit parameters (Nnucl,Naitken1,Naitken2,Naccum1,Naccum2,Dpg1,Dpg2,Dpg3,Dpg4,Dpg5,sg1,sg2,sg3,sg4,sg5)
	col 23-27: goodness of fitting (sse,rsquare,dfe,adjrsquare,rmse)           
	col 28-29: Ndata | Nfit (total concentration comparison)
rever: sample indexes that did not obtain a satisfactory fit  

Example input file: smps_example.csv
DateTime,index,10.2,10.6,10.9,11.3,11.8,12.2,12.6,13.1,13.6,14.1,14.6,15.1,15.7,16.3,16.8,17.5,18.1,18.8,19.5,20.2,20.9,21.7,22.5,23.3,24.1,25,25.9,26.9,27.9,28.9,30,31.1,32.2,33.4,34.6,35.9,37.2,38.5,40,41.4,42.9,44.5,46.1,47.8,49.6,51.4,53.3,55.2,57.3,59.4,61.5,63.8,66.1,68.5,71,73.7,76.4,79.1,82,85.1,88.2,91.4,94.7,98.2,101.8,105.5,109.4,113.4,117.6,121.9,126.3,131,135.8,140.7,145.9,151.2,156.8,162.5,168.5,174.7,181.1,187.7,194.6,201.7,209.1,216.7,224.7,232.9,241.4,250.3,259.5,269,278.8,289,299.6,310.6,322,333.8,346,358.7,371.8,385.4,399.5,414.2
03-May-2023 15:25:00,1,508.035,399.722,605.231,839.615,526.408,603.502,1036.964,892.617,1064.452,1170.194,1071.165,1425.32,1199.279,1157.385,978.889,1228.81,1585.92,1475.544,1931.611,1976.652,1953.999,2346.611,1750.419,2152.154,2917.031,2806.223,3442.247,3150.701,3221.422,3683.421,3589.016,3793.668,4191.613,4130.132,4294.928,4623.143,5046.053,4703.416,5256.329,5030.191,5190.476,5789.725,5873.153,5933.73,5982.692,6794.277,6830.97,7143.356,7163.835,7299.87,7880.72,7875.907,7589.922,8129.306,8400.353,8120.7,7996.731,7446.867,7603.501,8549.805,7624.743,7742.944,7547.694,6866.928,7178.879,7284.856,6087.694,6116.203,6261.514,5857.721,5319.979,5275.224,5320.995,4587.617,4801.232,3920.668,4105.797,4209.254,3827.644,3747.721,3363.826,3623.737,3515.688,3339.827,3108.362,2985.807,3123.145,2744.225,2756.175,2325.416,2794.835,2582.369,2375.429,2334.544,2078.949,2050.8,1826.369,1964.166,1514.322,1335.702,1095.643,1152.571,839.518,874.643
03-May-2023 15:30:00,2,964.1405,652.2465,549.738,623.6395,668.3835,588.272,623.3035,839.18,723.617,615.08,728.066,836.905,808.5925,862.494,1275.333,935.2,1136.393,1428.471,1683.3595,1664.3625,1807.0745,1882.508,2012.1605,2288.736,2385.6635,2568.55,2843.911,2965.9645,2872.112,3463.612,3594.3025,3600.7685,3747.9375,3872.797,4126.073,4188.2815,4466.9325,4393.889,4926.1885,5025.0215,5192.139,5738.877,5601.4835,5958.5195,6341.4655,6398.4025,6916.0045,7245.7255,6941.8875,7263.741,7467.1845,7718.7655,7695.7345,7913.064,7859.43,7682.857,7780.427,7643.5035,7631.5155,7553.556,7503.208,7387.364,7091.37,6796.338,6701.315,6479.2455,6216.996,5742.8035,5870.7175,5508.0385,5536.4495,5249.51,4744.809,4507.1005,4635.9675,4674.739,4197.814,3970.5395,3822.662,3750.646,3480.7025,3353.4555,3354.162,3176.1205,3229.0935,3104.709,2956.411,2750.89,2522.5405,2461.9905,2551.9015,2521.078,2280.202,2232.5335,2058.6185,1956.8965,1797.727,1491.2515,1406.4,1208.555,1107.7125,975.4565,996.829,764.2955
...
