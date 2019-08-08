#include <ConstField.h>
#include <Exception.h>
#include <FieldManager.h>
#include <KalmanFitterRefTrack.h>
#include <StateOnPlane.h>
#include <Track.h>
#include <TrackPoint.h>

#include <MaterialEffects.h>
#include <RKTrackRep.h>
#include <TGeoMaterialInterface.h>

#include <EventDisplay.h>

#include <PlanarMeasurement.h>

#include <TEveManager.h>
#include <TApplication.h>
#include <TGeoManager.h>
#include <TVector3.h>
#include <TRandom3.h>

#include <vector>
#include <iomanip>

#include "TDatabasePDG.h"
#include <TMath.h>


int main(int argc, char* argv[]) {


  //obtained with the following commands:
  // cd /wd/dune-it/enurec/analysis/kloe-simu
  // $$ .> myMinimalFit/ev_312_mu.txt
  // $$ tEvent->Scan("particles.pdg:Entry$:particles.tr.@digits.size():particles.tr.digits.x:particles.tr.digits.y:particles.tr.digits.z:particles.tr.digits.t","particles.pdg==13&&particles.tr.@digits.size()>200&&particles.primary==1&&z>10500&&Entry$==312")
  // $$ .>
  // x=$(tail -n +5 ev_312_mu.txt | head -n -5 | sed 's/*/ /g' | awk -F' ' '{print $6}'); for i in ${x[@]}; do echo -n "$i, "; done; echo ""
  // px = 751.35641
  // py = -669.0994
  // pz = 1104.7786
  // E  = 1497.9753
  // mass = 105.658
  // charge = -1
  
  // energy in GeV
  // distance in cm
  // magnetic field in kG
  
  // magnetic field is not read from geometry but should be provided explicitly
  // geometry should be in cm and g/cm^3 
  
  std::vector<double> x = {-1019.879, -1017.072, -1003.659, -1002.137, -1000.758, -999.2351, -985.0332, -982.2240, -968.6766, -965.7839, -950.1326, -947.3119, -934.1577, -932.6006, -929.7042, -915.0920, -912.2677, -898.5721, -895.6770, -879.9350, -877.0973, -863.4333, -860.5437, -858.9457, -844.6591, -841.8112, -827.1855, -824.2979, -809.2535, -806.3984, -791.3922, -788.5032, -773.7134, -770.8464, -769.1889, -756.1522, -753.3024, -738.2267, -736.5815, -733.7234, -720.1121, -717.2795, -701.0916, -698.9733, -684.9233, -682.0896, -680.4088, -665.4463, -662.5640, -648.6578, -645.8191, -644.1105, -629.5529, -626.6593, -611.1154, -608.2823, -593.5360, -591.8283, -590.6287, -588.9192, -575.7867, -574.0403, -571.2061, -555.6168, -552.7008, -537.4259, -534.5908, -519.2726, -516.3478, -501.3887, -498.5832, -496.8149, -481.1737, -478.3804, -464.3310, -461.5261, -459.7379, -444.6799, -441.7383, -425.9742, -423.1744, -408.0259, -406.2568, -403.3012, -388.2627, -385.4834, -369.4294, -366.4632, -351.2290, -348.4506, -347.4489, -330.5770, -327.5904, -316.1625, -310.0667, -308.1883, -293.3163, -291.4895, -288.4899, -273.2314, -270.4691, -254.0602, -251.0474, -249.1994, -234.3256, -231.5670, -214.5218, -211.4910, -196.0700, -194.1356, -191.3935, -174.8015, -171.7587, -156.6836, -154.7326, -151.9965, -134.8865, -131.8302, -116.1126, -114.1459, -111.4247, -94.79738, -91.72840, -76.46022, -74.47559, -71.77268, -54.48837, -51.40321, -35.68796, -33.68377, -31.21365, -13.91327, -10.82190, -8.863436, 6.3087737, 9.0081423, 11.025019, 26.761451, 28.746340, 31.859100, 47.406669, 49.442898, 50.093229, 52.130775, 69.719505, 72.846244, 89.494694, 92.174690, 94.230815, 110.92448, 114.07056, 116.10935, 132.61066, 135.28718, 154.48460, 157.65383, 174.64837, 177.31027, 179.40447, 196.25464, 198.16651, 201.53122, 217.38563, 219.49520, 222.14005, 240.35767, 243.55372, 245.66949, 260.67720, 262.80617, 265.44096, 284.75899, 287.97286, 304.51504, 306.65899, 309.27773, 311.42510, 329.48096, 332.71973, 350.86953, 353.48619, 355.66203, 374.65299, 377.91686, 395.32348, 397.51789, 400.12762, 420.18068, 423.46695, 425.73132, 441.81222, 444.02897, 446.63118, 466.09027, 468.38966, 471.70008, 489.92481, 492.51388, 494.75747, 514.76709, 518.10201, 536.75264, 539.01590, 541.59401, 561.52231, 563.90903, 567.2668, 585.77923, 588.06542, 590.62511, 611.09534, 613.53746, 616.92557, 635.47405, 637.78464, 640.32410, 642.63652, 661.22663, 663.71471, 667.12602, 687.34422, 689.86816, 692.19730, 713.77642, 717.76994, 736.97037, 739.31788, 741.82685, 744.17574};
  
  std::vector<double> y = {3943.5869, 3941.069, 3928.9983, 3927.6229, 3926.3758, 3924.9968, 3912.1187, 3909.5647, 3897.2023, 3894.5542, 3880.1955, 3877.6008, 3865.4629, 3864.0210, 3861.3371, 3847.7556, 3845.1166, 3832.2911, 3829.5740, 3814.7524, 3812.066, 3799.0956, 3796.3433, 3794.8198, 3781.1620, 3778.4286, 3764.3536, 3761.5668, 3747.0098, 3744.2395, 3729.6419, 3726.8219, 3712.3306, 3709.5072, 3707.8752, 3694.9329, 3692.0783, 3676.9264, 3675.2672, 3672.3828, 3658.5726, 3655.6787, 3639.0691, 3636.9055, 3622.3983, 3619.4693, 3617.7305, 3602.2008, 3599.2036, 3584.7152, 3581.7498, 3579.9629, 3564.7027, 3561.6587, 3545.2484, 3542.2494, 3526.6009, 3524.7836, 3523.5071, 3521.6879, 3507.69, 3505.8253, 3502.7972, 3486.1163, 3482.9884, 3466.5745, 3463.5214, 3446.9663, 3443.7898, 3427.4804, 3424.4012, 3422.4591, 3405.2355, 3402.4874, 3386.6113, 3383.5044, 3381.5213, 3364.7947, 3361.5215, 3343.9263, 3340.7864, 3323.7, 3321.6913, 3318.3338, 3301.1918, 3298.0064, 3279.5568, 3276.1356, 3258.5328, 3255.3164, 3254.2074, 3234.5820, 3231.1089, 3213.8342, 3210.6727, 3208.4743, 3191.0182, 3188.8699, 3185.3385, 3167.3006, 3164.0211, 3144.4499, 3140.8480, 3138.6371, 3120.8120, 3117.4983, 3096.9613, 3093.2968, 3074.5931, 3072.2377, 3068.8964, 3048.6274, 3044.9060, 3026.4247, 3024.0269, 3020.6629, 2999.5582, 2995.7741, 2976.2601, 2973.8121, 2970.4214, 2949.6230, 2945.7623, 2926.4921, 2923.9812, 2920.5591, 2898.6397, 2894.7150, 2874.6727, 2872.1128, 2868.9830, 2846.8630, 2842.9188, 2840.4178, 2821.0375, 2817.5857, 2815.0053, 2794.8273, 2792.2732, 2788.2655, 2768.2054, 2765.5678, 2764.7247, 2762.0842, 2739.2369, 2735.1674, 2713.4530, 2709.9518, 2707.2641, 2685.4101, 2681.2794, 2678.6001, 2656.8686, 2653.3439, 2628.0086, 2623.8217, 2601.2929, 2597.7494, 2594.9595, 2572.4827, 2570.8117, 2565.4273, 2544.1360, 2541.2907, 2537.7217, 2513.0749, 2508.7324, 2505.8554, 2485.4231, 2482.5268, 2478.9385, 2452.5577, 2448.1504, 2425.4375, 2422.4861, 2418.8798, 2415.9209, 2390.9930, 2386.5113, 2361.3614, 2357.7364, 2354.7201, 2328.3303, 2323.7919, 2299.5425, 2296.4815, 2292.8396, 2264.7820, 2260.1747, 2256.9982, 2234.4318, 2231.3177, 2227.6606, 2200.2424, 2196.9932, 2192.3111, 2166.5048, 2162.8314, 2159.6460, 2131.1880, 2126.4258, 2099.7579, 2096.5208, 2092.8321, 2064.2335, 2060.7969, 2055.9596, 2029.2348, 2025.9295, 2022.2259, 1992.5217, 1988.9754, 1984.0537, 1957.0293, 1953.6482, 1949.9306, 1946.5434, 1919.2044, 1915.5379, 1910.5088, 1880.6572, 1876.9321, 1873.4924, 1841.5402, 1835.6037, 1807.0685, 1803.5772, 1799.8439, 1796.3469};
  
  std::vector<double> z = {10792.080, 10796.201, 10815.877, 10818.109, 10820.130, 10822.363, 10843.185, 10847.304, 10867.156, 10871.392, 10894.253, 10898.363, 10917.512, 10919.777, 10923.991, 10945.225, 10949.332, 10969.194, 10973.382, 10996.117, 11000.213, 11019.917, 11024.079, 11026.380, 11046.925, 11051.012, 11071.982, 11076.119, 11097.640, 11101.721, 11123.125, 11127.236, 11148.256, 11152.328, 11154.681, 11173.283, 11177.364, 11198.921, 11201.270, 11205.349, 11224.787, 11228.835, 11251.928, 11254.923, 11274.930, 11278.946, 11281.327, 11302.479, 11306.538, 11326.085, 11330.067, 11332.463, 11352.843, 11356.891, 11378.616, 11382.568, 11403.108, 11405.479, 11407.144, 11409.515, 11427.721, 11430.139, 11434.062, 11455.608, 11459.637, 11480.689, 11484.587, 11505.633, 11509.654, 11530.200, 11534.068, 11536.505, 11558.017, 11561.819, 11581.084, 11584.921, 11587.367, 11607.927, 11611.933, 11633.361, 11637.163, 11657.701, 11660.092, 11664.085, 11684.378, 11688.126, 11709.737, 11713.720, 11734.115, 11737.822, 11739.033, 11761.597, 11765.560, 11784.514, 11788.763, 11791.245, 11810.885, 11813.293, 11817.244, 11837.312, 11840.937, 11862.443, 11866.381, 11868.794, 11888.158, 11891.734, 11913.792, 11917.712, 11937.621, 11940.116, 11943.653, 11965.007, 11968.915, 11988.239, 11990.736, 11994.237, 12016.073, 12019.966, 12039.957, 12042.456, 12045.914, 12067.002, 12070.882, 12090.148, 12092.648, 12096.052, 12117.753, 12121.616, 12141.247, 12143.746, 12146.789, 12168.291, 12172.146, 12174.586, 12193.437, 12196.785, 12199.283, 12218.733, 12221.179, 12225.013, 12244.104, 12246.600, 12247.397, 12249.892, 12271.392, 12275.208, 12295.471, 12298.722, 12301.215, 12321.396, 12325.189, 12327.646, 12347.486, 12350.688, 12373.595, 12377.363, 12397.536, 12400.695, 12403.177, 12423.097, 12425.478, 12429.314, 12447.980, 12450.456, 12453.558, 12474.856, 12478.590, 12481.060, 12498.529, 12500.998, 12504.050, 12526.370, 12530.082, 12549.096, 12551.554, 12554.555, 12557.013, 12577.618, 12581.299, 12601.806, 12604.750, 12607.195, 12628.472, 12632.119, 12651.501, 12653.937, 12656.831, 12678.998, 12682.613, 12685.102, 12702.719, 12705.143, 12707.985, 12729.156, 12731.648, 12735.232, 12754.851, 12757.628, 12760.031, 12781.381, 12784.925, 12804.666, 12807.055, 12809.773, 12830.708, 12833.205, 12836.715, 12855.988, 12858.355, 12861.003, 12882.093, 12884.592, 12888.054, 12906.934, 12909.274, 12911.843, 12914.181, 12932.942, 12935.442, 12938.865, 12959.072, 12961.586, 12963.901, 12985.279, 12989.218, 13008.080, 13010.374, 13012.824, 13015.117};
  
  std::vector<double> t = {1.0571961, 1.0760598, 1.1679143, 1.1827469, 1.1878996, 1.2001947, 1.2892225, 1.3080367, 1.3994331, 1.4191844, 1.5219752, 1.5407845, 1.6364621, 1.6413727, 1.6599336, 1.7554943, 1.7744030, 1.8666495, 1.8865678, 1.9898972, 2.0090058, 2.1022836, 2.1228283, 2.1349807, 2.2253611, 2.2447852, 2.3417562, 2.3609788, 2.4621326, 2.4820419, 2.5807631, 2.5998324, 2.7007079, 2.7216994, 2.7330108, 2.8168985, 2.8369912, 2.9417407, 2.9514838, 2.9690946, 3.0590212, 3.0784134, 3.1890366, 3.2114002, 3.2989581, 3.3203033, 3.3311757, 3.4285264, 3.4477217, 3.5428695, 3.5638069, 3.5771722, 3.6703554, 3.690503, 3.7960436, 3.8140965, 3.9160320, 3.9301608, 3.9387326, 3.9469126, 4.0358145, 4.0492360, 4.0654755, 4.1688241, 4.1880457, 4.2927326, 4.3105461, 4.4137240, 4.4345679, 4.5349133, 4.5558914, 4.5698769, 4.6725582, 4.6912246, 4.7862347, 4.8078071, 4.8196333, 4.9178491, 4.9383674, 5.0452294, 5.0636403, 5.1744994, 5.1801219, 5.1987766, 5.3010700, 5.3193995, 5.4284250, 5.4496418, 5.5531931, 5.5752078, 5.5867649, 5.6921819, 5.7120071, 5.8136339, 5.8359714, 5.8465455, 5.9493741, 5.9604156, 5.9786528, 6.0820125, 6.1004747, 6.2125696, 6.2351478, 6.2494643, 6.3463979, 6.3652641, 6.4802779, 6.5020414, 6.6091126, 6.6227242, 6.6385349, 6.7503601, 6.7717931, 6.8805022, 6.8892538, 6.9064741, 7.0221413, 7.0435904, 7.1523852, 7.1706660, 7.1840274, 7.2957402, 7.3175279, 7.4289260, 7.4368191, 7.4542510, 7.5717775, 7.5943676, 7.7063741, 7.7162122, 7.7333120, 7.8506469, 7.8752906, 7.8864367, 7.9883814, 8.0100698, 8.0224783, 8.1330561, 8.1428324, 8.1621328, 8.2697273, 8.2905312, 8.2948922, 8.3020347, 8.4203870, 8.4418335, 8.5569657, 8.5788030, 8.5909103, 8.7032839, 8.7281055, 8.7407825, 8.8521424, 8.8701120, 8.9998730, 9.0211081, 9.1385979, 9.1603865, 9.1725524, 9.2874869, 9.3046522, 9.3223288, 9.4318309, 9.4500244, 9.4635764, 9.5868331, 9.6104549, 9.6315501, 9.7302536, 9.7442387, 9.7602816, 9.8911321, 9.9134615, 10.033816, 10.043646, 10.061403, 10.082833, 10.19825, 10.220077, 10.346103, 10.367050, 10.380377, 10.507257, 10.529665, 10.652810, 10.666915, 10.683201, 10.818353, 10.842464, 10.864614, 10.969637, 10.985534, 11.000895, 11.134061, 11.153592, 11.171258, 11.296464, 11.317640, 11.330412, 11.463982, 11.487317, 11.61785, 11.633489, 11.649261, 11.788417, 11.801613, 11.822219, 11.9516, 11.969603, 11.983827, 12.123790, 12.143051, 12.161648, 12.290817, 12.308405, 12.323077, 12.346854, 12.466844, 12.484270, 12.504061, 12.644560, 12.66514, 12.678277, 12.827320, 12.850017, 12.988072, 12.999994, 13.019042, 13.033907};
  
  // init geometry and mag. field
  TGeoManager::Import("geo_v12.root");
  //TApplication* gApp = new TApplication("myMinimalFit",&argc,argv);
  //TEveManager* gEve = new TEveManager(1000, 500);
  
  const int dbgLvl = 0;

  const int pdg = 13;
  
  const double mm2cm = 0.1;
  const double kG2T = 0.1;
  const double m2cm = 100.;

  const double px = 0.7513564;
  const double py = -0.669099;
  const double pz = 1.1047786;

  TDatabasePDG pdgdb;
  pdgdb.ReadPDGTable();

  const double mass = pdgdb.GetParticle(pdg)->Mass();
  const int Z = pdgdb.GetParticle(pdg)->Charge()/3;
  const int sign = Z/TMath::Abs(Z);

  std::cout << "pdg   :" << std::setw(25) << pdg << std::endl;
  std::cout << "mass  :" << std::setw(25) << mass << std::endl;
  std::cout << "charge:" << std::setw(25) << Z << std::endl;

  const double x0 = x[0]*mm2cm;
  const double y0 = y[0]*mm2cm;
  const double z0 = z[0]*mm2cm;
  const double t0 = t[0];

  std::cout << "x0    :" << std::setw(25) << x0 << std::endl;
  std::cout << "y0    :" << std::setw(25) << y0 << std::endl;
  std::cout << "z0    :" << std::setw(25) << z0 << std::endl;
  std::cout << "t0    :" << std::setw(25) << t0 << std::endl;

  const double Bx = 6.; // kGauss

  std::cout << "Bx    :" << std::setw(25) << Bx << std::endl;

  std::cout << "px    :" << std::setw(25) << px << std::endl;
  std::cout << "py    :" << std::setw(25) << py << std::endl;
  std::cout << "pz    :" << std::setw(25) << pz << std::endl;

  const double pt = TMath::Sqrt(py*py+pz*pz);

  const double p = TMath::Sqrt(px*px+pt*pt);

  const double e = TMath::Sqrt(p*p + mass*mass);

  const double gamma = e/mass;

  const double beta = p/e;

  std::cout << "pt    :" << std::setw(25) << pt << std::endl;
  std::cout << "p     :" << std::setw(25) << p << std::endl;
  std::cout << "e     :" << std::setw(25) << e << std::endl;
  std::cout << "gamma :" << std::setw(25) << gamma << std::endl;
  std::cout << "beta  :" << std::setw(25) << beta << std::endl;

  const double bx = beta * px / p;
  const double by = beta * py / p;
  const double bz = beta * pz / p;

  const double bt = beta * pt / p;

  std::cout << "bx    :" << std::setw(25) << bx << std::endl;
  std::cout << "by    :" << std::setw(25) << by << std::endl;
  std::cout << "bz    :" << std::setw(25) << bz << std::endl;
  std::cout << "bt    :" << std::setw(25) << bt << std::endl;

  const double c = 299792458.;

  const double r = pt/(c/1E9*Bx*kG2T*TMath::Abs(Z))*m2cm;

  const double period = 2*TMath::Pi()*r/(bt*c);

  const double yC = y0 + pz/pt * r * sign;
  const double zC = z0 - py/pt * r * sign;

  std::cout << "r     :" << std::setw(25) << r << std::endl;
  std::cout << "T     :" << std::setw(25) << period << std::endl;
  std::cout << "yC    :" << std::setw(25) << yC << std::endl;
  std::cout << "zC    :" << std::setw(25) << zC << std::endl;

  std::cout << "q/|p| :" << std::setw(25) << Z/p << std::endl;
  std::cout << "px/pz :" << std::setw(25) << px/pz << std::endl;
  std::cout << "py/pz :" << std::setw(25) << py/pz << std::endl;
  std::cout << "x0    :" << std::setw(25) << x0 << std::endl;
  std::cout << "y0    :" << std::setw(25) << y0 << std::endl;

  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(Bx, 0., 0.)); // kGauss
  
  genfit::MaterialEffects::getInstance()->setDebugLvl(dbgLvl);


  // init event display
  genfit::EventDisplay* display = genfit::EventDisplay::getInstance();


  // init fitter
  genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();
  fitter->setDebugLvl(dbgLvl);

  // start values for the fit, e.g. from pattern recognition
  TRandom3 rand(0);
  // with a value around 2 the fit is not precise
  double xsigma = 0.1;
  double psigma = 0.1;
  
  // start values for the fit, e.g. from pattern recognition
  double x0_guess = x0+xsigma*rand.Gaus();
  double y0_guess = y0+xsigma*rand.Gaus();
  double z0_guess = z0+xsigma*rand.Gaus();
  double px_guess = px*(1+psigma*rand.Gaus());
  double py_guess = py*(1+psigma*rand.Gaus());
  double pz_guess = pz*(1+psigma*rand.Gaus());
  TVector3 pos(x0_guess, y0_guess, z0_guess);
  TVector3 mom(px_guess, py_guess, pz_guess);

  // trackrep
  genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);
  rep->setDebugLvl(dbgLvl);

  // create track
  genfit::Track fitTrack(rep, pos, mom);

  const int detId(0); // detector ID
  int planeId(0); // detector plane ID
  int hitId(0); // hit ID

  double detectorResolution(1.); // resolution of planar detectors
  TMatrixDSym hitCov(2);
  hitCov.UnitMatrix();
  hitCov *= detectorResolution*detectorResolution;


  // add some planar hits to track with coordinates I just made up
  TVectorD hitCoords(2);
  
  //int n = 5;
  int n = x.size();

  for(unsigned int j = 0; j < n; j++)
  { 
    hitCoords[0] = x[j]*mm2cm;
    hitCoords[1] = y[j]*mm2cm;
    genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, detId, ++hitId, nullptr);
    measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(0,0,z[j]*mm2cm), TVector3(1,0,0), TVector3(0,1,0))), ++planeId);
    fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));
  } 

  //check
  fitTrack.checkConsistency();

  // do the fit
  fitter->processTrack(&fitTrack);

  // print fit result
  double p_reco = fitTrack.getFittedState().getMom().Mag();
  double px_reco = fitTrack.getFittedState().getMom().X();
  double py_reco = fitTrack.getFittedState().getMom().Y();
  double pz_reco = fitTrack.getFittedState().getMom().Z();
  double pt_reco = TMath::Sqrt(px_reco*px_reco+py_reco*py_reco+pz_reco*pz_reco);
  double x0_reco = fitTrack.getFittedState().getPos().X();
  double y0_reco = fitTrack.getFittedState().getPos().Y();
  std::cout << "****** fit result *******" << std::endl;
  std::cout << "         " << std::setw(25) << "true" << std::setw(25) << "reco" << std::setw(25) << "residual(%)" << std::endl;
  std::cout << "q/|p|   :" << std::setw(25) << Z/p << std::setw(25) << fitTrack.getFittedState().getQop() << std::setw(25) << 100*TMath::Abs(1. - fitTrack.getFittedState().getQop()/(Z/p)) << std::endl;
  std::cout << "px/pz   :" << std::setw(25) << px/pz << std::setw(25) << px_reco/pz_reco << std::setw(25) << 100*TMath::Abs(1. - (px_reco/pz_reco)/(px/pz)) << std::endl;
  std::cout << "py/pz   :" << std::setw(25) << py/pz << std::setw(25) << py_reco/pz_reco << std::setw(25) << 100*TMath::Abs(1. - (py_reco/pz_reco)/(py/pz)) << std::endl;
  std::cout << "x0      :" << std::setw(25) << x0 << std::setw(25) << x0_reco << std::setw(25) << 100*TMath::Abs(1. - x0_reco/x0) << std::endl;
  std::cout << "y0      :" << std::setw(25) << y0 << std::setw(25) << y0_reco << std::setw(25) << 100*TMath::Abs(1. - y0_reco/y0) << std::endl;
  //fitTrack.getFittedState().Print();
  //fitTrack.Print();
  std::cout << "*************************" << std::endl;

  // print particle parameters
  std::cout << "****** particle momentum *******" << std::endl;
  std::cout << "         " << std::setw(25) << "true" << std::setw(25) << "reco" << std::setw(25) << "residual(%)" << std::endl;
  std::cout << "p       :" << std::setw(25) << p << std::setw(25) << p_reco << std::setw(25) << 100*TMath::Abs(1. - p_reco/p) << std::endl;
  std::cout << "px      :" << std::setw(25) << px << std::setw(25) << px_reco << std::setw(25) << 100*TMath::Abs(1. - px_reco/px) << std::endl;
  std::cout << "py      :" << std::setw(25) << py << std::setw(25) << py_reco << std::setw(25) << 100*TMath::Abs(1. - py_reco/py) << std::endl;
  std::cout << "pz      :" << std::setw(25) << pz << std::setw(25) << pz_reco << std::setw(25) << 100*TMath::Abs(1. - pz_reco/pz) << std::endl;
  std::cout << "*************************" << std::endl;

  // print track parameters
  double r_reco = -m2cm/(fitTrack.getFittedState().getQop()*c/1E9*Bx*kG2T);
  double ez = py_reco/pt_reco;
  double ey = -pz_reco/pt_reco;
  double yC_reco = y0_reco + ey * r_reco;
  double zC_reco = z0 + ez * r_reco;
  std::cout << "****** track parameters *******" << std::endl;
  std::cout << "         " << std::setw(25) << "true" << std::setw(25) << "reco" << std::setw(25) << "residual" << std::endl;
  std::cout << "R       :" << std::setw(25) << r << std::setw(25) << r_reco << std::setw(25) << 100*TMath::Abs(1. - r_reco/r) << std::endl;
  std::cout << "yC      :" << std::setw(25) << yC << std::setw(25) << yC_reco << std::setw(25) << 100*TMath::Abs(1. - yC_reco/yC) << std::endl;
  std::cout << "zC      :" << std::setw(25) << zC << std::setw(25) << zC_reco << std::setw(25) << 100*TMath::Abs(1. - zC_reco/zC) << std::endl;
  std::cout << "*************************" << std::endl;

  //check
  fitTrack.checkConsistency();


  display->addEvent(&fitTrack);


  delete fitter;

  // open event display
  display->open();

}


