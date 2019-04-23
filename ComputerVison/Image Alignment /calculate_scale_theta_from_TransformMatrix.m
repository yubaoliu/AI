# https://www.mathworks.com/help/vision/examples/find-image-rotation-and-scale-using-automated-feature-matching.html


T=[7.69676379e-01  -3.70207415e-01  1.50656113e+02;
3.71046317e-01  7.68901573e-01  -3.42393076e+01;
2.16076945e-06  6.55717989e-07  1.00000000e+00
];

Tinv = T.'
ss = Tinv(2,1);
sc = Tinv(1,1);
scaleRecovered = sqrt(ss*ss + sc*sc)
thetaRecovered = atan2(ss,sc)*180/pi