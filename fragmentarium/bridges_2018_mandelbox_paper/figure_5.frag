// Mandelbox4D ShapeInversion experiments with hypercubes and hypersphere/hypercube mixtures
//    created by Gregg Helt
// Allows replacing of (hyper)spherical inversions with shape inversion
//    using shapes that are hypercubes and unions, intersections, and blends of hyperspheres and hypercubes
// Also allows replacing of min radius with same shapes
//
// Uses techniques described in paper:
//      "Extending Mandelbox Fractals with Shape Inversions", Gregg Helt, July 2018, Bridges 2018 Conference Proceedings
//
// Includes presets to recreate figures 5a-d from paper.


#define providesInit
#include "DE-Raytracer.frag"
#include "MathUtils.frag"

#group Mandelbox4D_ShapeInversion

// Iteration params
uniform int Iterations;  slider[0,17,300]
uniform int ColorIterations;  slider[0,3,300]
uniform float Bailout;  slider[0,1000,2000]

// BoxFold params
uniform float BoxFoldLimit;  slider[-5,1,5]

// SphereFold params
// Shape types:
//    0 => 4D hypersphere (glome)
//    1 => 4D hypercube	(tesseract)		    
//    2 => 4D hypersphere/hypercube union
//    3 => 4D hypershpere/hypercube intersection
//    4 => 4D hypersphere/hypercube simple linear blend
uniform int InvType; slider[0,0,4]
uniform float InvScale;  slider[0.0001,1.0,4.0]
uniform vec4  InvCenter; slider[(-1,-1,-1,-1),(0,0,0,0),(1,1,1,1)]
uniform float InvParamA;  slider[-6.0,2.0,6.0]

uniform int MinType; slider[0,0,4]
uniform float MinScale;  slider[0.0001,0.25,4.0]
uniform vec4  MinCenter; slider[(-1,-1,-1,-1),(0,0,0,0),(1,1,1,1)]
uniform float MinParamA;  slider[-6.0,2.0,6.0]

uniform float Scale;  slider[-5.0,3.0,5.0]
uniform float MinScalingConstant;  slider[-5.0,1.0,5.0]
uniform float OffsetW; slider[-8.0,0.0,8.0]
uniform int InnerMode; slider[0,0,4]

const float pi =  3.14159265359;

struct Shape {
  int type;
  float scale;
  vec4 center;
  float a; // additional shape param (only used by some shapes)
};

Shape invShape;
Shape minShape;

float absScalem1 = abs(Scale - 1.0);
float AbsScaleRaisedTo1mIters = pow(abs(Scale), float(1-Iterations));

void init() {
  invShape = Shape(InvType,
		   InvScale,
		   InvCenter,
		   InvParamA);
  minShape = Shape(MinType,
		   MinScale, 
		   MinCenter,
		   MinParamA);
}

float getShape2CenterSqr(vec4 p, Shape s) {
  vec4 intersectP;
  float idistance2 = 1.0; // square of distance from origin to intersection of shape boundary along ray from origin through point P
  if (s.type == 0) { // hypersphere
    idistance2 =  s.scale * s.scale;
  }
  else if (s.type == 1) {  // hypercube
    if (p.x == 0.0 && p.y == 0.0 && p.z == 0.0 && p.w == 0.0) {
      intersectP.xyzw = p.xyzw;
      idistance2 = dot(intersectP, intersectP);
    }
    else {
      float maxdim = max(abs(p.x), max(abs(p.y), max(abs(p.z), abs(p.w))));
      intersectP = (p / maxdim) * s.scale;
      idistance2 =  dot(intersectP, intersectP);
    }
  }
  else if (s.type == 2 || s.type == 3) {  // union or intersection of hypercube and hypersphere
    if (p.x == 0.0 && p.y == 0.0 && p.z == 0.0 && p.w == 0.0) { idistance2 = 0.0; }
    else {
      float maxdim = max(abs(p.x), max(abs(p.y), max(abs(p.z), abs(p.w))));
      vec4 box_intersect = (p / maxdim) * s.scale;
      float boxd2 = dot(box_intersect,box_intersect);
      float sphered2 = s.a * s.a * s.scale * s.scale;

      if (s.type == 2) {  // union
	idistance2 = max(boxd2, sphered2);
      }
      else if (s.type == 3) { // intersection
        idistance2 = min(boxd2, sphered2);
      }
    }
  }
  else if (s.type == 4) {  // blend diff of cube and sphere inscribed in cube (sphere touches center of cube faces)
    // ParamA controls blend contribution from hypercube and hypersphere
    // blend = (a * distance(hypercube)) + ((1-a) * distance(hypersphere))
    // A = 0 ==> hypersphere
    // A = 1 ==> hypercube
    // 0 < A < 1 ==> blend in between hypersphere and hypercube surface
    // A < 0 || A > 1 ==> overdriven blends
    float maxdim = max(abs(p.x), max(abs(p.y), max(abs(p.z), abs(p.w))));
    vec4 box_intersect = (p / maxdim) * s.scale;
    float boxd2 = dot(box_intersect,box_intersect);
    float sphere_distance = s.scale;
    float box_distance = sqrt(boxd2);
    float sdistance = (1.0-s.a)*sphere_distance + s.a*box_distance;
    idistance2 = sdistance * sdistance;
  }
  else {  // out of range for ShapeTypes
    idistance2 =  s.scale;
  }
  return idistance2;
}


// modified inversion fold that can use shapes other than circles
void shapeFold(inout vec4 p, inout float der) {
  vec4 minp = p;

  // get point position _relative_ to shape center
  p -= invShape.center;
  minp -= minShape.center;

  // get distance from point P to center of invShape, minShape, maxShape
  float p2InvCenterSqr = dot(p,p);
  float p2MinCenterSqr = dot(minp, minp);

  // get point of intersection S with surface of shape (along a ray from shape center to point P)
  float s2InvCenterSqr = getShape2CenterSqr(p, invShape);
  float s2MinCenterSqr = getShape2CenterSqr(minp, minShape);

  // determine if point P is inside invShape, minShape
  bool insideMinShape = (p2MinCenterSqr < s2MinCenterSqr);
  bool insideInvShape = (p2InvCenterSqr < s2InvCenterSqr);

  float ratio = 1.0;
  if (insideMinShape) {
    if (InnerMode == 0) { // scaling varies based on intersection with inv and min boundaries
      //   do inversion _but_ as if point were at edge of min shape
      ratio = s2InvCenterSqr/s2MinCenterSqr;
    }
    else if (InnerMode == 1) { // scaling is constant based on invShape and minShape scales
      ratio = (invShape.scale * invShape.scale) / (minShape.scale * minShape.scale);
    }
    else if (InnerMode == 2) { // scaling varies based on invShape scale and minShape intersection
      ratio = (invShape.scale * invShape.scale) /s2MinCenterSqr;
    }
    else if (InnerMode == 3) { // scaling varies based on invShape intersection and minShape scale
      ratio = s2InvCenterSqr/ (minShape.scale * minShape.scale);
    }
    // if InnerMode == 4, just use MinScalingConstant
    // testing an additional scaling constant 
    ratio *= MinScalingConstant;
    der *= ratio;
    p *= ratio;
  }
  else if (insideInvShape) {
    // if point is outside minShape but within invShape, do shape inversion relative to invShape
    ratio = s2InvCenterSqr/p2InvCenterSqr;
    p *= ratio;
    der *= ratio;
  }
  // get absolute point position (reverse offset relative to shape center)
  p += invShape.center;
}

// standard Mandelbox box fold in 4D
void boxFold(inout vec4 p, inout float der) {
  // fold x
  if (p.x > BoxFoldLimit) { p.x = 2.0*BoxFoldLimit - p.x; } 
  else if (p.x < -BoxFoldLimit) { p.x = -2.0*BoxFoldLimit - p.x; }
  // fold y
  if (p.y > BoxFoldLimit) { p.y = 2.0*BoxFoldLimit - p.y; }
  else if (p.y < -BoxFoldLimit) { p.y = -2.0*BoxFoldLimit -p.y; }
  // fold z
  if (p.z > BoxFoldLimit) { p.z = 2.0*BoxFoldLimit - p.z; }
  else if (p.z < -BoxFoldLimit) { p.z = -2.0*BoxFoldLimit - p.z; }
  // fold w
  if (p.w > BoxFoldLimit) { p.w = 2.0*BoxFoldLimit - p.w; }
  else if (p.w < -BoxFoldLimit) { p.w = -2.0*BoxFoldLimit  - p.w; }
}


float DE(vec3 p3d) {
  vec4 p4d = vec4(p3d.x, p3d.y, p3d.z, OffsetW);
  vec4 offset = p4d;
  float dr = 1.0;
  float distance2 = 0.0;
  for (int i = 0; i < Iterations; i++) {
    boxFold(p4d,dr);      
    if (i < ColorIterations) {
      float d2 = dot(p4d.xyz, p4d.xyz);
      orbitTrap = min(orbitTrap, abs(vec4(p4d.xyz, d2)));
    }
    shapeFold(p4d,dr);  
    p4d = (p4d * Scale) + offset; 
    dr = (dr * abs(Scale)) + 1.0;
    distance2 = dot(p4d,p4d);
    if (distance2 > Bailout) { break; }
  }
  return ((sqrt(distance2) - absScalem1) / abs(dr)) - AbsScaleRaisedTo1mIters;
}


#preset Default
FOV = 0.2156863
Eye = -7.8e-6,0,-0.3701195
Target = -0.0002615,-1e-7,2.942352
Up = 0,1,0
EquiRectangular = false
Gamma = 2.028688
ToneMapping = 4
Exposure = 0.7
Brightness = 3.262712
Contrast = 1.477273
Saturation = 1.427061
GaussianWeight = 1
AntiAliasScale = 2
Detail = -2.72379
DetailAO = -0.28574
FudgeFactor = 1
Dither = 0.5
NormalBackStep = 1
AO = 0,0,0,0.4670782
Specular = 0.5942029
SpecularExp = 8
SpecularMax = 9.255079
SpotLight = 1,1,1,1
SpotLightDir = -1,1
CamLight = 1,1,1,0.6502242
CamLightMin = 0.4697
Glow = 1,1,1,0
GlowMax = 0
Fog = 0.4980545
HardShadow = 0
ShadowSoft = 2
Reflection = 0
DebugSun = false
BaseColor = 1,1,1
OrbitStrength = 0.803532
X = 0.5,0.6,0.6,1
Y = 1,0.6666667,0,-0.4786151
Z = 0.8,0.78,1,0.5804481
R = 0.4,0.7,1,0.3156823
BackgroundColor = 0,0,0
GradientBackground = 0.3
CycleColors = true
Cycles = 5.499591
EnableFloor = false
FloorNormal = 0,0,0
FloorHeight = 0
FloorColor = 1,1,1
Iterations = 21
ColorIterations = 3
Bailout = 1000
BoxFoldLimit = 0.9545455
InvScale = 2.080727
InvCenter = 0,0,0.1154684,0.0980392
InvParamA = 0.9487522
MinScale = 1.171011
MinCenter = 0,0,0,0
MinParamA = 0.1469933
Scale = 2.134021
MinScalingConstant = 1
OffsetW = -1.722512
FocalPlane = 1
Aperture = 0
MaxRaySteps = 337
InvType = 3
MinType = 0
InnerMode = 0
#endpreset


#preset Figure 5a
FOV = 0.4745098
Eye = 3.322239,-2.738474,2.510036
Target = 3.008111,-2.882487,12.50404
Up = 0.9994209,-0.0133151,0.0312216
EquiRectangular = false
Gamma = 1.659836
ToneMapping = 4
Exposure = 0.58125
Brightness = 1.207627
Contrast = 1.590909
Saturation = 1
GaussianWeight = 1
AntiAliasScale = 2
Detail = -3.627016
DetailAO = -0.3068894
FudgeFactor = 0.3928571
Dither = 0.7318008
NormalBackStep = 1
AO = 0,0,0,0.600823
Specular = 0.3581781
SpecularExp = 16
SpecularMax = 10
SpotLight = 1,1,1,1
SpotLightDir = -1,0
CamLight = 1,1,1,0.8952381
CamLightMin = 0
Glow = 1,1,1,0
GlowMax = 20
Fog = 2
HardShadow = 0
ShadowSoft = 2
Reflection = 0
DebugSun = false
BaseColor = 1,1,1
OrbitStrength = 1
X = 0.5,0.6,0.6,0.5885947
Y = 1,0.6,0,0.4786151
Z = 0.8,0.78,1,0.4582485
R = 0.4,0.7,1,0.0264766
BackgroundColor = 0,0,0
GradientBackground = 0.04
CycleColors = false
Cycles = 0.1
EnableFloor = false
FloorNormal = 0,0,1
FloorHeight = 0
FloorColor = 1,1,1
Iterations = 20
ColorIterations = 3
Bailout = 1000
BoxFoldLimit = 1
InvType = 2
InvScale = 2.528
InvCenter = 0,0,0,0.1895425
InvParamA = 0.5309734
MinType = 1
MinScale = 1.709459
MinCenter = 0,0,0,0
MinParamA = 2.899777
Scale = 1.948454
MinScalingConstant = 1
OffsetW = 0.6182213
FocalPlane = 1
Aperture = 0
MaxRaySteps = 397
InnerMode = 1
#endpreset

  
#preset Figure 5b
FOV = 0.4745098
Eye = 3.322206,-2.73849,2.511104
Target = 3.008071,-2.88302,12.5051
Up = 0.9994209,-0.0133151,0.0312216
EquiRectangular = false
FocalPlane = 1
Aperture = 0
Gamma = 1.659836
ToneMapping = 4
Exposure = 0.58125
Brightness = 1.207627
Contrast = 1.590909
Saturation = 1
GaussianWeight = 1
AntiAliasScale = 2
Detail = -3.810484
DetailAO = -0.3068894
FudgeFactor = 0.3928571
MaxRaySteps = 397
Dither = 0.7318008
NormalBackStep = 1
AO = 0,0,0,0.600823
Specular = 0.3581781
SpecularExp = 16
SpecularMax = 10
SpotLight = 1,1,1,1
SpotLightDir = -1,0
CamLight = 1,1,1,0.8952381
CamLightMin = 0
Glow = 1,1,1,0
GlowMax = 20
Fog = 2
HardShadow = 0
ShadowSoft = 2
Reflection = 0
DebugSun = false
BaseColor = 1,1,1
OrbitStrength = 1
X = 0.5,0.6,0.6,0.5885947
Y = 1,0.6,0,0.4786151
Z = 0.8,0.78,1,0.4582485
R = 0.4,0.7,1,0.0264766
BackgroundColor = 0,0,0
GradientBackground = 0.04
CycleColors = false
Cycles = 0.1
EnableFloor = false
FloorNormal = 0,0,1
FloorHeight = 0
FloorColor = 1,1,1
Iterations = 20
ColorIterations = 3
Bailout = 1000
BoxFoldLimit = 1
InvType = 2
InvScale = 2.528
InvCenter = 0,0,0,0.1895425
InvParamA = 1.17
MinType = 1
MinScale = 1.709459
MinCenter = 0,0,0,0
MinParamA = 2.899777
Scale = 1.948454
MinScalingConstant = 1
OffsetW = -2.049892
InnerMode = 1
#endpreset

#preset Figure 5c
FOV = 0.2156863
Eye = -7.8e-6,0,-0.3701195
Target = -0.0002615,-1e-7,2.942352
Up = 0,1,0
EquiRectangular = false
Gamma = 2.028688
ToneMapping = 4
Exposure = 0.7
Brightness = 3.262712
Contrast = 1.477273
Saturation = 1.427061
GaussianWeight = 1
AntiAliasScale = 2
Detail = -2.72379
DetailAO = -0.28574
FudgeFactor = 1
Dither = 0.5
NormalBackStep = 1
AO = 0,0,0,0.4670782
Specular = 0.5942029
SpecularExp = 8
SpecularMax = 9.255079
SpotLight = 1,1,1,1
SpotLightDir = -1,1
CamLight = 1,1,1,0.6502242
CamLightMin = 0.4697
Glow = 1,1,1,0
GlowMax = 0
Fog = 0.4980545
HardShadow = 0
ShadowSoft = 2
Reflection = 0
DebugSun = false
BaseColor = 1,1,1
OrbitStrength = 0.803532
X = 0.5,0.6,0.6,1
Y = 1,0.6666667,0,-0.4786151
Z = 0.8,0.78,1,0.5804481
R = 0.4,0.7,1,0.3156823
BackgroundColor = 0,0,0
GradientBackground = 0.3
CycleColors = true
Cycles = 5.499591
EnableFloor = false
FloorNormal = 0,0,0
FloorHeight = 0
FloorColor = 1,1,1
Iterations = 21
ColorIterations = 3
Bailout = 1000
BoxFoldLimit = 0.9545455
InvScale = 2.080727
InvCenter = 0,0,0.1154684,0.0980392
InvParamA = 0.9487522
MinScale = 1.171011
MinCenter = 0,0,0,0
MinParamA = 0.1469933
Scale = 2.134021
MinScalingConstant = 1
OffsetW = -1.722512
FocalPlane = 1
Aperture = 0
MaxRaySteps = 337
InvType = 3
MinType = 0
InnerMode = 0
#endpreset


#preset Figure 5d
FOV = 0.0509804
Eye = 0,0,-2.67657
Target = 0,0,-0.1562559
Up = 0,1,0
EquiRectangular = false
FocalPlane = 1
Aperture = 0
Gamma = 2
ToneMapping = 4
Exposure = 0.9375
Brightness = 2.817797
Contrast = 1.838843
Saturation = 1.828753
GaussianWeight = 1
AntiAliasScale = 2
Detail = -4.008064
DetailAO = -0.28574
FudgeFactor = 1
MaxRaySteps = 337
Dither = 0.5
NormalBackStep = 1
AO = 0,0,0,0.4670782
Specular = 0.5942029
SpecularExp = 8
SpecularMax = 9.255079
SpotLight = 1,1,1,1
SpotLightDir = -1,1
CamLight = 1,1,1,0.6502242
CamLightMin = 0.4697
Glow = 1,1,1,0
GlowMax = 0
Fog = 0.4980545
HardShadow = 0
ShadowSoft = 2
Reflection = 0
DebugSun = false
BaseColor = 1,1,1
OrbitStrength = 0.803532
X = 0.5,0.6,0.6,1
Y = 1,0.6666667,0,-0.4786151
Z = 0.8,0.78,1,0.5804481
R = 0.4,0.7,1,0.3156823
BackgroundColor = 0,0,0
GradientBackground = 0.3
CycleColors = true
Cycles = 5.302045
EnableFloor = false
FloorNormal = 0,0,0
FloorHeight = 0
FloorColor = 1,1,1
Iterations = 21
ColorIterations = 3
Bailout = 1000
BoxFoldLimit = 0.9545455
InvType = 4
InvScale = 2.080727
InvCenter = 0,0,-0.0108932,0.1198257
InvParamA = 0.2920354
MinType = 4
MinScale = 1.171011
MinCenter = 0,0,0,0
MinParamA = 0.013363
Scale = 2.072165
MinScalingConstant = 1
OffsetW = 0.33
InnerMode = 0
#endpreset

