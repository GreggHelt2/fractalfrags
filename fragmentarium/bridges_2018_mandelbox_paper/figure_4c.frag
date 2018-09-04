// Mandelbox3d ShapeInversion experiments with Gielis supershapes
//    created by Gregg Helt
// Allows replacing of spherical inversion with shape inversion where shape types are:
//    0 => sphere
//    1 => cube
//    2 => Gielis 3D supershape
//
// Also allows replacing of min radius with same shapes
//
// Uses techniques described in paper:
//      "Extending Mandelbox Fractals with Shape Inversions", Gregg Helt, July 2018, Bridges 2018 Conference Proceedings
//      http://archive.bridgesmathart.org/2018/bridges2018-547.html
//
// Includes preset to recreate figure 4c from paper
//
// For more info on supershapes, see paper:
//     "Superquadrics with rational and irrational symmetry", Gielis et al 2003:
//     https://dl.acm.org/citation.cfm?id=781647
// or Paul Bourke's overview: http://paulbourke.net/geometry/supershape/

#define providesInit
#include "DE-Raytracer.frag"
#include "MathUtils.frag"
#group Mandelbox3D_ShapeInversion

// Iteration params
uniform int Iterations;  slider[0,17,300]
uniform int ColorIterations;  slider[0,3,300]
uniform float Bailout;  slider[0,1000,2000]

// BoxFold params
uniform float BoxFoldLimit;  slider[-5,1,5]

// SphereFold params
uniform int InvType; slider[0,0,2]
uniform float InvScale;  slider[0.0001,1.0,4.0]
uniform vec3  InvCenter; slider[(-1,-1,-1),(0,0,0),(1,1,1)]
uniform vec3  InvRotation; slider[(-180,-180,-180),(0,0,0),(180,180,180)]

uniform int MinType; slider[0,0,2]
uniform float MinScale;  slider[0.0001,0.25,4.0]
uniform vec3  MinCenter; slider[(-1,-1,-1),(0,0,0),(1,1,1)]
uniform vec3  MinRotation; slider[(-180,-180,-180),(0,0,0),(180,180,180)]

uniform float ThetaPiOffset; slider[-2.0,0.0,2.0]
uniform float PhiPiOffset; slider[-2.0,0.0,2.0]

uniform float Scale;  slider[-5.0,3.0,5.0]
uniform float MinScalingConstant;  slider[-5.0,1.0,5.0]

// Rotation params
uniform float RotAngle; slider[0.00,0,180]
uniform vec3 RotVector; slider[(0,0,0),(1,1,1),(1,1,1)]

// supershape params (applies to either InvShape and/or MinShape if they are supershapes)
uniform float n1; slider[-1000.0,1.0,1000.0]
uniform float n2; slider[-1000.0,1.0,1000.0]
uniform float n3; slider[-1000.0,1.0,1000.0]
uniform float m1; slider[-20.0,0.0,20.0]

uniform float n4; slider[-1000.0,1.0,1000.0]
uniform float n5; slider[-1000.0,1.0,1000.0]
uniform float n6; slider[-1000.0,1.0,1000.0]
uniform float m2; slider[-20.0,0.0,20.0]


const float pi =  3.14159265359;

struct Shape {
  int type;
  float scale;
  vec3 center;
  mat3 rotMat;   // rotation matrix for shape rotation
};


// vars calculated in init()
mat3 rot;
Shape invShape;
Shape minShape;

float absScalem1 = abs(Scale - 1.0);
float AbsScaleRaisedTo1mIters = pow(abs(Scale), float(1-Iterations));

void init() {
  rot = rotationMatrix3(normalize(RotVector), RotAngle);
  invShape = Shape(InvType,
		   InvScale,
		   InvCenter,
		   rotationMatrixXYZ(InvRotation));
  minShape = Shape(MinType,
		   MinScale, 
		   MinCenter,
		   rotationMatrixXYZ(MinRotation));
}


vec3 getShapeIntersect(vec3 p, Shape s) {
  p *= s.rotMat; 
  float pdistance = length(p);
  vec3 intersectP = vec3(0.0,0.0,0.0);
  // GLSL acos return value is in range [0, +PI]
  //     (and value is undefined if |arg| > 1)
  float theta = acos(p.z/pdistance);
  // GLSL atan return value is in range [-PI, +PI]
  //     (and value is undefined if arg2 == 0)
  float phi;
  if (p.x == 0.0) {
    phi = atan(p.y,0.00001);
  }
  else {
    phi = atan(p.y,p.x);
  }

  theta -= pi/2.0;
  theta += ThetaPiOffset;
  phi += PhiPiOffset;

  if (s.type == 0) { // circle
    intersectP = s.scale * vec3(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
  }
  else if (s.type == 1) { // origin-centered cube  (1 for current, 6 is legacy number for cube mode
    if (p.x == 0.0 && p.y == 0.0 && p.z == 0.0) { intersectP.xyz = p.xyz; }
    else {
      float maxdim = max(abs(p.x), max(abs(p.y), abs(p.z)));
      intersectP = (p / maxdim) * s.scale;
    }
  }
  else if (s.type == 2) {
    // supershape formula expects phi range: [-pi/2,pi/2] and theta range [-pi,pi],
    // so switch phi and theta:
    float tmp  = phi;
    phi = theta;
    theta = tmp;

    float r1 = pow( (pow(abs(cos(m1*theta/4.0)),n2) + pow(abs(sin(m1*theta/4.0)),n3)), -1.0/n1);
    float r2 = pow(   (pow(abs(cos(m2*phi/4.0)),n5) +     pow(abs(sin(m2*phi/4.0)),n6)), -1.0/n4);
    intersectP.x = s.scale * r1 * cos(theta) * r2 * cos(phi);
    intersectP.y = s.scale * r1 * sin(theta) * r2 * cos(phi);
    intersectP.z = s.scale * r2 * sin(phi);
  }
  else {  // out of range for ShapeTypes
    intersectP =  vec3(1.0, 1.0, 1.0);
  }
  return intersectP;
}


// modified inversion that can use shapes other than circles
void shapeInversion(inout vec3 p, inout float der) {
  vec3 minp = p;

  // get point position _relative_ to shape center
  p -= invShape.center;
  minp -= minShape.center;

  // get distance from point P to center of invShape, minShape, maxShape
  float p2InvCenterSqr = dot(p,p);
  float p2MinCenterSqr = dot(minp, minp);

  // get point of intersection S with surface of shape (along a ray from shape center to point P)
  vec3 invIntersect = getShapeIntersect(p, invShape);
  vec3 minIntersect = getShapeIntersect(minp, minShape);

  // get distance from shape surface intersection point S to shape center
  float s2InvCenterSqr = dot(invIntersect, invIntersect);
  float s2MinCenterSqr = dot(minIntersect, minIntersect);

  // determine if point P is inside invShape, minShape
  bool insideMinShape = (p2MinCenterSqr < s2MinCenterSqr);
  bool insideInvShape = (p2InvCenterSqr < s2InvCenterSqr);

  if (insideMinShape) {
    float ratio = (invShape.scale * invShape.scale) / (minShape.scale * minShape.scale);
    der *= ratio;
    ratio *= MinScalingConstant;
    p *= ratio;
  }
  else if (insideInvShape) {
    // if point is outside minShape but within invShape, do shape inversion relative to invShape
    float ratio = s2InvCenterSqr/p2InvCenterSqr;
    p *= ratio;
    der *= ratio;
  }
  // get absolute point position (reverse offset relative to shape center)
  p += invShape.center;
}

// standard Mandelbox box fold in 2D
void boxFold(inout vec3 p, inout float der) {
  // fold x
  if (p.x > BoxFoldLimit) { p.x = 2.0*BoxFoldLimit - p.x; } 
  else if (p.x < -BoxFoldLimit) { p.x = -2.0*BoxFoldLimit - p.x; }
  // fold y
  if (p.y > BoxFoldLimit) { p.y = 2.0*BoxFoldLimit - p.y; }
  else if (p.y < -BoxFoldLimit) { p.y = -2.0*BoxFoldLimit -p.y; }
  // fold z
  if (p.z > BoxFoldLimit) { p.z = 2.0*BoxFoldLimit - p.z; }
  else if (p.z < -BoxFoldLimit) { p.z = -2.0*BoxFoldLimit -p.z; }
}


float DE(vec3 p) {
  vec3 offset = p;
  float dr = 1.0;
  for (int i = 0; i < Iterations; i++) {
    p *= rot;      // rotation (does not change derivative)
    boxFold(p,dr);       // conditional reflection and scale
    if (i < ColorIterations) {
      float d2 = dot(p,p);
      orbitTrap = min(orbitTrap, abs(vec4(p, d2)));
    }
    shapeInversion(p,dr);    // conditional shape inversion
    p = (p * Scale) + offset;  // scale and translate
    dr = (dr * abs(Scale)) + 1.0;
    float distance = dot(p,p);
    if (distance > Bailout) { break; }
  }
  return ((length(p) - absScalem1) / abs(dr)) - AbsScaleRaisedTo1mIters;
}


#preset Default
FOV = 0.4
Eye = 2.828954,2.325113,-3.060299
Target = 2.902178,7.246194,5.644637
Up = 0.0126848,-0.8704981,0.4920039
EquiRectangular = false
Gamma = 1.659836
ToneMapping = 4
Exposure = 0.5125
Brightness = 1.652542
Contrast = 1.466942
Saturation = 1
GaussianWeight = 1
AntiAliasScale = 2
Detail = -3.612903
DetailAO = -0.5
FudgeFactor = 0.373102
Dither = 0.7318008
NormalBackStep = 1
AO = 0,0,0,0.8482759
Specular = 0.4084507
SpecularExp = 16
SpecularMax = 10
SpotLight = 1,1,1,0.7714286
SpotLightDir = 1,0.6892779
CamLight = 1,1,1,0.8952381
CamLightMin = 0
Glow = 1,1,1,0
GlowMax = 20
HardShadow = 0
ShadowSoft = 2
Reflection = 0
DebugSun = false
BaseColor = 1,1,1
OrbitStrength = 1
X = 0.6666667,0.6666667,0.4980392,1
Y = 1,0.6,0,0.4175153
Z = 0.8,0.78,1,0.209776
R = 0.4,0.7,1,0.0509165
BackgroundColor = 0,0,0
GradientBackground = 0.04
CycleColors = true
Cycles = 18.14254
EnableFloor = false
FloorNormal = 0,0,1
FloorHeight = 0
FloorColor = 1,1,1
FocalPlane = 1
Aperture = 0
MaxRaySteps = 474
Fog = 2.7
Iterations = 17
ColorIterations = 3
Bailout = 1000.0
BoxFoldLimit = 1
InvType = 2
InvScale = 2.528
InvCenter = 0,0,0
InvRotation = 0,0,0
MinType = 0
MinScale = 1.502137
MinCenter = 0,0,0
MinRotation = 0,0,0
ThetaPiOffset = 0
PhiPiOffset = 0
Scale = 1.913828
MinScalingConstant = 1
RotAngle = 0
RotVector = 1,1,1
n1 = 5
n2 = 5
n3 = 5
m1 = 4
n4 = 5
n5 = 5
n6 = 5
m2 = 4
#endpreset

