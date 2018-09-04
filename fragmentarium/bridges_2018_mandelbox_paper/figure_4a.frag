// Mandelbox3d ShapeInversion experiments with octahedrons
//    created by Gregg Helt
// Allows replacing of spherical inversion with shape inversion where shape types are:
//    0 => sphere
//    1 => cube
//    2 => sphere/cube union
//    3 => sphere/cube intersection
//    4 => sphere/cube simple linear blend
//    5 => cube (but done via platonic intersection method, unlike 1)
//    6 => octahedron (via platonic intersection method)
//    (other platonic solids left out here, done instead in separate frags because GLSL compilation was getting slow)
//
// Also allows replacing of min radius with same shapes
//
// Uses techniques described in paper:
//      "Extending Mandelbox Fractals with Shape Inversions", Gregg Helt, July 2018, Bridges 2018 Conference Proceedings
//      http://archive.bridgesmathart.org/2018/bridges2018-547.html
// Includes preset to recreate figure 4a from paper

#define providesInit
#include "DE-Raytracer.frag"
#include "MathUtils.frag"
#group Mandelbox3D_ShapeInversion_Platonics1

// Iteration params
uniform int Iterations;  slider[0,17,300]
uniform int ColorIterations;  slider[0,3,300]
uniform float Bailout;  slider[0,1000,2000]

// BoxFold params
uniform float BoxFoldLimit;  slider[-5,1,5]
// uniform float BoxFoldScale;  slider[-5,2,5]

// SphereFold params
uniform int InvType; slider[0,0,6]
uniform float InvScale;  slider[0.0001,1.0,4.0]
uniform vec3  InvCenter; slider[(-1,-1,-1),(0,0,0),(1,1,1)]
uniform vec3  InvRotation; slider[(-180,-180,-180),(0,0,0),(180,180,180)]
uniform float InvParamA;  slider[-6.0,2.0,6.0]

uniform int MinType; slider[0,0,6]
uniform float MinScale;  slider[0.0001,0.25,4.0]
uniform vec3  MinCenter; slider[(-1,-1,-1),(0,0,0),(1,1,1)]
uniform vec3  MinRotation; slider[(-180,-180,-180),(0,0,0),(180,180,180)]
uniform float MinParamA;  slider[-6.0,2.0,6.0]

uniform float Scale;  slider[-5.0,3.0,5.0]
uniform float MinScalingConstant;  slider[-5.0,1.0,5.0]

// Rotation params
uniform float RotAngle; slider[0.00,0,180]
uniform vec3 RotVector; slider[(0,0,0),(1,1,1),(1,1,1)]

const float pi =  3.14159265359;

struct Shape {
  int type;
  float scale;
  vec3 center;
  mat3 rotMat;   // rotation matrix for shape rotation
  float a;
};

// vars calculated in init()
mat3 rot;
float minDistanceSqr;
Shape invShape;
Shape minShape;

float sqrt_i3 = 1.0 / sqrt(3.0);
vec3 zeros = vec3(0.0,0.0,0.0);
vec3 Nv_oct[20];
vec3 Nv_cube[20];
int sides_oct = 8;
int sides_cube = 6;

float absScalem1 = abs(Scale - 1.0);
float AbsScaleRaisedTo1mIters = pow(abs(Scale), float(1-Iterations));

void init() {

  Nv_oct[0] = vec3(sqrt_i3, sqrt_i3, -sqrt_i3);
  Nv_oct[1] = vec3(sqrt_i3, -sqrt_i3, sqrt_i3);
  Nv_oct[2] = vec3(-sqrt_i3, sqrt_i3, sqrt_i3);
  Nv_oct[3] = vec3(-sqrt_i3, -sqrt_i3, -sqrt_i3);
  Nv_oct[4] = vec3(sqrt_i3, sqrt_i3, sqrt_i3);
  Nv_oct[5] = vec3(-sqrt_i3, -sqrt_i3, sqrt_i3);
  Nv_oct[6] = vec3(-sqrt_i3, sqrt_i3, -sqrt_i3);
  Nv_oct[7] = vec3(sqrt_i3, -sqrt_i3, -sqrt_i3);
  for (int i=8; i<20; i++) { Nv_oct[i]=zeros; }

  Nv_cube[0] = vec3(1.0, 0, 0);
  Nv_cube[1] = vec3(-1.0, 0, 0);
  Nv_cube[2] = vec3(0, 1.0, 0);
  Nv_cube[3] = vec3(0, -1.0, 0);
  Nv_cube[4] = vec3(0, 0, 1.0);
  Nv_cube[5] = vec3(0, 0, -1.0);
  for (int i=6; i<20; i++) { Nv_cube[i]=zeros; }
  
  rot = rotationMatrix3(normalize(RotVector), RotAngle);
  // minDistanceSqr = MinRadius * MinRadius;
  invShape = Shape(InvType,
		   InvScale,
		   InvCenter,
		   rotationMatrixXYZ(InvRotation),
		   InvParamA);
  minShape = Shape(MinType,
		   MinScale, 
		   MinCenter,
		   rotationMatrixXYZ(MinRotation),
		   MinParamA);

}

vec3 platonic_intersect(vec3 p, Shape s, vec3 face_normals[20], int face_count) {
  // Find the closest cutting plane of the platonic solid (if any) that cuts the line between the origin and point p
  //
  // algorithm adapted from parts of fractower's approach for "Generalized Box Folds", see fractalforums thread for more info:
  //      http://www.fractalforums.com/amazing-box-amazing-surf-and-variations/generalized-box-fold/
  // more specifically, from Mandelbulber2 implementation of generalized box fold:
  //      https://github.com/buddhi1980/mandelbulber2/blob/master/mandelbulber2/formula/opencl/generalized_fold_box.cl
    vec3 L = p;
    vec3 Y = vec3(0.0,0.0,0.0); // Y is the origin in this case.
    // int side = -1;
    float a, b, c;
    vec3 shape_intersect;
    vec3 platsect;
    float sidot;
    float mindist = 1000.0;
    int i;
    
    for (i = 0; i < face_count; i++) {
      b = s.scale;
      c = dot(L, face_normals[i]);
      if (c != 0.0) {
	a = b / c;
        shape_intersect = p * a;
	sidot = dot(shape_intersect, shape_intersect);
	if (sidot < mindist) {
	  mindist = sidot;
	  platsect = shape_intersect;
	}
      }
    }
    if (mindist == 1000.0) {
      // no intersections found?
      // give up for now
      // intersectP = p;
      platsect = vec3(1.0,1.0,1.0);
    }
    return platsect;
}

float getShape2CenterSqr(vec3 p, Shape s) {
  p *= s.rotMat;
  float pdistance = length(p);
  vec3 intersectP = vec3(0.0,0.0,0.0);
  float sdistance2;
  
  if (s.type == 0) { // circle
    sdistance2 = s.scale * s.scale;
  }
  else if (s.type == 1) { // origin-centered cube 
    // equation for intersection of ray from origin to point p
    // based on comments from Henning Makholm in this thread:
    // https://math.stackexchange.com/questions/465723/mapping-from-ball-into-cube

    if (p.x == 0.0 && p.y == 0.0 && p.z == 0.0) { intersectP.xyz = p.xyz; }
    else {
      float maxdim = max(abs(p.x), max(abs(p.y), abs(p.z)));
      intersectP = (p / maxdim) * s.scale;
      sdistance2 = dot(intersectP, intersectP);
    }
  }
  // sphere/cube union
  else if (s.type == 2) { 
    if (p.x == 0.0 && p.y == 0.0 && p.z == 0.0) { intersectP.xyz = p.xyz; }
    else {
      float maxdim = max(abs(p.x), max(abs(p.y), abs(p.z)));
      vec3 box_intersect = (p / maxdim) * s.scale;
      float box_distance2 = dot(box_intersect, box_intersect);
      float sphere_distance2 = s.a * s.a * s.scale * s.scale;
      if (box_distance2 > sphere_distance2) {
	sdistance2 = box_distance2;
      }
      else {
	sdistance2 = sphere_distance2;
      }
    }
  }
  // sphere/cube intersection
  else if (s.type == 3) {
    float maxdim = max(abs(p.x), max(abs(p.y), abs(p.z)));
    vec3 box_intersect = (p / maxdim) * s.scale;
    float box_distance2 = dot(box_intersect, box_intersect);
    float sphere_distance2 = s.a * s.a * s.scale * s.scale;
      
    if (box_distance2 < sphere_distance2) {
      sdistance2 = box_distance2;
    }
    else {
      sdistance2 = sphere_distance2;
    }
  }
  else if (s.type == 4) {  // blend diff of cube and sphere inscribed in cube (sphere touches center of cube faces)
    // ParamA controls blend contribution from cube and sphere
    // A = 0 ==> only sphere
    // A = 1 ==> only cube
    // blend = (a * intersect(cube)) + ((1-a) * intersect(sphere))
    float maxdim = max(abs(p.x), max(abs(p.y), abs(p.z)));
    vec3 box_intersect = (p / maxdim) * s.scale;
    float box_distance2 = dot(box_intersect, box_intersect);
    float sphere_distance2 = s.scale * s.scale;
    float box_distance = sqrt(box_distance2);
    float sphere_distance = sqrt(sphere_distance2);
    float sdistance = (1.0-s.a)*sphere_distance + s.a*box_distance;
    sdistance2 = sdistance * sdistance;
  }
  else if (s.type == 5) { // cube based on platonic intersection algorithm
    intersectP = platonic_intersect(p, s, Nv_cube, sides_cube);
    sdistance2 = dot(intersectP, intersectP);
  }
  else if (s.type == 6) { // octahedron based on platonic intersection algorithm
    intersectP = platonic_intersect(p, s, Nv_oct, sides_oct);
    sdistance2 = dot(intersectP, intersectP);
  }
  else {  // out of range for ShapeTypes
    sdistance2 = 1.0;
  }
  return sdistance2;
}


// modified inversion that can use shapes other than circles
void shapeInversion(inout vec3 p, inout float der) {
  vec3 minp = p;

  // get point position _relative_ to shape center
  p -= invShape.center;
  minp -= minShape.center;

  // get distance from point P to center of invShape, minShape 
  float p2InvCenterSqr = dot(p,p);
  float p2MinCenterSqr = dot(minp, minp);

  // get distance from shape surface intersection point S to shape center
  float s2InvCenterSqr = getShape2CenterSqr(p, invShape);
  float s2MinCenterSqr = getShape2CenterSqr(p, minShape);;

  // determine if point P is inside invShape, minShape 
  bool insideMinShape = (p2MinCenterSqr < s2MinCenterSqr);
  bool insideInvShape = (p2InvCenterSqr < s2InvCenterSqr);

  if (insideMinShape) {
    // if point is within min shape do inversion _but_ as if point were at edge of min shape
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
Exposure = 0.7625
Brightness = 1.652542
Contrast = 1.136364
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
Fog = 3.116371
Iterations = 17
ColorIterations = 3
Bailout = 1000
BoxFoldLimit = 1
InvType = 6
InvScale = 2.528
InvCenter = 0,0,0
InvRotation = 45,-19.07834,45
InvParamA = 1.141593
MinType = 0
MinScale = 1.502137
MinCenter = 0,0,0
MinRotation = 0,0,0
MinParamA = 2
Scale = 1.913828
MinScalingConstant = 1
RotAngle = 0
RotVector = 1,1,1
#endpreset

