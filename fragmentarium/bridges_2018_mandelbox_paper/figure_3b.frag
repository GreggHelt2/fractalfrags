// Mandelbox2d ShapeInversion experiments
//     created by Gregg Helt
// 2D version of the Mandelbox that allows replacing of circle inversion with shape inversion where shape types are:
//   0 => circle
//   1 => ellipse
//   2 => polygon
//   3 => superellipse
//   4 => FG-squircle (NOT the same as superellipse squircle)
// in addition to parameters to set the inversion shape:
//   replaced the min radius (used to determine whether to apply inversion or scaling) with parameterized minShape
//   added new shape maxShape to determine outer bound of region where inversion is applied
//   invShape, minShape, maxShape can all be rotated and moved
//   added color overlays to show inversion shape, minshape, maxshape, boxfold params, and more

#define providesInit
#info Mandelbox2D_ShapeInversion_v8

#include "Progressive2DColoring.frag"

#group ParamOverlay
uniform bool ShowInversionOverlay; checkbox[false]
uniform bool ShowBoxfoldOverlay; checkbox[false]
uniform bool ShowTestPoint; checkbox[false]
uniform bool ShowShapeIntersect; checkbox[false]
uniform vec2 TestPoint; slider[(-10,-10),(0.3846,1.5),(10,10)]
uniform vec3 TestPointColor; color[1,0,1]

uniform vec3 InvShapeColor; color[1,1,0]
uniform vec3 MinShapeColor; color[1,0,0]
uniform vec3 MaxShapeColor; color[0,1,0]
uniform vec3 BoxColor; color[0,1,1]
uniform vec3 BoxShiftColor; color[1,0,1]
uniform vec3 ScaleColor; color[1,0.541176,0.215686]
uniform vec3 OffsetColor; color[0.745098,0.411765,0.203922]
uniform float OverlayLineWidth; slider[2,5,30]
uniform float OverlayPointSize; slider[5,10,60]

// TestPoint Mode == 0, test point is point before boxfold
// TestPoint Mode == 1, test point is point after boxfold (allows more direct testing of inversion effect)
uniform int TestPointMode; slider[0,0,1]    


#group InversionShapes
                       // shape types:
                       // 0 => circle
                       // 1 => ellipse
                       // 2 => polygon
                       // 3 => superellipse
		       // 4 => FG-squircle (NOT the same as superellipse squircle)
uniform int InvShapeType; slider[0,0,4]
uniform vec2 InvCenter; slider[(-5,-5),(0,0),(5,5)]
uniform float InvShapeScale; slider[-5,1,5]
uniform float InvShapeRotate; slider[-1,0,1]
uniform float InvShapeParamA; slider[-5,1,20]
uniform float InvShapeParamB; slider[-5,4,20]
uniform float InvShapeParamC; slider[0.0001,1.5,8]

uniform int MinShapeType; slider[0,0,4]
uniform vec2 MinShapeCenter; slider[(-5,-5),(0,0),(5,5)]  
uniform float MinShapeScale; slider[-5,0.3,5]
uniform float MinShapeRotate; slider[-1,0,1]
uniform float MinShapeParamA; slider[-5,1,20]
uniform float MinShapeParamB; slider[-5,4,20]
uniform float MinShapeParamC; slider[0.0001,1.5,8]

uniform int MaxShapeType; slider[0,0,4]
uniform vec2 MaxShapeCenter; slider[(-5,-5),(0,0),(5,5)]
uniform float MaxShapeScale; slider[-5,1,5]
uniform float MaxShapeRotate; slider[-1,0,1]
uniform float MaxShapeParamA; slider[-5,1,20]
uniform float MaxShapeParamB; slider[-5,4,20]
uniform float MaxShapeParamC; slider[0.0001,1.5,8]


#group Mandelboxen
uniform float Scale; slider[-5,2,5]
uniform float BoxFoldLimit; slider[-5,1,5]
uniform float BoxFoldShift; slider[-5,0,5]
uniform vec2 BoxCenter; slider[(-5,-5),(0,0),(5,5)]
uniform bool UseBoxFoldShift; checkbox[true]
uniform bool UseMinScalingConstant; checkbox[false]
uniform float MinScalingConstant; slider[-10,2,10]
uniform bool UseMaxScalingConstant; checkbox[false]
uniform float MaxScalingConstant; slider[-10,1,10]
uniform bool UseMaxShape; checkbox[false]
uniform float ZFactor1; slider[-5,0,5]
uniform float ZFactor2; slider[-5,0,5]
uniform bool MinMaxCenterRelative; checkbox[true]

const float PI =  3.14159265359;
const float PI2 = 6.28318530718;

struct Shape {
  int stype;
  vec2 origin;
  float scale;
  float rotation; // in radians
  vec3 drawColor;
  float a;
  float b;
  float c;
};

Shape minShape;
Shape maxShape;
Shape invShape;
float z1, z2;
float overlayLineCoordWidth;
float overlayPointCoordRadius;
float boxFoldLow;
float boxFoldHigh;
float boxMinLow;
float boxMaxHigh;

void init() {
  z1 = ZFactor1 * ZFactor1;
  if (ZFactor1 < 0.0) { z1 = -z1; }
  z2 = ZFactor2 * ZFactor2;
  if (ZFactor2 < 0.0) { z2 = -z2; }
  // OverlayLineWidth is pixels
  // trying to keep overlay constant pixel size
  overlayLineCoordWidth = OverlayLineWidth * pixelSize.x / zoomscale / 2.0;
  overlayPointCoordRadius = OverlayPointSize * pixelSize.x / zoomscale;

  invShape = Shape(InvShapeType, InvCenter,
                   InvShapeScale, (InvShapeRotate * PI), InvShapeColor,
                   InvShapeParamA, InvShapeParamB, InvShapeParamC);
  minShape = Shape(MinShapeType,
                   MinMaxCenterRelative ? (InvCenter + MinShapeCenter): MinShapeCenter,
                   MinShapeScale, (MinShapeRotate * PI), MinShapeColor,
                   MinShapeParamA, MinShapeParamB, MinShapeParamC);
  maxShape = Shape(MaxShapeType,
                   MinMaxCenterRelative ? (InvCenter + MaxShapeCenter): MaxShapeCenter,
                   MaxShapeScale, (MaxShapeRotate * PI), MaxShapeColor,
                   MaxShapeParamA, MaxShapeParamB, MaxShapeParamC);
  boxFoldHigh = BoxFoldLimit;
  boxFoldLow = -BoxFoldLimit;
  if (UseBoxFoldShift) {
    boxMaxHigh = boxFoldHigh + BoxFoldShift;
    boxMinLow = boxFoldLow - BoxFoldShift;
  }
  else {
    boxMaxHigh = boxFoldHigh;
    boxMinLow = boxFoldLow;
  }
}


void boxFold(inout vec2 p) {
  p -= BoxCenter;
  float distanceToFold;
  if (p.x > boxMaxHigh) {
    distanceToFold = abs(p.x - boxFoldHigh);
    if  (p.x > boxFoldHigh) { p.x = boxFoldHigh - distanceToFold; }
    else {  p.x = boxFoldHigh + distanceToFold; } // p.x <= boxFoldHigh
  }
  else if (p.x < boxMinLow) {
    distanceToFold = abs(p.x - boxFoldLow);
    if (p.x < boxFoldLow) { p.x = boxFoldLow + distanceToFold; }
    else { p.x = boxFoldLow - distanceToFold; } // p.x >= boxFoldLow
  }
  if (p.y > boxMaxHigh) {
    distanceToFold = abs(p.y - boxFoldHigh);
    if  (p.y > boxFoldHigh) { p.y = boxFoldHigh - distanceToFold; }
    else {  p.y = boxFoldHigh + distanceToFold; } // p.y <= boxFoldHigh
  }
  else if (p.y < boxMinLow) {
    distanceToFold = abs(p.y - boxFoldLow);
    if (p.y < boxFoldLow) { p.y = boxFoldLow + distanceToFold; }
    else { p.y = boxFoldLow - distanceToFold; } // p.y >= boxFoldLow
  }
  p += BoxCenter;
}


vec2 getShapeIntersect(vec2 p, Shape shp) {
  float t = atan(p.y, p.x);  // theta angle
  t += shp.rotation;
  vec2 intersectP = vec2(1.0,1.0);
  float r = 1.0;
  if (shp.stype == 0) { // circle
    r = shp.scale;
  }
  else if (shp.stype == 1) {  // ellipse
    float sinA2 = shp.a * shp.a * sin(t) * sin(t);
    float cosB2 = shp.b * shp.b * cos(t) * cos(t);
    r = shp.scale * (shp.a * shp.b) / sqrt(sinA2 + cosB2);
  }
  else if (shp.stype == 2) { // polygon
    // param A => number of sides to the polygon
    float n = floor(shp.a);
    r = cos(PI/n) / cos(mod(t,PI2/n) - PI/n);
    r *= shp.scale;;
  }
  else if (shp.stype == 3) { // superellipse
    float sinAPow = pow(abs(shp.a*sin(t)), shp.c);
    float cosBPow = pow(abs(shp.b*cos(t)), shp.c);
    r = (shp.a * shp.b) / pow((sinAPow + cosBPow), 1.0/shp.c);
    r *= shp.scale;
  }
  else if (shp.stype == 4) {  // Fernandez Guasti squircle
    // calcs based on Fong, 2016: https://arxiv.org/abs/1604.02174
    //    with r ==> shape.scale and s ==> shape.a
    // NOT the same as the superellipse squircle! (superellipse with a=b and c=4)
    // for distinction between two definitions, http://mathworld.wolfram.com/Squircle.html
    float s = shp.a;
    float sin2t = sin(2.0*t);
    if (s == 0.0 || sin2t == 0.0 || sin(t) == 0.0 || cos(t) == 0.0) { r = shp.scale; } // special-casing circle to avoid division by zero
    else {
      float x = (shp.scale * sign(cos(t)) / (s * sqrt(2.0) * abs(sin(t)))) * sqrt(1.0 - sqrt(1.0 - s*s*sin2t*sin2t));
      float y = (shp.scale * sign(sin(t)) / (s * sqrt(2.0) * abs(cos(t)))) * sqrt(1.0 - sqrt(1.0 - s*s*sin2t*sin2t));
      if (x == 0.0) { y = shp.scale; }
      r = sqrt(x*x + y*y);
    }
  }
  else {  // out of range for shape types, should never get here but if do, just return a constant
    r = 1.0;
  }
  t -= shp.rotation;
  intersectP = vec2(r*cos(t), r*sin(t));
  return intersectP;
}

// modified inversion that can use shapes other than circles
void shapeInversion(inout vec2 p) {
  vec2 minp = vec2(p.x, p.y);
  vec2 maxp = vec2(p.x, p.y);
  p -= invShape.origin;
  minp -= minShape.origin;
  maxp -= maxShape.origin;
  float p2InvOriginSqr = dot(p,p);
  float p2MinOriginSqr = dot(minp, minp);
  float p2MaxOriginSqr = dot(maxp, maxp);
  vec2 invIntersect = getShapeIntersect(p, invShape);
  float p2InvEdgeSqr = dot(invIntersect, invIntersect);
  
  vec2 minIntersect = getShapeIntersect(minp, minShape);
  float p2MinEdgeSqr = dot(minIntersect, minIntersect);
  vec2 maxIntersect;
  float p2MaxEdgeSqr;
  if (UseMaxShape) {
    maxIntersect = getShapeIntersect(maxp, maxShape);
    p2MaxEdgeSqr = dot(maxIntersect, maxIntersect);
  }
  else { // don't use maxShape, consider max same as invShape
    maxIntersect = invIntersect;
    p2MaxEdgeSqr = p2InvEdgeSqr;
  }
  p2InvOriginSqr += z1;
  p2InvEdgeSqr += z2;

  bool insideMinShape = (p2MinOriginSqr < p2MinEdgeSqr);
  bool insideInvShape = (p2InvOriginSqr < p2InvEdgeSqr);
  bool insideMaxShape = (p2MaxOriginSqr < p2MaxEdgeSqr);
  
  if (insideMinShape) {
    if (UseMinScalingConstant)  {
      // if point is within min shape and constant option on, scale by constant
      p *= MinScalingConstant;
    } else {
      // if point is within min shape and constant option off,
      //    do inversion but as if point were at edge of min shape
      p *= p2InvEdgeSqr / p2MinEdgeSqr;
    }
  }
  else if (insideMaxShape) {
    // if point is outside minShape but with maxShape, do shape inversion relative to invShape
    p *= p2InvEdgeSqr/p2InvOriginSqr;
  }
  else if (UseMaxScalingConstant) {
    // try doing a conditional scaling if point is outside of maxShape?
    p *= MaxScalingConstant;
  }
  p += invShape.origin;
}

// Mandelbox2D iteration function
//   (called iteratively in Progessive2DColoring.frag)
// for every iteration:
//   p' = (Scale(Inversion(Boxfold(p))) + c
vec2 formula(vec2 p, vec2 c) {
  boxFold(p);
  shapeInversion(p);
  p *= Scale;
  p += c;
  return p;
}


bool colorShape(vec2 p, Shape shp, inout vec3 color) {
  if (!ShowInversionOverlay) { return false; }
  vec2 pc = p;
  pc -= shp.origin;  
  vec2 shpIntersect = getShapeIntersect(pc, shp);
  if (distance(pc, shpIntersect) <= overlayLineCoordWidth) {
    color = shp.drawColor;
    return true;
  }
  else {
    return false;
  }
}

bool colorBox(vec2 p, inout vec3 color) {
  vec2 pc = p;
  pc -= BoxCenter;
  bool hit = false;
  if ((abs(pc.x - BoxFoldLimit) <= overlayLineCoordWidth) ||
      (abs(pc.y - BoxFoldLimit) <= overlayLineCoordWidth) ||
      (abs(pc.x + BoxFoldLimit) <= overlayLineCoordWidth) ||
      (abs(pc.y + BoxFoldLimit) <= overlayLineCoordWidth)) {
    color = BoxColor;
    hit = true;
  }
  else if (UseBoxFoldShift) {
    if ((abs(pc.x - boxMaxHigh) <= overlayLineCoordWidth) || 
        (abs(pc.y - boxMaxHigh) <= overlayLineCoordWidth) ||
        (abs(pc.x - boxMinLow) <= overlayLineCoordWidth) ||
        (abs(pc.y - boxMinLow) <= overlayLineCoordWidth)) {
      color = BoxShiftColor;
      hit = true;
    }
  }
  return hit;
}

bool colorTestPoint(vec2 p, inout vec3 color) {
  bool hit = false;
  vec2 offset = vec2(TestPoint.x, TestPoint.y);
  vec2 q = vec2(TestPoint.x, TestPoint.y);

  
  if (TestPointMode == 0) {
    if (distance(q, p) < overlayPointCoordRadius) {
      color = TestPointColor;
      hit = true;
    }
    if (!hit) {
      boxFold(q);
      if (distance(q, p) < overlayPointCoordRadius) {
        color = BoxColor;
        hit = true;
      }
    }
  }
  else {
    if (distance(q, p) < overlayPointCoordRadius) {
      color = BoxColor;
      hit = true;
    }
  }

  if (!hit && ShowShapeIntersect) {
    // test for nearness to intersection with ray origin->boxFold(q)
    q -= invShape.origin;
    p -= invShape.origin;
    vec2 ishape = getShapeIntersect(q, invShape);

    if (distance(ishape,p) < overlayPointCoordRadius) {
      color = invShape.drawColor;
      hit = true;
    }
    q += invShape.origin;
    p += invShape.origin;
  }



  if (!hit) {
    shapeInversion(q);
    if (distance(q, p) < overlayPointCoordRadius) {
      color = invShape.drawColor;
      hit = true;
    }
  }
  if (!hit) {
    q *= Scale;
    if (distance(q, p) < overlayPointCoordRadius) {
      color = ScaleColor;
      hit = true;
    }
  }
  if (!hit) {
    q += offset;
    if (distance(q, p) < overlayPointCoordRadius) {
      color = OffsetColor;
      hit = true;
    }
  }
  return hit;
}


bool applyColorOverlay(vec2 p, inout vec3 color) {
  bool hit = false;
  if (ShowTestPoint) {
    hit = colorTestPoint(p, color);
  }
  if (!hit && ShowInversionOverlay) {
    hit = colorShape(p, invShape, color);
    if (!hit) { hit = colorShape(p, minShape, color); }
    if (!hit && UseMaxShape) { hit = colorShape(p, maxShape, color); }
  }
  if (!hit && ShowBoxfoldOverlay) {
    hit = colorBox(p, color);
  }
  return hit;
}

#preset Default
Gamma = 1.2
Brightness = 1.184669
Contrast = 0.9966777
Saturation = 1.014799
Center = 0,-0.0041119
Zoom = 0.5007709
ToneMapping = 1
Exposure = 1.078224
AARange = 2
AAExp = 1
GaussianAA = true
Iterations = 7
PreIterations = 1
R = 1
G = 0.1518625
B = 0
C = 0.1485714
EscapeSize = 0.68
ColoringType = 1
ColorFactor = 0.15
FlattenEscapeColor = true
GrayScale = false
Fade = 1
EscapeColor = 0,0,0
ShowInversionOverlay = false
ShowBoxfoldOverlay = false
ShowTestPoint = false
ShowShapeIntersect = false
TestPoint = 1.581197,1.367521
TestPointColor = 1,0,1
InvShapeColor = 1,1,0
MinShapeColor = 1,0,0
MaxShapeColor = 0,1,0
BoxColor = 0,1,1
BoxShiftColor = 1,0,1
ScaleColor = 1,0.541176,0.215686
OffsetColor = 0.745098,0.411765,0.203922
OverlayLineWidth = 30
OverlayPointSize = 20
TestPointMode = 0
InvShapeType = 4
InvCenter = 0,0
InvShapeScale = 0.4883721
InvShapeRotate = 0
InvShapeParamA = 0.6
InvShapeParamB = 2.981221
InvShapeParamC = 4.120141
MinShapeType = 2
MinShapeCenter = 0,0
MinShapeScale = 0.4460094
MinShapeRotate = 0
MinShapeParamA = 6.458333
MinShapeParamB = 1
MinShapeParamC = 1
MaxShapeType = 0
MaxShapeCenter = 0,0
MaxShapeScale = 1.933638
MaxShapeRotate = -0.2494172
MaxShapeParamA = 0.9
MaxShapeParamB = 4
MaxShapeParamC = 1.5
Scale = -1.666667
BoxFoldLimit = 0.8550186
BoxFoldShift = 0
BoxCenter = 0,0
UseBoxFoldShift = false
UseMinScalingConstant = false
MinScalingConstant = 1.604938
UseMaxScalingConstant = false
MaxScalingConstant = 1
UseMaxShape = false
ZFactor1 = -0.3556485
ZFactor2 = 0.5439331
MinMaxCenterRelative = false
#endpreset
