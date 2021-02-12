#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <GL/glut.h>
#include <GL/freeglut_ext.h>

#include "curves.hh"

// ***********
// * Options *
// ***********

size_t degree = 3;
enum ApproximationType {
  INTERP_PIEGL_SIMPLE, INTERP_PIEGL_PARABOLA, INTERP_SKETCHES, APPROXIMATE,
  PROXIMITY, PROXIMITY_FIT, PROXIMITY_DISPLACEMENT
} approximation_type = INTERP_SKETCHES;
size_t const approximation_complexity = 5;
bool show_curvature = false;
bool centripetal_parameterization = false;
int const resolution = 800;	// number of points in the curve
int const width = 800;
int const height = 600;
int const tolerance = 8;	// largest allowed hor/ver deviation in pixels
int const curvature_sparsity = 1;
std::string message;

// IDs for the menu elements / keyboard shortcuts
enum MenuCommand { MENU_RESET, MENU_DELETE_LAST, MENU_CUBIC, MENU_QUARTIC, MENU_CURVATURE,
                   MENU_INC_CURVATURE, MENU_DEC_CURVATURE,
                   MENU_PARAM_ARC, MENU_PARAM_CENTRIPETAL,
                   MENU_INTERPOLATE_SIMPLE, MENU_INTERPOLATE_PARABOLA, MENU_INTERPOLATE_SKETCHES,
                   MENU_APPROXIMATE,
                   MENU_PROXIMITY, MENU_PROXIMITY_FIT, MENU_PROXIMITY_DISPLACEMENT,
                   MENU_INC_DEPTH, MENU_DEC_DEPTH,
                   MENU_INC_ALPHA, MENU_DEC_ALPHA, MENU_DEFAULT_ALPHA,
                   MENU_LOAD, MENU_SAVE, MENU_PRINT, MENU_QUIT };
enum LoadSave { NOTHING, LOADING, SAVING } loadsave = NOTHING;

// ********************
// * Global variables *
// ********************

BSplineCurve curve;
int dragging = -1;		// -1 means no dragging, otherwise the point#
PointVector points;
struct SavedPoints {
  bool save;
  DoubleVector params;
  PointVector points;
} saved_points;
double curvature_magnification = 0.1;
size_t depth = 0;
double alpha = 0;

// ***********
// * Display *
// ***********

Point getObjectCoordinates(int x, int y)
{
  double model[16], proj[16];
  int view[4];
  double rx, ry, rz;
  glGetDoublev(GL_MODELVIEW_MATRIX, model);
  glGetDoublev(GL_PROJECTION_MATRIX, proj);
  glGetIntegerv(GL_VIEWPORT, view);
  gluUnProject(x, view[3] - y, 0, model, proj, view, &rx, &ry, &rz);
  return Point(rx, ry);
}

Point getWindowCoordinates(Point p)
{
  double model[16], proj[16];
  int view[4];
  double rx, ry, rz;
  glGetDoublev(GL_MODELVIEW_MATRIX, model);
  glGetDoublev(GL_PROJECTION_MATRIX, proj);
  glGetIntegerv(GL_VIEWPORT, view);
  gluProject(p.x, p.y, p.z, model, proj, view, &rx, &ry, &rz);
  return Point(rx, view[3] - ry);
}

Point interpolate(Point a, double t, Point b)
{
  a.x += (b.x - a.x) * t;
  a.y += (b.y - a.y) * t;
  a.z += (b.z - a.z) * t;
  return a;
}

void drawCurve()
{
  // knot vector
  Point min(-0.9, -0.9, 0.0), max(0.9, -0.9, 0.0);
  glLineWidth(2.0);
  glColor3d(0.0, 0.0, 1.0);
  glBegin(GL_LINES);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(max.x, max.y, max.z);
  glEnd();
  glPointSize(8.0);
  glBegin(GL_POINTS);
  for(int i = 0, ie = curve.knots.size(); i != ie; ++i) {
    Point const tmp = interpolate(min, curve.knots[i], max);
    glColor3d(0.0, 0.4, 0.4);
    glVertex3d(tmp.x, tmp.y, tmp.z);
  }
  glEnd();

  // interpolated point boxes
  glColor3d(0.0, 0.4, 0.4);
  glPointSize(10.0);
  glBegin(GL_POINTS);
  for(PointVector::const_iterator i = points.begin(); i != points.end(); ++i)
    glVertex3d(i->x, i->y, i->z);
  glEnd();

  if(curve.cp.empty())
    return;

  // curvature
  if (show_curvature) {
    glLineWidth(1.0);
    glColor3d(1.0, 0.0, 1.0);
    glBegin(GL_LINES);
    for(int i = 0; i < resolution; i += curvature_sparsity) {
      double const t = (double)i / ((double)resolution - 1.0);
      if(curve.knots[degree] <= t && t <= curve.knots[curve.knots.size() - degree - 1]) {
        VectorVector der;
        Point p = curve.derivatives(t, 2, der);
        double const curvature = (der[1] ^ der[2]).norm() / std::pow(der[1].norm(), 3.0);
        Vector perp(der[1].y, -der[1].x);
        perp.normalize();
        if (perp * der[2] >= 0)
          perp *= -1.0;
        perp *= curvature * curvature_magnification;
        glVertex3d(p.x, p.y, p.z);
        glVertex3d(p.x + perp.x, p.y + perp.y, p.z + perp.z);
      }
    }
    glEnd();
  }

  // control polygon
  glLineWidth(2.0);
  glColor3d(0.8, 0.5, 0.0);
  glBegin(GL_LINE_STRIP);
  for(PointVector::const_iterator i = curve.cp.begin(); i != curve.cp.end(); ++i)
    glVertex3d(i->x, i->y, i->z);
  glEnd();

  // control point boxes
  glColor3d(1.0, 0.0, 0.0);
  glPointSize(10.0);
  glBegin(GL_POINTS);
  for(PointVector::const_iterator i = curve.cp.begin(); i != curve.cp.end(); ++i)
    glVertex3d(i->x, i->y, i->z);
  glEnd();

  // the curve itself
  glLineWidth(3.0);
  glColor3d(0.0, 1.0, 0.0);
  glBegin(GL_LINE_STRIP);
  for(int i = 0; i < resolution; ++i) {
    double const t = (double)i / ((double)resolution - 1.0);
    if(curve.knots[degree] <= t && t <= curve.knots[curve.knots.size() - degree - 1]) {
      Point p = curve.evaluate(t);
      if(saved_points.save) {
	saved_points.params.push_back(t);
	saved_points.points.push_back(p);
      }
      glVertex3d(p.x, p.y, p.z);
    }
  }
  glEnd();

  // Segments with inflections
  glColor3d(0.0, 0.5, 0.0);
  glLineWidth(3.0);
  std::list<BezierCurve> list = curve.convertToBezierCurves();
  for (std::list<BezierCurve>::const_iterator i = list.begin(), ie = list.end(); i != ie; ++i)
    if (i->hasInflections()) {
      glBegin(GL_LINE_STRIP);
      for(size_t j = 0, je = resolution / list.size(); j < je; ++j) {
        double const t = (double)j / (double)(je - 1);
        Point p = i->evaluate(t);
        glVertex3d(p.x, p.y, p.z);
      }
      glEnd();
    }
}

void drawInfo()
{
  std::string s_degree, s_param, s_interp;
  if (degree == 3)
    s_degree = "Deg: 3";
  else
    s_degree = "Deg: 4";
  if (centripetal_parameterization)
    s_param = "Par: cpt";
  else
    s_param = "Par: arc";
  switch (approximation_type) {
  case INTERP_PIEGL_SIMPLE: s_interp = "Int: nor"; break;
  case INTERP_PIEGL_PARABOLA: s_interp = "Int: par"; break;
  case INTERP_SKETCHES: s_interp = "Int: ske"; break;
  case APPROXIMATE: s_interp = "Approx."; break;
  case PROXIMITY:
    {
      std::ostringstream s;
      s << "Deg: " << curve.n;
      s_degree = s.str();
    }
    {
      std::ostringstream s;
      s << "Alp: ";
      if (alpha)
        s << alpha;
      else
        s << "1/3";
      s_param = s.str();
    }
    s_interp = "Pro: sub";
    break;
  case PROXIMITY_FIT:
    {
      std::ostringstream s;
      s << "Deg: " << curve.n;
      s_degree = s.str();
    }
    {
      std::ostringstream s;
      s << "Smo: " << alpha;
      s_param = s.str();
    }
    s_interp = "Pro: fit";
    break;
  case PROXIMITY_DISPLACEMENT:
    {
      std::ostringstream s;
      s << "Deg: " << curve.n;
      s_degree = s.str();
    }
    {
      std::ostringstream s;
      s << "Alp: " << alpha;
      s_param = s.str();
    }
    s_interp = "Pro: dis";
    break;
  }
  glColor3d(0.0, 0.0, 0.0);
  glRasterPos2f(-0.98, 0.92);
  glutBitmapString(GLUT_BITMAP_HELVETICA_18, (unsigned char *)s_degree.c_str());
  glRasterPos2f(-0.98, 0.85);
  glutBitmapString(GLUT_BITMAP_HELVETICA_18, (unsigned char *)s_param.c_str());
  glRasterPos2f(-0.98, 0.78);
  glutBitmapString(GLUT_BITMAP_HELVETICA_18, (unsigned char *)s_interp.c_str());
  glRasterPos2f(-0.98, 0.71);
  glutBitmapString(GLUT_BITMAP_HELVETICA_18, (unsigned char *)message.c_str());
  message = "";
}

void display()
{
  glMatrixMode(GL_MODELVIEW);
  glClearColor(1.0, 1.0, 1.0, 0.0);
  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  drawInfo();
  glLoadIdentity();
  drawCurve();
  glutSwapBuffers();
}

// ******************
// * User interface *
// ******************

void uniform_knots(int n)
{
  curve.knots.clear();
  for(int i = 0; i < n; ++i)
    curve.knots.push_back((double)i / (double)(n - 1.0));
}

void reset()
{
  saved_points.save = false;
  points.clear();
  curve.p = degree;
  curve.n = 0;
  curve.cp.clear();
  uniform_knots(2);
}

void printData()
{
  std::cout << "B-spline curve:" << std::endl;
  std::cout << "\tDegree: " << degree << std::endl;
  std::cout << "\tKnots: [";
  for(DoubleVector::const_iterator i = curve.knots.begin(); i != curve.knots.end(); ++i) {
    if (i != curve.knots.begin())
      std::cout << ", ";
    std::cout << *i;
  }
  std::cout << "]\n\tControls: [";
  for(PointVector::const_iterator i = curve.cp.begin(); i != curve.cp.end(); ++i) {
    if (i != curve.cp.begin())
      std::cout << ", ";
    std::cout << "[" << i->x << ", " << i->y << "]";
  }
  std::cout << "]" << std::endl;
}

void reconstructCurve()
{
  if (points.size() < degree)
    return;
  switch(approximation_type) {
  case INTERP_PIEGL_SIMPLE:
    curve = BSplineCurve::interpolate(points, std::min(degree, points.size() - 1),
                                      centripetal_parameterization);
    break;
  case INTERP_PIEGL_PARABOLA:
    curve = BSplineCurve::interpolateWithParabolaRule(points, degree, centripetal_parameterization);
    break;
  case INTERP_SKETCHES:
    curve = BSplineCurve::interpolateAsInSketches(points, degree, centripetal_parameterization);
    break;
  case APPROXIMATE:
    curve = BSplineCurve::approximate(points, std::min(degree, points.size() - 1),
                                      std::min(approximation_complexity, points.size() - 1),
                                      centripetal_parameterization);
    break;
  case PROXIMITY:
    curve = BSplineCurve::proximity(points, depth, alpha);
    break;
  case PROXIMITY_FIT:
    curve = BSplineCurve::proximityFit(points, depth, alpha);
    break;
  case PROXIMITY_DISPLACEMENT:
    curve = BSplineCurve::proximityDisplacement(points, depth, alpha);
    break;
  }
}

void executeCommand(int command)
{
  switch(static_cast<MenuCommand>(command)) {
  case MENU_RESET: reset(); break;
  case MENU_DELETE_LAST: if (!points.empty()) points.pop_back(); reconstructCurve(); break;
  case MENU_CUBIC: degree = 3; reconstructCurve(); break;
  case MENU_QUARTIC: degree = 4; reconstructCurve(); break;
  case MENU_CURVATURE: show_curvature = !show_curvature; break;
  case MENU_INC_CURVATURE: curvature_magnification *= 2; glutPostRedisplay(); break;
  case MENU_DEC_CURVATURE: curvature_magnification /= 2; glutPostRedisplay(); break;
  case MENU_PARAM_ARC: centripetal_parameterization = false; reconstructCurve(); break;
  case MENU_PARAM_CENTRIPETAL: centripetal_parameterization = true; reconstructCurve(); break;
  case MENU_INTERPOLATE_SIMPLE: approximation_type = INTERP_PIEGL_SIMPLE; reconstructCurve(); break;
  case MENU_INTERPOLATE_PARABOLA:
    approximation_type = INTERP_PIEGL_PARABOLA; reconstructCurve(); break;
  case MENU_INTERPOLATE_SKETCHES: approximation_type = INTERP_SKETCHES; reconstructCurve(); break;
  case MENU_APPROXIMATE: approximation_type = APPROXIMATE; reconstructCurve(); break;
  case MENU_PROXIMITY: approximation_type = PROXIMITY; depth = alpha = 0; reconstructCurve(); break;
  case MENU_PROXIMITY_FIT: approximation_type = PROXIMITY_FIT; depth = alpha = 0;
    reconstructCurve(); break;
  case MENU_PROXIMITY_DISPLACEMENT: approximation_type = PROXIMITY_DISPLACEMENT; depth = alpha = 0;
    reconstructCurve(); break;
  case MENU_INC_DEPTH: depth++; reconstructCurve(); break;
  case MENU_DEC_DEPTH: if (depth) depth--; reconstructCurve(); break;
  case MENU_DEC_ALPHA:
    alpha = (alpha > 0) ? alpha - 0.05 : 0.5;
    if (std::abs(alpha) < 1e-5)
      alpha = 0;
    reconstructCurve();
    break;
  case MENU_INC_ALPHA:
    alpha = (alpha < 1) ? alpha + 0.05 : 0.5;
    if (std::abs(alpha) < 1e-5)
      alpha = 0;
    reconstructCurve();
    break;
  case MENU_DEFAULT_ALPHA: alpha = alpha ? 0 : 1; reconstructCurve(); break;
  case MENU_LOAD: loadsave = LOADING; message = "Load from slot (0-9)"; glutPostRedisplay(); break;
  case MENU_SAVE: loadsave = SAVING; message = "Save in slot (0-9)"; glutPostRedisplay(); break;
  case MENU_PRINT: printData(); break;
  case MENU_QUIT: exit(0);
  }
  glutPostRedisplay();
}

void loadFile(const char *filename) {
  std::ifstream f(filename);
  if (!f.good()) {
    message = "Cannot open file.";
    glutPostRedisplay();
    return;
  }
  size_t n;
  f >> n;
  points.clear();
  points.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    double x, y;
    f >> x >> y;
    points.push_back(Point(x, y));
  }
  f.close();
  reconstructCurve();
  message = std::string("Loaded from slot ") + filename[13] + ".";
}

void saveFile(const char *filename) {
  std::ofstream f(filename);
  if (!f.good()) {
    message = "Cannot open file.";
    glutPostRedisplay();
  }
  f << points.size();
  for (PointVector::const_iterator i = points.begin(), ie = points.end(); i != ie; ++i)
    f << ' ' << i->x << ' ' << i->y;
  f.close();
  message = std::string("Saved in slot ") + filename[13] + ".";
}

void keyboard(unsigned char key, int x, int y)
{
  switch (loadsave) {
  case NOTHING:
    switch (key) {
    case 'r' : executeCommand(MENU_RESET); break;
    case 'd' : executeCommand(MENU_DELETE_LAST); break;
    case '3' : executeCommand(MENU_CUBIC); break;
    case '4' : executeCommand(MENU_QUARTIC); break;
    case 'C' : executeCommand(MENU_CURVATURE); break;
    case '*' : executeCommand(MENU_INC_CURVATURE); break;
    case '/' : executeCommand(MENU_DEC_CURVATURE); break;
    case 'a' : executeCommand(MENU_PARAM_ARC); break;
    case 'c' : executeCommand(MENU_PARAM_CENTRIPETAL); break;
    case 'n' : executeCommand(MENU_INTERPOLATE_SIMPLE); break;
    case 'p' : executeCommand(MENU_INTERPOLATE_PARABOLA); break;
    case 's' : executeCommand(MENU_INTERPOLATE_SKETCHES); break;
    case 'A' : executeCommand(MENU_APPROXIMATE); break;
    case 'b' : executeCommand(MENU_PROXIMITY); break;
    case 'f' : executeCommand(MENU_PROXIMITY_FIT); break;
    case 'D' : executeCommand(MENU_PROXIMITY_DISPLACEMENT); break;
    case '+' : executeCommand(MENU_INC_DEPTH); break;
    case '-' : executeCommand(MENU_DEC_DEPTH); break;
    case '8' : executeCommand(MENU_DEC_ALPHA); break;
    case '9' : executeCommand(MENU_INC_ALPHA); break;
    case '0' : executeCommand(MENU_DEFAULT_ALPHA); break;
    case 'L' : executeCommand(MENU_LOAD); break;
    case 'S' : executeCommand(MENU_SAVE); break;
    case 'P' : executeCommand(MENU_PRINT); break;
    case 'q' : executeCommand(MENU_QUIT); break;
    }
    break;
  default:
    if (key >= '0' && key <= '9') {
      char filename[] = "quartic-test-?.clubs";
      filename[13] = key;
      if (loadsave == LOADING)
        loadFile(filename);
      else
        saveFile(filename);
      loadsave = NOTHING;
      glutPostRedisplay();
    }
  }
}

int nearestPoint(PointVector v, int x, int y)
{
  // Returns the index of the point in v near to (x, y),
  // and -1 if no point is near enough
  int result = -1;

  for(size_t i = 0, ie = v.size(); i < ie; ++i) {
    Point p = getWindowCoordinates(v[i]);
    if(std::fabs(p.x - x) < tolerance && fabs(p.y - y) < tolerance)
      result = i;
  }

  return result;
}

void mouseButton(int button, int state, int x, int y)
{
  if(button == GLUT_LEFT_BUTTON && state == GLUT_UP)
    dragging = -1;
  else if(button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
    if((dragging = nearestPoint(points, x, y)) == -1) {
      // insert new point
      Point p = getObjectCoordinates(x, y);
      points.push_back(p);
      reconstructCurve();
      glutPostRedisplay();
    }
  } else if(button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN) {
    saved_points.save = true;
    display();
    saved_points.save = false;
    int nearest = nearestPoint(saved_points.points, x, y);
    if(nearest >= 0)
      std::cout << "Selected: " << saved_points.params[nearest] << std::endl;
  }
}

void mouseMotion(int x, int y)
{
  if(dragging < 0)
    return;
  points[dragging] = getObjectCoordinates(x, y);
  reconstructCurve();
  glutPostRedisplay();
}

int buildPopupMenu()
{
  int menu = glutCreateMenu(executeCommand);
  glutAddMenuEntry("Reset\t(r)", MENU_RESET);
  glutAddMenuEntry("Delete last point\t(d)", MENU_DELETE_LAST);
  glutAddMenuEntry("Toggle curvature\t(C)", MENU_CURVATURE);
  glutAddMenuEntry("Scale up curvature\t(*)", MENU_INC_CURVATURE);
  glutAddMenuEntry("Scale down curvature\t(/)", MENU_DEC_CURVATURE);
  glutAddMenuEntry("----", -1);
  glutAddMenuEntry("Use cubic curve\t(3)", MENU_CUBIC);
  glutAddMenuEntry("Use quartic curve\t(4)", MENU_QUARTIC);
  glutAddMenuEntry("----", -1);
  glutAddMenuEntry("Arc-length parameterization\t(a)", MENU_PARAM_ARC);
  glutAddMenuEntry("Centripetal parameterization\t(c)", MENU_PARAM_CENTRIPETAL);
  glutAddMenuEntry("----", -1);
  glutAddMenuEntry("Interpolation: Piegl/normal\t(n)", MENU_INTERPOLATE_SIMPLE);
  glutAddMenuEntry("Interpolation: Piegl/parabola\t(p)", MENU_INTERPOLATE_PARABOLA);
  glutAddMenuEntry("Interpolation: Sketches\t(s)", MENU_INTERPOLATE_SKETCHES);
  glutAddMenuEntry("----", -1);
  glutAddMenuEntry("Approximation\t(A)", MENU_APPROXIMATE);
  glutAddMenuEntry("----", -1);
  glutAddMenuEntry("Proximity Bezier curve\t(b)", MENU_PROXIMITY);
  glutAddMenuEntry("Proximity Bezier fit\t(f)", MENU_PROXIMITY_FIT);
  glutAddMenuEntry("Proximity w/displacement\t(D)", MENU_PROXIMITY_DISPLACEMENT);
  glutAddMenuEntry("Increase proximity\t(+)", MENU_INC_DEPTH);
  glutAddMenuEntry("Decrease proximity\t(-)", MENU_DEC_DEPTH);
  glutAddMenuEntry("Decrease prox. alpha\t(8)", MENU_DEC_ALPHA);
  glutAddMenuEntry("Increase prox. alpha\t(9)", MENU_INC_ALPHA);
  glutAddMenuEntry("Default prox. alpha\t(0)", MENU_DEFAULT_ALPHA);
  glutAddMenuEntry("----", -1);
  glutAddMenuEntry("Load clubs\t(L)", MENU_LOAD);
  glutAddMenuEntry("Save clubs\t(S)", MENU_SAVE);
  glutAddMenuEntry("Print curve data\t(P)", MENU_PRINT);
  glutAddMenuEntry("Quit\t(q)", MENU_QUIT);
  return menu;
}

// ******************
// * Initialization *
// ******************

int main(int argc, char *argv[])
{
  if(argc > 1)
    degree = atoi(argv[1]);

  reset();

  // glut window initialization
  glutInit(&argc, argv);
  glutInitWindowSize(width, height);
  glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
  glutCreateWindow("BSpline Interpolation");

  // register callbacks
  glutDisplayFunc(display);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouseButton);
  glutMotionFunc(mouseMotion);

  // create popup menu
  buildPopupMenu();
  glutAttachMenu(GLUT_RIGHT_BUTTON);

  // turn the flow of control over to glut
  glutMainLoop ();

  return 0;
}
