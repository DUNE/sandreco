#ifndef SANDDISPLAYUTILS_H
#define SANDDISPLAYUTILS_H 1

#include <TROOT.h>

// using namespace std;

class TBox;
class TEllipse;

class SANDDisplayUtils {
 public:

  static void DrawBox(TBox *box, Color_t color, Style_t style = 1001);
  static void DrawEllipse(TEllipse *ellipse, Color_t color, Style_t style = 1001);

};

 #endif
