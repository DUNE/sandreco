#include "SANDDisplayUtils.h"
#include <TBox.h>
#include <TEllipse.h>

//----------------------------------------------------------------------------
void SANDDisplayUtils::DrawEllipse(TEllipse *ellipse, Color_t color,
                                   Style_t style)
{
  ellipse->SetFillStyle(style);        // set ellipse fill style
  ellipse->SetLineColor(color);        // set ellipse color
  ellipse->SetFillColor(color);        // set ellipse color
  ellipse->SetBit(TBox::kCannotMove);  // avoid moving ellipse
  ellipse->SetBit(kCanDelete);         // delete during canvas clear
  ellipse->Draw();                     // draw box
}

//----------------------------------------------------------------------------
void SANDDisplayUtils::DrawBox(TBox *box, Color_t color, Style_t style)
{
  box->SetFillStyle(style);        // set box fill style
  box->SetLineColor(color);        // set box color
  box->SetFillColor(color);        // set box color
  box->SetBit(TBox::kCannotMove);  // avoid moving box
  box->SetBit(kCanDelete);         // delete during canvas clear
  box->Draw();                     // draw box
}
