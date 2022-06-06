#ifndef _examples_double_bump_multipagecanvas_h
#define _examples_double_bump_multipagecanvas_h

#include <string>

#include <TCanvas.h>

class MultiPageCanvas : public TCanvas {
private:
  std::string fName;
public:
  MultiPageCanvas(const std::string& name) : fName(name) {
    if (!fName.empty()) {
      TCanvas::Print((fName + "[").c_str());
    }
  }

  ~MultiPageCanvas() {
    if (!fName.empty()) {
      TCanvas::Print((fName + "]").c_str());
    }
  }

  void Print(bool isLogx = false, bool isLogy = false) {
    if (!fName.empty()) {
      TCanvas::SetLogx(isLogx);
      TCanvas::SetLogy(isLogy);
      TCanvas::Print(fName.c_str());
    }
  }
};

#endif