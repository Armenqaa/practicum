#include <math.h>
#include <stdio.h>
#include <fstream>
#include <iostream>

using namespace std;

double a21(1. / 5);
double a31(3. / 40);
double a32(9. / 40);
double a41(44. / 45);
double a42(-56. / 15);
double a43(32. / 9);
double a51(19372. / 6561);
double a52(-25360. / 2187);
double a53(64448. / 6561);
double a54(-212. / 729);
double a61(9017. / 3168);
double a62(-355. / 33);
double a63(46732. / 5247);
double a64(49. / 176);
double a65(-5103. / 18656);
double a71(35. / 384);
double a72(0.);
double a73(500. / 1113);
double a74(125. / 192);
double a75(-2187. / 6784);
double a76(11. / 84);
double c2(1. / 5);
double c3(3. / 10);
double c4(4. / 5);
double c5(8. / 9);
double c6(1.);
double c7(1.);
double b1(5179. / 57600);
double b2(0.);
double b3(7571. / 16695);
double b4(393. / 640);
double b5(-92097. / 339200);
double b6(187. / 2100);
double b7(1. / 40);
double lowB1(35. / 384);
double lowB2(0.);
double lowB3(500. / 1113);
double lowB4(125. / 192);
double lowB5(-2187. / 6784);
double lowB6(11. / 84);
double tol(1.e-12);
double fac(0.9);
double facmin(1.e-10);
double facmax(3.);

using namespace std;

double FuncY(double t, double y, double x) {
    return -x;
}

double FuncX(double t, double y, double x) {
    return y;
}

void FindK(double t, double &y, double &x, double h,  double& k1Y, double& k2Y, double& k3Y, double& k4Y, double& k5Y, double& k6Y, double& k1X, double& k2X, double& k3X, double& k4X, double& k5X, double& k6X) {
    k1X = FuncX(t, y, x);
    k1Y = FuncY(t, y, x);
    k2X = FuncX(t + c2 * h, y + h * a21 * k1Y, x + h * a21 * k1X);
    k2Y = FuncY(t + c2 * h, y + h * a21 * k1Y, x + h * a21 * k1X);
    k3X = FuncX(t + c3 * h, y + h * (a31 * k1Y + a32 * k2Y), x + h * (a31 * k1X + a32 * k2X));
    k3Y = FuncY(t + c3 * h, y + h * (a31 * k1Y + a32 * k2Y), x + h * (a31 * k1X + a32 * k2X));
    k4X = FuncX(t + c4 * h, y + h * (a41 * k1Y + a42 * k2Y + a43 * k3Y), x + h * (a41 * k1X + a42 * k2X + a43 * k3X));
    k4Y = FuncY(t + c4 * h, y + h * (a41 * k1Y + a42 * k2Y + a43 * k3Y), x + h * (a41 * k1X + a42 * k2X + a43 * k3X));
    k5X = FuncX(t + c5 * h, y + h * (a51 * k1Y + a52 * k2Y + a53 * k3Y + a54 * k4Y), x + h * (a51 * k1X + a52 * k2X + a53 * k3X + a54 * k4X));
    k5Y = FuncY(t + c5 * h, y + h * (a51 * k1Y + a52 * k2Y + a53 * k3Y + a54 * k4Y), x + h * (a51 * k1X + a52 * k2X + a53 * k3X + a54 * k4X));
    k6X = FuncX(t + c6 * h, y + h * (a61 * k1Y + a62 * k2Y + a63 * k3Y + a64 * k4Y + a65 * k5Y), x + h * (a61 * k1X + a62 * k2X + a63 * k3X + a64 * k4X + a65 * k5X));
    k6Y = FuncY(t + c6 * h, y + h * (a61 * k1Y + a62 * k2Y + a63 * k3Y + a64 * k4Y + a65 * k5Y), x + h * (a61 * k1X + a62 * k2X + a63 * k3X + a64 * k4X + a65 * k5X));
    y = y + h * (lowB1 * k1Y + lowB2 * k2Y + lowB3 * k3Y + lowB4 * k4Y + lowB5 * k5Y + lowB6 * k6Y);
    x = x + h * (lowB1 * k1X + lowB2 * k2X + lowB3 * k3X + lowB4 * k4X + lowB5 * k5X + lowB6 * k6X);
}

double StepChoice(double t, double &y, double &x, double h) {
    // cout << "-----------------------" << endl;
    for (;;) {
        // cout << h << endl;
        double locX = x, locY = y, wX = x, wY = y, k1Y, k2Y, k3Y, k4Y, k5Y, k6Y, k1X, k2X, k3X, k4X, k5X, k6X;
        FindK(t, locY, locX, h, k1Y, k2Y, k3Y, k4Y, k5Y, k6Y, k1X, k2X, k3X, k4X, k5X, k6X);
        double k7X = FuncX(t + c7 * h, y + h * (a71 * k1Y + a72 * k2Y + a73 * k3Y + a74 * k4Y + a75 * k5Y + a76 * k6Y), x + h * (a71 * k1X + a72 * k2X + a73 * k3X + a74 * k4X + a75 * k5X + a76 * k6X));
        double k7Y = FuncY(t + c7 * h, y + h * (a71 * k1Y + a72 * k2Y + a73 * k3Y + a74 * k4Y + a75 * k5Y + a76 * k6Y), x + h * (a71 * k1X + a72 * k2X + a73 * k3X + a74 * k4X + a75 * k5X + a76 * k6X));
        wX = x + h * (b1 * k1X + b2 * k2X + b3 * k3X + b4 * k4X + b5 * k5X + b6 * k6X + b7 * k7X);
        wY = y + h * (b1 * k1Y + b2 * k2Y + b3 * k3Y + b4 * k4Y + b5 * k5Y + b6 * k6Y + b7 * k7Y);
        double d1 = fabs(locX + (locX - wX) / 31.);
        double d2 = fabs(locY + (locY - wY) / 31.);
        double err = max(fabs(locY - wY) / d2, fabs(locX - wX) / d1) / 31.;
        h *= min(facmax, max(facmin, fac * pow((tol / err), 1. / 6)));
        if (err <= tol) {
            x = locX;
            y = locY;
            return h;
        }
    }
}

int main() {
    ofstream out;
    int draw = 0;
    double nextH = 0.1, x = 0., y = 1., currH = 0.1;
    out.open("auto_rk.txt");
    for (double T = 0.; T <= 100000 * 3.1415926535897932; T += currH) {
        currH = nextH;
        if (!(draw % 3000)) {
            cout << x - sin(T) << " " << y - cos(T) << endl;
        }
        nextH = StepChoice(T, y, x, currH);
        draw++;
    }

    out.close();
}
