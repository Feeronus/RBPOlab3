#include <iostream>
import Math;
using namespace std;
int main()
{
	double re, im;
	cout << "Input Real and Imag num:" << endl;
	cin >> re >> im;
	Complex c_num(re, im);

	double nom, denom;
	cout << "Input Nominator and Denominator num:" << endl;
	cin >> nom >> denom;
	Rational r_num(nom, denom);
	double d_num;
	cout << "Input double:" << endl;
	cin >> d_num;

	Rational answ = f(r_num);
	double a = answ.Nominator(), b = answ.Denominator();

	cout << "f(" << c_num << ") = " << f(c_num) << endl;
	cout << "f(" << r_num << ") = " << f(r_num) << " = " << a / b << endl;
	cout << "f(" << d_num << ") = " << f(d_num) << endl;

	return 0;
}