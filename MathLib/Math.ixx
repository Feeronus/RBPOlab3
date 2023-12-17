module;
#include <math.h>
#include <ostream>
export module Math;

export class Complex
{
private:
	double m_mod;
	double m_arg;
public:
	Complex(double _mod)
	{
		m_mod = _mod;
		m_arg = 0;
	}
	Complex()
	{
		m_mod = 0;
		m_arg = 0;
	}
	Complex(double _re, double _im)
	{
		m_mod = sqrt(_re * _re + _im * _im);
		m_arg = atan2(_im, _re);
	}

	static Complex FromExponentialForm(double _mod, double _arg)
	{
		Complex exp_obj;
		exp_obj.m_mod = _mod;
		if (_mod < 1E-5)
			exp_obj.m_mod = 0;
		exp_obj.m_arg = _arg;
		if (_arg < 1E-5)
			exp_obj.m_arg = 0;
		return exp_obj;
	}
	static Complex FromAlgebraicForm(double _re, double _im)
	{
		Complex alg_obj(_re, _im);
		return alg_obj;
	}
	double Re() const
	{
		return m_mod*cos(m_arg);
	}
	double Im() const
	{
		return m_mod*sin(m_arg);
	}
	double Mod() const
	{
		return m_mod;
	}
	double Arg() const
	{
		return m_arg;
	}

	explicit operator double() const {
		return m_mod*cos(m_arg);
	}
	Complex operator-()
	{
		Complex obj(*this);
		obj.m_mod = sqrt(pow(obj.Re()*(-1), 2) + pow(obj.Im()*(-1), 2));
		obj.m_arg = atan2(obj.Im() * (-1), obj.Re() * (-1));
		return obj;
	}
	Complex& operator++()
	{
		double _re = Re() + 1;
		double _im = Im();
		m_mod = sqrt(pow(_re,2)+pow(_im,2));
		m_arg = atan2(_im, _re);
		return (*this);
	}
	Complex operator++(int _post_inc)
	{
		Complex obj(*this);
		++(*this);
		return obj;
	}
	Complex& operator--()
	{
		double _re = Re() - 1;
		double _im = Im();
		m_mod = sqrt(pow(_re, 2) + pow(_im, 2));
		m_arg = atan2(_im, _re);
		return (*this);
	}
	Complex operator--(int _post_dec)
	{
		Complex obj(*this);
		--(*this);
		return obj;
	}
	Complex& operator+=(Complex _obj)
	{

		m_mod = sqrt(pow(Re()+_obj.Re(),2) + pow(Im() + _obj.Im(), 2));
		m_arg = atan2(Im() + _obj.Im(), Re() + _obj.Re());
		return (*this);
	}
	Complex& operator-=(Complex _obj) {
		double _re = Re() - _obj.Re();
		double _im = Im() - _obj.Im();
		m_mod = sqrt(pow(_re, 2) + pow(_im, 2));
		m_arg = atan2(_im, _re);
		return (*this);
	}
	Complex& operator*=(Complex _obj) {
		double temp_re = Re();
		double temp_im = Arg();
		double _re = temp_re * _obj.Re() - temp_im * _obj.Im();
		double _im = temp_re * _obj.Im() + temp_im * _obj.Re();
		m_mod = sqrt(pow(_re, 2) + pow(_im, 2));
		m_arg = atan2(_im, _re);
		return (*this);
	}
	Complex& operator/=(Complex _obj) {
		double temp_re1 = Re(), temp_im1 = Im();
		double temp_re2 = _obj.Re(), temp_im2 = _obj.Im();
		double _re = (temp_re1 * temp_re2 + temp_im1 * temp_im2) / (pow(temp_re2, 2) + pow(temp_im2, 2));
		double _im = (temp_re2 * temp_im1 - temp_re1 * temp_im2) / (pow(temp_re2, 2) + pow(temp_im2, 2));
		m_mod = sqrt(pow(_re, 2) + pow(_im, 2));
		m_arg = atan2(_im, _re);
		return (*this);
	}
	friend Complex operator+ (const Complex& _obj1, const Complex& _obj2);
	friend Complex operator- (const Complex& _obj1, const Complex& _obj2);
	friend Complex operator* (const Complex& _obj1, const Complex& _obj2);
	friend Complex operator/ (const Complex& _obj1, const Complex& _obj2);

	friend Complex operator ""i(long double _im);
	friend Complex operator ""i(unsigned long long _im);

	friend std::ostream& operator<<(std::ostream& stream, const Complex& _obj);
};
export Complex operator+(const Complex& _obj1, const Complex& _obj2)
{
	return Complex(_obj1.Re() + _obj2.Re(), _obj1.Im() + _obj2.Im());
}
export Complex operator-(const Complex& _obj1, const Complex& _obj2)
{
	double _re1 = _obj1.m_mod * cos(_obj1.m_arg);
	double _re2 = _obj2.m_mod * cos(_obj2.m_arg);
	double _im1 = _obj1.m_mod * sin(_obj1.m_arg);
	double _im2 = _obj2.m_mod * sin(_obj2.m_arg);
	double _re = round((_re1 - _re2) * 1E5) / 1E5;
	double _im = round((_im1 - _im2) * 1E5) / 1E5;
	return Complex(_re, _im);
}
export Complex operator*(const Complex& _obj1, const Complex& _obj2)
{
	return Complex((_obj1.Re() * _obj2.Re() - _obj1.Im() * _obj2.Im()),
		(_obj1.Re() * _obj2.Im() + _obj1.Im() * _obj2.Re()));
}
export Complex operator/(const Complex& _obj1, const Complex& _obj2)
{
	double _re = ((_obj1.Re() * _obj2.Re() + _obj1.Im() * _obj2.Im()) / (pow(_obj2.Re(), 2) + pow(_obj2.Im(), 2)));
	double _im = ((_obj2.Re() * _obj1.Im() - _obj1.Re() * _obj2.Im()) / (pow(_obj2.Re(), 2) + pow(_obj2.Im(), 2)));
	return Complex(_re,_im);
}
export Complex operator""i(long double _im)
{
	return Complex(0.0, static_cast<double>(_im));
}
export Complex operator""i(unsigned long long _im)
{
	return Complex(0.0, static_cast<double>(_im));
}
export std::ostream& operator<<(std::ostream& stream, const Complex& _obj)
{
	if (_obj.Im() < 0)
	{
		stream << _obj.Re() << " " << _obj.Im() << "i";
	}
	else
	{
		stream << _obj.Re() << " + " << _obj.Im() << "i";
	}
	return stream;
}

export int FindGreatestCommonDivisor(int a, int b)
{
	int r;
	if (a < 0)
		a *= -1;
	if (b < 0)
		b *= -1;
	while (true)
	{
		if (b == 0)
			return a;
		r = a % b;
		a = b;
		b = r;
	}
}
export int FindLeastCommonMultiple(int x, int y) {
	return abs(x * y) / FindGreatestCommonDivisor(x, y);
}

//========================================================================================

export class Rational {
	int m_nominator;
	int m_denominator;
	
public:
	void normalize()
	{
		int nod = FindGreatestCommonDivisor(m_nominator, m_denominator);
		m_nominator /= nod;
		m_denominator /= nod;
		if (m_denominator < 0) {
			m_denominator *= -1;
			m_nominator *= -1;
		}
	}

	Rational()
	{
		m_nominator = 0;
		m_denominator = 1;
	}
	Rational(int _nominator, int _denominator) {
		m_denominator = _denominator;
		m_nominator = _nominator;
		normalize();
	}
	Rational(int _nominator) {
		m_nominator = _nominator;
		m_denominator = 1;
	}
	int Nominator() const {
		return m_nominator;
	}
	int Denominator() const {
		return m_denominator;
	}
	explicit operator double() const {
		return double(m_nominator) / m_denominator;
	}
	Rational operator-() {
		Rational obj(*this);
		obj.m_nominator *= -1;
		return obj;
	}
	Rational& operator++ () {
		m_nominator += m_denominator;
		return (*this);
	}
	Rational operator++ (int _param) {
		Rational obj(*this);
		(*this).m_nominator += m_denominator;
		return obj;
	}
	Rational& operator-- () {
		m_nominator -= m_denominator;
		return (*this);
	}
	Rational operator-- (int _param) {
		Rational obj(*this);
		(*this).m_nominator -= m_denominator;
		return obj;
	}
	Rational& operator+=(Rational _obj) {
		int new_den = FindLeastCommonMultiple(m_denominator, _obj.m_denominator);
		m_nominator = new_den / m_denominator * m_nominator;
		m_nominator += new_den / _obj.m_denominator * _obj.m_nominator;
		m_denominator = new_den;
		normalize();
		return (*this);
	}
	Rational& operator-=(Rational _obj) {
		int new_d = FindGreatestCommonDivisor(m_denominator, _obj.m_denominator);
		m_nominator = new_d / m_denominator * m_nominator;
		m_nominator -= new_d / _obj.m_denominator * _obj.m_nominator;
		m_denominator = new_d;
		normalize();
		return (*this);
	}
	Rational& operator*=(Rational _obj) {
		m_denominator *= _obj.m_denominator;
		m_nominator *= _obj.m_nominator;
		normalize();
		return (*this);
	}
	Rational& operator/=(Rational _obj) {
		m_denominator *= _obj.m_nominator;
		m_nominator *= _obj.m_denominator;
		normalize();
		return (*this);
	}
	friend Rational operator+ (const Rational& _obj1, const Rational& _obj2);
	friend Rational operator- (const Rational& _obj1, const Rational& _obj2);
	friend Rational operator* (const Rational& _obj1, const Rational& _obj2);
	friend Rational operator/(const Rational& _obj1, const Rational& _obj2);

	friend bool operator==(const Rational& _obj1, const Rational& _obj2);
	friend bool operator>(const Rational& _obj1, const Rational& _obj2);
	friend bool operator<(const Rational& _obj1, const Rational& _obj2);
	friend bool operator>=(const Rational& _obj1, const Rational& _obj2);
	friend bool operator<=(const Rational& _obj1, const Rational& _obj2);

	friend std::ostream& operator<<(std::ostream& stream, const Rational& _obj);
};

export Rational operator+ (const Rational& _obj1, const Rational& _obj2) {
	int denominator = FindLeastCommonMultiple(_obj1.m_denominator, _obj2.m_denominator);
	int nominator = denominator / _obj1.m_denominator * _obj1.m_nominator;
	nominator += denominator / _obj2.m_denominator * _obj2.m_nominator;
	return Rational{ nominator, denominator };
}

export Rational operator-(const Rational& _obj1, const Rational& _obj2)
{
	int denominator = FindLeastCommonMultiple(_obj1.m_denominator, _obj2.m_denominator);
	int nominator = denominator / _obj1.m_denominator * _obj1.m_nominator;
	nominator -= denominator / _obj2.m_denominator * _obj2.m_nominator;
	return Rational{ nominator, denominator };
}

export Rational operator*(const Rational& _obj1, const Rational& _obj2)
{
	return Rational{ _obj1.m_nominator * _obj2.m_nominator, _obj2.m_denominator * _obj1.m_denominator };
}

export Rational operator/(const Rational& _obj1, const Rational& _obj2)
{
	return Rational{ _obj1.m_nominator * _obj2.m_denominator,_obj1.m_denominator * _obj2.m_nominator };
}

export bool operator==(const Rational& _obj1, const Rational& _obj2)
{
	return _obj1.m_nominator == _obj2.m_nominator && _obj1.m_denominator == _obj2.m_denominator;
}

export bool operator>(const Rational& _obj1, const Rational& _obj2)
{
	int den = FindLeastCommonMultiple(_obj1.m_denominator, _obj2.m_denominator);
	return den / _obj1.m_denominator * _obj1.m_nominator > den / _obj2.m_denominator * _obj2.m_nominator;
}
export bool operator<(const Rational& _obj1, const Rational& _obj2)
{
	int den = FindLeastCommonMultiple(_obj1.m_denominator, _obj2.m_denominator);
	return den / _obj1.m_denominator * _obj1.m_nominator < den / _obj2.m_denominator * _obj2.m_nominator;
}
export bool operator>=(const Rational& _obj1, const Rational& _obj2)
{
	int den = FindLeastCommonMultiple(_obj1.m_denominator, _obj2.m_denominator);
	return den / _obj1.m_denominator * _obj1.m_nominator >= den / _obj2.m_denominator * _obj2.m_nominator;
}
export bool operator<=(const Rational& _obj1, const Rational& _obj2)
{
	int den = FindLeastCommonMultiple(_obj1.m_denominator, _obj2.m_denominator);
	return den / _obj1.m_denominator * _obj1.m_nominator <= den / _obj2.m_denominator * _obj2.m_nominator;
}

export std::ostream& operator<<(std::ostream& stream, const Rational& _obj) {
	stream << _obj.m_nominator << "|" << _obj.m_denominator;
	return stream;
}

//-----------------------------------------------------------
export Complex f(const Complex& z)
{
	Complex a(1, 0);
	
	Complex result = a * z * z - cos(2*double(z));
	return result;
}

export Rational f(const Rational& r) 
{
	Rational a(1, 1);
	Rational cosinus(round(cos(2 * double(r)) * 10000), 10000);
	Rational result = a * r * r - cosinus;
	return result;
}

export double f(const double& d) 
{
	double a = 1.0;
	double result = a * d * d - cos(2 * d);
	return result;
}