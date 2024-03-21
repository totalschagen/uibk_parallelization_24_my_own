#ifndef LIMITER_HPP
#define LIMITER_HPP

class limiter_base {
public:
	limiter_base();
	// TBD by students -> use sensible parameter list
	virtual double compute(double first, double second, double third) = 0;
};

class limiter_minmod : public limiter_base {
public:
	limiter_minmod(double theta);
	// TBD by students -> use sensible parameter list
	double compute(double first, double second, double third);

private:
	double theta;
};

#endif