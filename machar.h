/*
 * machar.h
 *
 *  Created on: Oct 6, 2017
 *      Author: di
 */

#ifndef MACHAR_H_
#define MACHAR_H_
#include <iostream>
#include <limits>
#include <stdlib.h>

struct Machar {
	int ibeta, it, irnd, ngrd, machep, negep, iexp, minexp, maxexp;
	float eps, epsneg, xmin, xmax;
	Machar() {
		/*Determines machine-specific parameters affecting floating-point arithmetic, including ibeta,
		 the floating-point radix; it, the number of base-ibeta digits in the floating-point mantissa;
		 eps, the smallest positive number that, added to 1.0, is not equal to 1.0; epsneg, the
		 smallest positive number that, subtracted from 1.0, is not equal to 1.0; xmin, the smallest
		 representable positive number; and xmax, the largest representable positive number. See
		 text for description of other returned parameters. Change Doub to float to find single
		 precision parameters.*/
		int i, itemp, iz, j, k, mx, nxres;
		float a, b, beta, betah, betain, one, t, temp, temp1, tempa, two, y, z,
				zero;
		one = float(1);
		two = one + one;
		zero = one - one;
		a = one; //Determine ibeta and beta by the method of M.
		do { //Malcolm.
			a += a;
			temp = a + one;
			temp1 = temp - a;
		} while (temp1 - one == zero);
		b = one;
		do {
			b += b;
			temp = a + b;
			itemp = int(temp - a);
		} while (itemp == 0);
		ibeta = itemp;
		beta = float(ibeta);
		it = 0; //Determine it and irnd.
		b = one;
		do {
			++it;
			b *= beta;
			temp = b + one;
			temp1 = temp - b;
		} while (temp1 - one == zero);
		irnd = 0;
		betah = beta / two;
		temp = a + betah;
		if (temp - a != zero)
			irnd = 1;
		tempa = a + beta;
		temp = tempa + betah;
		if (irnd == 0 && temp - tempa != zero)
			irnd = 2;
		negep = it + 3; //Determine negep and epsneg.
		betain = one / beta;
		a = one;
		for (i = 1; i <= negep; i++)
			a *= betain;
		b = a;
		for (;;) {
			temp = one - a;
			if (temp - one != zero)
				break;
			a *= beta;
			--negep;
		}
		negep = -negep;
		epsneg = a;
		machep = -it - 3; //Determine machep and eps.
		a = b;
		for (;;) {
			temp = one + a;
			if (temp - one != zero)
				break;
			a *= beta;
			++machep;
		}
		eps = a;
		ngrd = 0; //Determine ngrd.
		temp = one + eps;
		if (irnd == 0 && temp * one - one != zero)
			ngrd = 1;
		i = 0; //Determine iexp.
		k = 1;
		z = betain;
		t = one + eps;
		nxres = 0;
		for (;;) { //Loop until an underflow occurs, then exit.
			y = z;
			z = y * y;
			a = z * one; //Check here for the underflow.
			temp = z * t;
			if (a + a == zero || abs(z) >= y)
				break;
			temp1 = temp * betain;
			if (temp1 * beta == z)
				break;
			++i;
			k += k;
		}
		if (ibeta != 10) {
			iexp = i + 1;
			mx = k + k;
		} else { //For decimal machines only.
			iexp = 2;
			iz = ibeta;
			while (k >= iz) {
				iz *= ibeta;
				++iexp;
			}
			mx = iz + iz - 1;
		}
		for (;;) { //To determine minexp and xmin, loop until an
			xmin = y; //underflow occurs, then exit.
			y *= betain;
			a = y * one; //Check here for the underflow.
			temp = y * t;
			if (a + a != zero && abs(y) < xmin) {
				++k;
				temp1 = temp * betain;
				if (temp1 * beta == y && temp != y) {
					nxres = 3;
					xmin = y;
					break;
				}
			} else
				break;
		}
		minexp = -k; //Determine maxexp, xmax.
		if (mx <= k + k - 3 && ibeta != 10) {
			mx += mx;
			++iexp;
		}
		maxexp = mx + minexp;
		irnd += nxres; //Adjust irnd to reflect partial underflow.
		if (irnd >= 2)
			maxexp -= 2; //Adjust for IEEE-style machines.
		i = maxexp + minexp;
		/*Adjust for machines with implicit leading bit in binary mantissa, and machines with
		 radix point at extreme right of mantissa.*/
		if (ibeta == 2 && !i)
			--maxexp;
		if (i > 20)
			--maxexp;
		if (a != y)
			maxexp -= 2;
		xmax = one - epsneg;
		if (xmax * one != xmax)
			xmax = one - beta * epsneg;
		xmax /= (xmin * beta * beta * beta);
		i = maxexp + minexp + 3;
		for (j = 1; j <= i; j++) {
			if (ibeta == 2)
				xmax += xmax;
			else
				xmax *= beta;
		}
	}
	void report() {
		std::cout << "quantity: numeric_limits<float> says (we calculate)" << std::endl;
		std::cout << "radix: " << std::numeric_limits<float>::radix << " (" << ibeta
				<< ")" << std::endl;
		std::cout << "mantissa digits: " << std::numeric_limits<float>::digits << " ("
				<< it << ")" << std::endl;
		std::cout << "round style: " << std::numeric_limits<float>::round_style << " ("
				<< irnd << ") [our 5 == IEEE 1]" << std::endl;
		std::cout << "guard digits: " << "[not in numeric_limits]" << " (" << ngrd
				<< ")" << std::endl;
		std::cout << "epsilon: " << std::numeric_limits<float>::epsilon() << " (" << eps
				<< ")" << std::endl;
		std::cout << "neg epsilon: " << "[not in numeric_limits]" << " (" << epsneg
				<< ")" << std::endl;
		std::cout << "epsilon power: " << "[not in numeric_limits]" << " (" << machep
				<< ")" << std::endl;
		std::cout << "neg epsilon power: " << "[not in numeric_limits]" << " ("
				<< negep << ")" << std::endl;
		std::cout << "exponent digits: " << "[not in numeric_limits]" << " (" << iexp
				<< ")" << std::endl;
		std::cout << "min exponent: " << std::numeric_limits<float>::min_exponent << " ("
				<< minexp << ")" << std::endl;
		std::cout << "max exponent: " << std::numeric_limits<float>::max_exponent << " ("
				<< maxexp << ")" << std::endl;
		std::cout << "minimum: " << std::numeric_limits<float>::min() << " (" << xmin
				<< ")" << std::endl;
		std::cout << "maximum: " << std::numeric_limits<float>::max() << " (" << xmax
				<< ")" << std::endl;
	}
};

#endif /* MACHAR_H_ */
