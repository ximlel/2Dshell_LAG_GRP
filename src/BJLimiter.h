/*********************************************************************
/* File name: BJLimiter.h
/* 
/* Abstract: Barth-Jespersen Limiter
/*
/* Date: Apr 24, 2018
/* Author: Jian Cheng @ IAPCM
***********************************************************************/
#ifndef __BJLIMITER__H__
#define __BJLIMITER__H__
#include <algorithm>
#include <cmath>

double BJLimiter(double W_max, double W_min, double W_c, double W_xp)
{
	double eps(1.0e-8), mu(0);

	if( W_xp - W_c > eps )
		mu = std::min( 1.0, (W_max-W_c)/(W_xp - W_c) );
	else if( W_xp - W_c < -eps )
		mu = std::min( 1.0, (W_min-W_c)/(W_xp - W_c) );
	else
		mu = 1.0;

	return mu;
}

#endif