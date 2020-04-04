#include "stdafx.h"



Angle::Angle()
{
	this->r=0;
	this->valid = false;
}


Angle::~Angle()
{
}

Angle Angle::operator+(Angle ang)
{
	Angle sum;
	bool leftOpConvex,rightOpConvex,flag=false;

	if ( ( this->getP_y().x() != ang.getP_x().x() ) || ( this->getP_y().y() != ang.getP_x().y() ) )
	{
		std::cout << "error adding angles\n";
		sum.lowValid();
		return sum;
	}

	sum.raiseValid();
	sum.setP_x(this->getP_x());
	sum.setX(this->getX());
	sum.setP_y(ang.getP_y());
	sum.setR(this->getR() + ang.getR());

	Triangle_2 tri1(this->getP_x(),this->getX(),this->getP_y());
	Triangle_2 tri2(ang.getP_x(),ang.getX(),ang.getP_y());
	Triangle_2 tri3(this->getP_x(),this->getX(),ang.getP_y());
	//Triangle_2 tri4(this->getP_x(),ang.getP_y(),this->getX());

	if ( tri1.orientation() <= 0 )
		leftOpConvex = true;
	else
		leftOpConvex = false;
	if ( tri2.orientation() <= 0 )
		rightOpConvex = true;
	else
		rightOpConvex = false;

	if ( !rightOpConvex && !leftOpConvex )
		flag = true;
	if ( !rightOpConvex && tri3.orientation() < 0 )
		flag = true;
	if ( !leftOpConvex && tri3.orientation() < 0 )
		flag = true;

	if ( flag )
		sum.setR( sum.getR() + 1 );

	return sum;
}

double Angle::getVal()
{/*
	double val=0,nv1,nv2;
	Vector_2 v1,v2;
	v1 = this->getP_x() - this->getX();
	nv1 = std::sqrt( v1.x()*v1.x() + v1.y()*v1.y() );
	v2 = this->getP_y() - this->getX();
	nv2 = std::sqrt( v2.x()*v2.x() + v2.y()*v2.y() );

	val = ( v1.x() * v2.x() + v1.y() * v2.y() ) / (nv1*nv2);
	//double dot_product = val;
	val = std::acos(val) * 180.0f / M_PI;

	double cross_product = ( v1.x()*v2.y() ) - ( v1.y()*v2.x() );
	if ( cross_product >= 0 )
		return (val);
	else
		return ( 360 - val );



		
		*/
	return 0;

}
