#pragma once
class Angle
{
private:
	Point_2 p_x,x,p_y;
	int r;	
	bool valid;
public:
	Angle();
	Angle( Point_2 p_x , Point_2 x ,Point_2 p_y, int r=0 ):p_x(p_x),x(x),p_y(p_y),r(r){this->valid = true;}
	~Angle();
	Point_2 getP_x() {return this->p_x;}
	Point_2 getP_y() {return this->p_y;}
	Point_2 getX() {return this->x;}
	int getR() {return this->r;}
	void setP_x( Point_2 p_x ){this->p_x = p_x;}
	void setP_y( Point_2 p_y ){this->p_y = p_y;}
	void setX( Point_2 x ){this->x = x;}
	void setR( int r ){this->r = r;}
	void raiseValid() {this->valid = true;}
	void lowValid(){this->valid = false;}
	bool isValid() {return this->valid;}
	Angle operator+(Angle ang);

	double getVal();
};

