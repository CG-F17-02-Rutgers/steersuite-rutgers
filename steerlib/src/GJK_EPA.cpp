#include "obstacles/GJK_EPA.h"
#include "util/Geometry.h"
SteerLib::GJK_EPA::GJK_EPA(){}

//Look at the GJK_EPA.h header file for documentation and instructions



//Everything starts here!
bool SteerLib::GJK_EPA::intersect(float& return_penetration_depth, 
								Util::Vector& return_penetration_vector, 
								const std::vector<Util::Vector>& _shapeA, 
								const std::vector<Util::Vector>& _shapeB){


	//Define an empty point set W. Contains the points in the current test simplex.
	std::vector<Util::Vector> simplexW;

	//Run GJK Algorithm. Returns true if there is a collison.
	bool collisionExist = GJK(_shapeA, _shapeB, simplexW);

	if (collisionExist == true){
		//Collision exists. Calculate the penetration depth and penetration vector.
		EPA(return_penetration_depth, return_penetration_vector, _shapeA, _shapeB, simplexW);
		
		return true;
	}
	return false;
	// There is no collision
}


bool SteerLib::GJK_EPA::GJK(const std::vector<Util::Vector>& _shapeA, 
		const std::vector<Util::Vector>& _shapeB, 
		std::vector<Util::Vector>& simplexW) {

	//First starting point
	Util::Vector d(1, 0, 0);

	//Calculation for first Simplex point (Related Lecture Slides: Page 38, 50, 53)
	Util::Vector w0 = SimplexPointw(_shapeA, _shapeB, d);
	simplexW.push_back(w0);

	d = -d;

	while (true) {
		//Get the next point w and then add it to the simplex W.
		Util::Vector nextVector = SimplexPointw(_shapeA, _shapeB, d);
		simplexW.push_back(nextVector);

		//Check if this next w passes through origin
		float dotty = dot(nextVector, d);
		if (dotty < 0) {
			return false;
		}
		else {
			//the Simplex W contains the origin. We are done. Collision exists.
			if (containsOrigin(simplexW, d) ){
				return true;
			}
			//if it doesn't, repeats.
		}
	}
	return false;

}
bool SteerLib::GJK_EPA::EPA(float& return_penetration_depth, Util::Vector& return_penetration_vector, const std::vector<Util::Vector>& _shapeA, const std::vector<Util::Vector>& _shapeB, std::vector<Util::Vector>& simplexW) {
	return false;
}

bool SteerLib::GJK_EPA::containsOrigin(std::vector<Util::Vector>& simplexW, Util::Vector& d) {
	//declare all variable needed
	Util::Vector a;
	Util::Vector a0;
	Util::Vector b;
	Util::Vector c;
	Util::Vector ab;
	Util::Vector ac;
	Util::Vector abPerp;
	Util::Vector acPerp;

	//initialized
	a = simplexW.back();
	a0 = a * -1; 

	//algorithm
	if(simplexW.size() == 3){
		b = simplexW.front();
		c = simplexW.at(1);
		ab = b-a;
		ac = c-a;
		abPerp = cross(cross(ac,ab),ab);
		acPerp = cross(cross(ab,ac),ac);
		if(dot(abPerp,a0) > 0){
			simplexW.erase(simplexW.begin()+1);
			d = abPerp;
		}
		else if(dot(acPerp,a0) > 0){
			simplexW.erase(simplexW.begin());
			d = acPerp;
		}
		else{
			printf("ffffff\n");
			return true;
		}
	}
	else{
		b = simplexW.front();
		ab = b-a;
		abPerp = cross(cross(ac,ab),ab);
		d = abPerp;
	}
	return false;
}


Util::Vector SteerLib::GJK_EPA::SimplexPointw(const std::vector<Util::Vector>& _shapeA,
												const std::vector<Util::Vector>& _shapeB, 
												Util::Vector& d) {
	//Calculate Support Point(Furthest Point) for shapeA (Related Lecture Slide: Page 38)

	float highest = -FLT_MAX;
	Util::Vector supportA(0, 0, 0);

	for (int i=0; i<_shapeA.size(); i++) {
		Util::Vector v = _shapeA[i];
		float dotty = dot(v, d);
		if (dotty > highest) {
			highest = dotty;
			supportA = v;
		}
	}

	//Calculate Support Point(Furthest Point) for shapeB (Related Lecture Slide: Page 38)
	highest = -FLT_MAX;
	Util::Vector supportB(0, 0, 0);

	for (int i=0; i<_shapeA.size(); i++) {
		Util::Vector v = _shapeA[i];
		float dotty = dot(v, -d);
		if (dotty > highest) {
			highest = dotty;
			supportB = v;
		}
	}

	return supportA - supportB;
}
