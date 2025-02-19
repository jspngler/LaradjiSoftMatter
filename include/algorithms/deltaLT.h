#ifndef MPD_DELTALT
#define MPD_DELTALT

#include <cmath>
#include <iostream>
#include "dataTypes.h"

namespace mpd {

template <typename T>
struct deltaLT {
	//The system timestep and this length step, deltaT and deltaL,
	// together decide how fast this converges to the end length, endL.
	// RelaxStep decides how often it runs. The difference between the 
	// initial size.s[dim] and endL decide how many steps it takes.
	//E.g.
	// size = {25 55 100}
	// deltaT = 0.02
	// deltaL = 0.01
	// endL = 50 
	// dim = 2 (maybe Z or z in mpd file)
	// relaxStep = 10
	//Would result in 10*(100-50)/0.01=50000 steps taken to get
	// from 100 to 50, or 50000*0.02=1000 tau.
	//This would likely increase the system temperature rapidly
	// if the standard for tau is about 19 ns. Volume would be 1/2
	// within 19 us. A more reasonable deltaL would likely be 0.0001, 
	// shrinking the volume within 1.9 ms (100000 tau).
	//Once the Z dimension reaches that size, the system will 
	// no longer change size.
	T deltaL,endL;
	int dim, relaxStep;
	
	//Just ensure they are all zero, we will use deltaL==0 to indicate inactive
	deltaLT():deltaL(0),endL(0),dim(0),relaxStep(0){}
	
	//Return a new size based on the activity
	constexpr threeVector<T> newSize(threeVector<T> oldSize)
	{
		T diff=oldSize.s[dim]-endL;
		//Check if it is within the last step
		if(std::abs(diff)>deltaL)
		{
			T direction=diff/std::abs(diff);
			oldSize.s[dim]-=deltaL*direction;
		}
		else
			//this will make a smaller step to the end length
			oldSize.s[dim]=endL;
		
		return oldSize;
	}
	
	//this exists somewhere else, or at least it should. 
	// I should probably make a canned version of this
	constexpr threeVector<T> scaleFactor(threeVector<T> oldSize)
	{
		auto nextSize=newSize(oldSize);
		
		threeVector<T> scale=1.0;
		
		scale.s[dim]=1.0+(nextSize.s[dim]-oldSize.s[dim])/oldSize.s[dim];
		
		return scale;
	}
	
	constexpr bool ready(int step)
	{
		return (deltaL!=0 && step%relaxStep==0);
	}
	
	constexpr bool active()
	{
		return deltaL!=0;
	}
	
	constexpr int nWords() const {return 4;}
	
	std::istream &inStep(std::istream &stream, int wStep)
	{
		switch(wStep)
		{
			case 0: stream >> deltaL; break;
			case 1: stream >> endL; break;
			case 2: stream >> dim; 
				if(dim<0 || dim>2)
				{
					std::cerr << "deltaLT dimension, " << dim << ", out of bounds!" << std::endl;
					throw 1;
				}
				break;
			case 3: stream >> relaxStep; break;
			default:
				std::cerr << "Forgot to reset input word counter!";
				throw 1;
				break;
		}
		return stream;
	}
	
	std::ostream &outStep(std::ostream &stream, int wStep)
	{
		switch(wStep)
		{
			case 0: stream << deltaL << ' '; break;
			case 1: stream << endL << ' '; break;
			case 2: stream << dim << ' '; break;
			case 3: stream << relaxStep << '\n'; break;
			default:
				std::cerr << "Forgot to reset output word counter!";
				throw 1;
				break;
		}
		return stream;
	}
};

}
#endif
