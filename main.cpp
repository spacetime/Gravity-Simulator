/* Code for Physics Engine 
 
 To be linked with Graphics Engine.
 
 
 Code Started: 24 December 2010
 
 Last updation: 21 Feb 2011
 
 */

/* MISSING FEATURES
 
 * Auto Console Logging
 
 
 
 
 
 
 */
#include <iostream>
#include <math.h>
#include <fstream>
//#define GRAV 0.00000000006673   // Gravitational Constant
//#define DT 0.00001				//Time frame dt

const double GRAV=6.67e-11;
bool debug_mode=1; //verbose
bool give_output=1; //displays results in console
//set to true only when above two are false!!!
bool formal_mode=1;

int time_limit=10000; //seconds

const double DT=0.00001;

//Global Variables Start
//Constant Variables to Store Net number of Stars. Restricted to capacity of 'int'and 'double' datatype.

int netstars=0;					//constant after input_net()
double netmass=0;				// Net mass. Constant after externmass();


////reusable for movement through time

double fx=0,fy=0,fz=0;			// forces along axes, will be set and reset by star.calc_force()
//double usablemass=0;			// set and reset by star.use_mass();
double COM=0,cx=0,cy=0,cz=0;	//COM = Centre of Mass of Universe; cx,cy,cz= components of COM along axes


//Time Control

long Global_Time=0;
//Global Variables End

//Needed Prototypes
double usable_mass(int);

/*
 ////////////Class: Star/////////////
 This class contains unique characteristics or the individual star.
 Operations will be done on a large array of such 'stars'. 
 */

class star{
private:
	
	
	
public:
	int SelfID;					// ID to identify the star
	double velx,vely,velz;		// represent velocities in respective axes (m/s)
	int radius;					// radius of the star (m)
	bool immovable;				// set as 1 to render a star fixed relative to absolute zero.
	double x,y,z;				// x,y,z are positions of star in these axes (m)
	double mass;				// mass of the star (in kg)
	int skipiter;				// number of iterations to skip. For reducing calculation.
	
	
	double use_mass(){
		//include provision for relativistic mass later
		return (netmass-mass);
		
	  }
	
	
	//sets the values of fx,fy,fz;

	void find_forces(){
	//step1: find value of FORCE
	//step2: find components along x,y,z;
		
	
		if (debug_mode==true) std::cout<<"\tCoordinates of COM2: "<<cx<<","<<cy<<","<<cz<<"\n";
		
		//COM2=sqrt((c2x*c2x)+(c2y*c2y)+(c2z*c2z)) //Distance of COM of universe-selected from origin
		
		//using F= G M m / r^2
		
		double dist_btw_selected_and_COM2;
		dist_btw_selected_and_COM2=sqrt((cx-x)*(cx-x) + (cy-y)*(cy-y) + (cy*y)*(cy*y));
		
		double Fnet=(GRAV*usable_mass(SelfID)*mass)/(dist_btw_selected_and_COM2*dist_btw_selected_and_COM2);
		if (debug_mode==true) std::cout<<"Distance btw selected star and COM2: "<<dist_btw_selected_and_COM2<<"\n";
		if (debug_mode==true) std::cout<<"Net Force (NEwtons): "<<Fnet<<"\n";
		//refer to annexure A
		
		double k=(Fnet)/(dist_btw_selected_and_COM2);
		
		
		fx=(cx-x)*k;
		fy=(cy-y)*k;
		fz=(cz-z)*k;
				
		if (debug_mode==true) std::cout<<"\tForces in Axes: "<<fx<<","<<fy<<","<<fz<<"\n";
	}
	
	void move_through_time(){
		double ax=fx/mass, ay=fy/mass, az=fz/mass;
		if (debug_mode==true) std::cout<<"Acceleration delta changes: "<<ax*DT<<","<<ay*DT<<","<<az*DT<<"\n";
		if (debug_mode==true) std::cout<<"Acceleration changed values: "<<ax<<","<<ay<<","<<az<<"\n";

		velx+=(ax*DT);
		vely+=(ay*DT);
		velz+=(az*DT);
		if (debug_mode==true) std::cout<<"Velocity delta changes: "<<ax*DT<<","<<ay*DT<<","<<az*DT<<"\n";
		if (debug_mode==true) std::cout<<"Velocity changed values: "<<velx<<","<<vely<<","<<velz<<"\n";

		
		x+=(velx*DT);
		y+=(vely*DT);
		z+=(velz*DT);
		if (debug_mode==true) std::cout<<"Position delta Changes: "<<velx*DT<<","<<vely*DT<<","<<velz*DT<<"\n";
		if (debug_mode==true) std::cout<<"Position Changed Values: "<<x<<","<<y<<","<<z<<"\n";
	}
	
	
	
	void detect_collision(){
		
		
		
		
		
	}
	
	
	
};

//Star Class End

//declare objects;
star cluster[2];





/*
 
 ****************INTERNAL GLOBAL FUNCTIONS FOR MAIN****************
 
 */
// Calculates the mass of the universe EXCLUDING calling star
// Parameters: Array ID of star
// Returns: double, mass of universe - mass of calling star

double usable_mass(int SelfID){
	double usablemass=0;
	for (int i = 0; i < netstars; i++ ) {
		if (i==SelfID) continue;
		usablemass+=cluster[i].mass;		
	}
	
	return usablemass;
}

// Calculate net distance btw two stars.
// Parameters: array ID of first star, array ID of second star
// Returns: distance btw their coordinates (as double)

double dist_rel(int origin,int other_star){
	
	double x,y,z;
	x=(cluster[origin].x-cluster[other_star].x);
	y=(cluster[origin].y-cluster[other_star].y);
	z=(cluster[origin].z-cluster[other_star].z);
	
	return sqrt((x*x)+(y*y)+(z*z));
}


// Input number of particles to take part in universe as an integer value
// Parameters: NULL
// Returns: NULL, sets value global variable 'netstars'

void input_netnumber(){
	
		
	//temporary fix, assigning 2
	
	netstars=2;
	
	
}



// Calculates the net mass of the universe
// Parameters: NULL
// Returns: NULL, sets value of global variable 'netmass'

void externmass(){
	
	for (int i = 0; i < netstars; i++ ) netmass+=cluster[i].mass;
	
}





// Inputs intial states and assigns them to each star of universe.
// Parameters: NULL
// Returns: NULL, sets parameters of each of 'netstars' number of star objects

void input_initial_state(){
//This is a sample symmetrical state
	cluster[0].SelfID=0;
	cluster[0].x= 7600000;
	cluster[0].y=0;
	cluster[0].z=0;
	cluster[0].velx=0;
	cluster[0].vely=30;
	cluster[0].velz=0;
	cluster[0].mass=1000;
	cluster[0].immovable=false;
	
	cluster[1].SelfID=1;
	cluster[1].x=10;
	cluster[1].y=0;
	cluster[1].z=0;
	cluster[1].velx=0;
	cluster[1].vely=0;
	cluster[1].velz=0;
	cluster[1].mass= 5.9742e24;
	cluster[1].immovable=false;
	
	
	
}






//Calculates forces in x,y,z dirction when called with the starID
// Parameters: starID, the array ID of the star for which force is to calculated.
// Returns:
/*
void calc_force(int starID){
	
	
	double force;
		
		
	
	
}

 */

// Calcualates Center of Mass of current universe state minus selected star when called WITH RESPECT TO 0,0,0
// Parameters: NULL
// Returns: NULL, sets value of global COM,cx,cy,cz - Centre of Mass coordinates of the system minus selected star (double) WRT 0,0,0

void calc_relativeCOM(int starID){
	//formula => sigma (mass.distance) / sigma (mass)
	//imp: distances need to be RELATIVE to selected star
	
	cx=cy=cz=COM=0;
	for (int i = 0; i < netstars; i++ ) {
		if (i==starID) continue;
		cx+=(cluster[i].mass*(cluster[i].x));					//component along x axis = sigma[ mass*x ];
	}
	

	cx/=usable_mass(starID);
	if (debug_mode==true)	std::cout<<"\n**New Iterarion: StarID: "<<starID<<"   usable mass:"<<usable_mass(starID)<<"  CX after division: "<<cx<<"\n";
	
	for (int i = 0; i < netstars; i++ ) {
		if (i==starID) continue;
		cy+=(cluster[i].mass*(cluster[i].y));					//component along y axis = sigma[ mass*y ];
		
	}
	cy/=usable_mass(starID);
	
	for (int i = 0; i < netstars; i++ ) {
		if (i==starID) continue;
		cz+=(cluster[i].mass*(cluster[i].z));					//component along z axis = sigma[ mass*z ];
	}
	cz/=usable_mass(starID);
	
	//std::cout<<netmass;

	if (debug_mode==true)	std::cout<<"Coordinates of universal COM: "<<cx<<","<<cy<<","<<cz<<"\n";
	for (int i = 0; i < netstars; i++ )	COM+=sqrt((cx*cx)+(cy*cy)+(cz*cz));
	
}


 
/*
 
 ****************INTERNAL GLOBAL FUNCTIONS FOR MAIN END****************
 
 */








////////////BEGIN MAIN PROGRAM/////////////////////


//To DO: find out about namespaces
using namespace std;

int main (int argc, char * const argv[]) {
	cout.precision(15);
	
	//Initializ and Open a CSV document to store output
	if (formal_mode==true) freopen("PEc.csv", "w", stdout);
	if (formal_mode==true) cout<<"ID,XAxis,YAxis,ZAxis,GlobalTime\n";
	//ofstream write_file;
	//write_file.open("PE.csv");
//	write_file<<"Hello";
	
	//Input Net number of stars
	input_netnumber();
	
	
	//Accept Parameters for each
	input_initial_state();
	
	// Calculate Mass of particle universe
	externmass();
	
	/* EXPLANATION:
	 
	 For every 'iteration' in time, the absolute centre of mass (i.e. wrt 0,0,0) is calculated. let this be (cx,cy,cz).
	 Angles for component of force can be found using cx,cy,cz.
	 Components of Centre of Mass wrt selected star are then => c'x=cx-m*x, c'y=cy-m*y, c'z=cz-m*z ; where m is mass of selected star
	 & x,y,z are it's coordinates.
	 therefore, Centre of Mass of universe wrt to selected star => sqrt( c'x^2 + c'y^2 + c'z^2);
	 
	 Force = G M m / r^2
	 G= GRAV => Gravitational Constant
	 M= Universal Mass - self mass
	 m= self mass
	 r= distance btw COM of (universe - self) coordinates and self.x,self.y,self,z coordinates
	 
	 Annexure A
	 
	 
	 */
	
	int iter_count=0;
	
	bool write_mode=1;
	//double a1=5+2.68823e-2;
	//double a2=5;
	//cout<<(a1==a2);
	
	
	//loop0: this loop controls the time depending upon dt?
	while (Global_Time<=time_limit){   // Set it to time limit ************
		//cout<<Global_Time;

		//cout << "LOOP0<>"; //debug
		
		
		if (iter_count==100){
			
		//switch Writing Mode ON for next iteration
			write_mode=true;			
			
		}
		
		
		
		//loop1: Run iterations from Star ID = 0 to Star ID = netstars-1;
		for (int i = 0; i < netstars; i++ ){
		//	cout << "loop1<>"; //debug
			//Proceed IF iter & immovable allows it
			if (cluster[i].immovable == true) {
				if (debug_mode==true) std::cout<<"**Star Immovable, Skipped**\n\nStar 1:\n";
				continue; //skip this star, no action needed.
				
				}
			
			else if (cluster[i].skipiter == 0){
				
				calc_relativeCOM(i);
				
				//Calculate Force of external Particles
				cluster[i].find_forces();
			//	cout<<"find_forces in loop1 Over\n";
				
				
				//Move through time dt (predefined)
				cluster[i].move_through_time();
			//	cout<<"move_through_time in loop1 Over\n";
				
				
				//collision detection
				cluster[i].detect_collision();
				//cout<<"detect collisions in loop1 Over\n";
				
				
				if (write_mode==true) {
					//cout<<"Write Mode Initialized";
				//	cout<<"ID "<<i<<","<<cluster[i].x<<","<<cluster[i].y<<","<<cluster[i].z<<endl;
				//	write_file.write("sup",100);
				}
				
			} //skipiter condition end
			
			// underneath else if reduces the value of skipiter till its zero, then computation resumes
			
			else if (cluster[i].skipiter>0){
				cout<<"SkipIter reduced";
				cluster[i].skipiter--;
			}
			
		}//loop1 end
		/* 
		 ***************Function of Iteration*******************
		 The iteration serves a purpose to reduce calculation. Whenever the velocity is so 
		 high such that 'dt' can be made large without significant loss in accuracy, 
		 calculation is done for calculated 'n' iterations beforehand. The calculation 
		 steps are then skipped for n iterations. 
		 
		 */
		
		
		
		
		/*	************ How iter_count works... ***********
		 
		 only after 1/DT iterations will ONE second pass. iter_count increments
		 from 0 to <1/DT. Then it is reset to 0 and Global_Time incremented by 1 second.
		 
		 
		 RECHECK for existance of OBOB error. 
		 
		 */
		
		
		
		
		if (iter_count<(1/DT)){
			iter_count++;
			//cout<<"Iteration Number: "<<iter_count<<endl;
		}
		
		else {
			
			iter_count=0;
			Global_Time+=1;
			
			//Sampling the output by second:P
			if (Global_Time%100==0){
				
				if (formal_mode==true){
					for (int i=0;i<1; i++){
						cout<<i<<","<<cluster[i].x<<","<<cluster[i].y<<","<<cluster[i].velz<<","<<Global_Time<<endl;
					}			
				}	
				
				
			}
			
			
			if (give_output==true){
			cout<<"***************\n";
				cout<<"At "<<Global_Time<<"Seconds, x,y,z of Star ID 1"<<": "<<cluster[0].x<<","<<cluster[0].y<<","<<cluster[0].z<<endl;
				cout<<"At "<<Global_Time<<"Seconds, x,y,z of Star ID 2"<<": "<<cluster[1].x<<","<<cluster[1].y<<","<<cluster[1].z<<endl;
			}
			
			
		}
			
			
			//reset write mode
		if (write_mode==true){
		//	cout<<"Write Mode Reset";
			write_mode=false;
		}
		//loop0 ends
	}
//	write_file.close();
    return 0;
}
