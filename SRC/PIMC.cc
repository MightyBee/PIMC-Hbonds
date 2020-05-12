#define _USE_MATH_DEFINES
#include <vector>
#include <array>
#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <memory>
#include <ctime>
#include <cmath>
#include "ConfigFile.tcc"
using namespace std;

double hbar(10.54571628); // IUPAC

std::mt19937 rng(time(0));


/*################################# NOTES #################################/
	- theoretically x_0, x_1, ... , x_(N_slices)       , (N_slices+1) points
	- we consider boundary conditions, x_0 = x_(N_slices)
	- hence we only consider x_0, x_1, ... , x_(N_slices-1) , N_slices points
	- to get the "full picture" simply add one more point, equal to x_0
*/

/*###### PLAN OF THE CODE #####//
	- Part A : headears
			- A.1 : function headers
			- A.2 : definitions of class "Potential_ext" and inherited classes
			- A.3 : definitions of class "Potential_ext" and inherited classes
			- A.4 : definitions of class "System"

	- Part B : main
			- B.1 : parameters acquisition and system initialization
			- B.2 : metropolis algorithm
			- B.3 : statistics writing

	- Part A : definitions
			- C.1 : definitions of the methods of class "Potential_ext" and
			        inherited classes
			- C.2 : definitions of the methods of class "Potential_ext" and
			        inherited classes
			- C.3 : definitions of the methods of class "System"
			- C.4 : function definitions


*/

//########################################################################//
//########################################################################//
//######################                            ######################//
//######################      PART A : HEADERS      ######################//
//######################                            ######################//
//########################################################################//
//########################################################################//


//########################################################################//
//######################## A.1 : FUNCTION HEADERS ########################//
//########################################################################//

// Generate a random (uniform) double between 'min' and 'max'
double randomDouble(const double& min=0.0, const double& max=1.0,
                    const bool& closed=true);

// Generate a random double from a normal Cauchy distribution
double CauchyDistribution();

// Generate a random double from one of the implemeneted distributions
double GenerateDist(const double& h);



//########################################################################//
//##################### A.2 : CLASSES Potential_ext ######################//
//########################################################################//

// Abstract class for external potential
class Potential_ext {
public:
	// pure virtual method => abstract class
	// return V at point x
	virtual double operator()(const double& x) const = 0;
	double e0_estimator(const double& x) const{return 0;};
};


// Class for a null potential
class PotExt_null: public Potential_ext {
public:
	double operator()(const double& x) const {return 0.0;}
};


// Class for a harmonic potential
class PotExt_harm: public Potential_ext {
public:
	PotExt_harm(const ConfigFile& configFile);
	double operator()(const double& x) const;
private:
	// mass and squared frequency
	double m, omega2;
};


// Class for a double well potential
class PotExt_double: public Potential_ext {
public:
	PotExt_double(const ConfigFile& configFile);
	double operator()(const double& x) const;
private:
	// barrier height and position of the wells
	double V0, x0;
};


// Class for a square potential (barrier for V0>0 and well for V0<0)
class PotExt_square: public Potential_ext {
public:
	PotExt_square(const ConfigFile& configFile);
	double operator()(const double& x) const;
private:
	// potential height, position of the centre and width of the square well
	double V0, x0, L;
};


// Class for a sinusoidal potential
class PotExt_sin: public Potential_ext {
public:
	PotExt_sin(const ConfigFile& configFile);
	double operator()(const double& x) const;
private:
	// potential height and period
	double V0, L;
};


// Class for a Lennard-Jones potential
class PotExt_LJ: public Potential_ext {
public:
	PotExt_LJ(const ConfigFile& configFile);
	double operator()(const double& x) const;
private:
	double V0, x0;
};

// Class for a H-bond potential
class PotExt_OHbonds: public Potential_ext {
public:
	PotExt_OHbonds(const ConfigFile& configFile);
	// Morse potential
	double Vmorse(const double& x) const;
	// First derivative of Morse potential
	double dVmorse(const double& x) const;
	// Estimator of zero-point energy
	double e0_estimator(const double& x) const;
	// Return the value of the potential
	double operator()(const double& x) const;
private:
	// Characteristic constants of the potential
	double D, a, r0, delta1, b, R1;
	// R : distance between the donor and accpetor
	double R, DELTA;
};



//########################################################################//
//##################### A.2 : CLASSES Potential_int ######################//
//########################################################################//

// Abstract class for internal potential
class Potential_int {
public:
	// pure virtual method => abstract class
	// return V for particle 1,2 at position x1,x2 resp.
	virtual double operator()(const double& x1, const double& x2) const = 0;
};


// Class for a null internal potential (no interactions between particles)
class PotInt_null: public Potential_int {
public:
	double operator()(const double& x1, const double& x2) const {return 0.0;}
};


// Class for a harmonic potential between two particles
class PotInt_harm: public Potential_int {
public:
	PotInt_harm(const ConfigFile& configFile);
	double operator()(const double& x1, const double& x2) const;
private:
	// stifness and rest length
	double k, l0;
};


// Class for a Lennard-Jones potential between two particles
class PotInt_LJ: public Potential_int {
public:
	PotInt_LJ(const ConfigFile& configFile);
	// standard Lennard-Jones potential
	double LJ(const double& r) const;
	// Lennard-Jones with parameters defined below
	double operator()(const double& x1, const double& x2) const;
private:
	double V0, x0, G;
};



//########################################################################//
//######################### A.3 : CLASSES System #########################//
//########################################################################//


// Class System : - contains all the physical properties of the system
//                  considered, as well as the number of slices
//                - it can generate moves but it's not this class that
//                  will porpose it
class System {
public:
	// constructor
	System(const ConfigFile& configFile);
	// initlaize the system : hot start
	void initialize(const double& pos_min, const double& pos_max);
	// write in an output file the external potential used
	void write_potExt(const string& output);
	// return the number of particles in the system
	size_t nb_part() const {return N_part;}
	// return the number of slices in the system
	size_t nb_slices() const {return N_slices;}
	// return the number of time each site is visited by MH algorithm
	vector<vector<int>> get_visits() const {return verif;}
	// write the poisitions of all particles
	ostream& write(ostream& output) const;
	// return the kinetic term of a particle between a bead and a neighbour
	double kinetic(const int& particle, const int& bead, const int& bead_pm,
                  const double& displacement=0.0) const;
	// compute the whole euclidean action directly
	double energy();
	// return the euclidean action measured along the simulation's progress
	double get_H(){return H;}

	// returns if a move is accepted or not
	bool metropolisAcceptance();

	// different moves possible
	bool localMove(const double& h);
	bool globalDisplacement(const double& h);
	bool bisection(const double& h, const double& sRel);
	bool swap();
	bool inverse();
	bool symmetryCM();

	// energy measures
	void measure_energy(double, double);
	void average_energy();


private:
	////////  physical parameters  ////////
	unsigned int N_part;
	unsigned int N_slices;
	double tau;
	double d_tau;
	int q;
	vector<double> mass;
	double omega;
	unique_ptr<Potential_ext> ptr_Vext;
	unique_ptr<Potential_int> ptr_Vint;
	// table of positions : first index indicates the particle
	//                      second index indicates the slice
	vector<vector<double>> table;

	////////  utilitary variables  ////////
	// mm : time slice randomly selected during each iteration,
	// mm_plu=mm+1, mm_min=mm-1 with boundary conditions
	unsigned int mm, mm_plu, mm_min;
	// particle randomly selected during each iteration
	unsigned int nn;
	// displacement proposed
	double dis;
	// part of the action that is changed by moves
	double s_old, s_new;
	// euclidean action
	double H;

	vector<double> energies_psi;
	vector<double> energies_h;
	vector<vector<int>> verif;
};

ostream& operator<<(ostream& output, const System& s);








//########################################################################//
//########################################################################//
//######################                            ######################//
//######################       PART B : MAIN        ######################//
//######################                            ######################//
//########################################################################//
//########################################################################//


int main(int argc, char* argv[]){


	//####################################################################//
	//################### B.1 : PARAMETERS ACQUISITION ###################//
	//####################################################################//

	//######################### VILLARD LIBRARY #########################//

	// Default input configuration file
	string inputPath("configuration.in");
	if(argc>1) // specified input file specified by user
		inputPath = argv[1];

	// Parameters are read et stocked in a "map" of strings
	ConfigFile configFile(inputPath);

	for(int i(2); i<argc; ++i) // complementary inputs
		configFile.process(argv[i]);


	//######################### READ PARAMETERS #########################//

	// number of Monte Carlo iterations (aka sweeps)
	unsigned int N_sweeps(configFile.get<unsigned int>("N_sweeps"));
	// number of thermalisation sweeps
	unsigned int N_thermalisation(configFile.get<unsigned int>("N_thermal"));
	// initial minimum position
	double pos_min(configFile.get<double>("pos_min"));
	// initial maximal position
	double pos_max(configFile.get<double>("pos_max"));
	// displacement parameter of a point in the path
	vector<double> h(3,configFile.get<double>("h"));
	// proportion of each move
	double p_loc(configFile.get<double>("p_local"));
	double p_dsp(configFile.get<double>("p_displacement"));
	double p_swap(configFile.get<double>("p_swap"));
	double p_inv(configFile.get<double>("p_inverse"));
	double p_sym(configFile.get<double>("p_symmetryCM"));
	double p_bis(configFile.get<double>("p_bisection"));
	// relative size (0 to 1) of bisection move
	double s_bis(configFile.get<double>("s_bisection"));
	// numbers of tries for each moves
	vector<unsigned int> NbTries(6,0);
	// acceptance rates for the three moves
	vector<double>			accrate(6,0.0);
	// "instantaneous" acceptance rate for local moves
	double tmp_accrate(0.0);
	// ideal acceptance rate for local moves
	double idrate(configFile.get<double>("idrate"));
	// output is written every n_stride iterations
	size_t n_stride(configFile.get<size_t>("n_stride"));


	//Output files
	string output(configFile.get<string>("output"));
	string output_pos(output+"_pos.out");
	ofstream fichier_output(output_pos.c_str());
	fichier_output.precision(15);		// Precision

	string output_energy(output+"_nrg.out");
	ofstream fichier_energy(output_energy.c_str());
	fichier_energy.precision(15);		// Precision

	string output_rate(output+"_rate.out");
	ofstream fichier_rate(output_rate.c_str());
	fichier_rate.precision(15);

	// initialization of the system
	System s(configFile);
	s.initialize(pos_min,pos_max);
	s.write_potExt(output);
	fichier_output << s << endl;

	//UNCOMMENT IF YOU WANT TO COMPARE WITH OUTPUT2_NRG.OUT
	//s.measure_energy();
	//double V0(configFile.get<double>("R"));
	//double V0(configFile.get<double>("V0"));
	//double x0(configFile.get<double>("x0"));

	double last_measured_time(time(0));


	//####################################################################//
	//#################### B.2 : METROPOLIS ALGORITHM ####################//
	//####################################################################//

	//For every sweep...
	for(size_t i(0); i < N_sweeps; i++){

		// show progress of the simulation
		if(time(0) - last_measured_time >= 5){
			last_measured_time = time(0);
			cout << floor((double)i/N_sweeps*100) << " %" << endl;
		}

		//For every particle...
		/*(we try one average each move for each particle one time every sweep
		  if all moves proportions are set to 1) */
		for(size_t j(0); j < s.nb_part(); j++){

			//////////  local move  //////////
			/*if proportion of tries compared to the number of sweeps is
			  smaller than the target, try a move*/
			if(NbTries[0]*1.0/((i*s.nb_part()+j+1)*s.nb_slices()) < p_loc){
				tmp_accrate=0.0;
				for(size_t k(0); k < s.nb_slices(); k++){
					NbTries[0]++;
					if(s.localMove(h[0])){
						accrate[0]++;
						tmp_accrate++;
					}
				}
				tmp_accrate/=s.nb_slices();
				// adjust the parametr h to reach the ideal accpetance rate
				h[0]*=tmp_accrate/idrate;
			}

			//////////  global displacement  //////////
			/*if proportion of tries compared to the number of sweeps is
			  smaller than the target, try a move*/
			if(NbTries[1]*1.0/(i*s.nb_part()+j+1) < p_dsp){
				NbTries[1]++;
				if(s.globalDisplacement(h[1])){
					accrate[1]++;
				}
			}

			//////////  bisection //////////
			if(NbTries[2]*1.0/(i*s.nb_part()+j+1) < p_bis){
				NbTries[2]++;
				if(s.bisection(h[2], s_bis)){
					accrate[2]++;
				}
			}

			//////////  swap  //////////
			if(NbTries[3]*1.0/(i*s.nb_part()+j+1) < p_swap){
				NbTries[3]++;
				if(s.swap()){
					accrate[3]++;
				}
			}

			//////////  inverse  //////////
			if(NbTries[4]*1.0/(i*s.nb_part()+j+1) < p_inv){
				NbTries[4]++;
				if(s.inverse()){
					accrate[4]++;
				}
			}

			//////////  symmetryCM  //////////
			if(NbTries[5]*1.0/(i*s.nb_part()+j+1) < p_sym){
				NbTries[5]++;
				if(s.symmetryCM()){
					accrate[5]++;
				}
			}
		}


		//######################## OUTPUT IN FILE ##########################//
		if((i%n_stride) == 0){
			fichier_output << s << endl;
			fichier_energy << s.energy() << " " << s.get_H() << endl;
			fichier_rate << tmp_accrate << endl;

			//Energy measurement
			if(i >= N_thermalisation){
				//fichier_energy << s.energy() << " " << s.get_H() << endl;
				//s.measure_energy(V0, x0);
			}
		}
	}
	fichier_output.close();
	fichier_energy.close();
	fichier_rate.close();


	//####################################################################//
	//##################### B.3 : STATISTICS WRITING #####################//
	//####################################################################//

	// statistics on the moves
	// number of tries and acceptance's rates
	string output_stat(output+"_stat.out");
	fichier_output.open(output_stat.c_str());
	fichier_output.precision(15);
	for(size_t i(0); i<accrate.size(); i++){
		accrate[i]/=NbTries[i];
		fichier_output << NbTries[i] << " " << accrate[i] << endl;
	}

	// statistics on the visits of patricles and slices
	for(const auto& part : s.get_visits()){
		for(const auto& v : part){
			fichier_output << v/(NbTries[0]*1.0/(s.nb_part()*s.nb_slices()))
			               << " ";
		}
		fichier_output << endl;
	}

	fichier_output.close();

	//Energy
	s.average_energy();

	//############################ END OF MAIN ############################//
	return 0;
}









//########################################################################//
//########################################################################//
//######################                            ######################//
//######################    PART C : DEFINITIONS    ######################//
//######################                            ######################//
//########################################################################//
//########################################################################//



//########################################################################//
//########### C.1 : CLASS 'Potential_ext' METHODS DEFINITIONS ############//
//########################################################################//


//############################  PotExt_harm  #############################//
///// constructor /////
PotExt_harm::PotExt_harm(const ConfigFile& configFile) :
	Potential_ext(),
	m(configFile.get<double>("mass")),
	omega2(pow(configFile.get<double>("omega"),2))
	{}

///// potential opertor /////
double PotExt_harm::operator()(const double& x) const {
	return 0.5 * m * pow(x, 2) * omega2;
}


//###########################  PotExt_double  ############################//
///// constructor /////
PotExt_double::PotExt_double(const ConfigFile& configFile) :
	Potential_ext(),
	V0(configFile.get<double>("V0")),
	x0(configFile.get<double>("x0"))
	{}

///// potential opertor /////
double PotExt_double::operator()(const double& x) const {
	return V0*pow(pow(x/x0,2)-1,2);
}


//###########################  PotExt_square  ############################//
///// constructor /////
PotExt_square::PotExt_square(const ConfigFile& configFile) :
	Potential_ext(),
	V0(configFile.get<double>("V0")),
	x0(configFile.get<double>("x0")),
	L(configFile.get<double>("L"))
	{}

///// potential opertor /////
double PotExt_square::operator()(const double& x) const {
	if(abs(x - x0) < L/2){
		return V0;
	}else{
		return 0;
	}
}


//#############################  PotExt_sin  #############################//
///// constructor /////
PotExt_sin::PotExt_sin(const ConfigFile& configFile) :
	Potential_ext(),
	V0(configFile.get<double>("V0")),
	L(configFile.get<double>("L"))
	{}

///// potential opertor /////
double PotExt_sin::operator()(const double& x) const {
	return 0.5*V0*(1-cos(2*M_PI/L*x));
}


//#############################  PotExt_LJ  ##############################//
///// constructor /////
PotExt_LJ::PotExt_LJ(const ConfigFile& configFile) :
	Potential_ext(),
	V0(configFile.get<double>("V0")),
	x0(configFile.get<double>("x0"))
	{}

///// potential opertor /////
double PotExt_LJ::operator()(const double& x) const {
	return 4 * V0 * (pow(x/x0,-12) - pow(x/x0,-6) );
}

//###########################  PotExt_OHbonds  ###########################//
///// constructor /////
PotExt_OHbonds::PotExt_OHbonds(const ConfigFile& configFile) :
	Potential_ext(),
	D(83.402), a(2.2), r0(0.96), delta1(0.4*D), b(2.2), R1(2*r0+1/a),
	R(configFile.get<double>("R")),
	DELTA(delta1*exp(-b*(R-R1)))
	{}

///// Morse potential /////
double PotExt_OHbonds::Vmorse(const double& x) const{
	return D*(exp(-2*a*(x-r0))-2*exp(-a*(x-r0)));
}

///// first derivative of Morse potential /////
double PotExt_OHbonds::dVmorse(const double& x) const{
	return 2*a*D*(exp(-a*(x-r0))-exp(-2*a*(x-r0)));
}

///// potential opertor /////
double PotExt_OHbonds::operator()(const double& x) const{
	if(abs(x) < 8){
		return 0.5*(Vmorse(R/2+x)+Vmorse(R/2-x) - sqrt(pow(Vmorse(R/2+x)-Vmorse(R/2-x),2)+4*DELTA*DELTA));
	}else{
		return 0.0;
	}
}

///// estimataor of zero-point energy /////
double PotExt_OHbonds::e0_estimator(const double& x) const{
	return x/4 * (dVmorse(R/2+x) - dVmorse(R/2-x) - (Vmorse(R/2+x) - Vmorse(R/2-x))*(dVmorse(R/2+x) + dVmorse(R/2-x))/sqrt(pow(Vmorse(R/2+x) - Vmorse(R/2-x), 2) + 4*DELTA*DELTA));
}




//########################################################################//
//########### C.2 : CLASS 'Potential_int' METHODS DEFINITIONS ############//
//########################################################################//

//############################  PotInt_harm  #############################//
///// constructor /////
PotInt_harm::PotInt_harm(const ConfigFile& configFile) :
	Potential_int(),
	k(configFile.get<double>("k")),
	l0(configFile.get<double>("l0"))
	{}

///// potential opertor /////
double PotInt_harm::operator()(const double& x1, const double& x2) const {
	return 0.5*k*pow(abs(x2-x1)-l0,2);
}


//#############################  PotInt_LJ  ##############################//
///// constructor /////
PotInt_LJ::PotInt_LJ(const ConfigFile& configFile) :
	Potential_int(),
	V0(configFile.get<double>("Vmin")),
	x0(configFile.get<double>("x0")),
	G(configFile.get<double>("G"))
	{}

///// standard Lennard-Jones potential /////
double PotInt_LJ::LJ(const double& r) const {
	if(r>0.35){
		return pow(r,-12)-2*pow(r,-6);
	}else{
		return pow(0.35,-12)-2*pow(0.35,-6);
	}
}

///// potential opertor /////
double PotInt_LJ::operator()(const double& x1, const double& x2) const {
	return V0*LJ(abs(x1-x2)/x0);
}



//########################################################################//
//############### C.3 : CLASS 'System' METHODS DEFINITION ################//
//########################################################################//


//############################  constructor  #############################//
System::System(const ConfigFile& configFile) :
	N_part(configFile.get<unsigned int>("N_part")),
	N_slices(configFile.get<unsigned int>("N_slices")),
	tau(configFile.get<double>("tau")),
	d_tau(tau/N_slices),
	mass(N_part,configFile.get<double>("mass")),
	omega(configFile.get<double>("omega")),
	table(N_part, vector<double>(N_slices, 0.0)),
	mm(0), mm_plu(0), mm_min(0), nn(0),
	dis(0.0), s_old(0.0), s_new(0.0), H(0.0),
	verif(N_part, vector<int>(N_slices, 0))
	{
		for(unsigned int i(0); i<N_part; i++){
			mass[i]=configFile.get<double>("m"+to_string(i+1));
		}
		// choosing the external potential
		string V_ext(configFile.get<string>("V_ext"));
		if(V_ext=="null") ptr_Vext = move(unique_ptr<Potential_ext>(new PotExt_null()));
		else if(V_ext=="harmonic") ptr_Vext = move(unique_ptr<Potential_ext>(new PotExt_harm(configFile)));
		else if(V_ext=="double") ptr_Vext = move(unique_ptr<Potential_ext>(new PotExt_double(configFile)));
		else if(V_ext=="square") ptr_Vext = move(unique_ptr<Potential_ext>(new PotExt_square(configFile)));
		else if(V_ext=="sin") ptr_Vext = move(unique_ptr<Potential_ext>(new PotExt_sin(configFile)));
		else if(V_ext=="LJ") ptr_Vext = move(unique_ptr<Potential_ext>(new PotExt_LJ(configFile)));
		else if(V_ext=="OHbonds") ptr_Vext = move(unique_ptr<Potential_ext>(new PotExt_OHbonds(configFile)));
		else{
			cerr << "Please choose a valid potential." << endl;
		}
		// choosing the internal potential
		string V_int(configFile.get<string>("V_int"));
		if(V_int=="null") ptr_Vint = move(unique_ptr<Potential_int>(new PotInt_null()));
		else if(V_int=="harmonic") ptr_Vint = move(unique_ptr<Potential_int>(new PotInt_harm(configFile)));
		else if(V_int=="LJ") ptr_Vint = move(unique_ptr<Potential_int>(new PotInt_LJ(configFile)));
		else{
			cerr << "Please choose a valid potential." << endl;
		}
	}


// initialize random paths for each particles
void System::initialize(const double& pos_min, const double& pos_max){
	for(auto& particle : table){
		for(auto& pos : particle){
			pos = randomDouble(pos_min, pos_max);
		}
	}
	H=energy();
}


void System::write_potExt(const string& output){
	string output_pot(output+"_pot.out");
	ofstream f_pot(output_pot.c_str());
	f_pot.precision(15);
	size_t N(10000);
	double x(0.0);
	for(size_t i(0); i<N; i++){
		double xi(-50.0), xf(50.0);
		x=xi+i*(xf-xi)/(N-1);
		f_pot << x << " " << (*ptr_Vext)(x) << endl;
	}
	f_pot.close();
}


// write the paths of all particles in one line
ostream& System::write(ostream& output) const{
	for(const auto& particle : table){
		for(const auto& pos : particle){
			output << pos << " ";
		}
	}
	return output;
}


double System::kinetic(const int& particle, const int& bead, const int& bead_pm, const double& displacement) const{
	return 0.5*mass[particle]*pow(((table[particle][bead]+displacement)-table[particle][bead_pm])/d_tau,2);
}

double System::energy(){
	double E(0.0);
	for(size_t part(0); part<N_part; part++){
		for(size_t bead(0); bead<N_slices; bead++){
			E+=kinetic(part,bead,(bead+1)%N_slices);
			E+=(*ptr_Vext)(table[part][bead]);
			for(size_t part2(part+1); part2<N_part; part2++){
				E+=(*ptr_Vint)(table[part][bead],table[part2][bead]);
			}
		}
	}
	return E;
}



bool System::metropolisAcceptance(){
	return ( randomDouble(0,1) <= exp(-(0.1*d_tau/hbar) * (s_new - s_old)) );
}



//###############################  moves  ###############################//


bool System::localMove(const double& h){
	// random integer between 0 and N_slices-1
	mm = rng()%N_slices;
	// mm-1 with periodic boundary condition
	mm_min = (mm + N_slices - 1)%N_slices;
	// mm+1 with periodic boundary condition
	mm_plu = (mm + 1)%N_slices;
	// random integer between 0 and N_part-1
	nn = rng()%N_part;

	verif[nn][mm]++;

	dis=GenerateDist(h);

	/* as we take the difference of new and old action S_new-S_old, we can //
	// consider only the part of the action that is affected by the        //
	// proposed new position                                               */
	s_old = kinetic(nn,mm,mm_plu) + kinetic(nn,mm,mm_min)
			+ (*ptr_Vext)(table[nn][mm]);
	s_new = kinetic(nn,mm,mm_plu,dis) + kinetic(nn,mm,mm_min,dis)
			+ (*ptr_Vext)(table[nn][mm]+dis);

	if(N_part>1){
		for(size_t i(0); i<N_part; i++){
			if(i!=nn){
				s_old+=(*ptr_Vint)(table[i][mm],table[nn][mm]);
				s_new+=(*ptr_Vint)(table[i][mm],table[nn][mm]+dis);
			}
		}
	}

	if(metropolisAcceptance()){ // metropolis acceptance
		table[nn][mm] += dis;    // update position with new one
		H += (s_new - s_old);    // update total action
		return true;
	}else{
		return false;
	}
}




bool System::globalDisplacement(const double& h){
	// random integer between 0 and N_part-1
	nn = rng()%N_part;

	dis=GenerateDist(h);

	/* no relative move between the time slices //
	// --> only the potential action changes    */
	s_old=0.0;
	s_new=0.0;
	for(size_t j(0); j<N_slices; j++){
		s_old+=(*ptr_Vext)(table[nn][j]);
		s_new+=(*ptr_Vext)(table[nn][j]+dis);
		if(N_part>1){
			for(size_t i(0); i<N_part; i++){
				if(i!=nn){
					s_old+=(*ptr_Vint)(table[i][j],table[nn][j]);
					s_new+=(*ptr_Vint)(table[i][j],table[nn][j]+dis);
				}
			}
		}
	}

	if(metropolisAcceptance()){ // metropolis acceptance
		for(auto& pos : table[nn]){
			pos+=dis;
		}
		H += (s_new - s_old);
		return true;
	}else{
		return false;
	}
}




bool System::bisection(const double& h, const double& sRel){
	// random integer between 0 and N_slices-1
	mm = rng()%N_slices;
	// mm-1 with periodic boundary condition
	mm_min = (mm + N_slices - 1)%N_slices;
	// random integer between 0 and N_part-1
	nn = rng()%N_part;

	dis=GenerateDist(h);
	size_t l(N_slices*sRel);

	s_old=0.0;
	s_new=0.0;
	int ind_j(0);
	for(size_t j(0); j<l; j++){
		ind_j=(mm+j)%N_slices;
		s_old+=(*ptr_Vext)(table[nn][ind_j]);
		s_new+=(*ptr_Vext)(table[nn][ind_j]+dis);
		if(N_part>1){
			for(size_t i(0); i<N_part; i++){
				if(i!=nn){
					s_old+=(*ptr_Vint)(table[i][ind_j],table[nn][ind_j]);
					s_new+=(*ptr_Vint)(table[i][ind_j],table[nn][ind_j]+dis);
				}
			}
		}
	}
	s_old += kinetic(nn,mm,mm_min)
	       + kinetic(nn,(mm+l-1)%N_slices,(mm+l)%N_slices);
	s_new += kinetic(nn,mm,mm_min,dis)
	       + kinetic(nn,(mm+l-1)%N_slices,(mm+l)%N_slices,dis);

	if(metropolisAcceptance()){ // metropolis acceptance
		for(size_t i(0); i<l; i++){
			table[nn][(mm+i)%N_slices]+=dis;
		}
		H += (s_new - s_old);
		return true;
	}else{
		return false;
	}
}




bool System::swap(){
	if (N_part>1){
		// random integer between 0 and N_part-1
		mm_min = rng()%N_part;
		// another, but different, random integer between 0 and N_part-1
		mm_plu = (mm_min+1+rng()%(N_part-1))%N_part;
		// random integer between 0 and N_slices-1 (bead where the swap starts)
		mm = rng()%N_slices;
		//length of the swap (nb of slices swapped)
		nn = rng()%(N_slices-1)+1;

		s_old=0.0;
		s_new=0.0;
		int ind_j(0), ind_j_pm(0);

		ind_j=mm;
		ind_j_pm=(mm+N_slices-1)%N_slices;
		s_old += kinetic(mm_min,ind_j,ind_j_pm)
		       + kinetic(mm_plu,ind_j,ind_j_pm);
		s_new += kinetic(mm_min,ind_j,ind_j_pm,table[mm_plu][mm]-table[mm_min][mm])
		       + kinetic(mm_plu,ind_j,ind_j_pm,table[mm_min][mm]-table[mm_plu][mm]);

		for(size_t j(0); j<nn; j++){
			ind_j=(mm+j)%N_slices;
			ind_j_pm=(ind_j+N_slices-1)%N_slices;
			if(mass[mm_min]!=mass[mm_plu]){
				for(size_t i(0); i<N_part; i++){
					if(i!=mm_min and i!=mm_plu){
						// change for particle(mm_min)
						s_old+=(*ptr_Vint)(table[i][ind_j],table[mm_min][ind_j]);
						s_new+=(*ptr_Vint)(table[i][ind_j],table[mm_plu][ind_j]);
						// change for particle(mm_plu)
						s_old+=(*ptr_Vint)(table[i][ind_j],table[mm_plu][ind_j]);
						s_new+=(*ptr_Vint)(table[i][ind_j],table[mm_min][ind_j]);
					}
				}
			}
			/* in swapped part of paths :                      //
			// K1_new = m1*(K2_old/m2), K2_new=m2/m1*K1_old    */
			if(j){
				s_old += kinetic(mm_min,ind_j,ind_j_pm)
				       + kinetic(mm_plu,ind_j,ind_j_pm);
				s_new += mass[mm_min]/mass[mm_plu]*kinetic(mm_plu,ind_j,ind_j_pm)
				       + mass[mm_plu]/mass[mm_min]*kinetic(mm_min,ind_j,ind_j_pm);
			}
		}
		ind_j=(mm+nn-1)%N_slices;
		ind_j_pm=(mm+nn)%N_slices;
		s_old += kinetic(mm_min,ind_j,ind_j_pm)
		       + kinetic(mm_plu,ind_j,ind_j_pm);
		s_new += kinetic(mm_min,ind_j,ind_j_pm,table[mm_plu][ind_j]-table[mm_min][ind_j])
		       + kinetic(mm_plu,ind_j,ind_j_pm,table[mm_min][ind_j]-table[mm_plu][ind_j]);


		if(metropolisAcceptance()){ // metropolis acceptance
			double tmp(0.0);
			for(size_t j(0); j<nn; j++){
				ind_j=(mm+j)%N_slices;
				tmp=table[mm_min][ind_j];
				table[mm_min][ind_j]=table[mm_plu][ind_j];
				table[mm_plu][ind_j]=tmp;
			}
			H += (s_new - s_old);
			return true;
		}
	}
	return false;
}




bool System::inverse(){
	nn = rng()%N_part; // random integer between 0 and N_part-1

	// no relative move between the time slices
	// --> only the potential action changes
	s_old=0.0;
	s_new=0.0;
	for(size_t j(0); j<N_slices; j++){
		s_old+=(*ptr_Vext)( table[nn][j]);
		s_new+=(*ptr_Vext)(-table[nn][j]);
		if(N_part>1){
			for(size_t i(0); i<N_part; i++){
				if(i!=nn){
					s_old+=(*ptr_Vint)(table[i][j], table[nn][j]);
					s_new+=(*ptr_Vint)(table[i][j],-table[nn][j]);
				}
			}
		}
	}

	if(metropolisAcceptance()){ // metropolis acceptance
		for(auto& pos : table[nn]){
			pos*=-1;
		}
		H += (s_new - s_old);
		return true;
	}else{
		return false;
	}
}




bool System::symmetryCM(){
	nn = rng()%N_part; // random integer between 0 and N_part-1
	dis= 0;
	for(const auto& pos : table[nn]){
		dis+=pos;
	}
	dis*=-2.0/table[nn].size();

	// no relative move between the time slices
	// --> only the potential action changes
	s_old=0.0;
	s_new=0.0;
	for(size_t j(0); j<N_slices; j++){
		s_old+=(*ptr_Vext)(table[nn][j]);
		s_new+=(*ptr_Vext)(table[nn][j]+dis);
		if(N_part>1){
			for(size_t i(0); i<N_part; i++){
				if(i!=nn){
					s_old+=(*ptr_Vint)(table[i][j],table[nn][j]);
					s_new+=(*ptr_Vint)(table[i][j],table[nn][j]+dis);
				}
			}
		}
	}

	if(metropolisAcceptance()){ // metropolis acceptance
		for(auto& pos : table[nn]){
			pos+=dis;
		}
		H += (s_new - s_old);
		return true;
	}else{
		return false;
	}
}




void System::measure_energy(double V0, double x0){
	double temp_energy_H(0), temp_energy_ETH(0);

	double R(V0);	//For H-bond, V0 is R
	double D(83.402), a(2.2), r0(0.96), delta1(0.4*D);
	double b(2.2), R1(2*r0+1/a), DELTA(delta1*exp(-b*(R-R1)));

	temp_energy_ETH += (*ptr_Vext)(table[0][0])
	                 + (*ptr_Vext).e0_estimator(table[0][0]);

	for(size_t i(1); i < table[0].size(); i++){
		temp_energy_ETH += (*ptr_Vext)(table[0][i])
		                 + (*ptr_Vext).e0_estimator(table[0][i]);
		temp_energy_H += mass[0]/2 * pow((table[0][i] - table[0][i-1])/d_tau, 2)
		               + (*ptr_Vext)(table[0][i]);
	}
	temp_energy_H += mass[0]/2 * pow((table[0][0] - table[0][N_slices-1])/d_tau, 2)
	               + (*ptr_Vext)(table[0][0]);

	energies_psi.push_back(temp_energy_ETH/N_slices);
	energies_h.push_back(temp_energy_H/N_slices);
}




void System::average_energy(){
	ofstream fichier_output;
	fichier_output.open("e0.out");
	fichier_output.precision(15);

	double temp_energy(0), temp_error(0);

	cout << "Finally, with d_tau = " << d_tau << endl;

	//PSI energies
	for(size_t i(0); i < energies_psi.size(); i++){
		temp_energy += energies_psi[i];
		temp_error += pow(energies_psi[i], 2);
	}
	temp_energy = temp_energy/energies_psi.size();
	temp_error = sqrt((temp_error/energies_psi.size() - pow(temp_energy, 2))
	                  /energies_psi.size());

	cout << "PSI: " << temp_energy << " +- " << temp_error << endl;
	fichier_output << "PSI: " << temp_energy << " +- " << temp_error << endl;

	//H energies
	temp_energy = 0;
	temp_error = 0;

	for(size_t i(0); i < energies_h.size(); i++){
		temp_energy += energies_h[i];
		temp_error += pow(energies_h[i], 2);
	}
	temp_energy = temp_energy/energies_h.size();
	temp_error = sqrt((temp_error/energies_h.size() - pow(temp_energy, 2))
	                  /energies_h.size());

	cout << "H:   " << temp_energy << " +- " << temp_error << endl;
	fichier_output << "H:   " << temp_energy << " +- " << temp_error << endl;

	fichier_output.close();
}




ostream& operator<<(ostream& output, const System& s){
	return s.write(output);
}



//########################################################################//
//##################### C.4 : FUNCTION DEFINITIONS #######################//
//########################################################################//

// Generate a random (uniform) double between 'min' and 'max'
double randomDouble(const double& min, const double& max,
                    const bool& closed){
	if(closed) return (min + (max-min) * (double)rng()/rng.max());
	else return (min + (max-min) * ((double)rng()+0.5)/(rng.max()+1.0));
}

// Generate a random double from a normal Cauchy distribution
double CauchyDistribution(){
	return tan(M_PI*(randomDouble(-0.5,0.5,false)));
}

// Generate a random double from one of the implemeneted distributions
double GenerateDist(const double& h){
	if(rng()%2){
		return h * randomDouble(-1.0,1.0); // proposed displacement
	}else{
		return h * CauchyDistribution(); // proposed displacement
	}
}
