/*
mainTSM_Code.cpp: main script to run the Tumour-Stroma model as described in the BMC Medical Informatics paper: 
"Multiscalar cellular automaton simulates in-vivo tumour-stroma patterns calibrated from
 in-vitro assay data" in 2016 - DOI: 10.1186/s12911-017-0461-1

There is unrestricted license to use this script and modify it as long as the above publication is correctly cited.

Created by: Juan A Delgado
juan.x.delgado@gsk.com
PhD Student at AstraZeneca, United Kingdom
2015
 
*/

#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>
#include <ctime>
#include <vector>
#include <string>
#include <typeinfo>
#include <time.h>
using namespace std;
 double round(double varin)
 {double cvar = ceil(varin);
	 if((cvar-varin)>.5)
	    cvar--;
  return cvar;}

vector<int> isCell(vector<char>C1, char C2)
{  vector<int> Out(0);
	int k = 0;
	for(unsigned int it = 0;it < C1.size();it++)
		if(C1[it] == C2) {Out.push_back(it);k++;};
	return Out;
}
void EuSpace(double n,double ind, double &x, double &y)
{y = floor(ind/n);
x = ind-n*y;}
vector<double> myProfile(const double Diff,double Po2max,int n, vector<int>Pos, vector<char>Cell)
{// Identify the cells
	vector<int> Tum = isCell(Cell,'T'), Hyp = isCell(Cell,'H'), Str = isCell(Cell,'S');
	double xt = 0,yt = 0,xs = 0,ys = 0,aa;
	vector <double> Prof(0);
	double f=1,ux=0,rn;
	// Measure Euclidean distances
for (unsigned int j = 0;j<Tum.size();j++){
	EuSpace(n, Pos[Tum[j]],xt,yt);
	for (unsigned int k = 0;k<Str.size();k++){
			EuSpace(n, Pos[Str[k]],xs,ys);
		// Calculate the inverse distance
			ux += 1/sqrt(pow(xt-xs,2)+pow(yt-ys,2));	
		}// end of k/stromal cells

	//ux = 1/pow(ux,1.0/Str.size());
			Prof.push_back(Po2max*exp(-20.0/3.0*Diff/ux));
			ux=0;f=1;
} // end of j epithelial cells
	return Prof;}

vector<int> createKernel(int n,int PosNow)
{	// primary kernel
	vector<int>kr2(0);
	signed int var[8] = {-1,1,-n,-n-1,-n+1,n,n+1,n-1};
	// find the kernel applies to Pos	
	for(int k = 0;k<8;k++)
		kr2.push_back(var[k]+PosNow);
    return kr2;
}
bool ArrScComp(vector<int> myVect, int myScalar)
// compare a vector and a scalar or list of integers
{	bool outVect = 0;
	for (unsigned int a = 0;a<myVect.size();a++)
			if(myVect[a] == myScalar)
			{outVect = 1;break;}
return outVect;
}
void Probabilities(const double par[10],vector<char>*Cell,vector<int>*Pos,vector<double>*Age,vector<double>Po2,int ns,double Dt,int counts[5],ofstream& fout)
 /*Probabilities function
		INPUTS:
			double par - list of parameters for the rule-based algorithm
			vector<char>Cell - vector indicating the type of cells in the system at that given time
			vector<int>Pos - vector of positions of each of the elements
			vector<double>Age - age each of the cells
			vector<double>Po2 - profiles of oxygen in all the grid
			double Dt - this is the time step
			int counts - this is the number of cell species at any given time - Array[4] - Tumour, Hypoxic, Necrotic, Stroma, Total
	OUTPUTS:
			None
*/
{
	// Defintion of the parameters
	signed int ku = -1;
	vector <int> kr(0);
	vector <int> krmod(0);
	bool myVect[8];
	double randN;
	int randNi;
	vector<int>FreeSpaces(0);
	int celldim = Cell->size();
	//ofstream fout("Profiles.csv");
	//ofstream fout("FailReport.txt");

	// The algorithm starts
	for(unsigned int ii = 0;ii<celldim;ii++){	
		// look at the surroundings
		// is there any space free?		    
			kr = createKernel(ns,Pos->at(ii));
			
			// Compare the kernel to the Pos
			for(unsigned int ll = 0;ll<kr.size();ll++)
			{myVect[ll] = ArrScComp(*Pos,kr[ll]);
			if(!myVect[ll])
				krmod.push_back(kr[ll]);
			}

		// Generate a random number
            randN = rand()/(double)(RAND_MAX + 1);
			randNi = round(rand()/(double)(RAND_MAX + 1)*(krmod.size()-1));
			//fout<<"\n Seq: "<<ii<<" Cell: " << Cell->at(ii)<<"Tum "<<counts[0]<<"Hyp "<<counts[1]<<"Nec "<<counts[2]<<"Str "<<counts[3]<<"Tot "<<counts[4];
		// Cell decitions
		switch(Cell->at(ii)){		
		   case 'T':
			{  ku++;fout<<","<<Po2[ku];
				if (Po2[ku]<par[6] && par[2]<Age->at(ii))
					{Cell->at(ii) = 'H';
					Age->at(ii) = 0;
					counts[0]--;counts[1]++;
					} // become hypoxic
				else if(par[1]<Age->at(ii) && krmod.size() != 0 && (1-exp(-(Po2[ku]/.46)/par[0]*(Age->at(ii)-par[1])))>randN)
					{Cell->push_back('T');
					Pos->push_back(krmod[randNi]);
					Age->push_back(0);
					Age->at(ii) = 0;
					counts[0]++;counts[4]++;
					krmod.erase(krmod.begin()+randNi);
					}
				else
					Age->at(ii) += Dt; 
				break;}
			case 'S':
				{   fout<<",";
					if(par[5]*Dt*krmod.size()>randN)
					{Pos->at(ii) = krmod[randNi]; // Move
				     krmod.erase(krmod.begin()+randNi);
				}	
				Age->at(ii) += Dt;
					break;} 
			case 'H':
				{fout<<",";
				if (Age->at(ii)>par[3])
				    {Cell->at(ii) = 'N';// Become necrotic
					 Age->at(ii) = 0;// Restart the clock
				     counts[2]++;counts[1]--;}
				else
					Age->at(ii) += Dt;
				break;}
			case 'N':
				{fout<<",";
				Age->at(ii) += Dt;}
		}; // end of switch
	// Update the free space
		for(int kk = 0;kk<krmod.size();kk++)
			FreeSpaces.push_back(krmod[kk]);
	// Re-establish krmod
	krmod.erase(krmod.begin(),krmod.end());
	// Save to disk
	fout<<","<<Cell->at(ii)<<","<<Pos->at(ii)<<","<<Age->at(ii)<<"\n";
	//if (Cell->at(ii)=='T')
	//	fout<<","<<Po2[ku];
	
	}; // end of for

// Calculation of the recruitment of stromal cells
	int oo = 0;
	while(round(par[4]*counts[0])!=oo)
	{oo++;
	randNi = round(rand()/(double)(RAND_MAX + 1)*(FreeSpaces.size()-1));
	Cell->push_back('S');
	Age->push_back(0);
	Pos->push_back(FreeSpaces[randNi]);
	counts[3]++;counts[4]++;
	fout<<",,S,"<<Pos->back()<<","<<0<<"\n";
	}; // recruit stroma
} // end of Probabilities fcn

int main(){
	// ****************** Declaration of variables ****************************
	double tEnd = 720, Dt = 5,// time in hours 1269 and 720 hours
		 par[8] = {9.23,18.17,139,184,0.0031,0.0013,0.064,7.54}, // MCF7 0 At, 1 Bt, 2 Bh, 3 Bn, 4 ks, 5 Mus, 6 hH, 7 kR
         //par[8] = {11.4,23.2,139,184,0.017,0.0013,0.064,185}, // Calu3 0 At, 1 Bt, 2 Bh, 3 Bn, 4 ks, 5 Mus, 6 hH, 7 kR
		 //par[8] = {8.7,17.2,139,184,0.0041,0.0013,0.064,185}, // Calu6 0 At, 1 Bt, 2 Bh, 3 Bn, 4 ks, 5 Mus, 6 hH, 7 kR
		A = 30,B = 0,Po2max = .46,Ckmax = .3,pi = 3.1415926535897, // paramters + initial conditions
		 ns = 50.0,//grid size
		 ff = .79,// Confluence of the injection (.79 for MCF7,Calu3\.27 for Calu6)
		 Droplet = .1, // cm3 droplet of the implant (equivalent to .1cm3)
		VWG = 2,//Volume of the whle grid
		CellVox= VWG/(pi/6.0*pow(0.003,3.0))/pow(ns,2.0);
	const double AR = pow(VWG/ns/ns,.5);// cm3/pixel Aspect ratio for a tumour of maximum VWG
	par[7] = par[7]*AR;// cm-1 --> voxel-1
	double t = 0,tTol = .001;

	// define other parameters of the program
	double randN = 0.0;
	int randNi = 0, refNum = 0, counts[5] = {0},
		//x0 = pow(0.1*ff*3.0/4.0*pow(3.14,.5),2.0/3.0)/pow(AR,2),
		x0 = Droplet/pow(AR,2.0),req = sqrt(x0/3.14);//voxel2/voxel
	vector<int>kr(8,0); 
	vector<bool>myVect(8,false);
	bool dummy = true;
		
	// Experiment design matrices
	int icount = 0;

	// Re-establish t
	t = 0;counts[0] =0;counts[1] =0;counts[2] =0;counts[3] =0;counts[4] =0;

	// get the time
	time_t now = time(0);  
	tm *ltm = localtime(&now);
	// convert to string
	string iist;ostringstream convert;convert << icount;iist = convert.str();
    string namefile1 = "50ResultsMCF72_21042016.csv"; 
	string namefile2 = "50ProfilesMCF72_21042016.csv"; 
	ofstream fout(namefile1.c_str());
	ofstream ffout(namefile2.c_str());
	//fout.open();
	fout << "Time, Tumour, Hypoxic, Necrotic, Stroma, Total, Elapsed Time";
	ffout << "Time, Oxygen, Cell, Positions, Age \n";

	// Initial values
   	vector <char>* cellType = new vector<char>(0);
	vector <double>* cellAge = new vector<double>(0); // cell time
	vector <int>* cellPosition = new vector<int>(0); // cell position
	vector<double>Po2(0);

			for	(int xx = 0;xx<ns;xx++)	// make sure x is fine
			{
				for(int yy = 0;yy<ns;yy++) // make sure y is fine
				{
				int rnow = sqrt(pow(yy-(ns-1)/2,2)+pow(xx-(ns-1)/2,2));
			
					// Tumour cells
					double myrand = rand()/(double)(RAND_MAX + 1);
					if(rnow<req&&myrand<ff) // is y lower than the circle? 
					{cellPosition->push_back(yy*ns+xx);
					cellAge->push_back(rand()/(double)(RAND_MAX + 1)*24);// Assign random age
					cellType->push_back('T');
					counts[0]++;counts[4]++;
					}
					else if(rnow>req&&rnow<2+req)
					{// Stromal cells
					cellPosition->push_back(yy*ns+xx);
					cellAge->push_back(rand()/(double)(RAND_MAX + 1)*24);// Assign random age
					cellType->push_back('S');
					counts[3]++;counts[4]++;
					}//if	
				}//yy
			}// xx

     //clock declaration
		clock_t start, end; 
	
   // ****************** The algorithm starts ****************************

	while(t < tEnd){  
    start = clock();
	t += Dt;
	ffout <<t;

// Calculate the profiles
	Po2 = myProfile(par[7],Po2max,ns,*cellPosition,*cellType);
// Probabilities - Decitions
	Probabilities(par,cellType,cellPosition,cellAge,Po2,ns,Dt,counts,ffout);

// calculate performance
	end = clock();
	
	fout<<"\n" << t <<","<< counts[0]<< ","<< counts[1]<< ","<<counts[2] <<"," << counts[3]<< ","<< counts[4]<< ","<<(double)(end-start)/CLOCKS_PER_SEC;

	}; // end of while
	ffout.close();
    fout.close();

return 0;}