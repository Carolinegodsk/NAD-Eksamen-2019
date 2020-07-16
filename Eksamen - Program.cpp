//============================================================================
// Name        : Program.cpp
// Author      : S128441
// Version     :
// Copyright   : Your copyright notice
// Description : NAD Eksamensprojekt 2019
//============================================================================

#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <chrono>
#include <ctime>
using namespace std;

#define NMAX 20
#define eps 0.0000000001

// Dele

void Brugerintro();
void Del_a(char Valgabcd, double A[NMAX][NMAX], double b[NMAX],
        int &n, int &m);
void Del_b(char Valgabcd, char ValgGemAbb, double A[NMAX][NMAX], double b[NMAX],
		double TM[NMAX][NMAX+1], int &n, int &m);
void Del_c(char Valgabcd, double A[NMAX][NMAX], double b[NMAX], int &n, int &m);
void Del_d(char Valgabcd, char &ValgK, double K1[NMAX][NMAX],
		double K2[NMAX][NMAX], int &n, int &m);
void Del_NL(char ValgNLelSVD, double A[NMAX][NMAX], double Mlr[NMAX][NMAX+1],
		double b[NMAX], double resx[NMAX], int &n, int &m);
void Del_SVD(char ValgNLelSVD, char ValgZ, char ValgK, char ValgAnvendSVD,
		char Valgabcd, double Z[NMAX][NMAX], double A[NMAX][NMAX], double b[NMAX],
		double K1[NMAX][NMAX], double K2[NMAX][NMAX], double An2[NMAX][NMAX],
		double d[NMAX], double bd[NMAX-1], double ev[NMAX], double Qv[NMAX][NMAX],
		int &n, int &m);

// Funktioner

// Opgave 1

void DanBilledMatrixK1(double K1[NMAX][NMAX], int n);
void DanBilledMatrixK2(double K2[NMAX][NMAX], int n);

// Opgave 2

double DanMaxAbsDif(double M1[NMAX][NMAX], double M2[NMAX][NMAX], int n, int m);
void IndtastRMatrix(double R[NMAX][NMAX], int &n, int &m);
void DantempcMatrix(double tempc[NMAX][NMAX], double kopiA[NMAX][NMAX], int n,
		int m, int m_gl);
void TjekSymmetri(double M1[NMAX][NMAX], int n, int m);

// Opgave 3

void FuldPivotering(double Mlr[NMAX][NMAX+1], int j, int xindex[NMAX], int n);
void BackwardsSubstitution(double TM[NMAX][NMAX+1], double x[NMAX], int n);
void Gauss(double TM[NMAX][NMAX+1], int n, int &bs);
void DelvisPivotering(double TM[NMAX][NMAX+1], int j, int n);
void Singularitet();

// Opgave 4

void KonstRk(double a[NMAX], double b[NMAX-1], double c[NMAX-1],
		double s[NMAX-1], int n);
void KonstAk(double a[NMAX], double b[NMAX-1], double c[NMAX-1],
		double s[NMAX-1], int n);
void UpdateQonIterk(double Q[NMAX][NMAX], double c[NMAX-1],
		double s[NMAX-1], int n);
void MaxVaerdi(double v[NMAX], double &maxb, int n);
void SamlMatrix(double A[NMAX][NMAX], double a[NMAX], int n);
void QR(double A[NMAX][NMAX], double a[NMAX], double b[NMAX-1],
		double c[NMAX-1], double s[NMAX-1], int n);
void Householder(double A[NMAX][NMAX], double Q0[NMAX][NMAX],
		double Ak[NMAX][NMAX], double a[NMAX], double b[NMAX-1], int n);
void Danaogb(double A[NMAX][NMAX], double a[NMAX], double b[NMAX-1], int n);
void QRWilkinson(double Q[NMAX][NMAX], double a[NMAX], double b[NMAX-1],
		double ev[NMAX], int n);

// Opgave 5

void Ombyt(double &x, double &y);
void Sorter(double LamVals[NMAX], double WMatUS[NMAX][NMAX], double
            LamValsSo[NMAX], double VMatSo[NMAX][NMAX], int n, int m);
int Rang(double LamVals[NMAX], int n);
void PseudoInvers(double U1[NMAX][NMAX], double D1M[NMAX][NMAX],
		double V1[NMAX][NMAX], double Aplus[NMAX][NMAX], int n, int r);
void SVD(double Z[NMAX][NMAX], double LamValsSo[NMAX], double
		VMatSo[NMAX][NMAX], double U1[NMAX][NMAX], double D1M[NMAX][NMAX],
		double V1[NMAX][NMAX], int n, int r);

// Hjaelpefunktioner

void IndhentMatrix(double A[NMAX][NMAX], int &n, int &m, string filnavn);
void IndhentTotalmatrix(double A[NMAX][NMAX], double y[NMAX], int &n, int &m,
		string filnavn);
void IndtastMatrix(double A[NMAX][NMAX], int &n, int &m);
void IndtastTotalmatrix(double A[NMAX][NMAX], double b[NMAX], int &n, int &m);
void IndtastVektor(double v[NMAX], int n);
void DanTotalmatrix(double TM[NMAX][NMAX+1],double A[NMAX][NMAX], double b[NMAX],
		int n, int m);
void UdskrivMatrix(double M[NMAX][NMAX], int n, int m);
void UdskrivTotalmatrix(double TM[NMAX][NMAX+1], int n, int m);
void UdskrivVektor(double v[NMAX], int n);
void GemMatrix(double A[NMAX][NMAX], int &n, int &m, string stdnavn);
void GemTotalmatrix(double TM[NMAX][NMAX+1], int &n, int &m, string stdnavn);
bool FilEksisterende(string filnavn);
void Enhedsmatrix(double Q[NMAX][NMAX], int n);
void Nulmatrix(double Q[NMAX][NMAX], int n, int m);
void KopierVektor(double v[NMAX], double w[NMAX], int n);
void KopierMatrix(double A[NMAX][NMAX], double B[NMAX][NMAX], int n, int m);
void NormerVektorer(double v[NMAX][NMAX], int n, int m);
void TransponerMatrix(double M[NMAX][NMAX], double MT[NMAX][NMAX], int n, int m);
void MatrixProdukt(double M1[NMAX][NMAX], double M2[NMAX][NMAX],
		double M3[NMAX][NMAX], int n, int m);
void MatrixVektorProdukt(double M[NMAX][NMAX], double v[NMAX],
		double vprod[NMAX],	int n, int m);
void VektorAddition(double v1[NMAX], double v2[NMAX], double res[NMAX], int n);
void VektorProdukt(double v1[NMAX], double v2[NMAX], double res[NMAX], int n);


int main() {
	int m, n;
	double A[NMAX][NMAX], b[NMAX], K1[NMAX][NMAX], K2[NMAX][NMAX],
	TM[NMAX][NMAX+1], resx[NMAX], d[NMAX], bd[NMAX-1],
	ev[NMAX], Qv[NMAX][NMAX],
	Z[NMAX][NMAX], An2[NMAX][NMAX];
	char Valgabcd, ValgGemAbb, ValgK,
	ValgNLelSVD, ValgZ, ValgAnvendSVD, ValgIgen;
	string filnavn, x, stdnavn;

	do {
		Brugerintro();

//		Del a/b/c/d ---------------------------------------------------------

		cout << "\nVil du: \na) Indlaes n, m, A og b fra en datafil \nb) "
			"Indtast n, m, A og b "
			"\nc) Dan A og b, saa Ax=b har en brugervalgt MKL x^* "
			"\nd) Fremstil en billedmatrix K(nxn) ud fra design 1: K=K_2 eller"
			" 2: K=K_2 ";
		do {
			cout << "\n\nVaelg a, b, c eller d: ";
			cin >> Valgabcd;
		}
		while (Valgabcd != 'a' && Valgabcd != 'b' && Valgabcd != 'c' &&
				Valgabcd != 'd' && Valgabcd != 'A' && Valgabcd != 'B' &&
				Valgabcd != 'C' && 	Valgabcd != 'D');
		switch (Valgabcd){
			case 'a':{
			case 'A':
				Del_a(Valgabcd,A,b,n,m);
				break;
			}
			case 'b':{
			case 'B':
				Del_b(Valgabcd,ValgGemAbb,A,b,TM,n,m);
				break;
			}
			case 'c':{
			case 'C':
				Del_c(Valgabcd,A,b,n,m);
				break;
			}
			case 'd':{
			case 'D':
				Del_d(Valgabcd,ValgK,K1,K2,n,m);
				break;
			}
		}

//		Del NL/SVD -----------------------------------------------------------

		cout << "\n-----------------------------------------------------------"
				"-------------------------------------\n";
		cout << "\nVil du: \n1) NL - Loese normalligninger og finde mindste "
			"kvadraters loesningen x^* "
			"\n2) SVD - Udfoere QR-metoden og SVD-faktorisering ";
		do {
			cout << "\n\nVaelg 1 eller 2: ";
			cin >> ValgNLelSVD;
		}
		while (ValgNLelSVD != '1' && ValgNLelSVD != '2');
		switch (ValgNLelSVD){
			case '1':
				Del_NL(ValgNLelSVD,A,Mlr,b,resx,n,m);
				break;
			case '2':
				Del_SVD(ValgNLelSVD,ValgZ,ValgK,ValgAnvendSVD,Valgabcd,Z,A,b,K1,
					K2,An2,d,bd,ev,Qv,n,m);

				break;
			}

		cout << "\n------------------------------------------------------------"
				"------------------------------------\n";
		cout << "\nVil du koere hele programmet forfra? Indtast j eller n: ";
		cin >> ValgIgen;
	}
	while (ValgIgen == 'j' or ValgIgen == 'J');
	if (ValgIgen == 'n' or 'N'){
		cout << "\nDu har valgt at afslutte programmet. Vi ses en anden gang!";
	}
	else {
		cout << "\nDu har ikke indtastet 'j' eller 'n', og programmet afsluttes"
			" derfor.";
	}
	return 0;
}

// Dele ------------------------------------------------------------------------

void Brugerintro(){
	cout << "\n----------------------------------------------------------------"
			"--------------------------------\n";
	cout << endl << "Dette program har foelgende indhold:" << endl << endl <<
		"1)	Find mindste kvadraters loesningen (MKL): x^* til det overbestemte "
			<< endl <<
			"	ligningssystem (OLS): Ax=b, hvor A(nxm) og b(nx1) med n>m." << endl
		<< "	x^* findes som loesning til normalligningerne (NL): A^TAx=A^Tb." <<
			endl << endl << "2)	Find egenloesningerne {λ_i,v_i} i=1,...,n for A(nxn),"
			" hvor A=A^T." << endl << endl <<
		"3)	Find SVD-faktorerne U_1, V_1 og D_1 for den smalle SVD-faktorisering "
			<< endl << "	af Z(nxm), hvor n>=m og Z=U_1D_1V_1^T." << endl <<
			"	-	Beregn x^*=A^+B med Z=A, hvor A^+ er pseudoinvers til A." << endl
		<< "	-	Beregn K_N, som er reduceret rang N tilnaermelse til en " << endl
		<< "		billedmatrix K, hvor Z=K." << endl << endl <<
		"Frembringelse af datagrundlag a), b), c) og d): " << endl <<
		"a)	Indlaes n, m, A og b fra en datafil." << endl <<
		"b)	Indtast n, m, A og b." << endl <<
		"c)	Dan A og b, saa Ax=b har en brugervalgt MKL x^*." << endl <<
		"d)	Fremstil en billedmatrix K(nxn) ud fra design 1: K=K_2 eller 2: K=K_2"
		<< endl;
	cout << "\n-----------------------------------------------------------------"
			"-------------------------------\n";
}

void Del_a(char Valgabcd, double A[NMAX][NMAX], double b[NMAX], int &n, int &m){
	cout << "\nDu har valgt: " << Valgabcd << endl
		<< "Indhold: Indlaes n, m, A og b fra en datafil" << endl;
	cout << "\nPaa dette sted i programmet udfoeres: "
		"	\nVaelg indfil" << endl;
	IndhentTotalmatrix(A,b,n,m,"x");
	cout << "\nPaa dette sted i programmet udfoeres: "
			"	\nLaes fra indfil: n, m, A, b" << endl;
	cout << "\nIndhentet matrix: \n";
	UdskrivMatrix(A,n,m);
	cout << "\nIndhentet vektor: \n";
	UdskrivVektor(b,n);
}

void Del_b(char Valgabcd, char ValgGemAbb, double A[NMAX][NMAX],
		double b[NMAX],	double TM[NMAX][NMAX+1], int &n, int &m){
	cout << "\nDu har valgt: " << Valgabcd << endl
		<< "Indhold: Indtast n, m, A og b" << endl;
	cout << "\nPaa dette sted i programmet udfoeres: "
		"	\nTast: n, m, A, b" << endl;
	char RetMat = 'j';
	do {
		IndtastTotalmatrix(A,b,n,m);
		DanTotalmatrix(TM,A,b,n,m);
		UdskrivTotalmatrix(TM,n,m);
		cout << "\nVil du rette matricen og/eller vektoren? \nIndtast j "
				"eller n: ";
		cin >> RetMat;
	}
	while (RetMat == 'j' or RetMat == 'J');
	cout << "\nVil du gemme indtastet data i en fil? " << endl;
	do {
		cout << "Indtast j eller n: ";
		cin >> ValgGemAbb;
	}
	while (ValgGemAbb != 'j' && ValgGemAbb != 'J' && ValgGemAbb != 'n'
		&& 	ValgGemAbb != 'N');
	switch (ValgGemAbb){
		case 'j':
		case 'J':
			cout << "\nDu har valgt: " << ValgGemAbb << endl <<
			"\nPaa dette sted i programmet udfoeres: "
				"	\nFoelgende udskrives til fil: n, m, A, b" << endl;
			GemTotalmatrix(TM,n,m,"Abb");
			break;
		case 'n':
		case 'N':
			cout << "\nDu har valgt: " << ValgGemAbb << endl <<
				"\nDu har valgt ikke at gemme data i filen." << endl;
			break;
	}
}

void Del_c(char Valgabcd, double A[NMAX][NMAX], double b[NMAX], int &n, int &m){
	double An[NMAX][NMAX], kopiA[NMAX][NMAX], TM[NMAX][NMAX+1], Q[NMAX][NMAX],
		R[NMAX][NMAX], xstj[NMAX], bp[NMAX], tempc[NMAX][NMAX], tstj[NMAX], e[NMAX];
	int m_gl;
	char Valgn, ValgGemAbc;
	cout << "\nDu har valgt: " << Valgabcd << endl
		<< "Indhold: Dan A og b, saa Ax=b har en brugervalgt MKL x^*" << endl;
	cout << "\nPaa dette sted i programmet udfoeres: "
		"	Fra 'A8.txt': Laes A8" << endl;
	IndhentMatrix(A,n,m,"A8.txt");
	cout << "\nVaelg dimension: n=4 eller n=8? " << endl;
	do {
		cout << "Indtast 4 eller 8: ";
		cin >> Valgn;
	}
	while (Valgn != '4' && Valgn != '8');
	switch (Valgn){
		case '4':
			cout << "\nDu har valgt: " << Valgn << endl <<
				"\nPaa dette sted i programmet udfoeres: "
				"	Foelgende dannes: A_n=A_4" << endl;
			n = 4;
			cout << "\nMatrix A(" << n << "x" << n << "): \n";
			KopierMatrix(A,An,n,m);
			UdskrivMatrix(An,n,n);
			break;
		case '8':
			cout << "\nDu har valgt: " << Valgn << endl <<
				"\nPaa dette sted i programmet udfoeres: "
				"	Foelgende dannes: A_n=A_8" << endl;
			n = 8;
			cout << "\nMatrix A(" << n << "x" << n << "): \n";
			KopierMatrix(A,An,n,m);
			UdskrivMatrix(An,n,n);
			break;
	}
	m = n;
	m_gl = m;						// Bruges til dannelse af ctemp
	KopierMatrix(A,kopiA,n,m);		// Bruges til dannelse af ctemp
	cout << "\nPaa dette sted i programmet udfoeres: "
		"	Foelgende dannes: A og b (vaelg x^*)" << endl;
	do {
		cout << "\nVaelg m, hvor m<n: ";
		while (!(cin >> m)){
			cout << "\nIndtast venligst et heltal";
			cin.clear();
            cin.ignore(100,'\n');   // Tilpasset fra https://books.google.dk/
                                    // books?id=kAJ-DwAAQBAJ&pg=PA221&lpg=PA221&dq=
                                    // cin.clear();+cin.ignore(100,%27%5Cn%27);&source
                                    // =bl&ots=AZQ0Mc3n4Q&sig=ACfU3U0NVO_EJVnCmWZEiIw_
                                    // MhjNf11-rg&hl=da&sa=X&ved=2ahUKEwiIm-Gbv7viAhW
                                    // InhQKHTs2C0oQ6AEwEnoECAoQAQ#v=onepage&q=cin.
                                    // clear()%3B%20cin.ignore(100%2C'%5Cn')%3B&f=false
		}
	}
	while (m>n or m==n);
	// Dannelse af ctemp
	DantempcMatrix(tempc,kopiA,n,m,m_gl);
	cout << "\nMatrix A(" << n << "x" << m << ") er nu som foelgende: \n";
	UdskrivMatrix(An,n,m);
	// Ortonormal Q dannes med m soejler
	KopierMatrix(A,Q,n,m);
	NormerVektorer(Q,n,m);
	cout << "\nMatrix Q(" << n << "x" << m << ") med m=" << m << " ortonormale soejler: \n";
	UdskrivMatrix(Q,n,m);
	// R(mxm) dannes
	cout << "\nNu indtastes R(mxm), som er en oevre trekantsmatrix \n";
	IndtastRMatrix(R,m,m);
	cout << "\nDen indtastede oevre trekantmatrix R(" << m << "x" << m << ") er: \n";
	UdskrivMatrix(R,m,m);
	// A=QR beregnes
	MatrixProdukt(Q,R,An,n,m);
	cout << "\nA=QR: \n";
	UdskrivMatrix(An,n,m);
	// x* indtastes og bp beregnes
	cout << "\nNu indtastes x^*(" << m << "x1" << "), som er den oenskede loesning \n";
	IndtastVektor(xstj,m);
	MatrixVektorProdukt(An,xstj,bp,n,m);
	cout << "\nbp=Ax^* (projektionen af vektor b): \n";
	UdskrivVektor(bp,n);
	// t* indtastes og e beregnes
	cout << "\nNu indtastes t^*(" << n-m << "x1" << ") \n";
	IndtastVektor(tstj,m);
	MatrixVektorProdukt(tempc,tstj,e,n,m);
	cout << "\nMatrix tempc: \n";
	UdskrivMatrix(tempc,n,m);
	cout << "\nVektor e=ct^* (fejlvektoren): \n";
	UdskrivVektor(e,n);
	// b beregnes
	cout << "\nVektor b faas nu ved b=bp+e:  \n";
	VektorAddition(bp,e,b,n);
	UdskrivVektor(b,n);
	cout << "\nVil du gemme indtastet data i en fil? Indtast j eller n: " << endl;
	do {
		cout << "Indtast j eller n: ";
		cin >> ValgGemAbc;
	}
	while (ValgGemAbc != 'j' && ValgGemAbc != 'J' && ValgGemAbc != 'n' && ValgGemAbc != 'N');
	switch (ValgGemAbc){
		case 'j':
		case 'J':
			cout << "\nDu har valgt: " << ValgGemAbc << endl <<
				"\nPaa dette sted i programmet udfoeres: "
				"	Foelgende udskrives til fil: n, m, A, b" << endl;
			DanTotalmatrix(TM,A,b,n,m);
			GemTotalmatrix(TM,n,m,"Abc");
			break;
		case 'n':
		case 'N':
			cout << "\nDu har valgt: " << ValgGemAbc << endl <<
				"\nDu har valgt ikke at gemme data i filen." << endl;
			break;
	}
}

void Del_d(char Valgabcd, char &ValgK, double K1[NMAX][NMAX],
	double K2[NMAX][NMAX], int &n, int &m){
	cout << "\nDu har valgt: " << Valgabcd << endl
	<< "Indhold: Fremstil en billedmatrix K(nxn) ud fra design "
		"1: K=K_2 eller 2: K=K_2" << endl;
	cout << "\nPaa dette sted i programmet udfoeres: "
		"	\nFoelgende indtastes: n" << endl;
	do {
		cout << "\nIndtast oensket dimension, n (max " << NMAX
				<< "): ";
		while (!(cin >> n)){
			cout << "\nIndtast venligst et heltal";
			cin.clear();
            cin.ignore(100,'\n');   // Tilpasset fra https://books.google.dk/
                                    // books?id=kAJ-DwAAQBAJ&pg=PA221&lpg=PA221&dq=
                                    // cin.clear();+cin.ignore(100,%27%5Cn%27);&source
                                    // =bl&ots=AZQ0Mc3n4Q&sig=ACfU3U0NVO_EJVnCmWZEiIw_
                                    // MhjNf11-rg&hl=da&sa=X&ved=2ahUKEwiIm-Gbv7viAhW
                                    // InhQKHTs2C0oQ6AEwEnoECAoQAQ#v=onepage&q=cin.
                                    // clear()%3B%20cin.ignore(100%2C'%5Cn')%3B&f=false
		}
	}
	while (n > NMAX);
	cout << "\nVaelg K: K=K_1 eller K=K_2? " << endl;
	do {
		cout << "Indtast 1 eller 2: ";
		cin >> ValgK;
	}
	while (ValgK != '1' && ValgK != '2');
	switch (ValgK){
		case '1':
			cout << "\nDu har valgt: " << ValgK << endl <<
			"\nPaa dette sted i programmet udfoeres: "
				"	\nFoelgende dannes: K=K_1" << endl;
			DanBilledMatrixK1(K1,n);
			cout << "\nBilledmatrix K1: \n";
			UdskrivMatrix(K1,n,n);
			break;
		case '2':
			cout << "\nDu har valgt: " << ValgK << endl
			<< "\nPaa dette sted i programmet udfoeres: "
				"	\nFoelgende dannes: K=K_2" << endl;
			DanBilledMatrixK2(K2,n);
			cout << "\nBilledmatrix K2: \n";
			UdskrivMatrix(K2,n,n);
			break;
	}
}

void Del_NL(char ValgNLelSVD, double A[NMAX][NMAX],
		double Mlr[NMAX][NMAX+1],
	double b[NMAX], double resx[NMAX], int &n, int &m){
	int bs;
	double AT[NMAX][NMAX], ATA[NMAX][NMAX], ATb[NMAX],
	ATAATb[NMAX][NMAX+1];
	cout << "\nDu har valgt: " << ValgNLelSVD << endl
		<< "Indhold: NL - Loese normalligninger og finde "
		"mindste kvadraters loesningen x^*" << endl;
	cout << "\nPaa dette sted i programmet udfoeres: "
		"	\nFoelgende dannes: NL: A^TAx=A^Tb" << endl;
	// A^TAx=A^Tb
	TransponerMatrix(A,AT,n,m);
	cout << "\nDen transponerede matrix A^T:\n";
	UdskrivMatrix(AT,m,n);
	MatrixProdukt(AT,A,ATA,n,m);
	cout << "\nMatrixprodukt A^T*A:\n";
	UdskrivMatrix(ATA,m,m);
	MatrixVektorProdukt(AT,b,ATb,n,m);
	cout << "\nA^T*b:\n";
	UdskrivVektor(ATb,m);
	DanTotalmatrix(ATAATb,ATA,ATb,m,m);

	// Loesning af ligningssystem
	cout << "\nPaa dette sted i programmet udfoeres: "
		"	\nNL loeses og MKL findes" << endl;
	Gauss(ATAATb,m,bs);
	cout << "\nDen Gauss−eliminerede matrix:\n";
	UdskrivTotalmatrix(ATAATb,m,m);
	if(bs==1){
		BackwardsSubstitution(ATAATb,resx,m);
		cout << "\nLoesningsvektor x:\n";
		UdskrivVektor(resx,m);
	}
	else{
		Singularitet();
	}
}

void Del_SVD(char ValgNLelSVD, char ValgZ, char ValgK, char ValgAnvendSVD,
		char Valgabcd, double Z[NMAX][NMAX], double A[NMAX][NMAX],
		double b[NMAX],	double K1[NMAX][NMAX], double K2[NMAX][NMAX], double
		An2[NMAX][NMAX], double d[NMAX], double bd[NMAX-1], double ev[NMAX],
		double Qv[NMAX][NMAX], int &n, int &m){
	double ZT[NMAX][NMAX], ZTZ[NMAX][NMAX], Q0[NMAX][NMAX], Q[NMAX][NMAX],
	evSo[NMAX], QvSo[NMAX][NMAX], U1[NMAX][NMAX],
	V1[NMAX][NMAX], U1D1[NMAX][NMAX], V1T[NMAX][NMAX], U1D1V1T[NMAX][NMAX],
	Aplus[NMAX][NMAX], Aplusb[NMAX], D1M[NMAX][NMAX], kva_ev[NMAX];
	int N, r;

	cout << "\nDu har valgt: " << ValgNLelSVD << endl
		<< "Indhold: SVD - Udfoere QR-algoritmen og SVD-faktorisering" << endl;
	cout << "\nDu har tidligere valgt del " << Valgabcd << endl;
	switch (Valgabcd){
		case 'A':
		case 'a':
		case 'B':
		case 'b':
		case 'C':
		case 'c':
			cout << "Foelgende dannes: Z=A\n";
			KopierMatrix(A,Z,n,n);
			cout << "\nMatrix Z: \n";
			UdskrivMatrix(Z,n,n);
			break;
		case 'D':
		case 'd':
			switch (ValgK){
				case '1':
					cout << "Foelgende dannes: Z=K1\n";
					KopierMatrix(K1,Z,n,n);
					cout << "\nMatrix Z: \n";
					UdskrivMatrix(Z,n,n);
					break;
				case '2':
					cout << "Foelgende dannes: Z=K2\n";
					KopierMatrix(K2,Z,n,n);
					cout << "\nMatrix Z: \n";
					UdskrivMatrix(Z,n,n);
					break;
			}
			break;
	}

	cout << "\nPaa dette sted i programmet udfoeres: "
		"	\nEgenloesninger {λ_i,v_i} findes for Z^TZ, idet "
		"Z^TZ er symmetrisk for alle Z\n";
	TransponerMatrix(Z,ZT,n,n);
	MatrixProdukt(ZT,Z,ZTZ,n,n);
	cout << "\nMatrix Z^T*Z: \n";
	UdskrivMatrix(ZTZ,n,n);
	cout << "\nResultat fra Householder: ";
	Householder(ZTZ,Q0,An2,d,bd,n);
	cout << "\nResultat fra QR med Givens rotationer og "
			"Wilkinson shifts:\n";
	QRWilkinson(Q,d,bd,ev,n);
	for (int i=0; i<n; i++){
		kva_ev[i] = sqrt(ev[i]);
	}
	cout << "\nKvadrerede (tilnærmede) egenværdier er: \n";
	UdskrivVektor(kva_ev,n);
	MatrixProdukt(Q0,Q,Qv,n,n);
	cout << "\nNormerede egenvektorer er: \n";
	UdskrivMatrix(Qv,n,n);

	cout << "\nSVD påbegyndes: \n";
	cout << "\nPaa dette sted i programmet udfoeres: "
		"	\nFoelgende dannes: U_1, D_1, V_1 for Z=U_1D_1V_1^T\n";
	Sorter(ev,Qv,evSo,QvSo,n,n);
	cout << "\nSorterede egenværdier for Z^TZ: \n";
	UdskrivVektor(evSo,n);
	cout << "\nSorterede normerede egenvektorer: \n";
	UdskrivMatrix(QvSo,n,n);
	r = Rang(evSo,n);
	SVD(Z,evSo,QvSo,U1,D1M,V1,n,r);
	cout << "\nMatrix U1: \n";
	UdskrivMatrix(U1,n,r);
	cout << "\nDiagonalmatrix D1 (med egenværdier for Z): \n";
	UdskrivMatrix(D1M,r,r);
	TransponerMatrix(V1,V1T,n,r);
	cout << "\nTransponeret matrix V1^T: \n";
	UdskrivMatrix(V1T,r,n);
	// Tjek af Z
	cout << "\nTjek af SVD påbegyndes: \n";
	MatrixProdukt(U1,D1M,U1D1,n,r);
	MatrixProdukt(U1D1,V1T,U1D1V1T,n,n);
	cout << "\nZ=U1*D1*V1^T: \n";
	UdskrivMatrix(U1D1V1T,n,n);
	cout << "\nZ (oprindelig): \n";
	UdskrivMatrix(Z,n,n);
	cout << "\nDen maksimale difference mellem Z og Z=U1*D1*V1^T = "
			<< DanMaxAbsDif(Z,U1D1V1T,r,r) << endl;
    
	cout << "\nVaelg naeste skridt: A^+ eller K_N? " << endl;
	do {
		cout << "Indtast 1 (for A^+) eller 2 (K_N): ";
		cin >> ValgAnvendSVD;
	}
	while (ValgAnvendSVD != '1' && ValgAnvendSVD != '2');
	switch (ValgAnvendSVD){
		case '1':
			cout << "\nDu har valgt: " << ValgAnvendSVD << endl <<
				"\nPaa dette sted i programmet udfoeres: "
				"	\nFoelgende findes: A^+ og x^*=A^+b" << endl;
			PseudoInvers(U1,D1M,V1,Aplus,n,n);
			cout << "\nMatrix A^+: \n";
			UdskrivMatrix(Aplus,n,n);
			cout << "\nVektor b: \n";
			UdskrivVektor(b,n);
			MatrixVektorProdukt(Aplus,b,Aplusb,n,n);
			cout << "\nLoesningsvektor x^*: \n";
			UdskrivVektor(Aplusb,n);
			break;
		case '2':
			cout << "\nDu har valgt: " << ValgAnvendSVD << endl
				<< "\nPaa dette sted i programmet udfoeres: "
				"	\nVaelg N. Herefter dannes K_N" << endl;
			do {
				cout << "\nIndtast oensket N, hvor 1<=N<=r (rang=" <<
						Rang(evSo,n) << "): ";
				while (!(cin >> N)){
					cout << "\nIndtast venligst et heltal";
					cin.clear();
                    cin.ignore(100,'\n');   // Tilpasset fra https://books.google.dk/
                                            // books?id=kAJ-DwAAQBAJ&pg=PA221&lpg=PA221&dq=
                                            // cin.clear();+cin.ignore(100,%27%5Cn%27);&source
                                            // =bl&ots=AZQ0Mc3n4Q&sig=ACfU3U0NVO_EJVnCmWZEiIw_
                                            // MhjNf11-rg&hl=da&sa=X&ved=2ahUKEwiIm-Gbv7viAhW
                                            // InhQKHTs2C0oQ6AEwEnoECAoQAQ#v=onepage&q=cin.
                                            // clear()%3B%20cin.ignore(100%2C'%5Cn')%3B&f=false
				}
			}
			while (N>Rang(evSo,n) or N<1);
			SVD(Z,evSo,QvSo,U1,D1M,V1,n,N);
			cout << "\nMatrix U1: \n";
			UdskrivMatrix(U1,n,N);
			cout << "\nDiagonalmatrix D1: \n";
			UdskrivMatrix(D1M,N,N);
			TransponerMatrix(V1,V1T,n,N);
			cout << "\nTransponeret matrix V1^T: \n";
			UdskrivMatrix(V1T,N,n);
			// Tjek af Z
			cout << "\nTjek af reduceret rang påbegyndes: \n";
			MatrixProdukt(U1,D1M,U1D1,n,N);
			MatrixProdukt(U1D1,V1T,U1D1V1T,n,n);
			cout << "\nZ=U1*D1*V1^T: \n";
			UdskrivMatrix(U1D1V1T,n,n);
			cout << "\nZ (oprindelig): \n";
			UdskrivMatrix(Z,n,n);
			cout << "\nDen maksimale difference mellem Z og Z=U1*D1*V1^T = "
					<< DanMaxAbsDif(Z,U1D1V1T,N,N) << endl;
			break;
	}
}

// Funktioner ---------------------------------------------------------------

// Opgave 1

void DanBilledMatrixK1(double K1[NMAX][NMAX], int n){
	int i;
	Nulmatrix(K1,n,n);
	for (i=0; i<n; i++){
		K1[i][i] = 1;
		K1[i+1][i] = 1;
		K1[i][i+1] = 1;
	}
	K1[n-1][n-1] = 1;
	for (i=0; i<n; i++){
		K1[n-i-1][i] = 1;
		K1[n-i-1][i+1] = 1;
		K1[n-i-2][i] = 1;
	}
	K1[0][n] = 1;
}

void DanBilledMatrixK2(double K2[NMAX][NMAX], int n){
	int i;
	Nulmatrix(K2,n,n);
	for (i=0; i<n; i++){
		K2[i][i] = 1;
		K2[i+1][i] = 1;
		K2[i+2][i] = 1;
		K2[i][i+1] = 1;
		K2[i][i+2] = 1;
	}
	K2[n-1][n-1] = 1;
	for (i=0; i<n; i++){
		K2[n-i][i] = 1;
		K2[n-i][i+1] = 1;
		K2[n-i-1][i] = 1;
		K2[n-i-2][i] = 1;
		K2[n-i-3][i] = 1;
	}
	K2[0][n] = 1;
}

// Opgave 2

double DanMaxAbsDif(double M1[NMAX][NMAX], double M2[NMAX][NMAX], int n,
		int m){
	double MaxAbsDif=0.0, M[NMAX][NMAX];
	int i, j;
	for (i=0; i<n; i++){
		for (j=0; j<m; j++){
			M[i][j] = M1[i][j] - M2[i][j];
		}
	}
	MaxAbsDif=0.0;
	for (i=0; i<n; i++){
		for (j=0; j<m; j++){
			if (fabs(M[i][j]) > MaxAbsDif){
				MaxAbsDif = M[i][j];
			}
		}
	}
	return MaxAbsDif;
}

void IndtastRMatrix(double R[NMAX][NMAX], int &n, int &m){
	int i, j;
	for(i=0; i<n; i++){
		for(j=0; j<m; j++){
			if (i>j){
				R[i][j] = 0;
			}
			else {
				cout << "Indtast element i raekke nr. " << i+1 << " "
						"og soejle nr. " << j+1 << " her: ";
				cin >> R[i][j];
			}
		}
	}
}

void DantempcMatrix(double tempc[NMAX][NMAX], double kopiA[NMAX][NMAX],
		int n, int m, int m_gl){
	int k, i, j;
	double tempv[NMAX];
	k = 0;
	for (i=0; i<n; i++){
		for (j=m; j<m_gl; j++){
			tempv[k] = kopiA[i][j];
			k = k+1;
		}
	}
	k = 0;
	for (i=0; i<n; i++){
		for (j=0; j<m; j++){
			tempc[i][j] = tempv[k];
			k = k+1;
		}
	}
}

void TjekSymmetri(double M1[NMAX][NMAX], int n, int m){
	double M2[NMAX][NMAX], I[NMAX], sum, MaxAbsDif=0;
	TransponerMatrix(M1,M2,n,n);
	cout << "\nMaxAbsDif er: " << DanMaxAbsDif(M1,M2,n,n) << endl;
	if (fabs(DanMaxAbsDif(M1,M2,n,n)) < eps){
		cout << "\nMatricen er symmetrisk" << endl;

	}
	else cout << "Matricen er ikke symmetrisk" << endl;
	// Test for ortonormalitet;
	for(int j=0; j<n; j++){
		sum = 0;
		for(int i=0; i<n; i++){
			sum = sum + M1[i][j] * M1[i][j];
		}
		I[j] = sum;
	}
	sum = 0;
	for (int j=0; j<n; j++){
		sum += I[j];
	}
	if (sum==n){
		cout << "\nMatricen er ortonormal" << endl;
	}
	else {
		cout << "\nMatricen er ikke ortonormal" << endl;
	}
}

// Opgave 3

void FuldPivotering(double Mlr[NMAX][NMAX+1], int j, int xindex[NMAX],
		int n){
	double max=0, temp=0, xindextemp=0;
	int i, k, maxposcol=0, maxposrow=0;
	j=0;
	max = Mlr[j][j];
	for (i=j; i<n; i++){
		for (k=j; k<n; k++){
			if (fabs(Mlr[i][k]) > fabs(max)){
				max = Mlr[i][k];
				maxposrow = i;
				maxposcol = k;
			}
		}
	}
	for (k=j; k<=n; k++) {
		temp = Mlr[j][k];
		Mlr[j][k] = Mlr[maxposrow][k];
		Mlr[maxposrow][k] = temp;
	}

	xindextemp = xindex[j];
	xindex[j] = xindex[maxposcol];
	xindex[maxposcol] = xindextemp;

	for (k=j; k<n; k++){
		for (i=j; i<j+1; i++){
			temp = Mlr[k][i];
			Mlr[k][i] = Mlr[k][maxposcol];
			Mlr[k][maxposcol] = temp;
		}
	}
}

void Gauss(double TM[NMAX][NMAX+1], int n, int &bs){
	double factor, max=0;
	int j, k, i, xindex[NMAX];
	bs=1;
	for(j=0;j<=n-2;j++){
		max = TM[j][j];
		for (i=j; i<n; i++){
			for (k=j; k<n; k++){
				if (fabs(TM[i][k]) > fabs(max)){
					max = TM[i][k];
				}
			}
		}
		if (fabs(max) >= eps){
			DelvisPivotering(TM,j,n);
		}
		else if (fabs(max) < eps){
			FuldPivotering(TM,j,xindex,n);
		}
		if(fabs(TM[j][j])<eps){
			bs=0;
			break;
		}
		for(int i=j+1; i<=n-1; i++){
			factor = -TM[i][j]/TM[j][j];
			TM[i][j]=0;
			for(k=j+1; k<=n; k++){
				TM[i][k] += factor*TM[j][k];
			}
		}
	}
	if(fabs(TM[n-1][n-1])<eps){
		bs=0;
	}
}

void BackwardsSubstitution(double TM[NMAX][NMAX+1], double x[NMAX],
		int n){
	double sum;
	x[n-1] = TM[n-1][n]/TM[n-1][n-1];
	for(int i=n-2; i>=0; i--){
		sum=0;
		for(int j=i+1; j<=n-1; j++){
			sum+=TM[i][j]*x[j];
		}
		x[i]=(TM[i][n]-sum)/TM[i][i];
	}
}

void DelvisPivotering(double TM[NMAX][NMAX+1], int j, int n){
	double max, temp;
	int maxpos, i, k;
	max = fabs(TM[j][j]);
	maxpos = j;
	for(i=j+1; i<=n-1; i++){
		if(fabs(TM[i][j]) > max){
			max = fabs(TM[i][j]);
			maxpos = i;
		}
	}
	for (k=j; k<=n; k++){
		temp = TM[j][k];
		TM[j][k] = TM[maxpos][k];
		TM[maxpos][k] = temp;
	}
}

void Singularitet(){
	cout << "\nUnder Gauss-eliminationen har det vist sig, at matricen"
			" formodes at vaere singulaer (ikke-invertibel). "
			"\nBackwards substitution kan derfor ikke udfoeres.";
}

// Opgave 4

void KonstRk(double a[NMAX], double b[NMAX-1], double c[NMAX-1],
		double s[NMAX-1], int n){
	double t = b[0], r;
	for (int i=0; i<n-1; i++){
		r = sqrt(pow(a[i],2)+pow(t,2));
		if (r==0){
			c[i] = 0;
			s[i] = 0;
		}
		else {
			c[i] = a[i]/r;
			s[i] = t/r;
		}
		a[i] = r;
		t = b[i];
		b[i] = t*c[i]+a[i+1]*s[i];
		a[i+1] = -t*s[i]+a[i+1]*c[i];
		t = b[i+1];
		b[i+1] = t*c[i];
	}
}

void KonstAk(double a[NMAX], double b[NMAX-1], double c[NMAX-1],
		double s[NMAX-1], int n){
	for (int i=0; i<n-1; i++){
		a[i] = a[i]*c[i]+b[i]*s[i];
		b[i] = a[i+1]*s[i];
		a[i+1] = a[i+1]*c[i];
	}
}

void UpdateQonIterk(double Q[NMAX][NMAX], double c[NMAX-1], double s[NMAX-1],
		int n){
	double v1, v2;
	int i,j;
	for(i=0; i<n-1; i++){
		for(j=0; j<n; j++){
			v1 = c[i]*Q[j][i]+s[i]*Q[j][i+1];
			v2 = (-s[i])*Q[j][i]+c[i]*Q[j][i+1];
			Q[j][i] = v1;
			Q[j][i+1] = v2;
		}
	}
}

void MaxVaerdi(double v[NMAX], double &maxb, int n){
	int tal=0, i;
	for (i=tal+1; i<n; i++){
		if (fabs(v[i]) < fabs(v[tal])){
			tal = i;
		}
	}
	maxb = v[tal];
}

void SamlMatrix(double A[NMAX][NMAX], double a[NMAX], int n){
	int i;
	Nulmatrix(A,n,n);
	for (i=0; i<n; i++){
		A[i][i] = a[i];
	}
}

void QR(double Q[NMAX][NMAX], double a[NMAX], double b[NMAX-1],
		double c[NMAX-1], double s[NMAX-1], int n){
	double maxb;
	int k=0, i;
	for (i=0; i<n; i++){
		s[i] = c[i] = 0;
	}
	do{
		KonstRk(a,b,c,s,n);		// Her beregnes a og b i Rk samt c og s i Gk
		KonstAk(a,b,c,s,n);		// Her beregnes a og b i Ak, givet Gk og Rk
		UpdateQonIterk(Q,c,s,n);// Opdaterer Q: Q = Q * G1T * ,…, Gn-1T
		MaxVaerdi(b,maxb,n-1);
		k++; 					// Opdaterer k
	}
	while (fabs(maxb) > eps);
}

void Householder(double A[NMAX][NMAX], double Q0[NMAX][NMAX],
		double Ak[NMAX][NMAX], double a[NMAX], double b[NMAX-1], int n){
	int k, i, j;
	double H[NMAX][NMAX], Hk[NMAX][NMAX], An2[NMAX][NMAX], x[NMAX],
	w[NMAX], v[NMAX], q[NMAX], s, sum, c, r;
	KopierMatrix(A,Ak,n,n);
	Enhedsmatrix(Q0,n);
	for (k=0; k<n-2; k++){
		sum = 0;
		c = 0;
		for (i=0; i<n; i++){
			x[i] = Ak[i][k];
		}
		for (i=k+1; i<n; i++){
			sum += pow(x[i],2);
		}
		s = sqrt(sum);
		r = sqrt(2*s*(s+fabs(x[k+1])));
		for (i=0; i<=k; i++){
			w[i]=0;
		}
		if (x[k+1] < 0){
			w[k+1] = -1*(fabs(x[k+1])+s);
		}
		else if (x[k+1] > 0){
			w[k+1] = (fabs(x[k+1])+s);
		}
		for (i=k+2; i<n; i++){
			w[i] = x[i];
		}
		for (i=0; i<n; i++){
			w[i] = (1/r)*w[i];
		}
		for (i=0; i<n; i++){
			for (j=0; j<n; j++){
				if (i==j){
					H[i][i] = 1-2*w[i]*w[j];
				}
				else {
					H[i][j] = -2*w[i]*w[j];
				}
			}
		}
		MatrixProdukt(Q0,H,Hk,n,n);
		KopierMatrix(Hk,Q0,n,n);
		for (i=0; i<n; i++){
			sum = 0;
			for (j=0; j<n; j++){
				sum += Ak[i][j]*w[j];
			}
			v[i] = sum;
		}
		for (i=0; i<n; i++){
			c +=v[i]*w[i];
		}
		for (i=0; i<n; i++){
			q[i] = v[i]-c*w[i];
		}
		for (i=0; i<n; i++){
			for (j=0; j<n; j++){
				Ak[i][j] = Ak[i][j]-2*w[i]*q[j]-(2*q[i]*w[j]);
			}
		}
		KopierMatrix(Ak,An2,n,n);
	}
	KopierMatrix(Ak,An2,n,n);
	Danaogb(An2,a,b,n);
	cout << "\nMatrix Zn2: \n";
	UdskrivMatrix(An2,n,n);
}

void QRWilkinson(double Q[NMAX][NMAX], double a[NMAX], double b[NMAX-1],
		double ev[NMAX], int n){
	double shift,sum=0, ev1, ev2, a1[NMAX], b1[NMAX-1], c[NMAX-1],
			s[NMAX-1];
	int k=0, i, m;
	for(i=0; i<n; i++){
			a1[i] = a[i];       // Kopi af hoveddiagonal
            b1[i] = b[i];       // Kopi af bidiagonal
	}
	cout << endl << "Startvaerdier for diagonalen: \n";
	UdskrivVektor(a1,n);
	cout << endl << "Startvaerdier for bidiagonalen: \n";
	UdskrivVektor(b1,n-1);
	Enhedsmatrix(Q,n);
	for(m=n-1; m>=0; m--){
		do {
			// Beregner egenvaerdier for M(2x2):
			ev1 = (a1[m]+a1[m-1])/2+pow((a1[m]+a1[m-1])*(a1[m]+a1[m-1])/4-
					((a1[m]*a1[m- 1])-(b1[m-1]*b1[m-1])),0.5);
			ev2 = (a1[m]+a1[m-1])/2-pow((a1[m]+a1[m-1])*(a1[m]+a1[m-1])/4-
					((a1[m]*a1[m- 1])-(b1[m-1]*b1[m-1])),0.5);
            // Shifts defineres
			if(fabs(a1[m]-ev1) <= fabs(a1[m]-ev2)){
				shift = ev1;
			}
			else{
				shift = ev2;
			}
			// Shifts akkumuleres
			sum += shift;
			// Shifts traekkes fra i diagonalen
			for(int i=0; i<n; i++){
				a1[i] -= shift;
			}
			// QR-algoritmen med Givens rotationer udfoeres
            QR(Q,a1,b1,c,s,m+1);
            k++;
		}
		while(fabs(b1[m-1]) >= eps);
		ev[m] = a1[m]+sum;
	}
	cout << "\nAntal iterationer: " << k << endl;
	cout << "\nEgenvaerdierne er: " << endl;
	UdskrivVektor(ev,n);
}

void Danaogb(double A[NMAX][NMAX], double a[NMAX], double b[NMAX-1],
		int n){
	int i;
	for (i=0; i<n; i++){
		a[i] = A[i][i];
	}
	for (i=0; i<n-1; i++){
		b[i] = A[i][i+1];
	}
}

// Opgave 5

void Ombyt(double &x, double &y){
	double temp;
	temp = x;
	x = y;
	y = temp;
}

void Sorter(double LamVals[NMAX], double WMatUS[NMAX][NMAX], double
		LamValsSo[NMAX], double VMatSo[NMAX][NMAX], int n, int m){
	int i, j, k;
	KopierVektor(LamVals,LamValsSo,n);
	KopierMatrix(WMatUS,VMatSo,n,m);
	for (i=n-1; i>=1; i--){
		for (j=0; j<=i-1; j++){
			if (fabs(LamValsSo[j]) < fabs(LamValsSo[j+1])){
				Ombyt(LamValsSo[j],LamValsSo[j+1]);
				for (k=0; k<n; k++){
					Ombyt(VMatSo[k][j],VMatSo[k][j+1]);
				}
			}
		}
	}
}

int Rang(double LamVals[NMAX], int n){
	int rang=0, i;
	for (i=0; i<n; i++){
		if (LamVals[i] > eps){
			rang += 1;
		}
	}
	return rang;
}

void PseudoInvers(double U1[NMAX][NMAX], double D1M[NMAX][NMAX], double
		V1[NMAX][NMAX], double Aplus[NMAX][NMAX], int n, int r){
	double D1inv[NMAX][NMAX], U1T[NMAX][NMAX], V1D1inv[NMAX][NMAX];
	int i, j;
	for (i=0; i<r; i++){    // D1inv
		for (j=0; j<r; j++){
			if (i==j){
				D1inv[i][i] = 1/D1M[i][i];
			}
			else {
				D1inv[i][j] = 0;
			}
		}
	}
	TransponerMatrix(U1,U1T,r,n);   // U1T
	MatrixProdukt(V1,D1inv,V1D1inv,n,r);    // V1*D1
	MatrixProdukt(V1D1inv,U1T,Aplus,n,r);    // Aplus=A^+
}

void SVD(double Z[NMAX][NMAX], double LamValsSo[NMAX], double
		VMatSo[NMAX][NMAX], double U1[NMAX][NMAX], double D1M[NMAX][NMAX],
		double V1[NMAX][NMAX], int n, int r){
	double D1[NMAX], D1inv[NMAX], V1D1inv[NMAX][NMAX], D1invM[NMAX][NMAX];
	int i, j;
	for (i=0; i<n; i++){    // V1
		for (j=0; j<r; j++){
			V1[i][j] = VMatSo[i][j];
		}
	}
	for (i=0; i<r; i++){    // D1
		D1[i] = sqrt(fabs(LamValsSo[i]));
	}
	for (i=0; i<r; i++){    // D1inv
		D1inv[i] = 1/D1[i];
	}
	// U1
	SamlMatrix(D1M,D1,n);
	SamlMatrix(D1invM,D1inv,n);
	MatrixProdukt(V1,D1invM,V1D1inv,n,r);
	MatrixProdukt(Z,V1D1inv,U1,n,r);
}

// Hjaelpefunktioner -------------------------------------------------------

void IndhentMatrix(double A[NMAX][NMAX], int &n, int &m, string filnavn){
	ifstream Fil;
	if (filnavn == "x" or filnavn == "x"){
		cout << "\nIndtast filnavn paa den fil, der oenskes indhentet: ";
		cin >> filnavn;
		filnavn = filnavn + ".txt";
	}
	do {
		if (FilEksisterende(filnavn) == false){
			cout << "\nDenne fil findes ikke. Indtast venligst et andet navn."
					" \n";
			cout << "\nIndtast filnavn paa den fil, der oenskes indhentet: ";
			cin >> filnavn;
			filnavn = filnavn + ".txt";
		}
	}
	while(FilEksisterende(filnavn) == false);
	Fil.open(filnavn);
	Fil >> n;
	Fil >> m;
	for (int i=0;i<n;i++){
		for (int j=0; j<m; j++) Fil >> A[i][j];
	}
	Fil.close();
}

void IndhentTotalmatrix(double A[NMAX][NMAX], double y[NMAX], int &n,
		int &m, string filnavn){
	ifstream Fil;
	if (filnavn == "x" or filnavn == "X"){
		cout << "\nIndtast filnavn paa den fil, der oenskes indhentet: ";
		cin >> filnavn;
		filnavn = filnavn + ".txt";
	}
	do {
		if (FilEksisterende(filnavn) == false){
			cout << "\nDenne fil findes ikke. Indtast venligst et andet navn."
					" \n";
			cout << "\nIndtast filnavn paa den fil, der oenskes indhentet: ";
			cin >> filnavn;
			filnavn = filnavn + ".txt";
		}
	}
	while(FilEksisterende(filnavn) == false);
	Fil.open(filnavn);
	Fil >> n;
	Fil >> m;
	for (int i=0;i<n;i++){
		for (int j=0; j<m; j++) Fil >> A[i][j];
		Fil>>y[i];
	}
	Fil.close();
}

void IndtastMatrix(double A[NMAX][NMAX], int &n, int &m){
	cout << "Indtast antal raekker, n: "; cin >> n;
	cout << "Indtast antal soejler, m: "; cin >> m;
	cout << "Indtast matricen A: \n";
		for(int i=0; i<n; i++){
			for(int j=0; j<m; j++){
				cout << "Indtast plads [" << i+1 << "][" << j+1 << "] her: ";
				cin >> A[i][j];
			}
		}
}

void IndtastVektor(double v[NMAX], int n){
	for(int k=0; k<n; k++){
		cout << "Indtast element nr. " << k+1 << " her: ";
		cin >> v[k];
	}
}

void IndtastTotalmatrix(double A[NMAX][NMAX], double b[NMAX], int &n, int &m){
	cout << "Indtast antal raekker, n: "; cin >> n;
	cout << "Indtast antal soejler, m: "; cin >> m;
	cout << "Foerst indtastes matricen A: \n";
		for(int i=0; i<n; i++){
			for(int j=0; j<m; j++){
				cout << "Indtast plads [" << i+1 << "][" << j+1 << "] her: ";
				cin >> A[i][j];
			}
		}
		cout << "Nu indtastes vektoren b: \n";
		for(int i=0; i<n; i++){
			cout << "Indtast plads [" << i+1 << "] her: ";
			cin >> b[i];
		}
}

void DanTotalmatrix(double TM[NMAX][NMAX+1],double A[NMAX][NMAX], double
		b[NMAX], int n, int m){
	for(int i=0; i<n; i++){
		for(int j=0; j<m; j++){
			TM[i][j] = A[i][j];
		}
	}
	for(int i=0; i<n; i++){
		TM[i][m] = b[i];
	}
}

void UdskrivMatrix(double M[NMAX][NMAX], int n, int m){
	int i, j;
	for(i=0; i<n; i++){
		for(j=0; j<m; j++){
			if (fabs(M[i][j]) < eps){
				cout << setw(20) << setprecision(14) << 0.0;
			}
			else {
				cout << setw(20) << setprecision(14) << M[i][j];
			}
		}
	cout << endl;
	}
}

void UdskrivTotalmatrix(double TM[NMAX][NMAX+1], int n, int m){
	int i, j;
	for(i=0; i<n; i++){
		for(j=0; j<m+1; j++){
			if (fabs(TM[i][j]) < eps){
				cout << setw(20) << setprecision(14) << 0.0;
			}
			else {
				cout << setw(20) << setprecision(14) << TM[i][j];
			}
		}
	cout << endl;
	}
}

void UdskrivVektor(double v[NMAX], int n){
	int i;
	for (i=0; i<n; i++){
		if (fabs(v[i]) < eps){
			cout << setw(20) << setprecision(14) << 0.0 << endl;
		}
		else {
			cout << setw(20) << setprecision(14) << v[i] << endl;
		}
	}
}

void GemMatrix(double A[NMAX][NMAX], int &n, int &m, string stdnavn){
	int svarNavn;
	ofstream Fil;
	string filnavn;
    auto now = std::chrono::system_clock::now();    // Tilpasset fra:
                                                    // https://stackoverflow.com/questions/
                                                    // 17223096/outputting-date-and-time-in
                                                    // c-using-stdchrono
	std::time_t now_time = std::chrono::system_clock::to_time_t(now);
	int svar;
	if (stdnavn != "x" or filnavn != "x"){
		filnavn = stdnavn;
	}
	cout << "\nVil du: \n1) Vaelge standardnavn " << stdnavn+".txt" << "\n2) "
			"Autogenere navn \n3) Indtaste navn manuelt \n";
	do {
		cout << "\nIndtast 1, 2 eller 3: ";
		cin >> svar;
	}
	while (svar != 1 && svar != 2 && svar != 3);
	switch (svar){
		case 1:
			cout << "\nStandardnavn " << stdnavn << " er valgt." << endl;
			break;
		case 2:
			cout << "\nNavn autogeneres." << endl;
			filnavn = std::ctime(&now_time);
			break;
		case 3:
			cout << "\nNavn indtastes manuelt." << endl;
			cout << "\nIndtast oensket navn for fil: ";
			cin >> filnavn;
			break;
	}
	filnavn = filnavn + ".txt";
	do {
		if (FilEksisterende(filnavn) == true){
			cout << "\nDer goeres opmaerksom paa, at dette navn allerede "
					"eksisterer og derfor vil overskrive den nuvaerende fil.";
			cout << "\nVil du: \n1) Bevare navn og overskrive nuvaerende fil "
					"\n2) Indtaste nyt navn \n";
			do {
				cout << "\nIndtast 1 eller 2: ";
				cin >> svarNavn;
			}
			while (svarNavn != 1 && svarNavn != 2);
			switch (svarNavn){
				case 1:
					cout << "\nNavn bevares og overskriver nuvaerende navn.\n";
					break;
				case 2:
					cout << "\nIndtast oensket navn for fil: ";
					cin >> filnavn;
				}
			}
		}
	while (FilEksisterende(filnavn) == true);
	Fil.open(filnavn);
	Fil << n << endl;
	Fil << m << endl;
	for (int i=0; i<n; i++){
		for (int j=0;j<m;j++) Fil << A[i][j] << " ";
		Fil << endl;
	}
	cout << "\nFilen er nu gemt med navnet: " << filnavn << endl;
	Fil.close();
}

void GemTotalmatrix(double TM[NMAX][NMAX+1], int &n, int &m, string stdnavn){
	int svarNavn;
	ofstream Fil;
	string filnavn;
    auto now = std::chrono::system_clock::now();    // Tilpasset fra:
                                                    // https://stackoverflow.com/questions/
                                                    // 17223096/outputting-date-and-time-in
                                                    // c-using-stdchrono
	std::time_t now_time = std::chrono::system_clock::to_time_t(now);
	int svar;
	if (stdnavn != "x" or filnavn != "x"){
		filnavn = stdnavn;
	}
	cout << "\nVil du: \n1) Vaelge standardnavn " << stdnavn+".txt" << "\n2) "
			"Autogenere navn \n3) Indtaste navn manuelt \n";
	do {
		cout << "\nIndtast 1, 2 eller 3: ";
		cin >> svar;
	}
	while (svar != 1 && svar != 2 && svar != 3);
	switch (svar){
		case 1:
			cout << "\nStandardnavn " << stdnavn << " er valgt." << endl;
			break;
		case 2:
			cout << "\nNavn autogeneres." << endl;
			filnavn = std::ctime(&now_time);
			break;
		case 3:
			cout << "\nNavn indtastes manuelt." << endl;
			cout << "\nIndtast oensket navn for fil: ";
			cin >> filnavn;
			break;
	}
	filnavn = filnavn + ".txt";
	do {
		if (FilEksisterende(filnavn) == true){
			cout << "\nDer goeres opmaerksom paa, at dette navn allerede "
					"eksisterer og derfor vil overskrive den nuvaerende fil.";
			cout << "\nVil du: \n1) Bevare navn og overskrive nuvaerende fil "
					"\n2) Indtaste nyt navn \n";
			do {
				cout << "\nIndtast 1 eller 2: ";
				cin >> svarNavn;
			}
			while (svarNavn != 1 && svarNavn != 2);
			switch (svarNavn){
				case 1:
					cout << "\nNavn bevares og overskriver nuvaerende navn.\n";
					break;
				case 2:
					cout << "\nIndtast oensket navn for fil: ";
					cin >> filnavn;
				}
			}
		}
	while (FilEksisterende(filnavn) == true);
	Fil.open(filnavn);
	Fil << n << endl;
	Fil << m << endl;
	for (int i=0; i<n; i++){
		for (int j=0;j<m+1;j++) Fil << TM[i][j] << " ";
		Fil << endl;
	}
	cout << "\nFilen er nu gemt med navnet: " << filnavn << endl;
	Fil.close();
}

bool FilEksisterende(string filnavn)
{
	ifstream file(filnavn);
	if(!file)
		return false;    // Filen var ikke fundet
	else
		return true;     // Filen var fundet
}

void Enhedsmatrix(double Q[NMAX][NMAX], int n){
	int i, j;
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			Q[i][j] = 0;
		}
	}
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			Q[i][i] = 1;
		}
	}
}

void Nulmatrix(double Q[NMAX][NMAX], int n, int m){
	int i, j;
	for(i=0; i<n; i++){
		for(j=0; j<m; j++){
			Q[i][j] = 0;
		}
	}
}

void KopierVektor(double v[NMAX], double w[NMAX], int n){
	int i;
	for (i=0; i<n; i++){
		w[i] = v[i];
	}
}

void KopierMatrix(double A[NMAX][NMAX], double B[NMAX][NMAX], int n, int m){
	int i, j;
	for(i=0; i<n; i++){
		for(j=0; j<m; j++){
			B[i][j] = A[i][j];
		}
	}
}

void NormerVektorer(double v[NMAX][NMAX], int n, int m){
	int i, j;
	double temp[NMAX];
	for(j=0 ; j<n ; j++){
		for(i=0 ; i<n ; i++) {
			temp[j]+=pow(v[i][j],2);
		}
		temp[j]=sqrt(temp[j]);
	}
	for(j=0 ; j<n ; j++) {
		for(i=0 ; i<n ; i++){
			v[i][j]=(1/temp[j]*v[i][j]);
	  }
	}
}

void TransponerMatrix(double M[NMAX][NMAX], double MT[NMAX][NMAX], int n,
		int m){
	for (int i=0; i<n; i++){
		for (int j=0; j<m; j++){
			MT[j][i] = M[i][j];
			if (fabs(MT[j][i]) < eps){
				MT[j][i] = 0.0;
			}
		}
	}
}

void MatrixProdukt(double M1[NMAX][NMAX], double M2[NMAX][NMAX],
		double M3[NMAX][NMAX], int n, int m){
	int i, j, k;
	int q = n;
	if (n<m){
		q = m;
	}
	for(i=0; i<n; i++){
		for(j=0; j<m; j++){
			M3[i][j] = 0;
			for(k=0; k<q; k++) {
				M3[i][j] += M1[i][k]*M2[k][j];
			}
		}
	}
}

void MatrixVektorProdukt(double M[NMAX][NMAX], double v[NMAX], double
		vprod[NMAX], int n, int m){
	int i, j;
	for (i=0; i<n; i++){
		vprod[i] = 0;
	}
    for (i=0; i<n; i++){
    	for (j=0; j<m; j++){
    		vprod[i] += M[i][j]*v[j];
    	}
    }
}

void VektorAddition(double v1[NMAX], double v2[NMAX], double res[NMAX],
		int n){
	for(int i=0; i<n; i++){
		res[i] = v1[i]+v2[i];
	}
}

void VektorProdukt(double v1[NMAX], double v2[NMAX], double res[NMAX],
		int n){
	for(int i=0; i<n; i++){
		res[i] = v1[i]*v2[i];
	}
}
