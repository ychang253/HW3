#pragma once
class Grid
{
public:
	Grid();
	Grid(int m);
	~Grid();

	void SolveByJacobi(int iLoopMax, double dRelax); // HW1
	void SolveByMultiGrid(int iLoopMax);             // HW2
    void OutputResultAsVTK();                        // HW1
    
	void SolveTwoDWave();
	void TwoDWaveStep1();
	void TwoDWaveStep2();
	void UpateTime();

	//Spectral Method
	int m_iChebyN;
	int m_iChebyInner2DGrids;
	double *m_pUInner, *m_pLaU, *m_pLaULaU;
	double *m_pUtInner, *m_pLaUt;
	double **m_pDNij, **m_pDN2ij;           //0~N , dimension N+1
	double **m_pDN2ij_Tu, **m_pIij_Tu;          //0~N-2, dimension N-1   
	double **m_pLNij;                           // 0~ (N-1)^2-1,  dimension (N-1)^2   
	double ChebyGridX(int j);
	void NewSpectralMemory();
	void SetChebyshevDifferMatrix();
	void SetChebyshevDNSquare();
	void SetChebyLaplacian();
	void UpdateUInnerForSpectral();
	void SolveTwoDWavebySpectral();
	void ComputeChebyLaplacianU();
	void ComputeChebyLaplacianUt();
	void TwoDWaveStep1bySpectral();
	void TwoDWaveStep2bySpectral();
	void OutputResultAsVTKSpectral();

	int m_iMeshM;
	int m_iNtotal;
	int m_iT;
	int m_iIterations;
	int m_iJacobiLoopMax;
	double m_dH;
	double m_dMaxError;
	double m_dt;

	double **m_pUij, **m_pU0ij , **m_pRij;
	double **m_pU00ij;
};

