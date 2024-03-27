// AF-ICP: Proposed algorithm in five iterations
//-----------------------------------------------------------//
#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>
#include	<float.h>
#include	<string.h>
#include	<unistd.h>
#include	<assert.h>
#include	<omp.h>

//-- Parallel processing --//
int	NumThreads;//Number of threads for use

//-- Ordinary part --//
int	N;//Number of objects (documents)
int	D;//Dimensionlaity: Number of terms (Full)
int	K;//Number of clusters (centroids, means)
int	*vNT;//Number of distinct terms for each document
double	SumNT;
int	**mTID;//Term ID for each 
double	**mX;//td-idf for each document
int	**mClstMember;//Cluster member ID
int	*vClstSize;//Cluster size: Number of cluster members
int	*vAssID;//ID of cluster to which j-th object is assigned
int	*vInOut;//Boolean flag for object changing its cluster: each iter.
int	*vNew;//Boolean flag whether cluster members change or not
struct	chobj{int out; int in;};
struct	chobj *vChObj;//Object ID that changes its assigned cluster
struct	optParams{int thterm; double thval;};

//-- Sparse expression --//
struct	tidDF{int tid; int df;};
struct	tidDF *vSortTID;//Sorted term ID in ascending order of df
int	NumUTD;//Number of used terms with a dense format

//-- Inverted file --//
int	*vNM;//Number of means (mf) that use j-th term with dense format
int	**mMID;//ID of h-th mean in means that use j-th term
double	**mMean;//Centorid positions
double	**mSim;//Distance to all centroids from j-th object
int	**mAppear;//Number of appeared centroids
int	*vKcount;//Counter for mClstMember
int	*vTcount;//Counter for term

double	*vTFIDF;//TFIDF vector for calculating a mean
double	*vAssSim;//Similarity of each object to its assigned cluster's mean
int	*vDense2Full;

int	ThTerm;
int	*vOrderTh;//Order at which each object's tid exceeds ThTerm.
int	PartTerms;//NumUTD-ThTermMin
double	**mMeanVal;

double	**mObjPL1norm;
double	*vObjPL1norm;
double	*vObjL1norm;
int	*vNMlarge, *vNMmove;
int	*vObjStatus;//Boolean flag for invariant or moving for each object
int	**mExactCal;
int	*vND;//Number of objects (df) that use j-th term with dense format

double	DFMF;//sum_tid=1^D (doc_freq*mean_freq)=Mult
int	ThTermMin;
double	ThVal, ThValMin, ThValMax, ThValStep1, ThValStep2;
double *vAvgObj;
double	*vAvgMean, *vAvgSim;

int	NumMvMeans;
int	*vMvMID;


void	readData_sparse(char *fn)
{//Original data --> dense format (termID) --> Preparing sorted-df structure
	FILE	*fp;
	int	i, j;
	int	id, tid, dense_tid;//tid: Full_termID
	double	v;//Current v: tf, but after that, tf-idf
	int	dummy, dummy0, dummy1;
	int	*vFull2Dense;//term ID transformed from that in full to dense expression
	int	num_terms;

	if((fp = fopen(fn, "r")) == NULL){
		printf("Unknown File = %s\n", fn); exit(1);}

	SumNT = 0.0;
	dummy = fscanf(fp, "%d %d %d", &N, &D, &dummy0);
	vFull2Dense = (int *) malloc(sizeof(int)*D);
	for(i = 0; i < D; ++i) vFull2Dense[i] = -1;

	vNT = (int *) malloc(sizeof(int)*N);
	mTID = (int **) malloc(sizeof(int *)*N);
	mX = (double **) malloc(sizeof(double *)*N);
	for(i = 0, dense_tid = 0; i < N; i++){  
		dummy = fscanf(fp, "%d", &vNT[i]);
		num_terms = vNT[i];
		mTID[i] = (int *) malloc(sizeof(int)*num_terms);
		mX[i] = (double *) malloc(sizeof(double)*num_terms);
		SumNT += (double) num_terms;
		for(j = 0; j < num_terms; j++){
			dummy = fscanf(fp, "%d:%lf", &id, &v);
			tid = id-1;
			mX[i][j] = v;
			if(vFull2Dense[tid] == -1){
				mTID[i][j] = dense_tid; 
				vFull2Dense[tid] = dense_tid++;
			}else{
				mTID[i][j] = vFull2Dense[tid];
			}
		}
		dummy = fscanf(fp, "%d %d", &dummy0, &dummy1);
	}
	if(dummy < 0) printf("ERROR: Somethong strange\n");
	fclose(fp);

	vDense2Full = (int *) malloc(sizeof(int)*dense_tid);
	for(i = 0, j = 0; i < D; ++i){
		if(vFull2Dense[i]	== -1) continue;
		vDense2Full[vFull2Dense[i]] = i; j++;
	}
	if(j != dense_tid){
		printf("ERROR: Mismatch of Full2Dense and Dense2Full.\n"); exit(1);}

	NumUTD = dense_tid;
	vSortTID = (struct tidDF *) malloc(sizeof(struct tidDF)*NumUTD);
	for(i = 0, j = 0; i < D; ++i){
		if((dense_tid = vFull2Dense[i]) == -1) continue;
		vSortTID[dense_tid].tid = dense_tid;
		vSortTID[dense_tid].df = -1; j++;//Initialization
	}
	free(vFull2Dense);

	printf("N=%d D=%d NumUTD=%d K=%d ThTermMin=%d\n",N,D,NumUTD,K,ThTermMin);
	fprintf(stderr,"N=%d D=%d NumUTD=%d K=%d ThTermMin=%d\n",N,D,NumUTD,K,ThTermMin);
}

//-- QuickSort: vSortTID --//
void	swapTID(int p, int q)
{
	struct tidDF	tmp;
	tmp = vSortTID[p]; vSortTID[p] = vSortTID[q]; vSortTID[q] = tmp;
}
int	SortTIDpartition(int low, int high)
{
	int i;
	int pvt, trgt;//pivot and target

	pvt = high;
	trgt = low;
	for(i = low; i < high; ++i){
		if(vSortTID[i].df < vSortTID[pvt].df){
			swapTID(i,trgt); trgt++;
		}else if(vSortTID[i].df == vSortTID[pvt].df
				&& vSortTID[i].tid < vSortTID[pvt].tid){
			swapTID(i,trgt); trgt++;
		}
	}
	swapTID(pvt,trgt);
	return trgt;
}
void	quickSortTID(int low, int high)
{
	int pvt;
	if(high > low){
		pvt = SortTIDpartition(low,high);
		quickSortTID(low,pvt-1);
		quickSortTID(pvt+1,high);
	}
}

void	swapObjElmt(int i, int p, int q)//Swap elements of object
{
	int	tmp_tid;
	double	tmp_val;
	tmp_tid = mTID[i][p]; mTID[i][p] = mTID[i][q]; mTID[i][q] = tmp_tid;
	tmp_val = mX[i][p]; mX[i][p] = mX[i][q]; mX[i][q] = tmp_val;
}
int	ObjElmtpartition(int i, int low, int high)
{
	int j;
	int pvt, trgt;//pivot and target

	pvt = high;
	trgt = low;
	for(j = low; j < high; ++j){
		if(mTID[i][j] < mTID[i][pvt]){
			swapObjElmt(i,j,trgt); trgt++;
		}
	}
	swapObjElmt(i,pvt,trgt);
	return trgt;
}
void	qsObjElmt(int i, int low, int high)
{
	int pvt;
	if(high > low){
		pvt = ObjElmtpartition(i,low,high);
		qsObjElmt(i,low,pvt-1);
		qsObjElmt(i,pvt+1,high);
	}
}

void	swapMean(int *mMID, double *mMean, int p, int q)
{
	int	tmpMID;
	double	tmpMean;
	tmpMID = mMID[p]; mMID[p] = mMID[q]; mMID[q] = tmpMID;
	tmpMean = mMean[p]; mMean[p] = mMean[q]; mMean[q] = tmpMean;
}
int	SortMeanpartition(int *mMID, double *mMean, int low, int high)
{
	int i;
	int pvt, trgt;//pivot and target

	pvt = high;
	trgt = low;
	for(i = low; i < high; ++i){
		if(mMean[i] > mMean[pvt]){
			swapMean(mMID,mMean,i,trgt); trgt++;
		}
	}
	swapMean(mMID,mMean,pvt,trgt);
	return trgt;
}
void	quickSortMean(int *mMID, double *mMean, int low, int high)
{
	int pvt;
	if(high > low){
		pvt = SortMeanpartition(mMID,mMean,low,high);
		quickSortMean(mMID,mMean,low,pvt-1);
		quickSortMean(mMID,mMean,pvt+1,high);
	}
}
//-- END: QuickSort --//

void	calTFIDF()
{// Making vSortedTID[*].{tid,df}
	int	i, j;
	double	v;
	double	*vIDF;//Inverted document frequency
	int	num_terms; // tmp_tid;
	int	*vInvMap;
	double	invN = (double) 1.0/N;

	vND = (int *) calloc(NumUTD,sizeof(int));
	for(i = 0; i < N; ++i){
		num_terms = vNT[i];
		for(j = 0; j < num_terms; ++j) vND[mTID[i][j]]++;
	}
#pragma omp parallel num_threads(NumThreads) private(i)
{
	#pragma omp for
	for(i = 0; i < NumUTD; ++i) vSortTID[i].df = vND[vSortTID[i].tid];
}//omp
	quickSortTID(0,NumUTD-1);//vSortTID[i] was arranged in ascending order of df.

	vInvMap = (int *) malloc(sizeof(int)*NumUTD);
#pragma omp parallel num_threads(NumThreads) private(i)
{
	#pragma omp for
	for(i = 0; i < NumUTD; ++i){
		vND[i] = vSortTID[i].df;
		vInvMap[vSortTID[i].tid] = i;//Mapping of old TID to new TID
	}
}//omp

	//-- Term ID in each object is sorted in ascending order of ID. --//
#pragma omp parallel num_threads(NumThreads) private(i,num_terms)
{
	#pragma omp for
	for(i = 0; i < N; ++i){
		num_terms = vNT[i];
		for(j = 0; j < num_terms; ++j){
			mTID[i][j] = vInvMap[mTID[i][j]];	
		}
		qsObjElmt(i,0,num_terms-1);
	}
}//omp

	vIDF = (double *) malloc(sizeof(double)*NumUTD);
	v = log(1.0*N);
#pragma omp parallel num_threads(NumThreads) private(i)
{
	#pragma omp for
	for(i = 0; i < NumUTD; i++) vIDF[i] = v;
	#pragma omp for
	for(i = 0; i < NumUTD; i++) vIDF[i] -= log(1.0*vND[i]);
}//omp

	//-- Average object-feature value --//
	vAvgObj = (double *) calloc(NumUTD,sizeof(double));
	for(i = 0; i < N; i++){
		num_terms = vNT[i];
		for(j = 0; j < num_terms; ++j){
			mX[i][j] *= vIDF[mTID[i][j]];
		}
		for(j = 0, v = 0.0; j < num_terms; ++j) v += mX[i][j]*mX[i][j];
		for(j = 0, v = 1.0/sqrt(v); j < num_terms; ++j){
			mX[i][j] *= v;
			vAvgObj[mTID[i][j]] += mX[i][j];
		}
	}
#pragma omp parallel num_threads(NumThreads) private(j)
{
	#pragma omp for
	for(j = 0; j < NumUTD; ++j) vAvgObj[j] *= invN;
}//omp

	free(vIDF); free(vInvMap);
}
void	initData(int seedN, int seedK)
{
	int	i, j;
	int	id, tid;//term_ID
	int	cid;//Cluster ID
	int	order;//ID for choosing init centroids
	int	num_terms;
	int	*vID, *vPnt;//ID and Point
	int	*vCount;

	srand(seedN);
	vAssID = (int *) malloc(sizeof(int)*N);
	vClstSize = (int *) calloc(K,sizeof(int));
	vChObj = (struct chobj *) malloc(sizeof(struct chobj)*N);
	vNM = (int *) calloc(NumUTD,sizeof(int));
	for(i = 0; i != N; ++i){
		cid = rand() % K; 
		vAssID[i] = cid; vClstSize[cid]++;
		vChObj[i].out = -1; vChObj[i].in = -1; 
	}

	//-- Choose initial centroids 
	//-- Means are already normalized regarding their lengths --//
	srand(seedK);
	vID = (int *) malloc(sizeof(int)*N);
	for(i = 0; i != N; ++i) vID[i] = i;
	vPnt = (int *) calloc(K,sizeof(int));
	for(i = 0, j = N; i != K; ++i){
		order = rand() % j;
		id = vID[order]; vPnt[i] = id;
		if(order != j-1){vID[order] = vID[j-1];}
		j--;
	}//K-loop

	vNM = (int *) calloc(NumUTD,sizeof(int));
	for(i = 0; i != K; ++i){
		id = vPnt[i];
		num_terms = vNT[id];
		for(j = 0; j < num_terms; ++j){
			tid = mTID[id][j]; vNM[tid]++;
		}
	}

	mMID = (int **) malloc(sizeof(int *)*NumUTD);
	mMean = (double **) malloc(sizeof(double *)*NumUTD);
	for(i = 0; i != NumUTD; ++i){
		mMID[i] = (int *) malloc(sizeof(int)*vNM[i]);
		mMean[i] = (double *) malloc(sizeof(double)*vNM[i]);
	}

	vCount = (int *) calloc(NumUTD,sizeof(int));
	for(i = 0; i != K; ++i){
		id = vPnt[i];
		num_terms = vNT[id];
		for(j = 0; j < num_terms; ++j){
			tid = mTID[id][j];
			mMID[tid][vCount[tid]] = i;
			mMean[tid][vCount[tid]++] = mX[id][j];
		}
	}
	free(vPnt); free(vID); free(vCount);

	mSim = (double **) malloc(sizeof(double *)*NumThreads);
	mAppear = (int **) malloc(sizeof(int *)*NumThreads);
	mExactCal = (int **) malloc(sizeof(int *)*NumThreads);
	mObjPL1norm = (double **) malloc(sizeof(double *)*NumThreads);
	for(i = 0; i != NumThreads; ++i){
		mSim[i] = (double *) malloc(sizeof(double)*K);
		mAppear[i] = (int *) malloc(sizeof(int)*K);
		mExactCal[i] = (int *) malloc(sizeof(int)*K);
		mObjPL1norm[i] = (double *) malloc(sizeof(double)*K);
	}

	vObjL1norm = (double *) calloc(N,sizeof(double));
#pragma omp parallel num_threads(NumThreads) private(i,j)
{
	#pragma omp for
	for(i = 0; i < N; i++){
		for(j = 0; j < vNT[i]; ++j) vObjL1norm[i] += mX[i][j];
	}
}//omp

	mClstMember = (int **) calloc(K,sizeof(int *));
	vKcount = (int *) calloc(K,sizeof(int));
	vTFIDF = (double *) calloc(NumUTD,sizeof(double));
	vTcount = (int *) calloc(NumUTD,sizeof(int));
	vAssSim = (double *) calloc(N,sizeof(double));
	vNMmove = (int *) malloc(sizeof(int)*NumUTD);
	vMvMID = (int *) malloc(sizeof(int)*K);
	vObjStatus = (int *) calloc(N,sizeof(int));

	vInOut = (int *) calloc(K,sizeof(int));
	vNew = (int *) malloc(sizeof(int)*K);
	for(i = 0; i != K; ++i) vNew[i] = 1;
}

	//-- Estimating lower bound on ThVal --//
double estThValLB(double thval_step)
{//Required: **mMean, *vNM, NumUTD, *vAssSim, *vObjL1norm
	int	j, h;
	int	th_id, tid, bin;
//	double	bin_width = 1e-4;
	double	bin_width = thval_step;
	double	scale = 1.0/bin_width;
	int num_bins;
	double	mult, min_mult, min_thv;
	int	min_bin;
	int	df;

	unsigned long	**mMult12, *vMult12;// #mult. in Regions 1,2
	int	*vUnPrune;
	double *vMult3;// #multiplications in Region3

	num_bins = ((int) floor(scale)) +1;
	vMult12 = (unsigned long *) calloc(num_bins,sizeof(unsigned long));
	mMult12 = (unsigned long **) malloc(sizeof(unsigned long *)*NumThreads);
	for(h = 0; h < NumThreads; h++)
		mMult12[h] = (unsigned long *) calloc(num_bins,sizeof(unsigned long));
	vUnPrune = (int *) calloc(num_bins,sizeof(int));
	vMult3 = (double *) malloc(sizeof(double)*num_bins);

	//-- #multiplications in Regions 1 and 2 --//
#pragma omp parallel num_threads(NumThreads) private(th_id,tid,h,bin,df)
{
	th_id = omp_get_thread_num();
	#pragma omp for
	for(tid = 0; tid < NumUTD; tid++){
		df = vND[tid];
		for(h = 0; h < vNM[tid]; ++h){
			bin = (int) floor(scale*mMean[tid][h]);
			mMult12[th_id][bin] += (unsigned long) df;
		}
	}//t-loop
}//omp
	for(j = 0; j < num_bins; j++)
		for(h = 0; h < NumThreads; ++h) vMult12[j] += mMult12[h][j];
	for(h = num_bins-1; h > 0; h--) vMult12[h-1] += vMult12[h];

	//-- #multiplications in Region 3 --//
	for(j = 0; j < N; j++){
		bin = (int) floor(scale*vAssSim[j]/vObjL1norm[j]);
		vUnPrune[bin] += vNT[j];
	}
	for(h = 1; h < num_bins; h++) vUnPrune[h] += vUnPrune[h-1];
	assert(vUnPrune[num_bins-1] == SumNT);

	min_mult = DBL_MAX; min_bin = num_bins;
	for(h = 0; h < num_bins; h++){
		vMult3[h] = (double) K*vUnPrune[h];
		mult = (double) vMult12[h] +vMult3[h];
		if(mult < min_mult){min_mult = mult; min_bin = h;}
	}
	fprintf(stderr,"min_bin=%d val=%.4f min_mult=%e\n",min_bin,min_bin*bin_width,min_mult);

	for(h = 0; h < NumThreads; h++) free(mMult12[h]);
	free(mMult12);
	if(vUnPrune != NULL) free(vUnPrune); 
	if(vMult3 != NULL) free(vMult3);

	return min_thv = (double) min_bin*bin_width;
}


	//-- Approximate #multiplications: Differentital --//
int	findUpperThVal(int i, int j, double v, double* thval)
{
	int	mid;
	if(v <= thval[i]) return(i); else if(v > thval[j]) return(j+1); 
	else if(i+1 == j) return(j);
	mid = (int)((i+j)/2);
	if(v <= thval[mid]) return(findUpperThVal(i+1,mid,v,thval));
	else return(findUpperThVal(mid,j-1,v,thval));
}
struct optParams	approxMult(
	int thterm_min, double thv_max, double thv_min, double thv_step)
{
	int	num_thv = (int) floor(thv_max/thv_step -thv_min/thv_step) +1;
	double	*vThVal;
	struct optParams params;

	int	tid, ptid;// term
	int j, h, s;// itr;
	double	invK = (double) 1.0/K;
	double	Ke = (double) (K)/exp(1.0);
	double	mult3, prev_mult3;// Number of multiplications	
	double	mult;
	int	min_tid;
	double	min_mult, min_thv;
	double	eps = 1e-4;// Magic number
	int	**mMFL;// Mean frequency lower than PreThVal
	double	*vLocalMin_mult;
	int	*vLocalMin_tid;
	int	partTerms = NumUTD -thterm_min;	

	double	avgDelta_ubmean, gamma;
	int	**mNTH;//Number of terms in each object whose termID >= Thterm
	double	**mPowGamma;
	double	**mXpart;
	int	**mDIDpart;;
	int	*vPTIDcount;

	//-- Sort mMean[tid][] --//
#pragma omp parallel num_threads(NumThreads) private(tid)
{
	#pragma omp for
	for(tid = thterm_min; tid < NumUTD; tid++){
		quickSortMean(mMID[tid],mMean[tid],0,vNM[tid]-1);
	}
}//omp
	// Averages of mean features
	vAvgMean = (double *) calloc(NumUTD,sizeof(double));
	vAvgSim = (double *) calloc(N,sizeof(double));
#pragma omp parallel num_threads(NumThreads) private(tid,h)
{
	#pragma omp for
	for(tid = 0; tid < NumUTD; tid++){
		for(h = 0; h < vNM[tid]; ++h) vAvgMean[tid] += mMean[tid][h];
		vAvgMean[tid] *= invK;
	}
}//omp
#pragma omp parallel num_threads(NumThreads) private(j,h)
{
	#pragma omp for
	for(j = 0; j < N; ++j){//vAvgSim[j] <-- initialized by calloc()
		for(h = 0; h < vNT[j]; ++h) vAvgSim[j] += mX[j][h]*vAvgMean[mTID[j][h]];
	}
	#pragma omp barrier
	#pragma omp for
	for(j = 0; j < N; ++j) vAvgSim[j] = vAssSim[j] -vAvgSim[j];
}//omp

	//-- Prepare ThVal_candidates --//
	vThVal = (double *) malloc(sizeof(double)*num_thv);
	for(j = 0; j < num_thv; ++j) vThVal[j] = thv_min +thv_step*j;
	assert(vThVal[num_thv-1] < thv_max+eps && vThVal[num_thv-1] > thv_max-eps);

	//-- Average of mean features lower than ThVal --//
	mMFL = (int **) malloc(sizeof(int *)*num_thv);
	for(j = 0; j < num_thv; ++j) mMFL[j] = (int *) calloc(partTerms,sizeof(int));
#pragma omp parallel num_threads(NumThreads) private(tid,ptid,h,s)
{
	#pragma omp for
	for(tid = thterm_min; tid < NumUTD; tid++){//Not parallel for findUpperThVal()
		ptid = tid -thterm_min;
		for(h = vNM[tid]-1; h >= 0; h--)
			if((s = findUpperThVal(0,num_thv-1,mMean[tid][h],vThVal)) < num_thv)
				mMFL[s][ptid] += 1;
	}
}//omp
#pragma omp parallel num_threads(NumThreads) private(tid,ptid,s)
{
	#pragma omp for
	for(tid = thterm_min; tid < NumUTD; tid++){
		ptid = tid -thterm_min;
		for(s = 1; s < num_thv; ++s) mMFL[s][ptid] += mMFL[s-1][ptid];
	}
}//omp

	//-- # multiplications with a naive TAAT --//
	DFMF = 0.0;
#pragma omp parallel num_threads(NumThreads) private(tid) reduction(+: DFMF)
{
	#pragma omp for
	for(tid = 0; tid < NumUTD; tid++) DFMF += (double) vND[tid] *vNM[tid];
}//omp

	//-- MEM alloc for document-inverted-index
	mPowGamma = (double **) malloc(sizeof(double *)*num_thv);
	for(s = 0; s < num_thv; s++){
		mPowGamma[s] = (double *) malloc(sizeof(double)*N);
		for(j = 0; j < N; j++) mPowGamma[s][j] = 1.0;
	}
	vPTIDcount = (int *) calloc(partTerms,sizeof(int));
	mDIDpart = (int **) malloc(sizeof(int *)*partTerms);
	mXpart = (double **) malloc(sizeof(double *)*partTerms);
	for(tid = thterm_min; tid < NumUTD; tid++){
		ptid = tid -thterm_min;
		mDIDpart[ptid] = (int *) malloc(sizeof(int)*vND[tid]);
		mXpart[ptid] = (double *) malloc(sizeof(double)*vND[tid]);
	}
	//-- Making document inverted-index in partial
	for(j = 0; j < N; j++){
		for(h = vNT[j]-1; h >= 0; h--){
			ptid = mTID[j][h] -thterm_min;
			if(ptid < 0) break;
			mDIDpart[ptid][vPTIDcount[ptid]] = j;
			mXpart[ptid][vPTIDcount[ptid]++] = mX[j][h];
		}
	}

	//-- Estimate #multiplications and determine two parameters --//
	mNTH = (int **) malloc(sizeof(int *)*num_thv);
	for(s = 0; s < num_thv; s++) mNTH[s] = (int *) calloc(N,sizeof(int));
	vLocalMin_mult = (double *) malloc(sizeof(double)*num_thv);
	vLocalMin_tid = (int *) malloc(sizeof(int)*num_thv);
#pragma omp parallel num_threads(NumThreads) private(s,tid,ptid,j,h,avgDelta_ubmean,gamma,mult3,prev_mult3,mult)
{
	#pragma omp for
	for(s = 0; s < num_thv; s++){
		prev_mult3 = 0.0; mult3 = 0.0;
		vLocalMin_mult[s] = DFMF; mult = DFMF;

		for(tid = NumUTD-1; tid >= thterm_min; tid--){
			ptid = tid -thterm_min;

			avgDelta_ubmean = 0.0;
			for(h = vNM[tid] -mMFL[s][ptid]; h < vNM[tid]; h++)
				avgDelta_ubmean += vThVal[s] -mMean[tid][h];
			avgDelta_ubmean += vThVal[s]*(K -vNM[tid]);//mean-feature=0
			avgDelta_ubmean *= invK;
			
			for(h = 0; h < vND[tid]; h++){
				j = mDIDpart[ptid][h];
				mult3 -= ((double) mNTH[s][j])*mPowGamma[s][j];

				gamma = avgDelta_ubmean*mXpart[ptid][h]/vAvgSim[j];//Denominator means Delta.
				//1/K can NOT be merged to mPowGamma[s][j] due to *=
				if((mPowGamma[s][j] *= pow(Ke,gamma)) > Ke){mPowGamma[s][j] = Ke;}//Not use K

				mNTH[s][j]++;
				mult3 += ((double) mNTH[s][j])*mPowGamma[s][j];
			}
			mult += mult3 -prev_mult3 - ((double) vND[tid])*mMFL[s][ptid];

			prev_mult3 = mult3;
			if(mult < vLocalMin_mult[s]){vLocalMin_mult[s] = mult; vLocalMin_tid[s] = tid;}

		}//tid-loop
	}//num_thv

}//omp

	min_mult = DFMF; min_tid = NumUTD-1; min_thv = 1.0;
	for(s = 0; s < num_thv; s++){
		printf("%f\t%d\t%e\n",vThVal[s],vLocalMin_tid[s],vLocalMin_mult[s]);
		fprintf(stderr,"%f\t%d\t%e\n",vThVal[s],vLocalMin_tid[s],vLocalMin_mult[s]);

		if(vLocalMin_mult[s] < min_mult){
			min_mult = vLocalMin_mult[s];
			min_tid = vLocalMin_tid[s]; min_thv = vThVal[s];
		}
	}
	params.thterm = min_tid;
	params.thval = min_thv;

	for(s = 0; s < num_thv; ++s){free(mMFL[s]); free(mPowGamma[s]); free(mNTH[s]);}
	free(mMFL); free(mPowGamma); free(mNTH);
	free(vLocalMin_mult); free(vLocalMin_tid); free(vThVal);

	free(vPTIDcount);
	for(ptid = 0; ptid < partTerms; ptid++){free(mDIDpart[ptid]); free(mXpart[ptid]);}
	free(mDIDpart); free(mXpart);

	return params;
}


int	calMeans_sparse_1st(int i, double all_distcalc)
{
	int	j, h, k, t;
	double	v, w;
	double	obj_funct;//Objective function value
	int	id, tid, ptid;//term_ID
	int	init_cid, final_cid;
	long num_obsim;//Number of multiplications for sim. calc.
	unsigned long num_mult;//Number of multiplications for sim. calc.
	int	num_chobj;//Number of objects whose assignments change
	int	num_unchclst;//Number of clusters whose memeber remains unchanged
	int	ctrM, ctrI;//Tentative counters for moving and invariant means 
	double	ass_sim;//Similarity from an object to its assigned centroid
	double	min_thv;//Lower bound on ThVal
	int	out_cid, in_cid;//Centroid ID of out-going and in-coming
	struct optParams params;
	int	th_id;//Thread ID
	double	dfmf, rate0, rate1, rate2;
	double	start_time, end_time;

	num_chobj = 0; num_mult = 0; obj_funct = 0.0; num_obsim = 0;

	start_time = omp_get_wtime();
	//-- Object assignment --//
#pragma omp parallel num_threads(NumThreads) private(th_id,k,t,h,init_cid,final_cid,ass_sim,tid,v) reduction(+:num_mult,num_chobj,obj_funct,num_obsim)
{
	th_id = omp_get_thread_num();
	#pragma omp for
	for(j = 0; j < N; j++){ 

		init_cid = vAssID[j];//j's assigned cluster (or centroid)
		final_cid = init_cid;

		//-- Similarity calculation --//
		for(k = 0; k < K; k++){mSim[th_id][k] = 0.0; mAppear[th_id][k] = 0;}
		for(t = 0; t < vNT[j]; t++){
			tid = mTID[j][t]; 
			for(h = 0; h < vNM[tid]; ++h){
				k = mMID[tid][h];
				mAppear[th_id][k] = 1;
				mSim[th_id][k] += mX[j][t]*mMean[tid][h];
				num_mult++;
			}
		}

		//-- Assignment --//
		for(k = 0, ass_sim = 0.0; k < K; k++){
			num_obsim += mAppear[th_id][k];
			if(mSim[th_id][k] > ass_sim){
				ass_sim = mSim[th_id][k]; final_cid = k;
			}
		}
		vAssSim[j] = ass_sim;

		if(final_cid != init_cid){
			num_chobj++; 
			vAssID[j] = final_cid;
			vChObj[j].out = init_cid; vChObj[j].in = final_cid;
		}
		obj_funct += ass_sim;
	}//j-loop: Object
}//omp

	for(j = 0; j < N; ++j){
		if(vChObj[j].out == -1) continue;
		out_cid = vChObj[j].out; in_cid = vChObj[j].in;
		vInOut[out_cid] = 1; vInOut[in_cid] = 1; 
		vClstSize[out_cid]--; vClstSize[in_cid]++;
		vChObj[j].out = -1; vChObj[j].in = -1;
	}

	num_unchclst = K;
	for(k = 0; k < K; k++){
		if(vClstSize[k] == 0){
			printf("ERROR: ClstSize[%d]=0\n",k); exit(1);}
		mClstMember[k] = (int *) calloc(vClstSize[k],sizeof(int));
		vKcount[k] = 0;
		num_unchclst -= vInOut[k];
	}
	NumMvMeans = K -num_unchclst;

	for(j = 0; j < N; ++j){
		k = vAssID[j];
		mClstMember[k][vKcount[k]++] = j;
	}

	dfmf = 0.0;
#pragma omp parallel num_threads(NumThreads) private(tid) reduction(+: dfmf)
{
	#pragma omp for
	for(tid = 0; tid < NumUTD; tid++) dfmf += (double) vND[tid]*vNM[tid];
}//omp

	rate0 = (double) num_chobj/N;
	rate1 = 1.0 -(num_obsim/all_distcalc);
	rate2 = (double) (num_mult)/dfmf;
	printf("KM%d %1.10e %d %e %ld %e %lu %e %d ",
		i+1,obj_funct,num_chobj,rate0,num_obsim,rate1,num_mult,rate2,num_unchclst); 
	fprintf(stderr,"KM%d %1.10e %d %e %ld %e %lu %e %d ",
		i+1,obj_funct,num_chobj,rate0,num_obsim,rate1,num_mult,rate2,num_unchclst); 
	
	end_time = omp_get_wtime();
	printf("\n"); fprintf(stderr,"\n");
	printf("Assignment(i=%d) time=%f\n",i+1,end_time-start_time);
	fprintf(stderr,"Assignment(i=%d) time=%f\n",i+1,end_time-start_time);
	
	if(num_chobj == 0)
		return 1;//All assignments remain unchanged.


	start_time = omp_get_wtime();
	//-- Mean update: Sparse expression --//
#pragma omp parallel num_threads(NumThreads)
{
	#pragma omp for
	for(j = 0; j < NumUTD; j++){
		vNM[j] = 0;
		free(mMID[j]); free(mMean[j]);
	}
}//omp

	//-- to avoid a large memory consumption. --//
	//-- NumUTD * (K || NumThreads) is forbidden. --//
	ctrM = 0;
	for(k = 0; k < K; k++) vMvMID[k] = -1;
	for(k = 0; k < K; k++){
		vNew[k] = vInOut[k]; vInOut[k] = 0;
		if(vNew[k] == 1) vMvMID[ctrM++] = k;

		for(j = 0; j < vClstSize[k]; j++){
			id = mClstMember[k][j];
			vObjStatus[id] = vNew[k];
			for(h = 0; h < vNT[id]; h++) vTcount[mTID[id][h]] = 1;
		}
		for(tid = 0; tid < NumUTD; ++tid){
			if(vTcount[tid] == 0) continue;
			vNM[tid]++; vTcount[tid] = 0;
		}
	}//for-K
	assert(ctrM == NumMvMeans);

	for(j = 0; j < NumUTD; ++j){
		mMID[j] = (int *) malloc(sizeof(int)*vNM[j]);
		mMean[j] = (double *) malloc(sizeof(double)*vNM[j]);
	}

	//-- K-loop --//
	for(k = 0; k < K; ++k){
		for(j = 0; j < vClstSize[k]; ++j){
			id = mClstMember[k][j];
			for(h = 0; h < vNT[id]; ++h){
				tid = mTID[id][h];
				vTFIDF[tid] += mX[id][h];
			}
		}//for-j

		w = 1.0/vClstSize[k];
		v = 0.0;
#pragma omp parallel num_threads(NumThreads) private(j) reduction(+:v)
{
		#pragma omp for
		for(j = 0; j < NumUTD; j++){
			if(vTFIDF[j] > 0.0){
				vTFIDF[j] *= w;
				v += vTFIDF[j]*vTFIDF[j];
			}
		}
}//omp
		if(v < DBL_MIN){
			printf("Mean[%d] vector length is too short or negative.\n",k);
			exit(1);
		}
		v = 1.0/sqrt(v);

#pragma omp parallel num_threads(NumThreads) private(j)
{
		#pragma omp for
		for(j = 0; j < NumUTD; j++){
			if(vTFIDF[j] > 0.0){
				vTFIDF[j] *= v;
				mMID[j][vTcount[j]] = k;
				mMean[j][vTcount[j]++] = vTFIDF[j];
			}
		}
}//omp
		
#pragma omp parallel num_threads(NumThreads) private(j,h,id,tid,ass_sim)
{
		#pragma omp for
		for(j = 0; j < vClstSize[k]; ++j){
			id = mClstMember[k][j];
			if(vObjStatus[id] == 0) continue;
			for(h = 0, ass_sim = 0.0; h < vNT[id]; ++h){
				tid = mTID[id][h];
				ass_sim += mX[id][h]*vTFIDF[tid];
			}
			if(ass_sim >= vAssSim[id]) vObjStatus[id] = 0;
			vAssSim[id] = ass_sim;
		}
}//omp

#pragma omp parallel num_threads(NumThreads)
{
		#pragma omp for
		for(j = 0; j < NumUTD; j++) vTFIDF[j] = 0.0;
}//omp

	}//for-k

	min_thv = estThValLB(ThValStep1);
	if(min_thv > 0) ThValMin = min_thv;
	//-- Estimate params of ThTerm and ThVal --//
	params = approxMult(ThTermMin,ThValMax,ThValMin,ThValStep1);
	ThTerm = params.thterm; ThVal = params.thval;
	printf("#-- ThTerm=%d ThVal=%f\n",ThTerm,ThVal);
	fprintf(stderr,"#-- ThTerm=%d ThVal=%f\n",ThTerm,ThVal);


	//-- Rewrite information related to ThTerm --//
	vOrderTh = (int *) malloc(sizeof(int)*N);
	vObjPL1norm = (double *) calloc(N,sizeof(double));

#pragma omp parallel num_threads(NumThreads) private(j,h)
{
	#pragma omp for
	for(j = 0; j < N; ++j){
		vOrderTh[j] = vNT[j];
		for(h = 0; h < vNT[j]; ++h)
			if(mTID[j][h] >= ThTerm){vOrderTh[j] = h; break;}
		for(h = vOrderTh[j]; h < vNT[j]; ++h)
			vObjPL1norm[j] += mX[j][h];
	}
}//omp

	//-- Make partial inverted-file 
	PartTerms = NumUTD -ThTerm;
	vNMlarge = (int *) malloc(sizeof(int)*PartTerms);

	mMeanVal = (double **) malloc(sizeof(double *)*PartTerms);
	for(ptid = 0; ptid < PartTerms; ++ptid) //including ThTerm ID
		mMeanVal[ptid] = (double *) malloc(sizeof(double)*K);

// -- mMean[][] have been sorted in descending order in approxMult().
#pragma omp parallel num_threads(NumThreads) private(tid,ptid,h)
{
	#pragma omp for
	for(tid = ThTerm; tid < NumUTD; tid++){
		ptid = tid -ThTerm;
		for(h = 0; h < vNM[tid]; ++h)
			if(mMean[tid][h] < ThVal) break;
		vNMlarge[ptid] = h;

		for(h = 0; h < K; h++) mMeanVal[ptid][h] = 0.0;
		for(h = vNMlarge[ptid]; h < vNM[tid]; ++h)
			mMeanVal[ptid][mMID[tid][h]] = mMean[tid][h];
	}
}//omp

	//-- Partition mean arrays into two blocks based on vNew[k] --//
#pragma omp parallel num_threads(NumThreads) private(j,h,k,v,tid,ptid,ctrI)
{
	#pragma omp for
	for(tid = ThTerm; tid < NumUTD; tid++){
		ptid = tid -ThTerm;
		ctrI = vNMlarge[ptid] -1;
		for(h = 0; h < vNMlarge[ptid]; ++h){
			if(h > ctrI) break;
			k = mMID[tid][h]; v = mMean[tid][h];
			if(vNew[k] == 0){
				for(j = ctrI; j > h; j--){
					if(vNew[mMID[tid][j]] == 1){
						mMID[tid][h] = mMID[tid][j]; mMean[tid][h] = mMean[tid][j];
						mMID[tid][j] = k; mMean[tid][j] = v;
						break;
					}
				}
				ctrI = j-1;
			}
		}//for-vNMlarge[ptid]
		vNMmove[tid] = ctrI +1; //Number of moving centroids
	}

	#pragma omp for
	for(tid = 0; tid < ThTerm; tid++){
		ctrI = vNM[tid] -1;
		for(h = 0; h < vNM[tid]; ++h){
			if(h > ctrI) break;
			k = mMID[tid][h]; v = mMean[tid][h];
			if(vNew[k] == 0){
				for(j = ctrI; j > h; j--){
					if(vNew[mMID[tid][j]] == 1){
						mMID[tid][h] = mMID[tid][j]; mMean[tid][h] = mMean[tid][j];
						mMID[tid][j] = k; mMean[tid][j] = v;
						break;
					}
				}
				ctrI = j-1;
			}
		}
		vNMmove[tid] = ctrI +1;
	}//for-tid
}//omp

#pragma omp parallel num_threads(NumThreads)
{
	#pragma omp for
	for(k = 0; k < K; k++) free(mClstMember[k]);
	#pragma omp for
	for(j = 0; j < NumUTD; j++) vTcount[j] = 0;
}//omp

	end_time = omp_get_wtime();
	printf("Update(i=%d) time=%f\n",i+1,end_time-start_time);
	fprintf(stderr,"Update(i=%d) time=%f\n",i+1,end_time-start_time);

	return 0;
}

int	calMeans_sparse_2nd(int i, double all_distcalc)
{
	int	j, h, k, t;
	double	v, w;
	double	obj_funct;//Objective function value
	int	id, tid, ptid;//term_ID
	int	init_cid, final_cid;
	long num_obsim;//Number of multiplications for sim. calc.
	unsigned long num_mult;//Number of multiplications for sim. calc.
	int	num_chobj;//Number of objects whose assignments change
	int	num_unchclst;//Number of clusters whose memeber remains unchanged
	double	ass_sim;//Similarity from an object to its assigned centroid
	double	min_thv;
	int	out_cid, in_cid;//Centroid ID of out-going and in-coming
	int	th_id;//Thread ID
	double	invThVal;//Inverse of ThVal to avoid devisions
	int	ctrM, ctrI;//Tentative counters for Moving and Invariant means
	double	dfmf, rate0, rate1, rate2;
	int num_means;
	double	start_time, end_time;
	struct optParams params;

	num_chobj = 0; num_mult = 0; obj_funct = 0.0;
	num_obsim = 0;

	start_time = omp_get_wtime();
	//-- Object assignment --//
#pragma omp parallel num_threads(NumThreads) private(j,th_id,k,t,h,init_cid,final_cid,ass_sim,tid,ptid,v,num_means) reduction(+:num_mult,num_chobj,obj_funct,num_obsim)
{
	th_id = omp_get_thread_num();
	#pragma omp for
	for(j = 0; j < N; j++){// Outermost loop
		init_cid = vAssID[j];//j's assigned cluster (or centroid)
		final_cid = init_cid;

		//-- Similarity calculation --//
		for(k = 0; k < K; k++){
			mSim[th_id][k] = 0.0; mAppear[th_id][k] = 0;
			mExactCal[th_id][k] = -1;
			mObjPL1norm[th_id][k] = vObjPL1norm[j];
		}

		if(vObjStatus[j] == 0){//Invariant means
			for(t = 0; t < vOrderTh[j]; t++){// middle loop
				tid = mTID[j][t];
				for(h = 0; h < vNMmove[tid]; h++){// innermost loop
					k = mMID[tid][h];
					mSim[th_id][k] += mX[j][t]*mMean[tid][h];
					num_mult++;
				}
			}
			for(t = vOrderTh[j]; t < vNT[j]; t++){// middle loop
				tid = mTID[j][t];
				for(h = 0; h < vNMmove[tid]; h++){// innermost loop
					k = mMID[tid][h];
					mSim[th_id][k] += mX[j][t]*mMean[tid][h];
					num_mult++;
					mObjPL1norm[th_id][k] -= mX[j][t];
				}
			}
			num_means = 0;
			for(h = 0; h < NumMvMeans; h++){//Pruning using L1norm in subspace
				k = vMvMID[h];
				v = mSim[th_id][k] +ThVal*mObjPL1norm[th_id][k];
				if(vAssSim[j] > v) continue;
	
				mSim[th_id][num_means] = mSim[th_id][k];
				mExactCal[th_id][num_means++] = k;
			}
		}else{
			for(t = 0; t < vOrderTh[j]; t++){// middle loop
				tid = mTID[j][t];
				for(h = 0; h < vNM[tid]; h++){// innermost loop
					k = mMID[tid][h];
					mSim[th_id][k] += mX[j][t]*mMean[tid][h];
					num_mult++;
				}
			}
			for(t = vOrderTh[j]; t < vNT[j]; t++){// middle loop
				tid = mTID[j][t]; ptid = tid -ThTerm;
				for(h = 0; h < vNMlarge[ptid]; h++){
					k = mMID[tid][h];
					mSim[th_id][k] += mX[j][t]*mMean[tid][h];
					num_mult++;
					mObjPL1norm[th_id][k] -= mX[j][t];
				}
			}
			num_means = 0;
			for(k = 0; k < K; k++){//Pruning using L1norm in subspace
				v = mSim[th_id][k] +ThVal*mObjPL1norm[th_id][k];
				if(vAssSim[j] > v) continue;
	
				mSim[th_id][num_means] = mSim[th_id][k];
				mExactCal[th_id][num_means++] = k;
			}
		}

		for(t = vOrderTh[j]; t < vNT[j]; t++){
			tid = mTID[j][t]; ptid = tid -ThTerm;
			for(h = 0; h < num_means; h++){
				k = mExactCal[th_id][h];
				mSim[th_id][h] += mX[j][t]*mMeanVal[ptid][k];
				mAppear[th_id][k] = 1;
				num_mult++;
			}
		}	

		for(h = 0, ass_sim = vAssSim[j]; h < num_means; h++){
			num_obsim += mAppear[th_id][h];
			if(mSim[th_id][h] > ass_sim){
				ass_sim = mSim[th_id][h];
				final_cid = mExactCal[th_id][h];
			}
		}
		vAssSim[j] = ass_sim;
		//-- END of sim. calc. --//

		if(final_cid != init_cid){
			num_chobj++; 
			vAssID[j] = final_cid;
			vChObj[j].out = init_cid; vChObj[j].in = final_cid;
		}
		obj_funct += ass_sim;

	}//j-loop: Object
}//omp

	for(j = 0; j < N; ++j){
		if(vChObj[j].out == -1) continue;
		out_cid = vChObj[j].out; in_cid = vChObj[j].in;
		vInOut[out_cid] = 1; vInOut[in_cid] = 1; 
		vClstSize[out_cid]--; vClstSize[in_cid]++;
		vChObj[j].out = -1; vChObj[j].in = -1;
	}

	num_unchclst = K;
	for(k = 0; k < K; k++){
		if(vClstSize[k] == 0){
			printf("ERROR: ClstSize[%d]=0\n",k); exit(1);}
		mClstMember[k] = (int *) calloc(vClstSize[k],sizeof(int));
		vKcount[k] = 0;
		num_unchclst -= vInOut[k];
	}
	NumMvMeans = K -num_unchclst;

	for(j = 0; j < N; ++j){
		k = vAssID[j];
		mClstMember[k][vKcount[k]++] = j;
	} 

	dfmf = 0.0;
#pragma omp parallel num_threads(NumThreads) private(tid) reduction(+: dfmf)
{
	#pragma omp for
	for(tid = 0; tid < NumUTD; tid++) dfmf += (double) vND[tid]*vNM[tid];
}//omp

	rate0 = (double) num_chobj/N;
	rate1 = 1.0 -(num_obsim/all_distcalc);
	rate2 = (double) (num_mult)/dfmf;
	printf("KM%d %1.10e %d %e %ld %e %lu %e %d ",
		i+1,obj_funct,num_chobj,rate0,num_obsim,rate1,num_mult,rate2,num_unchclst); 
	fprintf(stderr,"KM%d %1.10e %d %e %ld %e %lu %e %d ",
		i+1,obj_funct,num_chobj,rate0,num_obsim,rate1,num_mult,rate2,num_unchclst); 

	end_time = omp_get_wtime();
	printf("\n"); fprintf(stderr,"\n");
	printf("Assignment(i=%d) time=%f ThVal=%f\n",i+1,end_time-start_time,ThVal);
	fprintf(stderr,"Assignment(i=%d) time=%f ThVal=%f\n",i+1,end_time-start_time,ThVal);

	if(num_chobj == 0)
		return 1;//All assignments remain unchanged.

	//-- Mean update: Sparse expression --//
	start_time = omp_get_wtime();
#pragma omp parallel num_threads(NumThreads)
{
	#pragma omp for
	for(j = 0; j < NumUTD; j++){
		vNM[j] = 0;
		vNMmove[j] = 0;
		free(mMID[j]); free(mMean[j]);
	}
}//omp

	//-- The following loop was NOT parallelized --//
	//-- to avoid a large memory consumption. --//
	//-- NumUTD * (K || NumThreads) is forbidden. --//
	ctrM = 0;
	for(k = 0; k < K; k++) vMvMID[k] = -1;
	for(k = 0; k < K; k++){
		vNew[k] = vInOut[k]; vInOut[k] = 0;
		if(vNew[k] == 1) vMvMID[ctrM++] = k;

		for(j = 0; j < vClstSize[k]; j++){
			id = mClstMember[k][j];
			vObjStatus[id] = vNew[k];
			for(h = 0; h < vNT[id]; h++)
				vTcount[mTID[id][h]] = 1;
		}
		for(j = 0; j < NumUTD; ++j)
			if(vTcount[j] == 1){vNM[j]++; vTcount[j] = 0;}
	}//for-K
	assert(ctrM == NumMvMeans);


	for(j = 0; j < NumUTD; ++j){
		mMID[j] = (int *) malloc(sizeof(int)*vNM[j]);
		mMean[j] = (double *) malloc(sizeof(double)*vNM[j]);
	}
	for(k = 0; k < K; ++k){
		for(j = 0; j < vClstSize[k]; ++j){
			id = mClstMember[k][j];
			for(h = 0; h < vNT[id]; ++h){
				tid = mTID[id][h];
				vTFIDF[tid] += mX[id][h];
			}
		}//for-j

		w = 1.0/vClstSize[k];
		v = 0.0;
#pragma omp parallel num_threads(NumThreads) reduction(+:v)
{
		#pragma omp for
		for(j = 0; j < NumUTD; j++){
			if(vTFIDF[j] > 0.0){
				vTFIDF[j] *= w;
				v += vTFIDF[j]*vTFIDF[j];
			}
		}
}//omp
		if(v < DBL_MIN){
			printf("Mean[%d] vector length is too short or negative.\n",k);
			exit(1);
		}
		v = 1.0/sqrt(v);

#pragma omp parallel num_threads(NumThreads)
{
		#pragma omp for
		for(j = 0; j < NumUTD; j++){
			if(vTFIDF[j] > 0.0){
				vTFIDF[j] *= v;
				mMID[j][vTcount[j]] = k;
				mMean[j][vTcount[j]++] = vTFIDF[j];
			}
		}
}//omp

#pragma omp parallel num_threads(NumThreads) private(j,h,id,tid,ass_sim)
{
		#pragma omp for
		for(j = 0; j < vClstSize[k]; ++j){
			id = mClstMember[k][j];
			if(vObjStatus[id] == 0) continue;
			for(h = 0, ass_sim = 0.0; h < vNT[id]; ++h){
				tid = mTID[id][h];
				ass_sim += mX[id][h]*vTFIDF[tid];
			}
			if(ass_sim > vAssSim[id]) vObjStatus[id] = 0;
			vAssSim[id] = ass_sim;
		}
}//omp

#pragma omp parallel num_threads(NumThreads) private(j)
{
		#pragma omp for
		for(j = 0; j < NumUTD; j++) vTFIDF[j] = 0.0;
}//omp

	}//for-k

	for(ptid = 0; ptid < PartTerms; ptid++) free(mMeanVal[ptid]);
	free(mMeanVal); free(vNMlarge);

	min_thv = estThValLB(ThValStep2);
	if(min_thv > 0) ThValMin = min_thv;
	//-- Estimate params of ThTerm and ThVal --//
	params = approxMult(ThTermMin,ThValMax,ThValMin,ThValStep2);
	ThTerm = params.thterm; ThVal = params.thval;
	printf("#-- ThTerm=%d ThVal=%f\n",ThTerm,ThVal);
	fprintf(stderr,"#-- ThTerm=%d ThVal=%f\n",ThTerm,ThVal);

#pragma omp parallel num_threads(NumThreads) private(j,h)
{
	#pragma omp for
	for(j = 0; j < N; ++j) vObjPL1norm[j] = 0.0;
	#pragma omp barrier
	#pragma omp for
	for(j = 0; j < N; ++j){
		vOrderTh[j] = vNT[j];
		for(h = 0; h < vNT[j]; ++h)
			if(mTID[j][h] >= ThTerm){vOrderTh[j] = h; break;}
		for(h = vOrderTh[j]; h < vNT[j]; ++h)
			vObjPL1norm[j] += mX[j][h];
	}
}//omp

	//-- Make partial inverted-file 
	// Reset PartTerms
	PartTerms = NumUTD -ThTerm;
	vNMlarge = (int *) malloc(sizeof(int)*PartTerms);

	mMeanVal = (double **) malloc(sizeof(double *)*PartTerms);
	for(ptid = 0; ptid < PartTerms; ++ptid) //including ThTerm ID
		mMeanVal[ptid] = (double *) malloc(sizeof(double)*K);

//-- Scale object features by ThVal --//
	invThVal = 1.0/ThVal;
#pragma omp parallel num_threads(NumThreads) private(j,t,h)
{
	#pragma omp for
	for(j = 0; j < N; j++){
		vObjPL1norm[j] = 0.0;
		for(t = 0; t < vOrderTh[j]; ++t) mX[j][t] *= ThVal;
		for(t = vOrderTh[j]; t < vNT[j]; ++t){
			mX[j][t] *= ThVal;
			vObjPL1norm[j] += mX[j][t];
		}
	}
	#pragma omp barrier
	#pragma omp for
	for(t = 0; t < NumUTD; ++t)
		for(h = 0; h < vNM[t]; ++h){mMean[t][h] *= invThVal;}
}//omp

// -- mMean[][] have been sorted in descending order in approxMult().
#pragma omp parallel num_threads(NumThreads) private(tid,ptid,h)
{
	#pragma omp for
	for(tid = ThTerm; tid < NumUTD; tid++){
		ptid = tid -ThTerm;
		for(h = 0; h < vNM[tid]; ++h)
			//if(mMean[tid][h] < ThVal) break;
			if(mMean[tid][h] < 1.0) break;
		vNMlarge[ptid] = h;

		for(h = 0; h < K; h++) mMeanVal[ptid][h] = 0.0;
		for(h = vNMlarge[ptid]; h < vNM[tid]; ++h)
			mMeanVal[ptid][mMID[tid][h]] = mMean[tid][h];
	}
}//omp

#pragma omp parallel num_threads(NumThreads) private(j,h,k,v,tid,ptid,ctrI,num_means)
{
//	num_means = 0;
	#pragma omp for
	for(tid = ThTerm; tid < NumUTD; tid++){
		ptid = tid -ThTerm;
		num_means = vNMlarge[ptid];
		ctrI = num_means-1;
		for(h = 0; h < num_means; ++h){
			if(h > ctrI) break;
			k = mMID[tid][h]; v = mMean[tid][h];
			if(vNew[k] == 0){
				for(j = ctrI; j > h; j--){
					if(vNew[mMID[tid][j]] == 1){
						mMID[tid][h] = mMID[tid][j]; mMean[tid][h] = mMean[tid][j];
						mMID[tid][j] = k; mMean[tid][j] = v;
						break;
					}
				}
				ctrI = j-1;
			}
		}//for-num_means
		vNMmove[tid] = ctrI +1; //Number of moving centroids
	}

//	num_means = 0;
	#pragma omp for
	for(tid = 0; tid < ThTerm; tid++){
		num_means = vNM[tid];
		ctrI = num_means-1;
		for(h = 0; h < num_means; ++h){
			if(h > ctrI) break;
			k = mMID[tid][h]; v = mMean[tid][h];
			if(vNew[k] == 0){
				for(j = ctrI; j > h; j--){
					if(vNew[mMID[tid][j]] == 1){
						mMID[tid][h] = mMID[tid][j]; mMean[tid][h] = mMean[tid][j];
						mMID[tid][j] = k; mMean[tid][j] = v;
						break;
					}
				}
				ctrI = j-1;
			}
		}//for-num_means
		vNMmove[tid] = ctrI +1;
	}
}//omp

#pragma omp parallel num_threads(NumThreads)
{
	#pragma omp for
	for(k = 0; k < K; k++) free(mClstMember[k]);
	#pragma omp for
	for(j = 0; j < NumUTD; j++) vTcount[j] = 0;
}//omp

	end_time = omp_get_wtime();
	printf("Update(i=%d) time=%f\n",i+1,end_time-start_time);
	fprintf(stderr,"Update(i=%d) time=%f\n",i+1,end_time-start_time);

	return 0;
}

int	calMeans_sparse(int i, double all_distcalc)
{
	int	j, h, k, t;
	double	v, w;
	double	obj_funct;//Objective function value
	int	id, tid, ptid;//term_ID
	int	init_cid, final_cid;
	long num_obsim;//Number of multiplications for sim. calc.
	unsigned long num_mult;//Number of multiplications for sim. calc.
	int	num_chobj;//Number of objects whose assignments change
	int	num_unchclst;//Number of clusters whose memeber remains unchanged
	double	ass_sim;//Similarity from an object to its assigned centroid
	int	out_cid, in_cid;//Centroid ID of out-going and in-coming
	int	th_id;//Thread ID
	int	ctr_low;//Tentative counters for high and low w.r.t. ThVal
	int	ctrM, ctrI;//Tentative counters for Moving and Invariant means
	double	dfmf, rate0, rate1, rate2;
	int num_means;
	double	start_time, end_time;

	num_chobj = 0; obj_funct = 0.0; num_mult = 0; num_obsim = 0;

	start_time = omp_get_wtime();
	//-- Object assignment --//
#pragma omp parallel num_threads(NumThreads) private(j,th_id,k,t,h,init_cid,final_cid,ass_sim,tid,ptid,v,num_means) reduction(+:num_mult,num_chobj,obj_funct,num_obsim)
{
	th_id = omp_get_thread_num();
	#pragma omp for
	for(j = 0; j < N; j++){// Outermost loop
		init_cid = vAssID[j];//j's assigned cluster (or centroid)
		final_cid = init_cid;
		ass_sim = vAssSim[j];

		//-- Similarity calculation --//
		for(k = 0; k < K; k++){
			mSim[th_id][k] = 0.0; mAppear[th_id][k] = 0;
			mExactCal[th_id][k] = -1;
			mObjPL1norm[th_id][k] = vObjPL1norm[j];
		}

		if(vObjStatus[j] == 0){//Invariant means
			for(t = 0; t < vOrderTh[j]; t++){// middle loop
				tid = mTID[j][t];
				for(h = 0; h < vNMmove[tid]; h++){// innermost loop
					k = mMID[tid][h];
					mSim[th_id][k] += mX[j][t]*mMean[tid][h];
					num_mult++;
				}
			}
			for(t = vOrderTh[j]; t < vNT[j]; t++){// middle loop
				tid = mTID[j][t];
				for(h = 0; h < vNMmove[tid]; h++){// innermost loop
					k = mMID[tid][h];
					mSim[th_id][k] += mX[j][t]*mMean[tid][h];
					num_mult++;
					mObjPL1norm[th_id][k] -= mX[j][t];
				}
			}
			// Smaller than ThVal
			num_means = 0;
			for(h = 0; h < NumMvMeans; h++){//Pruning using L1norm in subspace
				k = vMvMID[h];
				v = mSim[th_id][k] +mObjPL1norm[th_id][k]; // num_mult++;
				if(vAssSim[j] > v) continue;
	
				mSim[th_id][num_means] = mSim[th_id][k];
				mExactCal[th_id][num_means++] = k;
			}
			for(t = vOrderTh[j]; t < vNT[j]; t++){
				tid = mTID[j][t]; ptid = tid -ThTerm;
				for(h = 0; h < num_means; h++){
					k = mExactCal[th_id][h];
					mSim[th_id][h] += mX[j][t]*mMeanVal[ptid][k];
					mAppear[th_id][k] = 1;
					num_mult++;
				}
			}
			for(h = 0; h < num_means; h++){
				num_obsim += mAppear[th_id][h];
				if(mSim[th_id][h] > ass_sim){
					ass_sim = mSim[th_id][h];
					final_cid = mExactCal[th_id][h];
				}
			}
			vAssSim[j] = ass_sim;

		}else{//vObjStatus[j] == 1
			for(t = 0; t < vOrderTh[j]; t++){// middle loop
				tid = mTID[j][t];
				for(h = 0; h < vNM[tid]; h++){// innermost loop
					k = mMID[tid][h];
					mSim[th_id][k] += mX[j][t]*mMean[tid][h];
					num_mult++;
				}
			}
			for(t = vOrderTh[j]; t < vNT[j]; t++){// middle loop
				tid = mTID[j][t]; ptid = tid -ThTerm;
				for(h = 0; h < vNMlarge[ptid]; h++){
					k = mMID[tid][h];
					mSim[th_id][k] += mX[j][t]*mMean[tid][h];
					num_mult++;
					mObjPL1norm[th_id][k] -= mX[j][t];
				}
			}
			// Smaller than ThVal
			num_means = 0;
			for(k = 0; k < K; k++){//Pruning using L1norm in subspace
				v = mSim[th_id][k] +mObjPL1norm[th_id][k]; // num_mult++;
				if(vAssSim[j] > v) continue;
	
				mSim[th_id][num_means] = mSim[th_id][k];
				mExactCal[th_id][num_means++] = k;
			}
			for(t = vOrderTh[j]; t < vNT[j]; t++){
				tid = mTID[j][t]; ptid = tid -ThTerm;
				for(h = 0; h < num_means; h++){
					k = mExactCal[th_id][h];
					mSim[th_id][h] += mX[j][t]*mMeanVal[ptid][k];
					mAppear[th_id][k] = 1;
					num_mult++;
				}
			}	
			for(h = 0; h < num_means; h++){
				num_obsim += mAppear[th_id][h];
				if(mSim[th_id][h] > ass_sim){
					ass_sim = mSim[th_id][h];
					final_cid = mExactCal[th_id][h];
				}
			}
			vAssSim[j] = ass_sim;
		}//if-else
		//-- END of sim. calc. --//

		if(final_cid != init_cid){
			num_chobj++; 
			vAssID[j] = final_cid;
			vChObj[j].out = init_cid; vChObj[j].in = final_cid;
		}
		obj_funct += ass_sim;//Sum of sim. to the confirmed assigned centroid 

	}//j-loop: Object
}//omp

	for(j = 0; j < N; ++j){
		if(vChObj[j].out == -1) continue;
		out_cid = vChObj[j].out; in_cid = vChObj[j].in;
		vInOut[out_cid] = 1; vInOut[in_cid] = 1; 
		vClstSize[out_cid]--; vClstSize[in_cid]++;
		vChObj[j].out = -1; vChObj[j].in = -1;
	}

	num_unchclst = K;
	for(k = 0; k < K; k++){
		if(vClstSize[k] == 0){
			printf("ERROR: ClstSize[%d]=0\n",k); exit(1);}
		mClstMember[k] = (int *) calloc(vClstSize[k],sizeof(int));
		vKcount[k] = 0;
		num_unchclst -= vInOut[k];
	}
	NumMvMeans = K -num_unchclst;

	for(j = 0; j < N; ++j){
		k = vAssID[j];
		mClstMember[k][vKcount[k]++] = j;
	} 

	dfmf = 0.0;
#pragma omp parallel num_threads(NumThreads) private(tid) reduction(+: dfmf)
{
	#pragma omp for
	for(tid = 0; tid < NumUTD; tid++) dfmf += (double) vND[tid]*vNM[tid];
}//omp

	rate0 = (double) num_chobj/N;
	rate1 = 1.0 -(num_obsim/all_distcalc);
	rate2 = (double) (num_mult)/dfmf;
	printf("KM%d %1.10e %d %e %ld %e %lu %e %d ",
		i+1,obj_funct,num_chobj,rate0,num_obsim,rate1,num_mult,rate2,num_unchclst); 
	fprintf(stderr,"KM%d %1.10e %d %e %ld %e %lu %e %d ",
		i+1,obj_funct,num_chobj,rate0,num_obsim,rate1,num_mult,rate2,num_unchclst); 

	end_time = omp_get_wtime();
	printf("\n"); fprintf(stderr,"\n");
	printf("Assignment(i=%d) time=%f ThVal=%f\n",i+1,end_time-start_time,ThVal);
	fprintf(stderr,"Assignment(i=%d) time=%f ThVal=%f\n",i+1,end_time-start_time,ThVal);

	if(num_chobj == 0)
		return 1;//All assignments remain unchanged.

	//-- Mean update: Sparse expression --//
	start_time = omp_get_wtime();
#pragma omp parallel num_threads(NumThreads)
{
	#pragma omp for
	for(j = 0; j < NumUTD; j++){
		vNM[j] = 0;
		vNMmove[j] = 0;
		free(mMID[j]); free(mMean[j]);
	}
}//omp

	//-- The following loop was NOT parallelized --//
	//-- to avoid a large memory consumption. --//
	//-- NumUTD * (K || NumThreads) is forbidden. --//
	ctrM = 0;
	for(k = 0; k < K; k++) vMvMID[k] = -1;
	for(k = 0; k < K; k++){
		vNew[k] = vInOut[k]; vInOut[k] = 0;
		if(vNew[k] == 1) vMvMID[ctrM++] = k;

		for(j = 0; j < vClstSize[k]; j++){
			id = mClstMember[k][j];
			vObjStatus[id] = vNew[k];
			for(h = 0; h < vNT[id]; h++){
				tid = mTID[id][h]; 
				vTcount[tid] = 1;
			}
		}
		for(j = 0; j < NumUTD; ++j)
			if(vTcount[j] == 1){vNM[j]++; vTcount[j] = 0;}
	}//for-K
	assert(ctrM == NumMvMeans);


	for(j = 0; j < NumUTD; ++j){
		mMID[j] = (int *) malloc(sizeof(int)*vNM[j]);
		mMean[j] = (double *) malloc(sizeof(double)*vNM[j]);
	}

	for(k = 0; k < K; ++k){
		for(j = 0; j < vClstSize[k]; ++j){
			id = mClstMember[k][j];
			for(h = 0; h < vNT[id]; ++h){
				tid = mTID[id][h];
				vTFIDF[tid] += mX[id][h];
			}
		}//for-j

		w = 1.0/vClstSize[k];
		v = 0.0;
#pragma omp parallel num_threads(NumThreads) reduction(+:v)
{
		#pragma omp for
		for(j = 0; j < NumUTD; j++){
			if(vTFIDF[j] > 0.0){
				vTFIDF[j] *= w;
				v += vTFIDF[j]*vTFIDF[j];
			}
		}
}//omp
		if(v < DBL_MIN){
			printf("Mean[%d] vector length is too short or negative.\n",k);
			exit(1);
		}
		v = 1.0/(ThVal*sqrt(v));

#pragma omp parallel num_threads(NumThreads)
{
		#pragma omp for
		for(j = 0; j < NumUTD; j++){
			if(vTFIDF[j] > 0.0){
				vTFIDF[j] *= v;
				mMID[j][vTcount[j]] = k;
				mMean[j][vTcount[j]++] = vTFIDF[j];
			}
		}
}//omp

		//-- Object-status judgement --//
#pragma omp parallel num_threads(NumThreads) private(j,h,id,tid,ass_sim)
{
		#pragma omp for
		for(j = 0; j < vClstSize[k]; ++j){
			id = mClstMember[k][j];
			if(vObjStatus[id] == 0) continue;
			for(h = 0, ass_sim = 0.0; h < vNT[id]; ++h){
				tid = mTID[id][h];
				ass_sim += mX[id][h]*vTFIDF[tid];
			}
			if(ass_sim > vAssSim[id]) vObjStatus[id] = 0;
			vAssSim[id] = ass_sim;
		}
}//omp

#pragma omp parallel num_threads(NumThreads) private(j)
{
		#pragma omp for
		for(j = 0; j < NumUTD; j++) vTFIDF[j] = 0.0;
}//omp

	}//for-k


#pragma omp parallel num_threads(NumThreads) private(j,h,k,v,tid,ptid,ctr_low)
{
	#pragma omp for
	for(tid = ThTerm; tid < NumUTD; tid++){
		ptid =  tid -ThTerm;
		ctr_low = vNM[tid] -1;
		for(h = 0; h < vNM[tid]; ++h){
			if(h > ctr_low) break;
			k = mMID[tid][h]; v = mMean[tid][h];
			if(mMean[tid][h] < 1.0){
				for(j = ctr_low; j > h; j--){
					if(mMean[tid][j] >= 1.0){
						mMID[tid][h] = mMID[tid][j]; mMean[tid][h] = mMean[tid][j];
						mMID[tid][j] = k; mMean[tid][j] = v;
						break;
					}
				}
				ctr_low = j-1;
			}
		}//for-num_means
		vNMlarge[ptid] = ctr_low +1;

		for(k = 0; k < K; ++k) mMeanVal[ptid][k] = 0.0;
		for(h = vNMlarge[ptid]; h < vNM[tid]; ++h)
			mMeanVal[ptid][mMID[tid][h]] = mMean[tid][h];
	}
}//omp

#pragma omp parallel num_threads(NumThreads) private(j,h,k,v,tid,ptid,ctrI,num_means)
{
	#pragma omp for
	for(tid = ThTerm; tid < NumUTD; tid++){
		ptid = tid -ThTerm;
		num_means = vNMlarge[ptid];
		ctrI = num_means-1;
		for(h = 0; h < num_means; ++h){
			if(h > ctrI) break;
			k = mMID[tid][h]; v = mMean[tid][h];
			if(vNew[k] == 0){
				for(j = ctrI; j > h; j--){
					if(vNew[mMID[tid][j]] == 1){
						mMID[tid][h] = mMID[tid][j]; mMean[tid][h] = mMean[tid][j];
						mMID[tid][j] = k; mMean[tid][j] = v;
						break;
					}
				}
				ctrI = j-1;
			}
		}//for-num_means
		vNMmove[tid] = ctrI +1; //Number of moving centroids
	}

	#pragma omp for
	for(tid = 0; tid < ThTerm; tid++){
		num_means = vNM[tid];
		ctrI = num_means-1;
		for(h = 0; h < num_means; ++h){
			if(h > ctrI) break;
			k = mMID[tid][h]; v = mMean[tid][h];
			if(vNew[k] == 0){
				for(j = ctrI; j > h; j--){
					if(vNew[mMID[tid][j]] == 1){
						mMID[tid][h] = mMID[tid][j]; mMean[tid][h] = mMean[tid][j];
						mMID[tid][j] = k; mMean[tid][j] = v;
						break;
					}
				}
				ctrI = j-1;
			}
		}//for-num_means
		vNMmove[tid] = ctrI +1;
	}
}//omp


#pragma omp parallel num_threads(NumThreads)
{
	#pragma omp for
	for(k = 0; k < K; k++) free(mClstMember[k]);
	#pragma omp for
	for(j = 0; j < NumUTD; j++) vTcount[j] = 0;
}//omp

	end_time = omp_get_wtime();
	printf("Update(i=%d) time=%f\n",i+1,end_time-start_time);
	fprintf(stderr,"Update(i=%d) time=%f\n",i+1,end_time-start_time);

	return 0;
}


void	prtAssign(char *fn)
{
	FILE	*fp;
	int	i;
	fp = fopen(fn, "w");
	fprintf(fp, "%d %d %d %d\n",N,K,D,NumUTD);
	for(i = 0; i < N; i++)
		fprintf(fp, "%d %d\n",i+1,vAssID[i]+1);
	fclose(fp);
}
void	prtProcInfo(char *fn)
{
	FILE	*fpr, *fpw;
	char	proc_file[256] = "/proc/";
	const char	fname[] = "/status";
	char	str_pid[128];
	int	c;
	pid_t	num_pid;
	
	num_pid = getpid();
	sprintf(str_pid, "%d", num_pid);
	strcat(proc_file,str_pid); strcat(proc_file,fname);
	fpw = fopen(fn,"w");
	if((fpr = fopen(proc_file,"r")) == NULL){
		printf("ERROR: Unknown File = %s\n", proc_file); exit(1);}
	while((c = getc(fpr)) != EOF) putc(c,fpw);
	fclose(fpr);
	fclose(fpw);
}
void	memFree()
{
	int	i;
	for(i = 0; i < N; ++i){free(mX[i]); free(mTID[i]);}
	free(mX); free(mTID); free(vNT);

	for(i = 0; i < NumUTD; ++i){
		if(mMean[i] != NULL) free(mMean[i]); 
		if(mMID[i] != NULL) free(mMID[i]);
	}
	free(mMean); free(mMID); free(vNM);

	for(i = 0; i < NumThreads; i++){
		if(mSim[i] != NULL) free(mSim[i]);
		if(mAppear[i] != NULL) free(mAppear[i]);
		if(mExactCal[i] != NULL) free(mExactCal[i]);
		if(mObjPL1norm[i] != NULL) free(mObjPL1norm[i]);
	}
	free(mSim); free(mAppear);
	free(mExactCal); free(mObjPL1norm);

	for(i = 0; i < PartTerms; ++i){
		if(mMeanVal[i] != NULL) free(mMeanVal[i]);
	}
	free(mMeanVal);

	free(vClstSize); free(mClstMember); free(vAssID); 
	free(vInOut); free(vNew);
	free(vChObj);
	free(vKcount); free(vTcount); free(vTFIDF);
	free(vSortTID); free(vAssSim);
	free(vDense2Full);
	free(vOrderTh);

	free(vNMlarge); free(vNMmove);
	free(vND); free(vAvgObj); free(vObjStatus);
	free(vAvgMean); free(vAvgSim);
	free(vObjL1norm); free(vMvMID);
}
int	main(int argc, char **argv)
{
	int	i;
	int	seedN, seedK;
	int	max_itr = 5;//Maximum repeat counts (iterations)
	double	all_distcalc;
	double	start_time, end_time;
	char	db_file[256];
	char	assign_file[256];
	char	info_file[256];

	NumThreads = atoi(argv[1]);

	ThTermMin = atoi(argv[2]);
	ThValMin = atof(argv[3]);
	ThValMax = atof(argv[4]);
	ThValStep1 = atof(argv[5]);
	ThValStep2 = atof(argv[6]);

	seedN = atoi(argv[7]);
	seedK = atoi(argv[8]);
	K = atoi(argv[9]);

	strcpy(db_file,argv[10]);
	strcpy(assign_file,argv[11]);
	strcpy(info_file,argv[12]); 

	readData_sparse(db_file);//Sparse coding
	calTFIDF();
	initData(seedN, seedK);
	all_distcalc = (double) N*K;

	start_time = omp_get_wtime();
	if(calMeans_sparse_1st(0,all_distcalc) == 1){
		end_time = omp_get_wtime();
		printf("%f\n",end_time-start_time);
		fprintf(stderr,"%f\n",end_time-start_time);
		prtAssign(assign_file);
		prtProcInfo(info_file);
		memFree();
		return 0;
	}
	end_time = omp_get_wtime();
	printf("%f\n",end_time-start_time);
	fprintf(stderr,"%f\n",end_time-start_time);
	
	start_time = omp_get_wtime();
	if(calMeans_sparse_2nd(1,all_distcalc) == 1){
		end_time = omp_get_wtime();
		printf("%f\n",end_time-start_time);
		fprintf(stderr,"%f\n",end_time-start_time);
		prtAssign(assign_file);
		prtProcInfo(info_file);
		memFree();
		return 0;
	}
	end_time = omp_get_wtime();
	printf("%f\n",end_time-start_time);
	fprintf(stderr,"%f\n",end_time-start_time);
	
	for(i = 2; i < max_itr; i++){
		start_time = omp_get_wtime();
		if(calMeans_sparse(i,all_distcalc) == 1) break;
		end_time = omp_get_wtime();
		printf("%f\n",end_time-start_time);
		fprintf(stderr,"%f\n",end_time-start_time);
	}
	end_time = omp_get_wtime();
	printf("%f\n",end_time-start_time);
	fprintf(stderr,"%f\n",end_time-start_time);

	prtAssign(assign_file);
	prtProcInfo(info_file);
	memFree();

	return 0;
}
