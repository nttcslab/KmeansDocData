// TA-ICP: Main filter == Threshold Algorithm
//		   Auxiliary filter == ICP
//-----------------------------------------------//
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
int	**mTID;//Term ID for each 
double	**mX;//td-idf for each document
int	**mClstMember;//Cluster member ID
int	*vClstSize;//Cluster size: Number of cluster members
int	*vAssID;//ID of cluster to which j-th object is assigned
int	*vInOut;//Boolean flag for object changing its cluster: each iter.
int	*vNew;//Boolean flag whether cluster members change or not
struct	chobj{int out; int in;};
struct	chobj *vChObj;//Object ID that changes its assigned cluster
int	*vDense2Full;

//-- Feature and Similarity to assigned cluster's mean --//
double	*vTFIDF;//TFIDF vector for calculating a mean
double	*vAssSim;//Similarity of each object to its assigned cluster's mean

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

//-- Moving-centroid inverted file --//
int	*vNMvM;//Number of moving means (mf) that change their positions
int	**mMvMID;
double	**mMvMean;
int	NumMvMeans;
int	*vObjStatus;//Flag for invariant or moving cluster

//-- Modified threshold algorithm --//
double	*vObjL1norm;//Object L1 norm
double	**mObjPL1norm;//Object partial L1 norm
double	**mMeanVal;//Full-expression mean-inverted-file
double	PartialDimRate;// = 0.90;
int	ThTerm;//Threshold on termID
int	PartTerms;
int	*vOrderTh;//Order at which each object's tid exceeds ThTerm

//-- Verification phase --//
int	**mExactCal;
int	*vND;//Number of objects (df) that use j-th term with dense format


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
	fclose(fp);
	if(dummy < 0) printf("ERROR: Somethong strange\n");

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

	printf("N=%d D=%d NumUTD=%d K=%d\n",N,D,NumUTD,K);
	fprintf(stderr,"N=%d D=%d NumUTD=%d K=%d\n",N,D,NumUTD,K);
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
	double	*vecIDF;//Inverted document frequency
	int	num_terms, tmp_tid;//=vNT[i]
	int	*vInvMap;

	vND = (int *) calloc(NumUTD,sizeof(int));
	for(i = 0; i < N; ++i){
		num_terms = vNT[i];
		for(j = 0; j < num_terms; ++j) vND[mTID[i][j]]++;
	}
#pragma omp parallel num_threads(NumThreads)
{
	#pragma omp for
	for(i = 0; i < NumUTD; ++i) vSortTID[i].df = vND[vSortTID[i].tid];
}//omp
	quickSortTID(0,NumUTD-1);//vSortTID[i] was arranged in ascending order of df.

	vInvMap = (int *) malloc(sizeof(int)*NumUTD);
#pragma omp parallel num_threads(NumThreads)
{
	#pragma omp for
	for(i = 0; i < NumUTD; ++i){
		vND[i] = vSortTID[i].df;
		vInvMap[vSortTID[i].tid] = i;//Mapping of old TID to new TID
	}
}//omp

	//-- Term ID in each object is sorted in ascending order of ID. --//
#pragma omp parallel num_threads(NumThreads)
{
	#pragma omp for
	for(i = 0; i < N; ++i){
		num_terms = vNT[i];
		for(j = 0; j < num_terms; ++j){
			tmp_tid = mTID[i][j];
			mTID[i][j] = vInvMap[tmp_tid];	
		}
		qsObjElmt(i,0,num_terms-1);
	}
}//omp

	vecIDF = (double *) malloc(sizeof(double)*NumUTD);
#pragma omp parallel num_threads(NumThreads) private(v)
{
	v = log(1.0*N);
	#pragma omp for
	for(i = 0; i < NumUTD; i++) vecIDF[i] = v;
	#pragma omp for
	for(i = 0; i < NumUTD; i++) vecIDF[i] -= log(1.0*vND[i]);
}//omp
	for(i = 0; i < N; i++){
		num_terms = vNT[i];
		for(j = 0; j < num_terms; ++j)
			mX[i][j] *= vecIDF[mTID[i][j]];
		for(j = 0, v = 0.0; j < num_terms; ++j) v += mX[i][j]*mX[i][j];
		for(j = 0, v = 1.0/sqrt(v); j < num_terms; ++j) mX[i][j] *= v;
	}

	free(vecIDF); free(vInvMap);
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
	vNMvM = (int *) malloc(sizeof(int)*NumUTD);
	for(i = 0; i != K; ++i){
		id = vPnt[i];
		num_terms = vNT[id];
		for(j = 0; j < num_terms; ++j){
			tid = mTID[id][j]; vNM[tid]++;
		}
	}

	mMID = (int **) malloc(sizeof(int *)*NumUTD);
	mMean = (double **) malloc(sizeof(double *)*NumUTD);
	mMvMID = (int **) malloc(sizeof(int *)*NumUTD);
	mMvMean = (double **) malloc(sizeof(double *)*NumUTD);
	for(i = 0; i != NumUTD; ++i){
		mMID[i] = (int *) malloc(sizeof(int)*vNM[i]);
		mMean[i] = (double *) malloc(sizeof(double)*vNM[i]);
	}

	ThTerm = (int) (PartialDimRate*NumUTD);
	PartTerms = NumUTD -ThTerm;
	mMeanVal = (double **) malloc(sizeof(double *)*PartTerms);
	for(i = 0; i < PartTerms; ++i)
		mMeanVal[i] = (double *) malloc(sizeof(double)*K);

	vOrderTh = (int *) malloc(sizeof(int)*N);
	vObjL1norm = (double *) calloc(N,sizeof(double));
#pragma omp parallel num_threads(NumThreads) private(i,j)
{
	#pragma omp for
	for(i = 0; i < N; ++i){
		vOrderTh[i] = vNT[i];
		for(j = 0; j < vNT[i]; ++j)
			if(mTID[i][j] >= ThTerm){vOrderTh[i] = j; break;}	
		for(j = vOrderTh[i]; j < vNT[i]; ++j)
			vObjL1norm[i] += mX[i][j];
	}
}//omp

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

	mClstMember = (int **) calloc(K,sizeof(int *));
	vKcount = (int *) calloc(K,sizeof(int));
	vTFIDF = (double *) calloc(NumUTD,sizeof(double));
	vTcount = (int *) calloc(NumUTD,sizeof(int));
	vAssSim = (double *) calloc(N,sizeof(double));
	vObjStatus = (int *) calloc(N,sizeof(int));

	vInOut = (int *) calloc(K,sizeof(int));
	vNew = (int *) malloc(sizeof(int)*K);
	for(i = 0; i != K; ++i) vNew[i] = 1;
}


int	calMeans_sparse_1st(int i, double all_distcalc)
{
	int	j, h, k, t;
	double	v, w;
	double	obj_funct;//Objective function value
	int	id, tid, ptid;//Object_ID, Term_ID, Partial Term_ID
	int	init_cid, final_cid;
	long num_obsim;//Number of multiplications for sim. calc.
	unsigned long num_mult;//Number of multiplications for sim. calc.
	unsigned long num_mult_ub;//Number of multiplications for UB calc.
	int	num_chobj;//Number of objects whose assignments change
	int	num_unchclst;//Number of clusters whose memeber remains unchanged
	double	ass_sim;//Similarity from an object to its assigned centroid
	int	out_cid, in_cid;//Centroid ID of out-going and in-coming
	int	th_id;//Thread ID
	double	dfmf, rate0, rate1, rate2;
	double	start_time, end_time;

	num_mult = 0; num_obsim = 0; num_chobj = 0; obj_funct = 0.0;
	num_mult_ub = 0;

	start_time = omp_get_wtime();
	//-- Object assignment --//
#pragma omp parallel num_threads(NumThreads) private(th_id,k,t,h,init_cid,final_cid,ass_sim,tid,v) reduction(+:num_mult,num_chobj,obj_funct,num_obsim,num_mult_ub)
{
	th_id = omp_get_thread_num();
	#pragma omp for
	for(j = 0; j < N; j++){// Outermost loop 

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

#pragma omp parallel num_threads(NumThreads) private(k)
{
	#pragma omp for
	for(k = 0; k < K; k++){
		if(vClstSize[k] == 0){
			printf("ERROR: ClstSize[%d]=0\n",k); exit(1);}
		vKcount[k] = 0;
	}
}//omp

	num_unchclst = K;
	for(k = 0; k < K; k++){
		mClstMember[k] = (int *) calloc(vClstSize[k],sizeof(int));
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
	printf("KM%d %1.10e %d %e %ld %e %lu %lu %e %d ",
		i+1,obj_funct,num_chobj,rate0,num_obsim,rate1,num_mult_ub,num_mult,rate2,num_unchclst); 
	fprintf(stderr,"KM%d %1.10e %d %e %ld %e %lu %lu %e %d ",
		i+1,obj_funct,num_chobj,rate0,num_obsim,rate1,num_mult_ub,num_mult,rate2,num_unchclst); 
	
	end_time = omp_get_wtime();
	printf("\n"); fprintf(stderr,"\n");
	printf("Assignment(i=%d) time=%f\n",i+1,end_time-start_time);
	fprintf(stderr,"Assignment(i=%d) time=%f\n",i+1,end_time-start_time);
	
	if(num_chobj == 0)
		return 1;//All assignments remain unchanged.


	start_time = omp_get_wtime();
	//-- Mean update: Sparse expression --//
#pragma omp parallel num_threads(NumThreads) private(j)
{
	#pragma omp for
	for(j = 0; j < NumUTD; j++){
		vNM[j] = 0; vNMvM[j] = 0;
		free(mMID[j]); free(mMean[j]);
	}
}//omp

	//-- to avoid a large memory consumption. --//
	//-- NumUTD * (K || NumThreads) is forbidden. --//
	for(k = 0; k < K; k++){
		vNew[k] = vInOut[k]; vInOut[k] = 0;
		for(j = 0; j < NumUTD; ++j) vTcount[j] = 0;

		for(j = 0; j < vClstSize[k]; ++j){
			id = mClstMember[k][j];
			vObjStatus[id] = vNew[k];
			for(h = 0; h < vNT[id]; ++h)
				vTcount[mTID[id][h]] = vNew[k] +1;
		}
		for(tid = 0; tid < NumUTD; ++tid){
			if(vTcount[tid] == 2){
				vNMvM[tid]++; vNM[tid]++; vTcount[tid] = 0;
			}else if(vTcount[tid] == 1){
				vNM[tid]++; vTcount[tid] = 0;
			}
		}
	}//for-K

	for(j = 0; j < NumUTD; ++j){
		mMID[j] = (int *) malloc(sizeof(int)*vNM[j]);
		mMean[j] = (double *) malloc(sizeof(double)*vNM[j]);
		mMvMID[j] = (int *) malloc(sizeof(int)*vNMvM[j]);
		mMvMean[j] = (double *) malloc(sizeof(double)*vNMvM[j]);
	}

	//-- for-k loop --//
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
			for(h = 0, ass_sim = 0.0; h < vNT[id]; ++h){
				tid = mTID[id][h];
				ass_sim += mX[id][h]*vTFIDF[tid];
			}
			vAssSim[id] = ass_sim;
		}
}//omp

#pragma omp parallel num_threads(NumThreads) private(j)
{
		#pragma omp for
		for(j = 0; j < NumUTD; j++) vTFIDF[j] = 0.0;
}//omp

	}//for-k

#pragma omp parallel num_threads(NumThreads) private(tid)
{
	#pragma omp for
	for(tid = 0; tid < NumUTD; tid++) 
		quickSortMean(mMID[tid],mMean[tid],0,vNM[tid]-1);
}//omp

	//-- Partition mean arrays into two blocks based on vNew[k] --//
#pragma omp parallel num_threads(NumThreads) private(tid,ptid,h,j,k,v)
{
	#pragma omp for
	for(tid = 0; tid < ThTerm; tid++){
		for(h = 0, j = 0; h < vNM[tid]; ++h){
			k = mMID[tid][h]; v = mMean[tid][h];
			if(vNew[k] == 1){
				mMvMID[tid][j] = k;
				mMvMean[tid][j++] = v;
			}
		}
		assert(vNMvM[tid] == j);	
	}
	#pragma omp for
	for(tid = ThTerm; tid < NumUTD; tid++){
		ptid = tid -ThTerm;
		for(h = 0, j = 0; h < vNM[tid]; ++h){
			k = mMID[tid][h]; v = mMean[tid][h];
			mMeanVal[ptid][k] = v;
			if(vNew[k] == 1){
				mMvMID[tid][j] = k;
				mMvMean[tid][j++] = v;
			}
		}
		assert(vNMvM[tid] == j);	
	}
}//omp

#pragma omp parallel num_threads(NumThreads) private(j,k)
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
	int	id, tid, ptid;//Object_ID, Term_ID, Partial Term_ID
	int	init_cid, final_cid;
	long num_obsim;//Number of multiplications for sim. calc.
	unsigned long num_mult;//Number of multiplications for sim. calc.
	unsigned long num_mult_ub;//Number of multiplications for UB calc.
	int	num_chobj;//Number of objects whose assignments change
	int	num_unchclst;//Number of clusters whose memeber remains unchanged
	double	ass_sim;//Similarity from an object to its assigned centroid
	double	thv;//Threshold on feature value for each object
	int	out_cid, in_cid;//Centroid ID of out-going and in-coming
	int	th_id;//Thread ID
	double	dfmf, rate0, rate1, rate2;
	int	candidates;
	double	start_time, end_time;

	num_mult = 0; num_obsim = 0; num_chobj = 0; obj_funct = 0.0;
	num_mult_ub = 0;

	start_time = omp_get_wtime();
	//-- Object assignment --//
#pragma omp parallel num_threads(NumThreads) private(th_id,j,k,t,h,init_cid,final_cid,ass_sim,tid,ptid,v,candidates,thv) reduction(+:num_mult,num_chobj,obj_funct,num_obsim,num_mult_ub)
{
	th_id = omp_get_thread_num();
	#pragma omp for
	for(j = 0; j < N; j++){// outermost loop
		init_cid = vAssID[j];//j's assigned cluster (or centroid)
		final_cid = init_cid;
		ass_sim = vAssSim[j];
		thv = ass_sim/vObjL1norm[j];

		//-- Similarity calculation --//
		for(k = 0; k < K; k++){
			mSim[th_id][k] = 0.0; mAppear[th_id][k] = 0;
			mExactCal[th_id][k] = -1;
			mObjPL1norm[th_id][k] = vObjL1norm[j];
		}

		if(vObjStatus[j] == 0){//Invariant or near-invariant objects
			for(t = 0; t < vOrderTh[j]; ++t){// middle loop
				tid = mTID[j][t];
				for(h = 0; h < vNMvM[tid]; ++h){
					k = mMvMID[tid][h]; v = mMvMean[tid][h];
					mSim[th_id][k] += mX[j][t]*v; num_mult++;
				}
			}
			for(t = vOrderTh[j]; t < vNT[j]; ++t){// middle loop
				tid = mTID[j][t];
				for(h = 0; h < vNMvM[tid]; ++h){
					k = mMvMID[tid][h]; v = mMvMean[tid][h];
					if(v < thv) break;

					mSim[th_id][k] += mX[j][t]*v; num_mult++;
					mObjPL1norm[th_id][k] -= mX[j][t];
				}
			}
		}else{
			for(t = 0; t < vOrderTh[j]; ++t){// middle loop
				tid = mTID[j][t];
				for(h = 0; h < vNM[tid]; ++h){
					k = mMID[tid][h]; v = mMean[tid][h];
					mSim[th_id][k] += mX[j][t]*v; num_mult++;
				}
			}
			for(t = vOrderTh[j]; t < vNT[j]; ++t){// middle loop
				tid = mTID[j][t];
				for(h = 0; h < vNM[tid]; ++h){// innermost loop
					k = mMID[tid][h]; v = mMean[tid][h];
					if(v < thv) break;

					mSim[th_id][k] += mX[j][t]*v; num_mult++;
					mObjPL1norm[th_id][k] -= mX[j][t];
				}
			}
		}//if-else

		// Gathering phase
		for(k = 0, candidates = 0; k < K; ++k){
			if(mSim[th_id][k] == 0) continue;

			v = mSim[th_id][k] +thv*mObjPL1norm[th_id][k]; num_mult_ub++;
			if(ass_sim >= v) continue;

			mSim[th_id][candidates] = mSim[th_id][k];
			mExactCal[th_id][candidates++] = k;
		}

		// Varification phase
		for(t = vOrderTh[j]; t < vNT[j]; ++t){// middle loop
			tid = mTID[j][t]; ptid = tid -ThTerm;
			for(h = 0; h < candidates; h++){// innermost loop
				k = mExactCal[th_id][h]; mAppear[th_id][h] = 1;
				v = mMeanVal[ptid][k];
				if(v >= thv) continue;

				mSim[th_id][h] += mX[j][t]*v; num_mult++;
			}
		}	
		for(h = 0; h < candidates; h++){
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

#pragma omp parallel num_threads(NumThreads) private(k)
{
	#pragma omp for
	for(k = 0; k < K; k++){
		if(vClstSize[k] == 0){
			printf("ERROR: ClstSize[%d]=0\n",k); exit(1);}
		vKcount[k] = 0;
	}
}//omp

	num_unchclst = K;
	for(k = 0; k < K; k++){
		mClstMember[k] = (int *) calloc(vClstSize[k],sizeof(int));
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
	printf("KM%d %1.10e %d %e %ld %e %lu %lu %e %d ",
		i+1,obj_funct,num_chobj,rate0,num_obsim,rate1,num_mult_ub,num_mult,rate2,num_unchclst); 
	fprintf(stderr,"KM%d %1.10e %d %e %ld %e %lu %lu %e %d ",
		i+1,obj_funct,num_chobj,rate0,num_obsim,rate1,num_mult_ub,num_mult,rate2,num_unchclst); 

	end_time = omp_get_wtime();
	printf("\n"); fprintf(stderr,"\n");
	printf("Assignment(i=%d) time=%f\n",i+1,end_time-start_time);
	fprintf(stderr,"Assignment(i=%d) time=%f\n",i+1,end_time-start_time);

	if(num_chobj == 0)
		return 1;//All assignments remain unchanged.

	//-- Mean update: Sparse expression --//
	start_time = omp_get_wtime();
#pragma omp parallel num_threads(NumThreads) private(j)
{
	#pragma omp for
	for(j = 0; j < NumUTD; j++){
		vNM[j] = 0; vNMvM[j] = 0;
		free(mMID[j]); free(mMean[j]);
		free(mMvMID[j]); free(mMvMean[j]);
	}
}//omp

	//-- to avoid a large memory consumption. --//
	//-- NumUTD * (K || NumThreads) is forbidden. --//
	for(k = 0; k < K; k++){
		vNew[k] = vInOut[k]; vInOut[k] = 0;
		for(j = 0; j < NumUTD; ++j) vTcount[j] = 0;

		for(j = 0; j < vClstSize[k]; j++){
			id = mClstMember[k][j];
			vObjStatus[id] = vNew[k];
			for(h = 0; h < vNT[id]; h++)
				vTcount[mTID[id][h]] = vNew[k] +1;
		}
		for(tid = 0; tid < NumUTD; ++tid){
			if(vTcount[tid] == 2){
				vNMvM[tid]++; vNM[tid]++; vTcount[tid] = 0;
			}else if(vTcount[tid] == 1){
				vNM[tid]++; vTcount[tid] = 0;
			}
		}
	}//for-K

	for(j = 0; j < NumUTD; ++j){
		mMID[j] = (int *) malloc(sizeof(int)*vNM[j]);
		mMean[j] = (double *) malloc(sizeof(double)*vNM[j]);
		mMvMID[j] = (int *) malloc(sizeof(int)*vNMvM[j]);
		mMvMean[j] = (double *) malloc(sizeof(double)*vNMvM[j]);
	}

	//-- for-k loop --//
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

#pragma omp parallel num_threads(NumThreads) private(tid)
{
	#pragma omp for
	for(tid = 0; tid < NumUTD; tid++) 
		quickSortMean(mMID[tid],mMean[tid],0,vNM[tid]-1);
}//omp

	//-- Partition mean arrays into two blocks based on vNew[k] --//
#pragma omp parallel num_threads(NumThreads) private(tid,ptid,h,j,k,v)
{
	#pragma omp for
	for(tid = 0; tid < ThTerm; tid++){
		for(h = 0, j = 0; h < vNM[tid]; ++h){
			k = mMID[tid][h]; v = mMean[tid][h];
			if(vNew[k] == 1){
				mMvMID[tid][j] = k;
				mMvMean[tid][j++] = v;//Already sorted
			}
		}
		assert(vNMvM[tid] == j);	
	}
	#pragma omp for
	for(tid = ThTerm; tid < NumUTD; tid++){
		ptid = tid -ThTerm;
		for(k = 0; k < K; ++k) mMeanVal[ptid][k] = 0.0;
		for(h = 0, j = 0; h < vNM[tid]; ++h){
			k = mMID[tid][h]; v = mMean[tid][h];
			mMeanVal[ptid][k] = v;
			if(vNew[k] == 1){
				mMvMID[tid][j] = k;
				mMvMean[tid][j++] = v;
			}
		}
		assert(vNMvM[tid] == j);	
	}
}//omp

#pragma omp parallel num_threads(NumThreads) private(j,k)
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
	free(mX); free(mTID); free(vNT); free(vObjL1norm);
	free(vOrderTh);

	for(i = 0; i < NumUTD; ++i){
		if(mMean[i] != NULL) free(mMean[i]); 
		if(mMID[i] != NULL) free(mMID[i]);
		if(mMvMean[i] != NULL) free(mMvMean[i]);
		if(mMvMID[i] != NULL) free(mMvMID[i]);
	}
	free(mMean); free(mMID); free(vNM);
	free(mMvMean); free(mMvMID); free(vNMvM);

	for(i = 0; i < PartTerms; ++i)
		if(mMeanVal[i] != NULL) free(mMeanVal[i]);
	free(mMeanVal);

	for(i = 0; i < NumThreads; i++){
		if(mSim[i] != NULL) free(mSim[i]);
		if(mAppear[i] != NULL) free(mAppear[i]);
		if(mExactCal[i] != NULL) free(mExactCal[i]);
		if(mObjPL1norm[i] != NULL) free(mObjPL1norm[i]);
	}
	free(mSim); 
	free(mAppear); free(mExactCal); free(mObjPL1norm);

	free(vClstSize); free(mClstMember); free(vAssID); 
	free(vInOut); free(vNew);
	free(vChObj);
	free(vKcount); free(vTcount); free(vTFIDF);
	free(vSortTID); free(vAssSim);
	free(vDense2Full);

	free(vND); free(vObjStatus);
}
int	main(int argc, char **argv)
{
	int	i;
	int	seedN, seedK;
	int	max_itr = 10000;//Maximum repeat counts (iterations)
	double	all_distcalc;
	double	start_time, end_time;
	char	db_file[256];
	char	assign_file[256];
	char	info_file[256];

	NumThreads = atoi(argv[1]);

	seedN = atoi(argv[2]);
	seedK = atoi(argv[3]);
	K = atoi(argv[4]);

	PartialDimRate = atof(argv[5]);

	strcpy(db_file,argv[6]);
	strcpy(assign_file,argv[7]);
	strcpy(info_file,argv[8]); 

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

	for(i = 1; i < max_itr; i++){
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
