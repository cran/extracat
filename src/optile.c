#include <R.h>
#include <Rinternals.h>
#include <math.h>






SEXP hammdist(SEXP dset){
	
	int n;
	int m;
	int i, i2;
	int j,k,d,d2;
	
	int dist;
	
	SEXP sDim = getAttrib(dset, R_DimSymbol);
	Rprintf("test");
	
	if (isInteger(sDim) && LENGTH(sDim) == 2) {
		int *dim = INTEGER(sDim);
		Rprintf("%d x %d matrix\n", dim[0], dim[1]);
		n = dim[0];
		m = dim[1];
	} else error("invalid dimensions");
	
	SEXP dm = allocVector(INTSXP,n*(n-1)/2);
	k = 0;
	for( i=0; i < n-1; i++ ){
		for( i2 = i+1; i2 < n; i2++ ){
			
			dist = 0;
			
			
			for( j=0; j < m; j++ ){
				d = n*j + i;
				d2 = n*j + i2;
				dist = (INTEGER(dset)[d] == INTEGER(dset)[d2]) ? dist : dist+1;
			}
			INTEGER(dm)[k] = dist;
			k = k+1;
			
		}
	}
	
	
	return dm;
}


float critperm3(int *m0, int n, int m, int nt){
	int i;
	int j,k,p,h;
	int cv2;
	float cv;
	float dist;
	float loss = 0.0;
	//Rprintf("critty3\n");
	int (*tab)[m*nt] = (int (*)[m*nt])m0;
	for( k = 0; k < nt; k++ ){  // nt subtables
		for( i=0; i < n; i++ ){
			for( j= 0+k*m; j < m +k*m; j++ ){
				
				
				cv = tab[i][j];
				//	Rprintf("current %d\n", tab[i][j]);
				for( h=0; h < i+1; h++ ){
					for( p=j; p < m +k*m; p++ ){
						dist =  pow( (i-h)*(i-h) + (p-j)*(p-j), 0.5 );
						cv2 = tab[h][p];
						loss = loss + cv*cv2*dist;
					}
				}
				
			}
		}
	}
	return loss;
}

float critperm4(int *m0, int n, int m, int nt){
	//measure of effectiveness
	int i,j,k;
	
	float loss = 0.0;
	
	int (*tab)[m*nt] = (int (*)[m*nt])m0;
	for( k = 0; k < nt; k++ ){  // nt subtables
		for( i=1; i < n; i++ ){
			//loss = loss + tab[i][j]*tab[i-1][j+1];
			for( j= 0+k*m; j < m +k*m - 1; j++ ){
				loss = loss + 10000/tab[i][j]*(tab[i-1][j]+tab[i][j+1]);
			}
		}
	}
	return loss;
}

float critperm1(int *m0, int n, int m, int nt){
	int i;
	int j,k;
	
	float loss = 0.0;

	int (*tab)[m*nt] = (int (*)[m*nt])m0;
	long M[n][m];
	
	
	//set the first row and the last column
	for( k = 0; k < nt; k++ ){
		for( i = 0; i < n; i++ ){ 
			M[i][m+k*m-1] = tab[i][m+k*m-1];
		}
	}
	for( j = 1; j < m*nt; j++ ){ 
		M[0][j] = tab[0][j];
	}
	
	for( k = 0; k < nt; k++ ){
		for( j= m +k*m - 2; j > 0 ; j-- ){
			for( i=0; i < n-1; i++ ){
				M[i][j] = tab[i][j] + M[i][j+1];
			}
		}
		for( i=1; i < n-1; i++ ){
			for( j = 1 + k*m; j < m + k*m ; j++ ){
				M[i][j] = M[i][j] + M[i-1][j];
				// add all rows but the first
				loss = loss + tab[i+1][j-1]*M[i][j];
			}
		}
		
		// add the first row terms
		for( j = 1 + k*m; j < m + k*m ; j++ ){
				loss = loss + tab[1][j-1]*M[0][j];
		}
		
	}
	return loss;
}




float critperm2(int *m0, int n, int m, int nt){
	int i;
	int j,k;
	
	float loss = 0.0;
	
	int (*tab)[m*nt] = (int (*)[m*nt])m0;
	long M2[n][m]; // cumsum by column
	long M1[n][m]; // cumsum by row
	
	
	//set the first row and the last column(s)
	for( k = 0; k < nt; k++ ){
		for( i = 0; i < n; i++ ){ 
			M1[i][m+k*m-2] = tab[i][m+k*m-1];
			M1[i][m+k*m-1] = 0;
		}
	}
	for( j = 0; j < m*nt; j++ ){ 
		M2[0][j] = 0;
		M2[1][j] = tab[0][j];
	}
	
	for( k = 0; k < nt; k++ ){
		
		// compute M1
		if(m > 2){ // CHECK ME
		for( j= m +k*m - 3; j >= 0 ; j-- ){
			for( i=0; i < n; i++ ){
				M1[i][j] = tab[i][j+1] + M1[i][j+1];
			}
		}
		}
		for( j= m +k*m - 3; j >= 0 ; j-- ){
			for( i=0; i < n; i++ ){
				M1[i][j] = M1[i][j] + M1[i][j+1];
			}
		}

		for( i=1; i < n; i++ ){
			for( j= m +k*m - 2; j >= 0 ; j-- ){
				M1[i][j] = M1[i][j] + M1[i-1][j];
			}
		}
			
		
		// compute M2	
		if(n > 2){
		for( i=2; i < n; i++ ){
			for( j = 0 + k*m; j < m + k*m ; j++ ){
				M2[i][j] = tab[i-1][j] + M2[i-1][j];
			}
		}
		}
		for( i=2; i < n; i++ ){
			for( j = 0 + k*m; j < m + k*m ; j++ ){
				M2[i][j] = M2[i][j] + M2[i-1][j];
			}
		}

		for( j = m + k*m - 2; j >= 0 ; j-- ){
			for( i=1; i < n; i++ ){
				M2[i][j] = M2[i][j] + M2[i][j+1];
			}
		}
		
		// compute the criterion
		for( i=0; i < n; i++ ){
			for( j = 0 + k*m; j < m + k*m ; j++ ){
				loss = loss + (M1[i][j]+M2[i][j])*tab[i][j];
			}
		}

	}
	return loss;
}





int bincoef(int n, int k){
	int s = 1;

	for(int i = (k+1); i < n+1; i++){
		s = s*i/(i-k);
	}
	
	return s;
}

// this function will perform permutations for the categories until a local optimum is reached
SEXP optisaic(SEXP M, SEXP dims, SEXP method, SEXP criterion, SEXP pv ){
	int pm = INTEGER(method)[0]; // 0,1 or 2 permutations
	int cm = INTEGER(criterion)[0]; // id number of critfun
	int n = INTEGER(dims)[0];
	int m = INTEGER(dims)[1];
	int nt = INTEGER(dims)[2]; // no. of tables
	int px = INTEGER(pv)[1];
	int py = INTEGER(pv)[0];
	//Rprintf("px=%d, py=%d, pm=%d",px,py,pm);
	if( px+py < pm ){
		pm = px+py;
	}
	//Rprintf("\nstart C code: pm = %d, cm = %d, n = %d, m = %d, nt = %d\n",pm,cm,n,m,nt);
	int s,si,i,j,z, tmp,i1,i2,j1,j2, bi;
	if(pm == 2){
		s = (bincoef(n,2)+1)*(bincoef(m,2)+1);

	}
	if(pm == 1){
		s = py*bincoef(n,2)+px*bincoef(m,2)+1;
	}
	if( pm == 0 ){
		s = 1;
	}
	
	float currmin;
	float best; // req check	
	float crtv[s];

	
	// write SEXP M in MX and MY
	int MX[n][m*nt];
	int MY[n][m*nt];
	for(i = 0; i < n; i++){
		for(j = 0; j < m*nt; j++){
			MX[i][j] = INTEGER(M)[n*j + i]; // the current state of the matrix
			MY[i][j] = MX[i][j];			// the temporary state of the matrix
		}
	}
	//Rprintf("MY = :\n");
	//for( i = 0; i < n; i++ ){
	//	for( j = 0; j < m-1; j++){
		//	Rprintf("%d\t",MY[i][j]);
			
		//}
		//Rprintf("%d\n",MY[i][m-1]);
	//} 

	// matrix with all permutations:
    int k = 0;
    int P[4][s];
	//Rprintf("S = %d\n",s);
if(pm == 2){
	for(j1 = 0; j1 < m-1; j1++){
		for(j2 = j1; j2 < m; j2++){
			if( j1 != j2 || (j1 == 0 && j2 == 0)  ){
				for(i1 = 0; i1 < n-1; i1++){
					for(i2 = i1; i2 < n; i2++){
						if( i1 != i2 || (i1 == 0 && i2 == 0)  ){
                                 P[0][k] = i1;
                                 P[1][k] = i2;
                                 P[2][k] = j1;
                                 P[3][k] = j2;
                             k++;
								 //Rprintf("k = %d\t",k);
                             }       
                         }
                      }
               }             
		}
	}
}

if(pm == 1){
	P[0][k] = 0;
	P[1][k] = 0;
	P[2][k] = 0;
	P[3][k] = 0;
	k++;
	if (px == 1) {
	for(j1 = 0; j1 < m-1; j1++){
		for(j2 = j1; j2 < m; j2++){
			if( j1 != j2  ){
                    P[0][k] = 0;
                    P[1][k] = 0;
                    P[2][k] = j1;
                    P[3][k] = j2;
                    k++;
				//Rprintf("(%d - %d %d %d %d)\t",k,P[0][k],P[1][k],P[2][k],P[3][k]);
			   }             
		}
	}
	}
	if(py == 1){
	for(i1 = 0; i1 < n-1; i1++){
		for(i2 = i1; i2 < n; i2++){
			if( i1 != i2   ){
                   P[0][k] = i1;
                   P[1][k] = i2;
                   P[2][k] = 0;
                   P[3][k] = 0;
                   k++;
               }             
		}
	}
	}
}
if(pm == 0 ){
         P[0][k] = 0;
         P[1][k] = 0;
         P[2][k] = 0;
         P[3][k] = 0;    
}

	//Rprintf("P = :\n");
	//for( i = 0; i < 4; i++ ){
		//for( j = 0; j < s-1; j++){
			//Rprintf("%d\t",P[i][j]);
			
		//}
		//Rprintf("%d\n",P[i][s-1]);
	//} 
	
	int* mv = &MY[0][0]; // a pointer to MY for the crit function to read from
    int opt = 0;
   	int ordx[m];
	int ordy[n];
	for( i = 0; i < n; i++ ){
         ordy[i] = i;
         }
 for( j = 0; j < m; j++ ){
         ordx[j] = j;
         }

	//int iter = 0;
// SCALING -----------------------------------------------------------------------------------o	
// ?!?
// SCALING -----------------------------------------------------------------------------------o	
	
   	while( opt == 0 ){
	
		for (si = 0; si < s; si++) { //goes through all permutations in P
			//Rprintf("start permutation %d",si);
			
				
				
				i1 = P[0][si];
				i2 = P[1][si];
				j1 = P[2][si];
				j2 = P[3][si];
				//Rprintf("Perm(%d, %d, %d, %d)",i1,i2,j1,j2);
				
							// permute the clomns in MY rowwise
							for(z = 0; z < n; z++){
								for( k = 0; k < nt; k++ ){
									MY[z][j1+k*m] = MX[z][j2+k*m];
									MY[z][j2+k*m] = MX[z][j1+k*m];
								}
							}
						
							// compute the rows of MY columnwise
							for(z = 0; z < (m*nt); z++){
									tmp = MY[i1][z];
									MY[i1][z] = MY[i2][z];
									MY[i2][z] = tmp;
							}
							//Rprintf("start crit number %d\n", cm);
							// compute the overall criterion for MY
										if(cm == 1){
											crtv[si] = critperm1(mv, n, m, nt);
											//Rprintf("val = %d",crtv[si]);
										}
										if(cm == 2){
											crtv[si] = critperm2(mv, n, m, nt);
										}
										if(cm == 3){
											crtv[si] = critperm3(mv, n, m, nt);
										}
			
			
			if(cm == 4){
				crtv[si] = critperm4(mv, n, m, nt);
			}
						
						// undo the current permutation
						for(z = 0; z < (m*nt); z++){
							tmp = MY[i1][z];
							MY[i1][z] = MY[i2][z];
							MY[i2][z] = tmp;
						}			
		
						// undo the current permutation #2
						for(z = 0; z < n; z++){
								for( k = 0; k < nt; k++ ){
									MY[z][j1+k*m] = MX[z][j1+k*m];
									MY[z][j2+k*m] = MX[z][j2+k*m];
								}
						}
						
	
    } // end si
		currmin = crtv[0];
		best = crtv[0];
		bi = 0;
		
		for( i = 0; i < s; i++ ){
			if(crtv[i] < currmin){
				bi = i;
				currmin = crtv[i];
			}
		}
		if( currmin < best ){
			best = currmin;
		}else {
			opt = 1;
		}
	
		// optimal permutation applied to MX and MY as well as ordx and ordy
			i1 = P[0][bi];
			i2 = P[1][bi];
			j1 = P[2][bi];
			j2 = P[3][bi];
			//Rprintf("best permutation is %d - %d,%d,%d,%d with value: %f\n",bi,i1,i2,j1,j2,best);
		//iter++;
		//if(iter>3){
		//	break;
		//}
			// permute the clomns in MY rowwise
			for(z = 0; z < n; z++){
				for( k = 0; k < nt; k++ ){
					MY[z][j1+k*m] = MX[z][j2+k*m];
					MY[z][j2+k*m] = MX[z][j1+k*m];
					MX[z][j1+k*m] = MY[z][j1+k*m];
					MX[z][j2+k*m] = MY[z][j2+k*m];
				}
			}
			// compute the rows of MY columnwise
			for(z = 0; z < (m*nt); z++){
				tmp = MY[i1][z];
				MY[i1][z] = MX[i2][z];
				MY[i2][z] = MX[i1][z];
				MX[i1][z] = MY[i1][z];
				MX[i2][z] = MY[i2][z];
			}
			tmp = ordx[j1]; 
			ordx[j1] = ordx[j2];
			ordx[j2] = tmp;
		
			tmp = ordy[i1]; 
			ordy[i1] = ordy[i2];
			ordy[i2] = tmp;

		
		
		
	} // end while
	SEXP out = allocVector(REALSXP,n+m+1);
	for( i = 0; i < n;i++ ){
		REAL(out)[i] = ordy[i]+1;
	}
	for( j = 0; j < m;j++ ){
		REAL(out)[n+j] = ordx[j]+1;
	}

	REAL(out)[n+m] = best;
	return out;
}

SEXP simplecrit(SEXP M, SEXP dims, SEXP criterion){
	int cm = INTEGER(criterion)[0]; // id number of critfun
	int n = INTEGER(dims)[0];
	int m = INTEGER(dims)[1];
	int nt = INTEGER(dims)[2]; // no. of tables

	int i,j;
	float val;
	int MX[n][m*nt];
	for(i = 0; i < n; i++){
		for(j = 0; j < m*nt; j++){
			MX[i][j] = INTEGER(M)[n*j + i]; 
		}
	}
	int* mv = &MX[0][0];
	if(cm == 1){
		val = critperm1(mv, n, m, nt);
	}
	if(cm == 2){
		val = critperm2(mv, n, m, nt);
	}
	if(cm == 3){
		val = critperm3(mv, n, m, nt);
	}
	
	if(cm == 4){
		val = critperm4(mv, n, m, nt);
	}

	SEXP out = allocVector(REALSXP,1);
	REAL(out)[0] = val;
	return out;

	
}

int mymax(int a, int b){
	if(a < b){
		return b;
	}else{
		return a;
	}
}
int mymin(int a, int b){
	if(a < b){
		return a;
	}else{
		return b;
	}
}



SEXP classcrit(SEXP M, SEXP dims, SEXP rord, SEXP cord){
	
	int n = INTEGER(dims)[0];
	int m = INTEGER(dims)[1];
	int i,j;
	
	int MX[n][m];
	for(i = 0; i < n; i++){
		for(j = 0; j < m; j++){
			MX[i][j] = INTEGER(M)[n*j + i]; 
		}
	}
	
	
	float loss = 0.0;
	
	float MY[n][m];
	
	//set the first row and the last column
	
	for( i = 0; i < n; i++ ){ 
		MY[i][m-1] = MX[ INTEGER(rord)[i] ][ INTEGER(cord)[m-1] ];
	}
	
	for( j = 1; j < m; j++ ){ 
		MY[0][j] = MX[ INTEGER(rord)[0] ][ INTEGER(cord)[j] ];
	}
	
	
		for( j = m - 2; j > 0 ; j-- ){
			for( i=0; i < n-1; i++ ){
				MY[i][j] = MX[ INTEGER(rord)[i] ][ INTEGER(cord)[j] ] + MY[i][j+1];
			}
		}
		for( i=1; i < n-1; i++ ){
			for( j = 1 ; j < m  ; j++ ){
				MY[i][j] = MY[i][j] + MY[i-1][j];
				// add all rows but the first
				loss = loss + MX[ INTEGER(rord)[i+1] ][ INTEGER(cord)[j-1] ] * MY[i][j];
			}
		}
		
		// add the first row terms
		for( j = 1 ; j < m  ; j++ ){
			loss = loss + MX[ INTEGER(rord)[1] ][ INTEGER(cord)[j-1] ]*MY[0][j];
		}
	
	SEXP out = allocVector(REALSXP,1);
	REAL(out)[0] = loss;
	
	return out;
}



// __________________________________________________________________________________________________________________________________________________________________________________________________________________________//
// __________________________________________________________________________________________________________________________________________________________________________________________________________________________//
// __________________________________________________________________________________________________________________________________________________________________________________________________________________________//



SEXP quicktile(SEXP M, SEXP dims, SEXP pv){
	
	// M is byrow
	int n = INTEGER(dims)[0];
	int m = INTEGER(dims)[1];
	
	int px = INTEGER(pv)[1];
	int py = INTEGER(pv)[0];
	
	
	
	//Rprintf("%d\t%d\t%d\t%d\t%d\t",n,m,px,py);
	
	int i,j,i2,j2,k,tmp;
	
	
	// index vectors with initial values 0..n-1 and 0..m-1 containing the actual orders
	int RI[n];
	int CI[m];
	for (i = 0; i < n; i++) {
		RI[i] = i;
	}
	for (j = 0; j < m; j++) {
		CI[j] = j;
	}
	
	// Delta matrices for rows and columns
	// crit( i -> i2 ) - crit( i2 -> i ) if i < i2 and -(...) else
	// i.e. this is 'lost' when putting i to i2
	
	
	
	
	// Sum of delta matrices 
	
		float DC2[m][m];
		float DR2[n][n];
		//int DC[m][m];
		//int DR[n][n];
		int MS[n][m];
		int MT[n][m];
	float DC[m][m];
	float DR[n][n];
	//double MS[n][m];
	//double MT[n][m];
	//double DC2[m][m];
	//double DR2[n][n];
	//Rprintf("check01");
	if (py == 1) {
		for (i = 0; i < n; i++) {
			DR2[i][i] = 0;
			for (i2 = 0; i2 < n; i2++) {
				DR[i][i2] = 0;
			}
		}
	}
	if (px == 1) {
		for (j = 0; j < m; j++) {
			DC2[j][j] = 0;
			
			for (j2 = 0; j2 < m; j2++) {
				DC[j][j2] = 0;
			}
		}
	}
	
	//Rprintf("check02");
	// the criterion
	float crtv = 0;
	float td;
	
	//an marray pointer to M
	
	 //int* mv = &INTEGER(M)[0];
	 //int (*MT)[m] = (int (*)[m])mv;
	//for (i = 0; i < n; i++) {
	//	for (j = 0; j < m; j++) {
	//	Rprintf("%d\t", MT[i][j]);
	//}
	//Rprintf("\n");
	//}

	for (j = 0; j < m; j++) {
		for (i = 0; i < n; i++) {
			MT[i][j] = INTEGER(M)[i+j*n];
		}
	}

	
	
// --------------------- compute the initial delta values --------------------- //
	
		// cumsum by columns
		for (j = 0; j < m; j++) {
			MS[0][j] = MT[0][j];
		}
		for (i = 1; i < n; i++) {
			for (j = 0; j < m; j++) {
				MS[i][j] = MT[i][j]+MS[i-1][j];
			}
		}
		
		// delta C
		for (j = 0; j < m-1; j++) {
			for (j2 = j+1; j2 < m; j2++) {
				for (i = 1; i < n; i++) {
					td = MT[i][j]*MS[i-1][j2];
					crtv += td;
					//DC[j][j2] += ( MT[i][j]*MS[i-1][j2] - MT[i][j2]*MS[i-1][j] ); // loss
					DC[j][j2] += ( td - MT[i][j2]*MS[i-1][j] ); // loss
				}
				DC[j2][j] = -DC[j][j2];
			}
		}
		// also computes the initial criterion
	
	if (py == 1) {
		// cumsum by rows
		for (i = 0; i < n; i++) {
			MS[i][0] = MT[i][0];
		}
		for (i = 0; i < n; i++) {
			for (j = 1; j < m; j++) {
				MS[i][j] = MT[i][j]+MS[i][j-1];
			}
		}
		
		// delta R
		for (i = 0; i < n-1; i++) {
			for (i2 = i+1; i2 < n; i2++) {
				for (j = 1; j < m; j++) {
					DR[i][i2] += ( MT[i][j]*MS[i2][j-1] - MT[i2][j]*MS[i][j-1] ); // loss
				}
				DR[i2][i] = -DR[i][i2];
			}
		}
	}
	

	
	//Rprintf("\n\nstart crtv = %f\n\n",crtv);
	
// --------------------- start the main algorithm --------------------- //	
	
	
	// look for the minimum among all cumulative sums in DC and DR starting from the diagonal (put i to i2)
	// 'best' is the decrease in crtv, r1, r2, c1, c2 are the corresp. indices
	float bestc = 0;
	float bestr = 0;
	float best = 0;
	int r1 = 0;
	int r2 = 0;
	int c1 = 0;
	int c2 = 0;
	int opt = 0;
	int steps = 0;

	
	while( opt == 0 ){	
		
		
		//look for the minimum among all cumulative sums in DC and DR starting from the diagonal (put i to k)
		
		if (py == 1) {
			for (i = 0; i < n-1; i++) {
				for (i2 = i+1; i2 < n; i2++) {
					DR2[i][i2] = DR2[i][i2-1] + DR[ RI[i] ][ RI[i2] ];
					
					if (DR2[i][i2] > bestr) {
						
						bestr = DR2[i][i2];
						r1 = i;
						r2 = i2;
					}
				}
				
			}
			
			for (i = 1; i < n; i++) {
				for (i2 = i-1; i2 >= 0; i2--) {
					DR2[i][i2] = DR2[i][i2+1] +  DR[ RI[i2] ][ RI[i] ];//i, i2?
					
					if (DR2[i][i2] > bestr) {
						
						bestr = DR2[i][i2];
						r1 = i;
						r2 = i2;
					}
				}
				
			}
		}
		
		
		if (px == 1) {
			for (j = 0; j < m-1; j++) {
				for (j2 = j+1; j2 < m; j2++) {
					DC2[j][j2] = DC2[j][j2-1] + DC[ CI[j] ][ CI[j2] ];
					
					if (DC2[j][j2] > bestc) {
						
						bestc = DC2[j][j2];
						c1 = j;
						c2 = j2;
					}
				}
				
			}
			
			for (j = 1; j < m; j++) {
				for (j2 = j-1; j2 >= 0; j2--) {
					DC2[j][j2] = DC2[j][j2+1] + DC[ CI[j2] ][ CI[j] ]; // j, j2??
					
					if (DC2[j][j2] > bestc) {
						
						bestc = DC2[j][j2];
						c1 = j;
						c2 = j2;
					}
				}
				
			}
		}
	
	// two columns OR two rows are exchanged. The DR/DC values will change whilst DC/DR remains unchanged.
	// the original data matrix and the delta matrices DC and DR and MS remain unchanged with repect to their row/column orders.
		// BUT the values change dependent on the permutations! Hence the sums (DR2 and DC2) use CI[j] and RI[i]
		
		//Rprintf("\n\n r1 = %d, r2 = %d, bestr = %f\n",r1,r2, bestr);
		//Rprintf("\n\n c1 = %d, c2 = %d, bestc = %f\n",c1,c2, bestc);
		
		//bestc = bestr+1;
		
		
		// reset bestc and bestr
		if( bestc > 0 || bestr > 0 ){	
			
		if (bestc > bestr) {
			best = bestc;
			
			if (c1 < c2) {
				// cumulative partial rowsums
				for (i = 0; i < n; i++) {
					MS[i][ CI[c2] ] = MT[i][ CI[c2] ];
				}
				if (c1+1 < c2) { // not neighboring
					for (k = c2-1; k > c1; k--) {
						for (i = 0; i < n; i++) {
							MS[i][ CI[k] ] = MT[i][ CI[k] ]+MS[i][ CI[k+1] ];
						}
					}
				}
				if (py == 1) {
					// changes to DR
					for (i = 0; i < n-1; i++) {
						for (i2 = i+1; i2 < n; i2++) {
							DR[i][i2] -= 2*(  MT[i2][ CI[c1] ]*MS[i][ CI[c1+1] ] -  MT[i][ CI[c1] ]*MS[i2][ CI[c1+1] ] );
							DR[i2][i] = -DR[i][i2];
						}
					}
				}
				
				// index vector CI changes
				tmp = CI[c1];
				for (i = c1; i < c2; i++) {
					CI[i] = CI[i+1];
				}
				CI[c2] = tmp;
				//Rprintf("------------M-------------\n");
				
				//for (i = 0; i < n; i++) {
				//	for (j = 0; j < m; j++) {
				//		Rprintf("%d\t",MT[ RI[i] ][ CI[j] ]);
				//	}
				//	Rprintf("\n");
				//}
				
				//Rprintf("----------------------------\n");
				
			}else {
				// cumulative partial rowsums
				for (i = 0; i < n; i++) {
					MS[i][ CI[c2] ] = MT[i][ CI[c2] ];
				}
				if (c2+1 < c1) { // not neighboring
					for (k = c2+1; k < c1; k++) {
						for (i = 0; i < n; i++) {
							MS[i][ CI[k] ] = MT[i][ CI[k] ]+MS[i][ CI[k-1] ];
						}
					}
				}
				if (py == 1) {
					// changes to DR
					for (i = 0; i < n-1; i++) {
						for (i2 = i+1; i2 < n; i2++) {
							DR[i][i2] -= 2*(  MT[i][ CI[c1] ]*MS[i2][ CI[c1-1] ] - MT[i2][ CI[c1] ]*MS[i][ CI[c1-1] ] ); // neg. of c1 < c2
							DR[i2][i] = -DR[i][i2];
						}
					}
				}
				
				// index vector CI changes
				tmp = CI[c1];
				for (i = c1; i > c2; i--) {
					CI[i] = CI[i-1];
				}
				CI[c2] = tmp;
			}

		}else{ // bestr >= bestc
		
			
			best = bestr;
			if (r1 < r2) {
				// cumulative partial colsums
				for (j = 0; j < m; j++) {
					MS[ RI[r2] ][j] = MT[ RI[r2] ][j];
				}
				if (r1+1 < r2) { // not neighboring
					for (k = r2-1; k > r1; k--) {
						for (j = 0; j < m; j++) {
							MS[ RI[k] ][j] = MT[ RI[k] ][j]+MS[ RI[k+1] ][j];
						}
					}
				}
			
				
				if (px == 1) {
					// changes to DC
					for (j = 0; j < m-1; j++) {
						for (j2 = j+1; j2 < m; j2++) {
							DC[j][j2] -= 2*( MT[ RI[r1] ][j2]*MS[ RI[r1+1] ][j] - MT[ RI[r1] ][j]*MS[ RI[r1+1] ][j2]);
							DC[j2][j] = -DC[j][j2];
						}
					}
				}
				
				// index vector RI changes
				tmp = RI[r1];
				for (i = r1; i < r2; i++) {
					RI[i] = RI[i+1];
				}
				RI[r2] = tmp;
			
			}else {
				// cumulative partial colsums
				for (j = 0; j < m; j++) {
					MS[ RI[r2] ][j] = MT[ RI[r2] ][j];
				}
				if (r2+1 < r1) { // not neighboring
					for (k = r2+1; k < r1; k++) {
						for (j = 0; j < m; j++) {
							MS[ RI[k] ][j] = MT[ RI[k] ][j]+MS[ RI[k-1] ][j];
						}
					}
				}
				if (px == 1) {
					// changes to DC
					for (j = 0; j < m-1; j++) {
						for (j2 = j+1; j2 < m; j2++) {
							DC[j][j2] -= 2*(  MT[ RI[r1] ][j]*MS[ RI[r1-1] ][j2] - MT[ RI[r1] ][j2]*MS[ RI[r1-1] ][j] );
							DC[j2][j] = -DC[j][j2];
						}
					}
				}
				
				// index vector RI changes
				tmp = RI[r1];
				for (i = r1; i > r2; i--) {
					RI[i] = RI[i-1];
				}
				RI[r2] = tmp;
				
				
				
			}

			
		}
		
		// reset
		bestc = 0;
		bestr = 0;
	
		
			crtv -= best;
			//Rprintf("step %d, crtv = %f\n",steps,crtv);
			//Rprintf("\n\nnew cv = %f\n\n",crtv);
			best = 0;
		}else { // best <= 0
			opt = 1;
		}
		//if( steps > 16 ){
		//	opt = 1;
		//}
		steps++;
	
}
	
	
	
	
	
	SEXP out = allocVector(REALSXP,n+m+1);
	for( i = 0; i < n;i++ ){
		REAL(out)[i] = RI[i]; //+1
	}
	for( j = 0; j < m;j++ ){
		REAL(out)[n+j] = CI[j]; //+1
	}
	
	REAL(out)[n+m] = crtv;
	//Rprintf("\nOptimization ended after %d steps with crtv=%f\n",steps,crtv);
	return out;
}




// __________________________________________________________________________________________________________________________________________________________________________________________________________________________//
// __________________________________________________________________________________________________________________________________________________________________________________________________________________________//
// __________________________________________________________________________________________________________________________________________________________________________________________________________________________//



float indcrit1(int *m0, int *rord, int *cord, int n, int m, int rev){
	

	int i,j;
	int (*MX)[m] = (int (*)[m])m0;
	
		
	float loss = 0.0;
	
	float MY[n][m];
	
	//set the first row and the last column
	if (rev == 1) {
		for( i = 0; i < n; i++ ){ 
			MY[i][0] = MX[ rord[i] ][ cord[0] ];
		}
		
		for( j = 1; j < m; j++ ){ 
			MY[0][j] = MX[ rord[0] ][ cord[j] ];
		}
		
		
		for( j = 1; j < m-1 ; j++ ){
			for( i=0; i < n-1; i++ ){
				MY[i][j] = MX[ rord[i] ][ cord[j] ] + MY[i][j-1];
			}
		}
		for( i=1; i < n-2; i++ ){
			for( j = 0 ; j < m-1  ; j++ ){
				MY[i][j] = MY[i][j] + MY[i-1][j];
				// add all rows but the first
				loss = loss + MX[ rord[i+1] ][ cord[j+1] ] * MY[i][j];
			}
		}
		
		// add the first row terms
		for( j = 1 ; j < m  ; j++ ){
			loss = loss + MX[ rord[1] ][ cord[j-1] ]*MY[0][j];
		}
		
	}else{
		for( i = 0; i < n; i++ ){ 
			MY[i][m-1] = MX[ rord[i] ][ cord[m-1] ];
		}
		
		for( j = 1; j < m; j++ ){ 
			MY[0][j] = MX[ rord[0] ][ cord[j] ];
		}
		
		
		for( j = m - 2; j > 0 ; j-- ){
			for( i=0; i < n-1; i++ ){
				MY[i][j] = MX[ rord[i] ][ cord[j] ] + MY[i][j+1];
			}
		}
		for( i=1; i < n-1; i++ ){
			for( j = 1 ; j < m  ; j++ ){
				MY[i][j] = MY[i][j] + MY[i-1][j];
				// add all rows but the first
				loss = loss + MX[ rord[i+1] ][ cord[j-1] ] * MY[i][j];
			}
		}
		
		// add the first row terms
		for( j = 1 ; j < m  ; j++ ){
			loss = loss + MX[ rord[1] ][ cord[j-1] ]*MY[0][j];
		}
	}
	
	
	
		
	return loss;
}


float indcrit2(int *m0, int *rord, int *cord, int n, int m, int rev){
// hamming distance criterion	
	
	int i,j;
	int (*MX)[m] = (int (*)[m])m0;
	
	
	float loss = 0.0;
	
	//float MY[n][m];
	long M2[n][m]; // cumsum by column
	long M1[n][m]; // cumsum by row
	
	//set the first row and the last column(s)
	for( i = 0; i < n; i++ ){ 
		M1[i][m-2] = MX[ rord[i] ][ cord[m-1] ];
		M1[i][m-1] = 0;
	}
	for( j = 0; j < m; j++ ){ 
		M2[0][j] = 0;
		M2[1][j] = MX[0][ cord[j] ];
	}
	
	
	// compute M1
	if(m > 2){ 
		for( j = m - 3; j >= 0 ; j-- ){
			for( i=0; i < n; i++ ){
				M1[i][j] = MX[ rord[i] ][ cord[j+1] ] + M1[i][j+1];
			}
		}
	}
	for( j = m - 3; j >= 0 ; j-- ){
		for( i = 0; i < n; i++ ){
			M1[i][j] = M1[i][j] + M1[i][j+1];
		}
	}
	
	for( i = 1; i < n; i++ ){
		for( j= m - 2; j >= 0 ; j-- ){
			M1[i][j] = M1[i][j] + M1[i-1][j];
		}
	}
	
	// compute M2	
	if(n > 2){
		for( i=2; i < n; i++ ){
			for( j = 0; j < m; j++ ){
				M2[i][j] = MX[ rord[i-1] ][ cord[j] ] + M2[i-1][j];
			}
		}
	}
	for( i=2; i < n; i++ ){
		for( j = 0; j < m; j++ ){
			M2[i][j] = M2[i][j] + M2[i-1][j];
		}
	}
	
	for( j = m - 2; j >= 0 ; j-- ){
		for( i = 1; i < n; i++ ){
			M2[i][j] = M2[i][j] + M2[i][j+1];
		}
	}
	
	// compute the criterion
	for( i = 0; i < n; i++ ){
		for( j = 0; j < m; j++ ){
			loss = loss + (M1[i][j] + M2[i][j]) * MX[ rord[i] ][ cord[j] ];
		}
	}
	return loss;
}




SEXP getclust(SEXP M, SEXP dims, SEXP tau0){
	
	int i,j,i2,j2, ki, kj, rd, cd;
	
	int n = INTEGER(dims)[0];
	int m = INTEGER(dims)[1];
	
	int MX[n][m];
	for (j = 0; j < m; j++) {
		for (i = 0; i < n; i++) {
			MX[i][j] = INTEGER(M)[i+j*n];
		}
	}
	
	float currtau;
	float currbest;
	long  a, b, c, d;
	float nn, zz;
	int k = 0;
	int ncl = 1;
	int colcuts[m];
	int rowcuts[n];
	colcuts[0] = 0;
	rowcuts[0] = 0;
	colcuts[1] = m;
	rowcuts[1] = n;
	
	int cord[m];
	int rord[n];
	for (i=0; i<n; i++) {
		rord[i] = i;
	}
	for (j=0; j<m; j++) {
		cord[j]=j;
	}
	
	int NE[n][m];
	int NW[n][m];
	int SE[n][m];
	int SW[n][m];
	
	for (j = 0; j < m; j++) {
		NE[0][j] = MX[0][j];
		SE[n-1][j] = MX[n-1][j];
		NW[0][j] = MX[0][j];
		SW[n-1][j] = MX[n-1][j];
	}
	//for (i = 0; i < n; i++) {
		//NE[i][m-1] = MX[i][m-1];
		//SE[i][m-1] = MX[i][m-1];
		//NW[i][0] = MX[i][0];
		//SW[i][0] = MX[i][0];
	//}
	for (i = 1; i < n; i++) {
		for (j = 0; j < m; j++) {
			NW[i][j] = NW[i-1][j] + MX[i][j];
			NE[i][j] = NW[i-1][j] + MX[i][j];
			SE[n-i-1][j] = SE[n-i][j] + MX[n-i-1][j];
			SW[n-i-1][j] = SW[n-i][j] + MX[n-i-1][j];
		}
	}
	for (i = 0; i < n; i++) {
		for (j = 1; j < m; j++) {
			NW[i][j] = NW[i][j-1] + NW[i][j];
			NE[i][m-j-1] = NE[i][m-j] + NE[i][m-j-1];
			SE[i][m-j-1] = SE[i][m-j] + SE[i][m-j-1];
			SW[i][j] = SW[i][j-1] + SW[i][j];
		}
	}
	//for (i = 0; i < n; i++) {
	//	for (j = 0; j < m; j++) {
	//		Rprintf("%d\t",NW[i][j]);
	//	}
	//	Rprintf("\n");
	//}
	
	
	///////////
	// a | b //
	//-------//
	// c | d //
	///////////
	
	while (k < ncl) {
		currbest = -0.001;
		
		rd = rowcuts[k+1]-rowcuts[k];
		cd = colcuts[k+1]-colcuts[k];
		
		if( (rd > 1)  && ( cd > 1)   ){
			
			for (i2 = rowcuts[k]+1; i2 < rowcuts[k+1]; i2++) { // i2 and j2 are first in the second part
				for (j2 = colcuts[k]+1; j2 < colcuts[k+1]; j2++) {
					
					a = NW[i2-1][j2-1];
					b = NE[i2-1][j2];
					c = SW[i2][j2-1];
					d = SE[i2][j2];
					if( k > 0 ){//if (colcuts[k] > 0) { //left+top
						a = NW[i2-1][j2-1] - NW[i2-1][ colcuts[k]-1 ] - NW[ rowcuts[k]-1 ][j2-1] + NW[ rowcuts[k]-1 ][ colcuts[k]-1 ];
						c = c - SW[i2][ colcuts[k]-1 ];
						b = b - NE[ rowcuts[k]-1 ][ j2 ];
					}
					if( k+1 < ncl ){//if (colcuts[k+1] < m-1) { //right+bottom
						d = SE[i2][j2] - SE[i2][ colcuts[k+1] ] - SE[ rowcuts[k+1] ][j2] + SE[ rowcuts[k+1] ][ colcuts[k+1] ];
						c = c - SW[ rowcuts[k+1] ][j2-1]; 
						b = b - NE[ i2-1 ][ colcuts[k+1] ];
					}
					if(k > 0 && k+1 < ncl){//if (colcuts[k+1] < m-1 && colcuts[k] > 0) { //all
						c = c + SW[ rowcuts[k+1] ][ colcuts[k]-1 ];
						b = b + NE[ rowcuts[k]-1 ][ colcuts[k+1] ];
					}
					//Rprintf(" a = %d, b = %d, c = %d, d = %d\n",a,b,c,d);
					
					zz = a*d - b*c;
					nn = pow((a+b)*(c+d),0.5)*pow((a+c)*(b+d),0.5);
					//Rprintf(" zz = %f, nn = %f\n",zz,nn);
					
					currtau = zz/nn;
					
					
					if(currtau > currbest){
						currbest = currtau;
						ki = i2;
						kj = j2;
						
					}
					
				}
			}
			//Rprintf("best crit = %f at x = %d and y = %d\n",currbest,kj,ki);
			if(currbest >= REAL(tau0)[0]){
				ncl++;
				for (i = ncl; i > k+1; i--) {
					colcuts[i] = colcuts[i-1];
				}
				colcuts[k+1] = kj;
				for (i = ncl; i > k+1; i--) {
					rowcuts[i] = rowcuts[i-1];
				}
				rowcuts[k+1] = ki;
			}else {
				k++;
			}
			
		}else {
			k++;
		}

		
				
		
	}
	
	SEXP out = allocVector(REALSXP,2*k);
	for( i = 0; i < k;i++ ){
		REAL(out)[i] = rowcuts[i+1];
		REAL(out)[i+k] = colcuts[i+1];
	}
	return out;
	
}




int getindex(int dims[], int ind[], int nd){

	int i,s;
	//int nd2 = sizeof(dims)/2;
	for (i=0; i<nd; i++) {
		if (ind[i] > dims[i]-1) {
			Rprintf("invalid index!\n");
				//for (int i2=0; i2<nd; i2++) {
				//	Rprintf("%d vs. %d,\t",ind[i2],dims[i2]);
				//}
			return -1;
		}
	}
	int cpx[nd];
	cpx[0] = dims[0];
	
	for (i = 1; i < nd; i++) {
		cpx[i] = cpx[i-1]*dims[i];
	}
	
	
	int k = ind[0]; 
	
	for (s = 1; s < nd; s++) {
		//Rprintf(" %d + %d * %d\n",k,ind[s],cpx[s-1]);
		k = k + ind[s]*cpx[s-1];
	}
	return k;
} 



void diagprod(int *mv1, int *mv0, float *cv, int *IM, int dims[], int index[], int step, int nd){
	int i,j,k, s;
	
	//int nd = sizeof(dims)/2;
	int val1, val2;
	
	int (*IX)[nd] = (int (*)[nd])IM;
	float tmp, tmp2;
	//Rprintf("\n welcome to diagprod\n");
	
	//Rprintf("step = %d\n",step);
	//Rprintf("nd = %d\n",nd);
	int indS[nd];
	for (i = 0; i< nd; i++) {
		indS[i] = index[i];
	}
	
	
	if( step == nd-1 ){ // the recursion in in the last dimension
		int ind2[nd];
		//Rprintf("\n IX in diagprod:\n");
		//for (int r = 0; r < nd; r++) {
		//	for (int t = 0; t < dims[r]; t++) {
		//		Rprintf("%d\t",IX[t][r]);
		//	}
		//	Rprintf("\n");
		//}
		for (j = 0; j < nd-1; j++) {
			//Rprintf("ind %d is %d\n",j,indS[j]);
			ind2[j] = IX[indS[j]-1][j];
			indS[j] = IX[indS[j]][j];
			//Rprintf("and changes to %d\t",indS[j]);
		}
		for (s = 1; s < dims[nd-1]; s++) {
				indS[nd-1] = IX[s][nd-1];
				ind2[nd-1] = IX[s-1][nd-1];
				//Rprintf("final ind:\n");
				//for (i=0; i<nd; i++) {
				//	Rprintf("%d\t",indS[i]);
				//}
				//Rprintf("\n");
				//for (i=0; i<nd; i++) {
				//	Rprintf("%d\t",ind2[i]);
				//}
				//Rprintf("\n ind1 = %d",getindex(dims,indS,nd));
				//Rprintf("\n ind2 = %d",getindex(dims,ind2,nd));
				//Rprintf("\n val1 = %d",mv0[ getindex(dims,indS,nd) ]);
				//Rprintf("\n val2 = %d",mv1[ getindex(dims,ind2,nd) ]);
				val1 = mv0[ getindex(dims,indS,nd) ];
				val2 = (int) mv1[ getindex(dims,ind2,nd) ];
				tmp = (float) val1*val2;
				tmp2 = cv[0];
				//Rprintf("\n tmp = %f",tmp);
				cv[0] = tmp + tmp2;//((float*) val1) * ((float*) val2);
				//Rprintf("\n cv = %f",cv[0]);
		}
	}else {
		
		for (s = 1; s < dims[step]; s++) {
				indS[step] = s;
			//Rprintf("temp ind:\n");
			//for (i=0; i<nd; i++) {
			//	Rprintf("%d\t",indS[i]);
			//}
			//Rprintf("\n");
				diagprod(mv1,mv0, cv, IM, dims, indS, step+1,  nd);
		}
	}
	
}

float mvclasscrit2(int *m0, int *IM, int dims[], int nd){
		
	
	int i,ix,j, j1,j2,s, s1, s2, k, rind ,lind;
	//int nd = sizeof(dims)/2;
	//Rprintf("welcome to classcrit\n");
	//Rprintf("dim = %d\n",nd);
	
	int (*IX)[nd] = (int (*)[nd])IM;
	
	
	int cp[nd];
	int cpr[nd];
	cp[0] = 1;
	cpr[nd-1] = 1;
	
	//for (i=0; i<nd; i++) {
	//	Rprintf("%d\t",dims[i]);
	//}
	for (i = 1; i < nd; i++) {
		cp[i] = cp[i-1]*dims[i-1];
		cpr[nd-1-i] = cpr[nd-i]*dims[nd-i];//CHECK
	}
	//for (i = 0; i < nd; i++) {
	//	Rprintf("%d\t",cp[i]);
	//}
	//for (i = 0; i < nd; i++) {
	//	Rprintf("%d\t",cpr[i]);
	//}
	// m1 array for the cumsums
	int ml = cp[nd-1]*dims[nd-1];
	int m1[ ml ];
	//int m2[ ml ];
	for (i = 0; i < ml; i++) {
		m1[i] = m0[i];
		//m2[i] = m0[i];
		//Rprintf("%d\t",m1[i]);
	}
	//Rprintf("m1 raw:\n");
	//for (i=0; i<ml; i++) {
	//	Rprintf("%d\t",m1[i]);
	//}
	//go through m1 and m0 with stepsize s acc. to dimension k
	// dims is in reverse order
	for (k = 0; k < nd; k++) {
		s1 = cp[k];
		s2 = cpr[k];
		//Rprintf("s1 = %d\n",s1);
		//Rprintf("s2 = %d\n",s2);
		//if (dims[k] > 2) {
			for (i = 1; i < dims[k]; i++) {
				for (j1 = 0; j1 < s2; j1++) {
						for (j2 = 0; j2 < s1; j2++) {
							rind = IX[i][k]*s1 + j1*s1*dims[k] + j2;
							lind = IX[i-1][k]*s1 + j1*s1*dims[k] + j2;
						//Rprintf("%d\t",rind );
						//Rprintf("%d\n",lind);
						m1[ rind ] += m1[ lind ];
					}
				}
			}
			//m1[  j1*s1*dims[k] + j2 - 1] = m0[ j1*s1*dims[k] + j2 - 1];
		//} 
		
	}

	// multiplication of m0 with m1[ ind -1 ] via diagprod
	float cv = 0.0;
	int *m1p = &m1[0];
	float *cvp = &cv;
	int init[nd];
	for (i=0; i<nd; i++) {
		init[i] = 0;
	}
	diagprod(m1p,m0,cvp,IM,dims,init,0, nd);
	//Rprintf("crit = %f\n",cv);
	return cv;
}



float mvclasscrit(int *m0, int *IM, int dims[], int nd){
	
	
	int i,ix,j, j1,j2,s, s1, s2, k, rind ,lind, d;

	int (*IX)[nd] = (int (*)[nd])IM;
	
	
	int cp[nd];
	int cpr[nd];
	cp[0] = 1;
	cpr[nd-1] = 1;

	for (i = 1; i < nd; i++) {
		cp[i] = cp[i-1]*dims[i-1];
		cpr[nd-1-i] = cpr[nd-i]*dims[nd-i];//CHECK
	}

	// m1 array for the cumsums
	int ml = cp[nd-1]*dims[nd-1];
	int m1[ ml ];
	int m2[ ml ];
	//int m2[ ml ];
	for (i = 0; i < ml; i++) {
		m1[i] = m0[i];
		m2[i] = 0;
	}
	
	//go through m1 and m0 with stepsize s acc. to dimension k
	// dims is in reverse order
	
	for (d = 0; d < nd; d++) {
		
		s1 = cp[d];
		s2 = cpr[d];
		
		// first k-vector is zero
		
		for (j1 = 0; j1 < s2; j1++) {
			for (j2 = 0; j2 < s1; j2++) {
				rind = IX[0][d]*s1 + j1*s1*dims[d] + j2;
				m1[ rind ] = 0;
			}
		}
		
		// cumsum over k = d 1st
		for (i = 1; i < dims[d]; i++) {
			for (j1 = 0; j1 < s2; j1++) {
				for (j2 = 0; j2 < s1; j2++) {
					rind = IX[i][d]*s1 + j1*s1*dims[d] + j2;
					lind = IX[i-1][d]*s1 + j1*s1*dims[d] + j2;
					m1[ rind ] = m1[ lind ] + m0[ lind ];
				}
			}
		}
		
		// cumsum over other dimensions once
		// add to m2 and reset m1
		for (k = 0; k < nd; k++) {
			s1 = cp[k];
			s2 = cpr[k];
			if (d != k) {
				for (i = 1; i < dims[k]; i++) {
					for (j1 = 0; j1 < s2; j1++) {
						for (j2 = 0; j2 < s1; j2++) {
								rind = IX[i][k]*s1 + j1*s1*dims[k] + j2;
								lind = IX[i-1][k]*s1 + j1*s1*dims[k] + j2;
								m1[ rind ] += m1[ lind ];
								m2[ rind ] += m1[ lind ]; // not rind
							}
							
						}
					}
				}
			} 
		
		for (i = 0; i < ml; i++) {
			m1[i] = m0[i];
		}
	}
	float cv = 0.0;
	for (i = 0; i < ml; i++) {
		cv = cv + m2[i] * m0[i];
	}
	//Rprintf("cv = %f\n",cv);
	return cv;
}



float mvhammcrit(int *m0, int *IM, int dims[], int nd){
	
	
	int i,ix,j, j1,j2,s, s1, s2, k, d, rind ,lind;

	int (*IX)[nd] = (int (*)[nd])IM;
	
	
	int cp[nd];
	int cpr[nd];
	cp[0] = 1;
	cpr[nd-1] = 1;

	for (i = 1; i < nd; i++) {
		cp[i] = cp[i-1]*dims[i-1];
		cpr[nd-1-i] = cpr[nd-i]*dims[nd-i];//CHECK
	}
	// m1 array for the cumsums
	int ml = cp[nd-1]*dims[nd-1];
	int m1[ ml ];
	int m2[ ml ];
	//int m2[ ml ];
	for (i = 0; i < ml; i++) {
		m1[i] = m0[i];
		m2[i] = 0;
	}
	//Rprintf("\n\n ml = %d\n\n",ml);

	//go through m1 and m0 with stepsize s acc. to dimension k
	// dims is in reverse order
	
	for (d = 0; d < nd; d++) {
		
			s1 = cp[d];
			s2 = cpr[d];
			
			// first k-vector is zero
			for (j1 = 0; j1 < s2; j1++) {
				for (j2 = 0; j2 < s1; j2++) {
					rind = IX[0][d]*s1 + j1*s1*dims[d] + j2;
					m1[ rind ] = 0;
				}
			}
			
			// cumsum over k = d 1st
			for (i = 1; i < dims[d]; i++) {
				for (j1 = 0; j1 < s2; j1++) {
					for (j2 = 0; j2 < s1; j2++) {
						rind = IX[i][d]*s1 + j1*s1*dims[d] + j2;
						lind = IX[i-1][d]*s1 + j1*s1*dims[d] + j2;
						m1[ rind ] = m1[ lind ] + m0[ lind ];
					}
				}
			}
		
			// cumsum over k = d 2nd
			for (i = 1; i < dims[d]; i++) {
				for (j1 = 0; j1 < s2; j1++) {
					for (j2 = 0; j2 < s1; j2++) {
						rind = IX[i][d]*s1 + j1*s1*dims[d] + j2;
						lind = IX[i-1][d]*s1 + j1*s1*dims[d] + j2;
						m1[ rind ] += m1[ lind ];
					}
				}
			}
		
		// cumsum over other dimensions once
		// add to m2 and reset m1
		for (k = 0; k < nd; k++) {
			s1 = cp[k];
			s2 = cpr[k];
			if (d != k) {
				for (j1 = 0; j1 < s2; j1++) {
					for (j2 = 0; j2 < s1; j2++) {
						rind = IX[0][k]*s1 + j1*s1*dims[k] + j2;
						m2[ rind ] += m1[ rind ];
						for (i = 1; i < dims[k]; i++) {
							rind = IX[i][k]*s1 + j1*s1*dims[k] + j2;
							lind = IX[i-1][k]*s1 + j1*s1*dims[k] + j2;
							m1[ rind ] += m1[ lind ];
							m2[ rind ] += m1[ rind ];
						}
					}
				}
			} 
		}
		for (i = 0; i < ml; i++) {
			m1[i] = m0[i];
		}
	}
	float cv = 0.0;
	for (i = 0; i < ml; i++) {
		cv += m2[i] * m0[i];
	}
	return cv;
}





// this function will perform permutations for the categories until a local optimum is reached
SEXP mvclass(SEXP M, SEXP dims, SEXP pv , SEXP vs){
	
	
	int r,s,t,i,j,k,ii,jj,kk, tmp;
	int nd = LENGTH(dims);
	int dimv[nd];
	for (i=0; i<nd; i++) {
		dimv[i] = INTEGER(dims)[i];
	}
	int mt = INTEGER(vs)[0];
	
	int ml = dimv[0];
	int biglen = ml;
	for(i = 1; i < nd; i++){
		biglen = biglen + dimv[i];
		if( ml < dimv[i] )
			ml = dimv[i];
	}
	int IM[ml][nd];
	int TIM[ml][nd];
	
	//Rprintf("init done");
	
	for (k = 0; k < nd; k++) {
		for (i = 0; i < dimv[k]; i++) {
			IM[i][k] = i;
			TIM[i][k] = i;
			//Rprintf("%d\t",TIM[i][k]);
		}
		//Rprintf("\n");
	}
	
	
	
	int* mv = &INTEGER(M)[0];
	int* TIMp = &TIM[0][0];
	
	
	int bestmove[3];
	float bestcrit = 0;//FLT_MAX;
	float tempcrit;
	int better;
	int opt = 0;
	//Rprintf("bestcrit = %f\n",bestcrit);
	
	while (opt == 0) {
		
		better = 0;
		for (k = 0; k < nd; k++) {
			
			
			if (INTEGER(pv)[k] == 1) {
				
				for (i = 0; i < dimv[k]; i++) { // move i to j in dim k

					for (j = 0; j < dimv[k]; j++) {
						//Rprintf("k = %d\t i = %d \tj = %d\t",k,i,j);
						if (i < j) {
							for (s = i; s < j; s++) {
								TIM[s][k] = IM[s+1][k];
							}
							TIM[j][k] = IM[i][k];
						}
						if (i > j) {
							for (s = i; s > j; s--) {
								TIM[s][k] = IM[s-1][k];
							}
							TIM[j][k] = IM[i][k];
						}		
					
						if( mt == 0 ){
							tempcrit =  mvclasscrit(mv, TIMp, dimv, nd);
						}
						if( mt == 1 ){
							tempcrit =  mvhammcrit(mv, TIMp, dimv, nd);
						}
												
						//Rprintf("tempcrit = %f\n", tempcrit);
						if (tempcrit > bestcrit) {
							bestcrit = tempcrit;
							bestmove[0] = k;
							bestmove[1] = i;
							bestmove[2] = j;
							better = 1;
						//	Rprintf("better = %f -----at %d  %d  %d------------ >>><<<\n", bestcrit,k,i,j);
						}
						//reset TIM
						if (i < j) {
							for (s = j; s >= i; s--) {
								TIM[s][k] = IM[s][k];
							}
						}
						if (i > j) {
							for (s = i; s >= j; s--) {
								TIM[s][k] = IM[s][k];
							}
						}

						
						
						
					} // j
				} // i
			} // if
		} // k
		
		if(better > 0){
			// changes to IM
			ii = bestmove[1];
			jj = bestmove[2];
			kk = bestmove[0];
			//Rprintf("permutation dim = %d, %d to %d\n ------->>><<<----->>><<<----#+###+#",kk,ii,jj);
			tmp = IM[ii][kk];
			if (ii < jj) {
				for (s = ii; s < jj; s++) {
					IM[s][kk] = IM[s+1][kk];
					TIM[s][kk] = IM[s][kk];
				}
				IM[jj][kk] = tmp;
				TIM[jj][kk] = IM[jj][kk];
			}
			if (ii > jj) {
				for (s = ii; s > jj; s--) {
					IM[s][kk] = IM[s-1][kk];
					TIM[s][kk] = IM[s][kk];
				}
				IM[jj][kk] = tmp;
				TIM[jj][kk] = IM[jj][kk];
			}
			
			//Rprintf("\n The new orders are:\n");
			//for (r = 0; r < nd; r++) {
			//	for (t = 0; t < dimv[r]; t++) {
			//		Rprintf("%d, ",IM[t][r]);
			//	}
			//	Rprintf("\n");
			//}
			
			
		}else{
			opt = 1;	
		}	
	}
	
	SEXP out = allocVector(REALSXP, biglen + 1);
	
	k = 0;
	for(i = 0; i < nd; i++){
		for (j = 0; j < dimv[i]; j++) {
			REAL(out)[k] = IM[j][i];
			k++;
		}
	}
	REAL(out)[biglen] = bestcrit;
	
	return out;
}








// ----------------------------------------------------------------------------------

static void quickSort (float *a, int *order, int lo, int hi)
{
	//  lo is the lower index, hi is the upper index
	//  of the region of array a that is to be sorted
    int i=lo, j=hi, ho, hi2,lo2;
	float atmp;
	
    float x=a[(lo+hi)/2];
	//	Rprintf("from %d to %d with pivot a[%d] = %f\n\n",i,j,(lo+hi)/2,x);
    //  partition
    do
    {    
        while (a[i]<x) i++; 
        while (a[j]>x) j--;
        if (i<=j)
        {
            atmp=a[i]; a[i]=a[j]; a[j]=atmp;
			ho = order[i]; order[i]=order[j]; order[j]=ho;
            i++; j--;
			//Rprintf("ex %d -- %d\n\n",i,j);
        }
    } while (i<=j);
	hi2 = j;
	lo2 = i;
    //  recursion
    if (lo<j) quickSort(a, order, lo, hi2);
    if (i<hi) quickSort(a, order, lo2, hi);
	
}




// this function will perform permutations for the categories until a local optimum is reached
SEXP preclass(SEXP M, SEXP dims, SEXP pv ){
	
	
	int r,s,t,i,j,k,ii,jj,kk, tmp, d;
	int nd = LENGTH(dims);
	int dimv[nd];
	for (i=0; i<nd; i++) {
		dimv[i] = INTEGER(dims)[i];
	}
	
	
	int ml = dimv[0];
	int biglen = ml;
	for(i = 1; i < nd; i++){
		biglen = biglen + dimv[i];
		if( ml < dimv[i] )
			ml = dimv[i];
	}
	int IM[ml][nd];
	int CIM[ml][nd]; // current IM
	int TIM[ml][nd]; // a dummy order matrix
	
	int* CIMp = &CIM[0][0];
	int* TIMp = &TIM[0][0];
	int* IX = &IM[0][0];
	
	float cvv[ml];
	float *cvp = &cvv[0];
	
	//Rprintf("init done");
	
	for (k = 0; k < nd; k++) {
		for (i = 0; i < dimv[k]; i++) {
			IM[i][k] = i;
			CIM[i][k] = i;
		}
	}
	
	for (k = 0; k < nd; k++) {
		for (i = 0; i < ml; i++) {
			TIM[i][k] = i;
		}
	}
	
	int* mv = &INTEGER(M)[0];
		
	float bestcrit = 0;//FLT_MAX;
	float newcrit, oldcrit;
	
	int opt = 0;
	//Rprintf("bestcrit = %f\n",bestcrit);
	
	
	
	int cp[nd];
	int cpr[nd];
	cp[0] = 1;
	cpr[nd-1] = 1;
	
	for (i = 1; i < nd; i++) {
		cp[i] = cp[i-1]*dimv[i-1];
		cpr[nd-1-i] = cpr[nd-i]*dimv[nd-i];//CHECK
	}
	// m1 array for the cumsums
	int mln = cp[nd-1]*dimv[nd-1];
	int m1[ mln ];
	int m1i[ mln ];
	
	int m0[ mln ];
	
	for (i = 0; i < mln; i++) {
		m0[i] = INTEGER(M)[i];
	}
	
	for (i = 0; i < mln; i++) {
		m1[i] = m0[i];
		m1i[i] = 0;
	}
	
	int* m1ip = &m1i[0];
	
	
	
	int ds, s1, s2, s1t, s2t, j1, j2, rind, lind;
	int tmpdimv[nd];
	int iweight;
	int ord[ml];
	int *ordp = &ord[0];
	
	while (opt == 0) {
		
		for (d = 0; d < nd; d++) {
			
			
			if (INTEGER(pv)[d] == 1) {
					
					s1 = cp[d];
					s2 = cpr[d];
					
					// first k-vector is zero
					for (j1 = 0; j1 < s2; j1++) {
						for (j2 = 0; j2 < s1; j2++) {
							rind = IM[0][d]*s1 + j1*s1*dimv[d] + j2;
							m1[ rind ] = m0[ rind ];
						}
					}
								
				
					// sum over k = d => the full sums are in i = dims[d]-1
					for (i = 1; i < dimv[d]; i++) {
						for (j1 = 0; j1 < s2; j1++) {
							for (j2 = 0; j2 < s1; j2++) {
								rind = IM[i][d]*s1 + j1*s1*dimv[d] + j2;
								lind = IM[i-1][d]*s1 + j1*s1*dimv[d] + j2;
								m1[ rind ] = m1[ lind ] + m0[ rind ];
							}
						}
					}
				
				
				//writing the sums to m1i
				rind = 0;
				for (j1 = 0; j1 < s2; j1++) {
					for (j2 = 0; j2 < s1; j2++) {
						//rind = 0*s1 + j1*s1*dimv[d] + j2;
						lind = IM[dimv[d]-1][d]*s1 + j1*s1*dimv[d] + j2;
						m1i[ rind*2 ] = m1[ lind ];
						//Rprintf("val %d from %d to %d\n", m1[ lind ],lind,rind);
						rind++;
					}
				}
				
				//temporary dimv
				ds = 1;
				tmpdimv[0] = 2;
				for (s = 0; s < 2; s++) {
					TIM[s][0] = s;
				}
				for (k = 0; k < nd; k++) {
						if (k != d) {
							for (s = 0; s < dimv[k]; s++) {
								TIM[s][ds] = IM[s][k];
							}
							tmpdimv[ds] = dimv[k];
							ds++;
						}
				}
							
				s1 = cp[d];
				s2 = cpr[d];
				
				for (i = 0; i < dimv[d]; i++) {
					rind = 0;
					// writing category i to m1i
					iweight = 0;
					for (j2 = 0; j2 < s1; j2++) {
						for (j1 = 0; j1 < s2; j1++) {
							//rind = 1*s1t + j1*s1t*dimv[d] + j2;
							lind = IM[i][d]*s1 + j1*s1*dimv[d] + j2;
							m1i[ rind*2+1 ] = m0[ lind ];
							iweight += m0[ lind ];
						//	Rprintf("val %d from %d to %d\n", m0[ lind ],lind,rind);
							rind++;
						}
					}
					//Rprintf("iweight = %d\n",iweight);
					//Rprintf("sums:\n");
					//for (s = 0; s < dimv[1-d]; s++) {
					//	for (j = 0; j < 2; j++) {
					//		Rprintf("%d\t\t", m1i[j+s*2]);
					//	}
					//	Rprintf("\n");
					//}
										
					cvv[i] = mvclasscrit(m1ip, TIMp, tmpdimv, nd) / iweight;
					//Rprintf("cvv[%d] = %f\n",i,cvv[i]);
				}	
				
				
			}
			for (s = 0; s < dimv[d]; s++) {
				ord[s] = CIM[s][d];
			}
			
			quickSort(cvp,ordp,0,dimv[d]-1);	
			
			//for (s = 0; s < dimv[d]; s++) {
			//	Rprintf("%d\t",ord[s]);
			//}
			//Rprintf("check06\n\n");
		
			//for (s = 0; s < dimv[d]; s++) {
			//	Rprintf("%f\t",cvv[s]);
			//}
			for (s = 0; s < dimv[d]; s++) {
				 CIM[s][d] = ord[s];
			}
			
			//Rprintf("\n\n");
			//for (s = 0; s < dimv[0]; s++) {
			//	for (j = 0; j < dimv[1]; j++) {
			//		Rprintf("%d \t",mv[  IM[s][0] + IM[j][1]*dimv[0]]);
			//	}
			//	Rprintf("\n");
			//}
			
			//Rprintf("\n\n");
			//for (s = 0; s < dimv[0]; s++) {
			//	for (j = 0; j < dimv[1]; j++) {
			//		Rprintf("%d \t",mv[  CIM[s][0] + CIM[j][1]*dimv[0]]);
			//	}
			//	Rprintf("\n");
			//}
			
			oldcrit = mvclasscrit(mv, IX, dimv, nd);
			//Rprintf("old crit = %f", oldcrit);
			newcrit = mvclasscrit(mv, CIMp, dimv, nd);
			if (newcrit > bestcrit) {
				bestcrit = newcrit;
				//Rprintf("new best crit = %f", bestcrit);
				for (s = 0; s < dimv[d]; s++) {
					IM[s][d] = CIM[s][d];
				}
			}else {
				//Rprintf("newcrit = %f stop", newcrit);
				opt = 1;
			}
		}
	}
				
	SEXP out = allocVector(REALSXP, biglen + 1);
	
	k = 0;
	for(i = 0; i < nd; i++){
		for (j = 0; j < dimv[i]; j++) {
			REAL(out)[k] = IM[j][i];
			k++;
		}
	}
	REAL(out)[biglen] = bestcrit;
	
	return out;
}



////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////



