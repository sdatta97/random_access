/******************************************************************************

Welcome to GDB Online.
GDB online is an online compiler and debugger tool for C, C++, Python, PHP, Ruby, 
C#, VB, Perl, Swift, Prolog, Javascript, Pascal, HTML, CSS, JS
Code, Compile, Run and Debug online from anywhere in world.

*******************************************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<math.h>
int main() 
	{
		double  V0[2], V1[2][193], V2[2][295][295], V11[2][193][193],V12[2][193][193], Vsub1[2][193][193], V122[2][193][193], V22[2][193][193];					// Initilization of value arrays
		double alph[193],beta[295][295]; 								// Initialization of decision arrays

		int NUM_ITERATIONS = 100;
             
		double vmax, tmax1, tex, t, tmax, vn, vn1, vn2, tdif, tdif1, aux, xmin, xmax, x,dif,dif1;
		double b1max, bxdif, bxdif1, exy, bmax, bmax1, xy, ex, ex2, exb, exyb, axdif, axdif1, a, ea, ea2, eax, eax2, amax1, amax, adif, adif1, b, y, ey, eb, eyb, bdif, bdif1;
		
		int i, j, k, l, m, n, it, ib, iby, lb, ia, iax, la,k1, kb, ibxy, ibx, kl;
		
		int m1 = 1 ;												// Precision parameters
		int m2 = m1 + 1;
		//double del = 0.01;
		double del = 0.005;
		double delf = 0.0005;
		double del1 = del*m1;
		double delf1 = delf*m1;        
        int flag = 0;												// Flag parameters
		int type = 0;

		V1[0][0] = 0;												// Array initializations
		V0[1] = 0;
		
		for(j=0; j<193; j++) V1[1][j] = 0;
			
		for(j=0; j<295; j++) {
			for(k=0; k<295; k++) {
				V2[1][j][k] = 0;
			}
		}
		for(j=0; j<193; j++) {
			for(k=0; k<193; k++) {
				V12[1][j][k] = 0;
				V11[1][j][k] = 0;
				Vsub1[1][j][k] = 0;	
				V122[1][j][k] = 0;
				V22[1][j][k] = 0;
			}
		}
	
		for(n=1; n<= NUM_ITERATIONS; n++) { 						// Main iteration loop 
			
			i = 1 - (n-(n/2)*2);									// Toggle between current and previous value estimates
			j = 0 + n - (n/2)*2;
			//printf("(%d,%d)\n", i, j);
			//Calculation of V0 and T
		
			vmax = 0;
			tmax1 = 0.99;
			//tmax1 = 0.099;
			//Coarse Search for T
			
			for(l=1; l<=94; l++) {
				
				t = tmax1 + l*del;
				tex = exp(-t);
				it  =  200 + l;
				
				vn = t*tex + (1 + t)*tex*V0[j] + V2[j][it-1][0]*(1- (1 + t)*tex);			// Value expression for V0
			
				if(vn>vmax){
					vmax = vn;										// Update values leding to maximum value
					tmax = t;
				}
			}
			printf("Course S0: %f, %f\n", vmax, tmax);
			// Fine search for T
            
			for(l=1; l<=19; l++) {
				
				//t = tmax - 0.01 + l*delf;
				t = tmax - del+ l*delf;
				
				if (t<0) continue;
				tex = exp(-t);
				aux = t/del;										// Generating values for linear interpolation
				it = aux;
				dif = aux - it;
				dif1 = 1 - dif;
				it = it + 1;
				//printf("%f, %d, %f, %f\n", aux,it,dif,dif1);
				vn = t*tex + (1 + t)*tex*V0[j] + (V2[j][it-1][0]*dif1 + V2[j][it][0]*dif)*(1 - (1 + t)*tex); // Value expression for V0 with interpolation
			
				if(vn>vmax) {
					tmax = t;										// Update values leding to maximum value
					vmax = vn;
				}
			}
		
			V0[i] = vmax;											// Update V0 for next iteration
			t = tmax;
			printf("Fine S0: %f, %f\n", vmax,tmax);
            
			x = V0[i] - V0[j];										// Update change in value for state
			xmin = x;
			xmax = x;
			printf("Value change S0: %f\n",x);
			//Calculate V2[0][y] and B[0][y], B>=0
		
			for(l=m2; l<=193; l=l+m1){
			
				y = (l-1)*del;
				ey = exp(-y);
				vn = 1 + V1[j][l-1];
				b = 0;
				bmax = 0;
				vmax = vn;

				// Fine search for B[0][y]
			    
				for(b=delf1;b<y;b=b+delf1){
				
					eb = exp(-b);
					eyb = ey/eb;
					aux = b/del1;									// Generating values for linear interpolation
					ib = aux;
					dif = aux - ib;
					dif1 = 1 - dif;
					iby = l - (ib + 1)*m1;
					ib = ib*m1 + 1;
					vn = (eb-ey)*(1 + V1[j][iby-1]*dif + V1[j][iby-1+m1]*dif1) + (1 - eb)*V2[j][0][ib-1+m1]*dif; // From value expression of VS2(0,y) with linear interpolation
					if(ib==1) {
						vn = (vn + (1-eb)*(1+V1[i][0])*bdif1)/(1-ey);				// Substituting value expression lim y -> 0 VS2(0,y) with linear interpolation
					}
				
					if(ib>1) {
						vn = (vn + (1-eb)*V2[j][0][ib-1]*bdif1)/(1-ey);				// From value expression of VS2(0,y) with linear interpolation
					}
				
					if(vn>vmax) {
						vmax = vn;									// Update values leding to maximum value
						bmax = b;
					}
				}
	
				beta[0][l-1] = bmax;
				V2[i][0][l-1] = vmax;								// Update V2 for next iteration
			    
				x = V2[i][0][l-1] - V2[j][0][l-1];					// Update change in value for state
			
				if(x<xmin) 	xmin = x;
				if(x>xmax) 	xmax = x;
				if (xmin < 0) {
    			    printf("V2[0]: Neg xmin: %f\n", xmin);
    			}
				
			}
		

			//Calculate V1[0] and alph[0]
		
			vmax = 0;
			//vn = 1 + V0[j];
			//vmax = vn;
			//Coarse Search for a
		
			for(l=1;l<=193;l=l+m1) {
				
				a = (l-1)*del;
			    if (a<0) continue;
				//if (a>0) break;
				ea = exp(-a);
				vn = ea*(1 + V0[j]) + Vsub1[j][0][l-1];			// From value expression VS1(0)
				if(vn>vmax) {										// Update values leding to maximum value
					amax = a;
					vmax = vn;
				}			
			}
			// (n+1) -> i, n->j
			//Fine Search for alpha
		    
			for(l=1; l<=19; l++) {

				a = amax - del1 + l*delf1;
			
				if(a<0) continue;
				//if (a>0) break;	
				ea = exp(-a);										// Generating values for linear interpolation
				aux = a/del1;
				ia = aux;
				dif = aux-ia;
				dif1 = 1 - dif;
				ia = ia*m1 + 1;
			
				vn = ea*(1+V0[j]) + Vsub1[j][0][ia + m1 - 1]*adif;		// From value expression VS1(0) and interpolation
				if(ia==1) {
					//vn = vn + (1-ea)*(1 + V1[j][0])*adif1;					// Substituting value expression lim y -> 0 VS2(0,y) with linear interpolation 
					vn = vn + 0;
				}
			
				if(ia>1) {
					vn = vn + Vsub1[j][0][ia-1]*adif1;					// From value expression VS1(0) and interpolation
				}
			
				if(vn>vmax) {
					amax = a;										// Update values leding to maximum value
					vmax = vn;
				}
			}
		
			alph[0] = amax;
			V1[i][0] = vmax;										// Update V1 for next iteration
			
			x = V1[i][0] - V1[j][0];								// Update change in value for state
			
			if(x<xmin) 	xmin = x;
			if(x>xmax) 	xmax = x;
			if (xmin < 0) {
			    printf("V1[0]: Neg xmin: %f\n", xmin);
			}

			//Calculate Vsub1[0][y], V122[0][y/2][y], V22[0][y], V11[0][y] and V12[0][y], y>0
			for(la = m2; la<=193; la = la+m1) {
				a =(la-1)*del;
				ea = exp(-a);
				ea2 = exp(-0.5*a);
				V122[i][0][la-1] = 1+V2[j][(la-1)/2][(la-1)/2];
				V22[i][0][la-1] = 0.5*(1+V12[j][0][la-1]+V22[j][0][la-1]);
				V12[i][0][la-1] = (ea2*(1-(1+a/2)*ea2)*(1+V2[j][(la-1)/2][0]) + (1-ea2-(a/2)*ea)*V122[j][0][la-1])/(1-(1+a)*ea);						// From value expression VS12(x)
				V11[i][0][la-1] = (a*ea*(1+V1[j][0]) + (1-(1+a)*ea)*V12[j][0][la-1])/(1-ea);
				Vsub1[i][0][la-1] = a*ea*(1+V1[j][0]) + (1-(1+a)*ea)*V12[j][0][la-1];
			}

			//Calculate Vsub1[x][0], V122[x][0][0], V22[x][0], V11[x][0] and V12[x][0], x>0
			for(l = m2; la<=193; l = l+m1) {
				x =(l-1)*del;
				ex = exp(-x);
				ex2 = exp(-x/2);
				V122[i][l-1][0] = (0.75*x*ex*(1+V2[j][0][0])+(1-(1+x)*ex)*V22[j][l-1][0])/(0.75*x*ex+(1-(1+x)*ex));
				V22[i][l-1][0] = (ex2*(1-(1+x/2)*ex2)*V22[j][(l-1)/2][0]+(x/2)*ex2*(1-ex2)*(1+V12[j][(l-1)/2][0])+(1-(1+x/2)*ex2)*V22[j][(l-1)/2][0])/(1-(1+x)*ex);
				Vsub1[i][l-1][0] = (1-(1+x)*ex)*V2[j][l-1][0]; 
				V12[i][l-1][0] = (0.25*x*ex*(1+V2[j][0][0])+(0.75*x*ex+1-(1+x)*ex)*V122[j][l-1][0])/(1-ex);						// From value expression VS12(x)
				V11[i][l-1][0] = 1+V1[j][l-1];
			}

			V122[i][0][0] = 1+V2[j][0][0];
			V22[i][0][0] = 0.5*(1+V12[j][0][0]+V22[j][0][0]);
			V12[i][0][0] = 0.25*(1+V2[j][0][0])+0.75*V122[j][0][0];						// From value expression VS12(0)
			V11[i][0][0] = 1+V1[j][0];	
			Vsub1[i][0][0] = 0; //1+V1[j][0]; 	
			//Calculate V1[x] and alph[x]
		
			for(l=m2; l<=193; l=l+m1) {
			
				x = (l-1)*del;
				ex = exp(-x);
				vmax = 0;
		
				//Coarse Search for alph[x]
			
				for(la = m2; la<=193; la = la+m1) {
				
					a =(la-1)*del;
					if (a<0) continue;
					//if (a>x) break;
					ea = exp(-a);
                    eax = ea/ex;
					if(a<=x) {
						vn = (ea-ex)*V1[j][l-la] + a*ea*(1+V0[j]) + (1-(1+a)*ea)*V2[j][la-1][0];	// From value expression VS1(x)
					}
				
					if(a>x) {
						vn = x*ea*(1 + V0[j]) + Vsub1[j][l-1][la-l];	
					}
				
					if(vn>vmax) {
						amax = a;
						vmax = vn;
					}
				}
			
				//Fine search for alph[x]
			
				for(la = 1; la<=19; la++) {
					
					a = amax - del1 + la*delf1;
					if (a<0) continue;
					//if (a>x) break;
					ea = exp(-a);
				
					if(a>x) {
						aux = (a-x)/del1;								// Generating values for linear interpolation
						iax = aux;
						dif = aux-iax;
						dif1 = 1 - dif;
						iax = iax*m1 + 1;
						vn = x*ea*(1+V0[j]) + Vsub1[j][l-1][iax-1]*dif1 + Vsub1[j][l-1][iax+m1-1]*dif;		// From value expression VS1(x) and interpolation
					}
				
					if(a<=x) {
					
						aux = a/del1;
						ia = aux;
						dif = aux- ia;
						dif1 = 1-dif;
						iax = l - (ia+1)*m1;
						ia = ia*m1 + 1;
						vn = (ea-ex)*(V1[j][iax-1]*dif + V1[j][iax+m1-1]*dif1) + a*ea*(1+V0[j]) + (1-(1+a)*ea)*(V2[j][ia-1][0]*dif1 + V2[j][ia+m1-1][0]*dif);	// From value expression VS1(x) and interpolation
					}
				
					if(vn>vmax) {
						amax = a;										// Update values leding to maximum value
						vmax = vn;
					}
				}
			
				alph[l-1] = amax;
				V1[i][l-1] = vmax/(1-ex);								// Update V1 for next iteration
			
				x = V1[i][l-1] - V1[j][l-1];							// Update change in value for state
			
				if(x<xmin) 	xmin = x;
				if(x>xmax) 	xmax = x;
				if (xmin < 0) {
    			    printf("Neg xmin: %f\n", xmin);
    			}
			}	

			//Calculate Vsub1[x][y], V122[x][y/2][y], V22[x][y], V11[x][y] and V12[x][y]
			for(l=m2; l<=193; l=l+m1) {		
				x = (l-1)*del;
				ex = exp(-x);
				ex2 = exp(-0.5*x);
				for(la = m2; la<=193-l; la = la+m1) {
					a =(l+la-2)*del;
					ea = exp(-a);
                    eax = ea/ex;
					eax2 = exp(-0.5*(a-x));
					V122[i][l-1][la-1] = (x*ex*(1-eax2-0.5*(a-x)*eax)*(1+V2[j][(la-1)/2][(la-1)/2])+(1-(1+x)*ex)*(1-(1+a-x)*eax)*V22[j][l-1][la-1])/(x*ex*(1-eax2-0.5*(a-x)*eax)+(1-(1+x)*ex)*(1-(1+a-x)*eax));
					V22[i][l-1][la-1] = (ex2*(1-(1+x/2)*ex2)*V22[j][(l-1)/2][la-1]+0.5*x*ex2*(1-ex2)*(1+V12[j][(l-1)/2][la-1])+(1-(1+0.5*x)*ex2)*V22[j][(l-1)/2][la-1])/(1-(1+x)*ex);
					V12[i][l-1][la-1] = (x*ex*eax2*(1-(1+0.5*(a-x))*eax2)*(1+V2[j][(la-1)/2][0])+ (x*ex*(1-eax2-0.5*(a-x)*eax)+(1-(1+x)*ex)*(1-(1+(a-x))*eax))*V122[j][l-1][la-1])/((1-ex)*(1-(1+a-x)*eax));						// From value expression VS12(x)
					V11[i][l-1][la-1] = ((a-x)*eax*(1+V1[j][l-1]) + (1-(1+a-x)*eax)*V12[j][l-1][la-1])/(1-eax);
					Vsub1[i][l-1][la-1] = (1-(1+x)*ex)*eax*V2[j][l-1][0] + (1-ex)*(a-x)*eax*(1+V1[j][l-1]) + (1-(1+a-x)*eax)*(1-ex)*V12[j][l-1][la-1];
				}
			}
			//Calculate V2[0][0]
		
			V2[i][0][0] = (V2[j][0][0] + V1[j][0] + 1)/2;				// Update V2 for next iteration and limiting expression lim x-> 0 VS2(x,0)
			
			x = V2[i][0][0] - V2[j][0][0];								// Update change in value for state
			
			if(x<xmin) 	xmin = x;
			if(x>xmax) 	xmax = x;
			if (xmin < 0) {
			    printf("Neg xmin: %f\n", xmin);
			}
			
			//Calculate V2[x][y] and B[x][y], x>0
		
			for(k=m2; k<=193; k = k+m1) {
			
				x = (k-1)*del;	
				ex = exp(-x);

				for(l=1;l<=193-k;l=l+m1) {
				
					y = (l-1)*del;
					ey = exp(-y);
					exy = ex*ey;
				
					b= 0;
					vn = 0;
					bmax = 0;
					b1max = 0;
					vmax = vn;
					flag = 0;
					type = 0;
				
					for(b=b+delf1; b<(x+y); b=b+delf1) {
					
						if(b>=x) {
							flag = 1;
							break;
						}
					
						eb = exp(-b);// Generating values for linear interpolation
						exb = ex/eb;
						exyb = ex*ey/eb;
						aux = b/del1;
						ib = aux;
						dif = aux - ib;
						dif1 = 1 - dif;
						ib = ib*m1 + 1;
						ibxy = k + l - ib - m1;
						ibx = k - ib - m1 + 1;
					
						// Value expression for V2(x,y)
						vn = (eb - ex - (x-b)*exy)*(V2[j][ibx-1][l-1]*dif + V2[j][ibx + m1 - 1][l-1]*dif1) + b*(eb - exy)*(1 + V1[j][ibxy-1]*dif + V1[j][ibxy + m1 -1]*dif1) + (1-(1+b)*eb)*(V2[j][ib-1][0]*dif1 + V2[j][ib + m1 - 1][0]*dif);
						if(vn>vmax) {									// Update values leding to maximum value
							vmax = vn;
							bmax = b;
						}
					}	
				
					// Now considering the case when B> x

					if(flag==1) {
						// Coarse search for B
						kl = k + l - m2;
						b1max = x;
					
						for(m = k; m<=kl; m = m + m1) {
						
							b = (m-1)*del;
							eb = exp(-b);
							exyb = ex*ey/eb;
							ibx = m - k + 1;
							ibxy = k + l - m;
						    
							vn = x*(eb-exy)*(1 + V1[j][ibxy-1]) + (1-ex - x*ex)*V2[j][k-1][ibx-1];	// Value expression for V2(x,y)
							if(vn>vmax) {								// Update values leding to maximum value
								b1max = b;
								vmax = vn;
								type = 1;
							}
						}
						// Fine search for B
						
						for(m = 1; m<=19; m++) {
							
							b = b1max - del1 + m*delf1;
						
							if(b<x) continue;
							
							eb = exp(-b);
						    exyb = ex*ey/eb;
							aux = b/del1;								// Generating values for linear interpolation
							ib = aux;
							bdif = aux - ib;
							bdif1 = 1 - bdif;
							ib = ib*m1 + 1;
							ibxy = k + l - ib - m1;
							ibx = ib - k + 1;
							// Value expression for V2(x,y) with interpolation
							vn = x*(eb - exy)*(1 + V1[j][ibxy-1]*bdif + V1[j][ibxy + m1 -1]*bdif1) + (1-ex-x*eb)*(V2[j][k-1][ibx-1]*bdif1 + V2[j][k-1][ibx + m1 -1]*bdif);
							if(vn>vmax) {								
								b1max = b;								// Update values leding to maximum value
								vmax = vn;
								type = 1;
							}
						}
						
					}
					
					if(type==0) b=bmax;
						else b=b1max;
								
					beta[k-1][l-1] = b;
					V2[i][k-1][l-1] = vmax/(1 - ex - x*exy);				// Update V2 for next iteration
					
					xy = V2[i][k-1][l-1] - V2[j][k-1][l-1];					// Update change in value for state
					if(xy<xmin) xmin = xy;
					if(xy>xmax) xmax = xy;
					if (xmin < 0) {
					    printf("curr V2 = %f, prev V2 = %f, k = %d, l = %d\n", V2[i][k-1][l-1], V2[j][k-1][l-1], k, l);
        			    printf("Neg xmin: %f\n", xmin);
        			}
				}
			}
				
			//Calculate v2[x][0] and B[x], x>0.96
		    
			for(l=194; l<=295; l++) {
				
				x = del*(l-1);
				ex = exp(-x);
				vn = 0;
				vmax = vn;

				for(b=b;b<=x;b=b+delf1) {
		
					eb = exp(-b);	// Generating values for linear interpolation
					exb = ex/eb;
					aux = b/del1;
					ib = aux;
					bdif = aux - ib;
					bdif1 = 1 - bdif;
					ib = ib*m1 + 1;
					aux = (x-b)/del1;
					ibx = aux;
					bxdif = aux - ibx;
					bxdif1 = 1 - bxdif;
					ibx = ibx*m1 + 1;
					
					// Value expression for V2(x,y) with interpolation
					vn = (eb - (1 + x- b)*ex)*(V2[j][ibx -1][0]*bxdif1 + V2[j][ibx+m1-1][0]*bxdif) + b*(eb-ex)*(1+V1[j][ibx-1]*bxdif1 + V1[j][ibx+m1-1]*bxdif) + (1-(1+b)*eb)*(V2[j][ib-1][0]*bdif1 + V2[j][ib+m1-1][0]*bdif);
					if(vn>vmax) {
						vmax = vn;											// Update values leding to maximum value
						bmax = b;	
					}
				}
			
				b = bmax;
				beta[l-1][0] = b;
				V2[i][l-1][0] = vmax/(1-(1+x)*ex);							// Update V2 for next iteration
			
				x = V2[i][l-1][0] - V2[j][l-1][0];							// Update change in value for state
				
				if(x<xmin) xmin = x;										// Update change in value for state
				if(x>xmax) xmax = x;
				if (xmin < 0) {
				    printf("curr V2 = %f, prev V2 = %f, l = %d\n", V2[i][l-1][0], V2[j][l-1][0], l);
    			    printf("ibx=%d, ib =%d\n", ibx, ib);
    			    printf("Neg xmin: %f\n", xmin);
    			}
			}
		}	

	    FILE *out_file  = fopen("out_nfcfs.txt", "w+"); 
        // test for files not existing. 
   	    if (out_file == NULL) 
        {   
            printf("Error! Could not open file\n"); 
            exit(-1); // must include stdlib.h 
        } 
		for(i=0;i<=192;i++) {
			fprintf(out_file,"%lf, %lf, %lf\n", i*del, V1[1][i], alph[i]);
		}
		
		printf("%.10lf-%.10lf=%.10lf\n", xmax, xmin, xmax-xmin);
	}