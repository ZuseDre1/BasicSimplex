import java.io.*;
import java.util.*;
import java.text.ParseException;
import org.apache.commons.lang3.ArrayUtils;
import Jama.Matrix;
import java.math.*;
import org.apache.commons.math3.fraction.*;


public class fracLPSolver {
	
	public fracLP currentLP;
	private BigFraction[] ibfs;
	private int[]bIndices;
	
	private fracSimplexRoutines SR;
	
	
	
	public fracLPSolver(fracLP lp_to_solve){
		
		this.currentLP=lp_to_solve;
		SR=new fracSimplexRoutines();
		
	}
	
	public fracLPSolver(fracLP lp_to_solve,BigFraction[] ibfs,int[]bIndices ){
		
		this.currentLP=lp_to_solve;
		this.bIndices=bIndices;
		this.ibfs=ibfs;
		SR=new fracSimplexRoutines();
		
	}




	
	
	
	
	
	public void solveLP(){
	
		if (this.ibfs==null){
			phase1();
		}
		else
			phase2();
	}
	
	
	
	
	
	
	
	
	
	
	
	/**
	Phase 1 of the revised simplex algorithm brings the LP into
	standard form and attemps to find an initial 
	basic feasible solution. Stores the standard-form problem in 
	the converted LPD in LP.
	*/
	public void phase1(){
		
		int slacks=SR.countSlacks(this.currentLP);
		int splits=SR.countUnbounded(this.currentLP);
		
		BigFraction[][] A= this.currentLP.oLPD.A;
		BigFraction[] c_orig = this.currentLP.oLPD.c;
		BigFraction[] b=this.currentLP.oLPD.b;
		int[] A_sense=this.currentLP.oLPD.A_sense;
		
		int new_n=A[0].length+slacks+splits;
		int new_m=A.length;
		
		BigFraction[][] A_new=new BigFraction[new_m][new_n];
		BigFraction[] c_new=new BigFraction[new_n];
		BigFraction[] b_new=b;
		
		
		int filled=0;
		int splitter=-1;
		int factor;
		
		int idx=0;
		for(int i=0;i<c_orig.length;i++){
			c_new[idx]=c_orig[i];
			if (this.currentLP.lowerBoundVector[i]<0){
				idx+=1;
				c_new[idx]=c_orig[i].multiply(splitter);
			}
			idx+=1;
		}
		
		
		
		
		
		
		filled=0;
		
		//fill in new matrix
		for(int i=0;i<A_new.length;i++){
			if (b[i].doubleValue()<0){
				factor=-1;
				b_new[i]=b_new[i].multiply(factor);
			}
			else
				factor=1;
			BigFraction[] old_coeffs=A[i];
			idx=0;
			for(int j=0;j<A[0].length;j++){
				
				A_new[i][idx]=old_coeffs[j].multiply(factor);
				
				if (this.currentLP.lowerBoundVector[j]<0){
					idx+=1;
					A_new[i][idx]=old_coeffs[j].multiply(factor).multiply(splitter);
				}
				idx+=1;
			}
			//eventually add slack
			if(A_sense[i]==1){
				A_new[i][idx+filled]=new BigFraction(-1);
				filled+=1;
			}
			if(A_sense[i]==-1){
				A_new[i][idx+filled]=new BigFraction(1);
				filled+=1;
			}
				
		}
		
		
		
		//save Problem in its standard form to cLPD in LP
		this.currentLP.cLPD=new fracLPDescriptor(A_new,b_new,c_new,new_n,new_m);
		BigFraction[][] A_a= this.currentLP.cLPD.A;
		BigFraction[] c_std = this.currentLP.cLPD.c;
		BigFraction[] b_a=this.currentLP.cLPD.b;
		int n=this.currentLP.cLPD.n;
		int m=this.currentLP.cLPD.m;
		
		
		
		
		
		
		
		
		
		//create auxiliary problem
		idx=new_n;
		BigFraction[] aux_b=b_new;
		BigFraction[][] aux_A=new BigFraction[new_m][new_n+new_m];
		BigFraction[] aux_c=new BigFraction[c_new.length+new_m];
		BigFraction[] aux_ibfs=new BigFraction[aux_c.length];
		int[] aux_bIndices=new int[new_m];
		
		
		
		for(int i=0;i<new_m;i++){
			aux_c[c_new.length+i]=new BigFraction(1);
			aux_ibfs[new_n+i]=aux_b[i];
		}
		
		for(int i=0;i<A_new.length;i++){
			for(int j=0;j<A_new[0].length;j++){
				aux_A[i][j]=A_new[i][j];
			}
			aux_A[i][idx]=new BigFraction(1);
			idx+=1;
		}
		for(int i=0;i<aux_bIndices.length;i++){
			aux_bIndices[i]=aux_c.length-new_m+i;
		}
		
		/*
		System.out.println("Aux_A"+Arrays.deepToString(aux_A));
		System.out.println("Aux_c"+Arrays.toString(aux_c));
		System.out.println("Aux_ibfs"+Arrays.toString(aux_ibfs));
		System.out.println("Aux_bIndices"+Arrays.toString(aux_bIndices));
		*/
		
		fracLP auxLP=new fracLP(aux_A,aux_b,aux_c,new_n+new_m,new_m);
		auxLP.cLPD=new fracLPDescriptor(aux_A,aux_b,aux_c,new_n+new_m,new_m);
		fracLPSolver auxSolver=new fracLPSolver(auxLP,aux_ibfs,aux_bIndices);
		auxSolver.solveLP(); 
		
		
		
		
		
		
		
		
		
		
		
		
		if (auxLP.optimal_cost.doubleValue()!=0)
			this.currentLP.optimal_cost=new BigFraction(Double.POSITIVE_INFINITY);
		
		//check Artificial Variables
		int[] aux_basis=auxLP.final_basis;
		System.out.println(Arrays.toString(A[0]));
		BigFraction[][] aux_B=this.SR.getSubMatrix(A_a,aux_basis);
		BigFraction[][] aux_B_inv=this.SR.invert(aux_B);
		BigFraction[][] B_invA=this.SR.mm_product(aux_B_inv,A_a);
		
		while(this.SR.checkAV(aux_basis,n)){
			
			aux_B_inv=this.SR.invert(aux_B);
			B_invA=this.SR.mm_product(aux_B_inv,A_a);
			
			//find artificial variable in basis
			int l=0;
			for(int i=0;i<aux_basis.length;i++){
				if (aux_basis[i]>=n){
					l=i;
					break;
				}
			} //AV in basis found
			
			//remove l from basis
			boolean check=true;
			for(int i=0;i<A_a[0].length;i++){
				if (A_a[l][i].doubleValue()!=0)
						check=false;
			}
			
			//case 1: check=true--> remove constraint l from A
			if(check){
				A_a=this.SR.RemoveRow(A_a, l);
				b_a=this.SR.RemoveComponent(b_a,l);
				break;
			}
			//case 2: check=false--> change basis vectors & update B
			int nbv=0;
			if(!check){
				for(int i=0;i<A_a[0].length;i++){
					
					if (A_a[l][i].doubleValue()!=0)
						nbv=i;
				}
				aux_basis[l]=nbv;
				aux_B=this.SR.getSubMatrix(A_a,aux_basis);
				
			}
		}
		
		aux_B_inv=this.SR.invert(aux_B);
		BigFraction[] x_b=this.SR.mv_product(aux_B_inv,b_a);
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		this.bIndices=auxLP.final_basis;
		this.ibfs=auxLP.x_values;
		this.solveLP();
		
		System.out.println("this"+this.currentLP.optimal_cost);
		
}
		
	
	
	
	
	
	
	
	
	
	
	
	
	
	/**
	Phase 2 of the revised simplex algorithm is started after an initial 
	basic feasible solution to the LP was found during Phase 1. In addition
	Phase 2 helps finding a solution to the auxiliary problem created
	during Phase 1.
	*/
	public double phase2(){
		
		BigFraction[][] A_a= this.currentLP.cLPD.A;
		BigFraction[] c_std = this.currentLP.cLPD.c;
		BigFraction[] b_a=this.currentLP.cLPD.b;
		int[] A_sense=this.currentLP.cLPD.A_sense;
		int n=this.currentLP.cLPD.n;
		int m=this.currentLP.cLPD.m;
		
		
		/*
		double[][] c_help= new double[1][n];
		c_help[0]=c_std;
		double[][] b_help= new double[1][n];
		Matrix A=new Matrix(A_a);
		Matrix c=new Matrix(c_help);
		Matrix b=new Matrix(b_help);
		
		System.out.println(Arrays.toString(this.range(0,currentLP.noOfConstraints-1)));
		double[][] dummy=new double[m][m];
		Matrix B=new Matrix(dummy);//=A.getMatrix(this.range(0,currentLP.noOfConstraints-1),this.bIndices);
		for(int l=0;l<this.bIndices.length;l++){
			int[] helpl={l};
			int[] helpl2={this.bIndices[l]};
			B.setMatrix(this.range(0,m-1),helpl, A.getMatrix(this.range(0,m-1),helpl2 ));	
		}
		*/
		BigFraction[][] B=SR.getSubMatrix(A_a,this.bIndices);
		BigFraction[][] B_inv=SR.invert(B);
		int[] helpi={0};
		
		BigFraction[] c_b=SR.getSubVector(c_std,this.bIndices);
		BigFraction[] p= SR.vm_product(c_b, B_inv);
		BigFraction[] ptA=SR.vm_product(p,A_a);
		BigFraction[] red_c = SR.vectorSubtract(c_std,ptA);
		
		BigFraction[] x=this.ibfs;
		boolean is_optimal=false;
		int u_min_index=10000;
		double u_pos_min=Double.POSITIVE_INFINITY;
		int[] helpj={0};
		//int[] helpr=this.range(0,m-1);
		BigFraction theta=new BigFraction(0);
		BigFraction[] u=new BigFraction[m];

		//as long as current BFS can be improved
		int counter=0;
		while (!is_optimal){
			counter++;
			//is set to false again as soon as something to be improved is found
			is_optimal=true;
			//check current nonbasic-variables for cost improvement
			int help=red_c.length;
			for(int j=0;j< red_c.length;j++){
				//if a cost-improving nonbasic variable j is found 
				if (red_c[j].doubleValue()<0){
					is_optimal=false;
					helpj[0]=j;
					//helpr=this.range(0,m-1);
					
					/**
					 Screen what basic variable determines how far we can go in jth basic direction without making any of them nonnegative.
					 This chosen basic variable i will be the one which is decreased to 0.
					 */
					BigFraction[] A_j=SR.getColoumn(A_a,j);
					u=SR.mv_product(B_inv, A_j);
					
					boolean u_all_neg=true;
					
					for (int i=0;i<this.bIndices.length;i++){
						
						if (u[i].doubleValue()>0){
							
							u_all_neg=false;
							BigFraction u_ratio=x[this.bIndices[i]].divide(u[i]);
				
							if (u_ratio.doubleValue()<u_pos_min){
								u_pos_min=u_ratio.doubleValue();
								u_min_index=i;
							}
						}
					}
					
					
					
					//here if existent the minimal positive u-ratio has been found 
					if (u_all_neg)
						return Double.NEGATIVE_INFINITY;
					//perform basis change
					theta=x[this.bIndices[u_min_index]].divide(u[u_min_index]);
					System.out.println(theta);
					x[j]=theta;
					for(int s=0;s<this.bIndices.length;s++){
						x[this.bIndices[s]]=x[this.bIndices[s]].subtract(theta.multiply(u[s]));
					}
					this.bIndices[u_min_index]=j;
					
					//update B
					B=SR.getSubMatrix(A_a,this.bIndices);
					B_inv=SR.invert(B);
					
					
					c_b=SR.getSubVector(c_std,this.bIndices);
					p= SR.vm_product(c_b, B_inv);
					ptA=SR.vm_product(p,A_a);
					red_c = SR.vectorSubtract(c_std,ptA);
					u_min_index=10000;
					u_pos_min=Double.POSITIVE_INFINITY;
					
					break;
					}
			}
			
			if (is_optimal)
				this.currentLP.optimal_cost=SR.vv_dot(c_std, x);     //to change!!!
				this.currentLP.x_values=x;
				this.currentLP.final_basis=this.bIndices;
				
			
		}
		//here loop has stopped and basis has not changed anymore
		double d=4;
		System.out.println("FINAL"+x.toString());
		return d;
		
		}			
					
					
					
					
					
					
				
	
		
		
			
		
			
		
	
	
	
	
	
	
	
	
	
	
	
	
		
}
