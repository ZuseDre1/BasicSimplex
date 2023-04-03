import java.io.*;
import java.util.*;
import java.text.ParseException;
import Jama.Matrix;
import java.math.BigInteger;
import org.apache.commons.math3.fraction.BigFraction;



		
		
	
public class fracSimplexRoutines{
	
	
	
	
	//////////////////////////////Routines for Phase 1/////////////////////////////////////////////
	
	/**
	 * counts the number of slack variables necessary to have only equality constraints
	 * @param constraints
	 * @return
	 */
	public int countSlacks(fracLP thisLP)
	{	
		int count_extra=0;
		int[] const_type=thisLP.oLPD.A_sense;
		
		for(int i=0;i<const_type.length;i++){
			
			if(const_type[i]!=0)
				count_extra+=1;
		}
		System.out.println("count_extra"+count_extra);
		return count_extra;
	}
	
	
	
	
	
	
	public int countUnbounded(fracLP thisLP){
		
		int count_extra=0;
		double[] l_bounds=thisLP.lowerBoundVector;
		double[] u_bounds=thisLP.upperBoundVector;
		
			
			for(int i=0;i<l_bounds.length;i++){
				
				if(l_bounds[i]<0)
					count_extra+=1;
				/*
				if(l_bounds[i]==0 && u_bounds[i]==Double.POSITIVE_INFINITY)
					count_extra+=1;
				*/
			}
		System.out.println("count_extra"+count_extra);
		return count_extra;
	}
	
	
	
	
	
	
	public boolean checkAV(int[] basis,int thresh){
		
		boolean check=false;
		for(int i=0;i<basis.length;i++){
			if (basis[i]>=thresh)
				check=true;
		}
		return check;
	}
					
		
	
	
	
	
	
	
	
	
	public BigFraction[][] RemoveRow( BigFraction[][] A,int row){
		
		 BigFraction[][] A_new=new  BigFraction[A.length-1][A[0].length];
		for(int i=0;i<row;i++){
			for(int j=0;j<A[0].length;j++){
				A_new[i][j]=A[i][j];
			}
		}
		for(int i=row+1;i<A.length;i++){
			for(int j=0;j<A[0].length;j++){
				A_new[i-1][j]=A[i][j];
			}
		}

		return A_new ;
	}

	
	
	
	

	
	
	
	
	
	
	
public BigFraction[] RemoveComponent(BigFraction[]c,int comp){
		
		BigFraction[] c_new =new BigFraction[c.length-1]; 
		for(int i=0;i<comp;i++){
			c_new[i]=c[i];
		}
		for(int i=comp+1;i<c.length;i++){
			c_new[i-1]=c[i];
		}
		return c_new;
	}
	



	
	
	
	
	
	
///////////////////////////////////////Routines for Phase 2/////////////////////////////////////
	
	
	
public BigFraction vv_dot(BigFraction[] v1,BigFraction[] v2){
	
	BigFraction p=v1[0].multiply(v2[0]);
	for(int i=1;i<v1.length;i++){
		p=p.add(v1[i].multiply(v2[i]));
	}
	return p;
}









public BigFraction[] sv_product(BigFraction[] c,BigFraction s){
	
	for(int i=0;i<c.length;i++){
		c[i]=c[i].multiply(s);
	}
	return c;
}








public BigFraction[] vm_product(BigFraction[]c, BigFraction[][]B){
	
	BigFraction[] p=new BigFraction[B[0].length];
	for(int i=0;i<B[0].length;i++){
		for(int j=0;j<c.length;j++){
			p[i]=p[i].add(c[j].multiply(B[j][i]));
			             
		}
	}
	return p;
	
}



public BigFraction[] mv_product(BigFraction[][]B,BigFraction[]c){
	
	BigFraction[] p=new BigFraction[B.length];
	for(int i=0;i<B.length;i++){
		for(int j=0;j<c.length;j++){
			p[i]=p[i].add(B[i][j].multiply(c[j]));
		}
	}
	return p;
	
}






public BigFraction[][] mm_product(BigFraction[][]B,BigFraction[][]A){
	
	BigFraction[][] BA=new BigFraction[B.length][A[0].length];
	for(int i=0;i<B.length;i++){
		for(int j=0;j<A[0].length;j++){
			for(int k=0;k<A.length;k++){
				BA[i][j].add(B[i][k].multiply(A[k][j]));
			}
		}
	}
	return BA;
}























public BigFraction[][] getSubMatrix(BigFraction[][] A,int[] cols){
	
	BigFraction[][] B=new BigFraction[A.length][cols.length];
	for(int i=0;i<A.length;i++){
		for(int j=0;j<cols.length;j++){
			B[i][j]=A[i][cols[j]];
		}
	}
	return B;
}
	

public BigFraction[] getColoumn(BigFraction[][] A,int col){
	
	BigFraction[] B=new BigFraction[A.length];
	for(int i=0;i<A.length;i++){
		B[i]=A[i][col];
		}
	return B;
}
	


public BigFraction[] getSubVector(BigFraction[] c,int[] idx){
	
	BigFraction[] c_b=new BigFraction[idx.length];
	for(int i=0;i<idx.length;i++){
		c_b[i]=c[idx[i]];
	}
	return c_b;
}


public BigFraction[] vectorSubtract(BigFraction[] c1, BigFraction[] c2){
	
	BigFraction[] c=new BigFraction[c1.length];
	for(int i=0;i<c1.length;i++){
		c[i]=c1[i].subtract(c2[i]);
	}
	return c;
}


public BigFraction[][] invert(BigFraction[][] B){
	/*
	 * should be changed to calculation on fractions without conversion
	 */
	FractionRoutines FR=new FractionRoutines();
	double[][] B_double=FR.getDoubleMatrix(B);
	Matrix B_help=new Matrix(B_double);
	Matrix B_inv_help=B_help.inverse();
	double[][] B_inv_double=B_inv_help.getArray();
	BigFraction[][] B_inv=FR.getFracMatrix(B_inv_double);
	
	return B_inv;
	
}


	
	
	
	
	
}
