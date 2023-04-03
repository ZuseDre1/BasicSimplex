import java.io.*;
import java.util.*;
import java.text.ParseException;
import Jama.Matrix;
import java.math.BigInteger;
import org.apache.commons.math3.fraction.BigFraction;




public class FractionRoutines {
	
	
	public BigFraction getFracValue(double s){
		
		BigFraction val=new BigFraction(s);
		return val;
	}
	
	
	
	
	
	
	
	
	public BigFraction[] getFracVector(double [] v){
		
		BigFraction[] v_frac=new BigFraction[v.length];
		
		for(int i=0;i<v.length;i++){
			v_frac[i]=new BigFraction(v[i]);
		}
		return v_frac;
	}
	
	
	
	
	
	
	
	
	public BigFraction[][] getFracMatrix(double [][] A){
			
			BigFraction[][] A_frac=new BigFraction[A.length][A[0].length];
			
			for(int i=0;i<A.length;i++){
				for(int j=0;j<A[0].length;j++){
					A_frac[i][j]=new BigFraction(A[i][j]);
				}
			}
			return A_frac;
	}
	
	
	
	
	
	public double[][] getDoubleMatrix(BigFraction[][] M){
		
		double[][] A=new double[M.length][M[0].length];
		for(int i=0;i<M.length;i++){
			for(int j=0;j<M[0].length;j++){
				A[i][j]=M[i][j].doubleValue();
			}
		}
		return A;
	}
	
	
}
