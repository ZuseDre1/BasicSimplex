import java.io.*;
import org.apache.commons.math3.fraction.BigFraction;
import java.util.*;
import java.text.ParseException;


public class fracLP {

	public FractionRoutines FR=new FractionRoutines();
	public LPReader LPR;
	public fracLPSolver LPS;
	public fracLPDescriptor oLPD;
	public fracLPDescriptor cLPD;
	public double[] lowerBoundVector;
	public double[] upperBoundVector;
	
	
	public BigFraction optimal_cost;
	public BigFraction[] x_values;
	public int[] final_basis;
	
	
	
	
	
	public fracLP(String filename){

		this.LPR= new LPReader(filename);
	
	

		try {
			LPR.readLP();
	    	} catch (ParseException p) {
		System.out.println(p);
	    	}
		catch (FileNotFoundException f) {
		System.out.println(f);
	    	}
		catch (IOException i) {
		System.out.println(i);
	    	}
		
		BigFraction[][] A=FR.getFracMatrix(this.LPR.constraintsMatrix());
		BigFraction[] b=FR.getFracVector(this.LPR.rhsVector());
		BigFraction[] c=FR.getFracVector(this.LPR.objectiveVector());
		
		
		this.oLPD=new fracLPDescriptor(A,this.LPR.senseVector(),b,c,this.LPR.objectiveSense(),this.LPR.noOfVariables(),this.LPR.noOfConstraints());
		
		this.lowerBoundVector = this.LPR.lowerBoundVector();
		this.upperBoundVector = this.LPR.upperBoundVector();
		
		
		

	
  	}
	
	
	
	
	
	
	
	
	/**
	 * Second Constructor to create auxiliary problem
	 */
	public fracLP(BigFraction[][] A,BigFraction[] b,BigFraction[] c,int n,int m){
		this.oLPD=new fracLPDescriptor(A,b,c,n,m);
		
	}
	
	
	
	
	
}
