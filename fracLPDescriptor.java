import java.math.BigInteger;
import org.apache.commons.math3.fraction.BigFraction;

public class fracLPDescriptor {
/**
 * Class represents basic set of properties
 * of an LP (no bounds, bounds are property of LP directly)
 */
	
	public BigFraction[][] A;
	public int[] A_sense;
	public BigFraction[] b;
	public BigFraction[] c;
	public int c_sense;
	public int n;
	public int m;
	
	
	public fracLPDescriptor( BigFraction[][] A,int[] A_sense,BigFraction[] b,BigFraction[] c,int c_sense,int n, int m){
		this.A=A;
		this.A_sense=A_sense;
		this.b=b;
		this.c=c;
		this.c_sense=c_sense;
		this.n=n;
		this.m=m;
	}
	
	public fracLPDescriptor( BigFraction[][] A,BigFraction[] b,BigFraction[] c,int n, int m){
		this.A=A;
		this.b=b;
		this.c=c;
		this.n=n;
		this.m=m;
	}
	
		

}
	

