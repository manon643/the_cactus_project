import java.io.*;
import Jcg.geometry.*;
import Jcg.polyhedron.*;
import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import Jama.SingularValueDecomposition;

import java.util.*;

public class Sorkine extends SurfaceModeling {
	
	public HashMap<Vertex<Point_3>, Point_3> originalPoints;
	public HashMap<Halfedge<Point_3>,Double> weights;
	public HashMap<Halfedge<Point_3>,Double> originalLengths;
	public HashMap<Vertex<Point_3>, Rotation_3> rotations;
	public LinkedList<Vertex<Point_3>> fixed;
	public Matrix invL;
	
	//To compute iteration by iteration
	private double prevEnergy;
	private PrintWriter writer;
	
	
	public Sorkine(Polyhedron_3<Point_3> polyhedron){
		super(polyhedron);
		
		this.originalPoints = new HashMap<Vertex<Point_3>, Point_3>();
		for (Vertex<Point_3> v: polyhedron.vertices){
			Point_3 p = new Point_3(v.getPoint());
			this.originalPoints.put(v, p);
		}
		
		weights = new HashMap<Halfedge<Point_3>,Double>();
		//Using wij=wji
		for (Halfedge<Point_3> h : polyhedron.halfedges){
			if (!weights.containsKey(h)) {
				//double w2 = ComputeWeight2(h);
				double w = ComputeWeight(h);
				//if (Math.abs(w-w2)>0.01)
				//	System.out.println(w+" != "+w2);
				this.weights.put(h, w);
				this.weights.put(h.opposite, w);
			}	
		}
		this.computeOriginalLengths();
		this.rotations = new HashMap<Vertex<Point_3>, Rotation_3>();	
	}

	public void add_fixed(LinkedList<Vertex<Point_3>> f) {
		this.fixed = f;
		computeInvL();
		return;
	}

	private void computeOriginalLengths() {
		this.originalLengths = new HashMap<Halfedge<Point_3>, Double>();
		for (Halfedge<Point_3> h : this.P_prime.halfedges) {
			Point_3 A = h.getVertex().getPoint();
			Point_3 B = h.getOpposite().getVertex().getPoint();
			originalLengths.put(h, A.distanceFrom(B).doubleValue());
		}
	}
	
	private double compareLenghts() {
		double r = 0;
		for (Halfedge<Point_3> h : this.P_prime.halfedges) {
			Point_3 A = h.getVertex().getPoint();
			Point_3 B = h.getOpposite().getVertex().getPoint();
			double l = A.distanceFrom(B).doubleValue();
			r+= Math.pow((l-originalLengths.get(h))/originalLengths.get(h), 2);
		}
		return Math.sqrt(r/this.P_prime.halfedges.size());
	}
	public void startSorkine() {
		this.prevEnergy = 0;
		this.writer = null;
		try {
			this.writer = new PrintWriter(new FileWriter("./Energy_graph.txt", true));
		} catch (IOException e) {
			e.printStackTrace();
		}
		writer.println("new");
	}
 	public double ComputeSorkineIteration(int count){
		//Step 1 : compute the rotations
		this.rotations = new HashMap<Vertex<Point_3>,Rotation_3>();
		for (Vertex<Point_3> v: this.P_prime.vertices){
			this.rotations.put(v, this.ComputeRotation(v));
		}

		
		//Step 2 : compute the points
		Point_3[] pts = this.ComputePoints(fixed);
		
		// Step 3 : check energy and difference
		double temp = this.prevEnergy;
		this.prevEnergy = this.LocalRigidityEnergy(pts);
		double diff;
		if (count>0)
			diff = temp-this.prevEnergy;
		else 
			diff = this.prevEnergy;
		
		writer.println(count +";"+ this.prevEnergy+";"+ diff);

		//Step 4 : move points
		this.movePts(pts);
		
		return Math.abs(diff);
	}
	
 	public void endSorkine(int count) {
		System.out.println("RMS is "+this.compareLenghts());
		writer.close() ;
		System.out.println("END : "+(count)+" iterations, energy is "+ this.prevEnergy);
 	}
 	
 	
	public void ComputeSorkineUntilThreshold(double epsilon){
		startSorkine();
		int count = 0;
		double abs_diff = 100;
		while (abs_diff<epsilon || count<=40){
			abs_diff = ComputeSorkineIteration(count);
			count++;
		}
	}
	
	
	private void movePts(Point_3[] pts) {
		for (Vertex<Point_3> v: this.P_prime.vertices) {
			v.setPoint(pts[v.index]);
		}
	}

	public double LocalRigidityEnergy(Point_3[] pts){
		double res=0;
		
		for (Vertex<Point_3> v: this.P_prime.vertices){
			int i = v.index;
			Point_3 pi = v.getPoint();
			Halfedge<Point_3> h= v.getHalfedge();
			Halfedge<Point_3> g = h.getOpposite();
			
			int count = 0;
			while (count == 0 || !g.opposite.equals(h)){
				count++;
				
				Vertex<Point_3> vj = g.getVertex();
				int j = vj.index;
				Point_3 pj = vj.getPoint();
				
				double[] m = new double[]{pi.x-pj.x,pi.y-pj.y,pi.z-pj.z};
				Matrix Eij = new Matrix(m,1);
				Eij = this.rotations.get(v).m.times(Eij.transpose());
				
				Point_3 pip = pts[i];
				Point_3 pjp = pts[j];
				
				double[] mp = new double[]{pip.x-pjp.x,pip.y-pjp.y,pip.z-pjp.z};
				Matrix Eijp = new Matrix(mp,1);
				
				Matrix temp = new Matrix(3,1);
				temp= (Eijp.transpose().minus(Eij)).times(this.weights.get(h));
				
				res+= Math.pow(temp.norm2(),2);
				
				g = g.opposite.next;
			}
		}
		return res;
	}
	
	private void computeInvL() {
		int n = this.P_prime.vertices.size();
		Matrix L = new Matrix(n, n);
		for (Vertex<Point_3> vi : this.P_prime.vertices){
			int i = vi.index;
			
			//if fixed we shouldn't compute the point
			if (fixed.contains(vi)) {
				L.set(i, i, 1);
			}
			else {
				Halfedge<Point_3> h = vi.getHalfedge();
				Halfedge<Point_3> g = h.getOpposite();
				Vertex<Point_3> vj;
				
				double wi = 0;
				int count = 0;
				while (count == 0 || !g.opposite.equals(h)){
					count++;
					double wij = this.weights.get(g);
					
					vj = g.getVertex();
					int j = vj.index;
					
					L.set(i, j, -wij);
					wi += wij;
					
					g = g.opposite.next;
				}
				L.set(i, i, wi);
			}
		}
		//System.out.println(Arrays.deepToString(L.getArray()));
		this.invL = L.inverse();
	}

	public Point_3[] ComputePoints(LinkedList<Vertex<Point_3>> fixed){
		int n = this.P_prime.vertices.size();
		Matrix b = new Matrix(n, 3);
		
		//construct matrix b
		
		for (Vertex<Point_3> vi : this.P_prime.vertices){
			int i = vi.index;
			Point_3 pip = vi.getPoint();
			Point_3 pi = this.originalPoints.get(vi);
			
			//if fixed we shouldn't compute the point
			if (fixed.contains(vi)) {
				b.set(i, 0, pip.x);
				b.set(i, 1, pip.y);
				b.set(i, 2, pip.z);
			}
			else {
				Halfedge<Point_3> h = vi.getHalfedge();
				Halfedge<Point_3> g = h.getOpposite();
				Vertex<Point_3> vj;
				
				int count = 0;
				Matrix row = new Matrix(1, 3);
				while (count == 0 || !g.opposite.equals(h)){
					count++;
					double wij = this.weights.get(g);
					
					vj = g.getVertex();
					
					Point_3 pj = this.originalPoints.get(vj);
					double[] m = new double[]{pi.x-pj.x,pi.y-pj.y,pi.z-pj.z};
					Matrix Eij = new Matrix(m,1);
					Matrix temp = new Matrix(3,1);
					temp = (this.rotations.get(vi).m.plus(this.rotations.get(vj).m)).times(Eij.transpose()).times(0.5*wij);
					row = row.plus(temp.transpose());

					g = g.opposite.next;
				}
				b.set(i, 0, row.get(0, 0));
				b.set(i, 1, row.get(0, 1));
				b.set(i, 2, row.get(0, 2));
			}
			
			
		}
		//System.out.println(Arrays.deepToString(L.getArrayCopy()));
		Matrix Pp = invL.times(b);
		
		//Transforming to array
		Point_3[] pts = new Point_3[n];
		for (Vertex<Point_3> vi : this.P_prime.vertices){
			int i = vi.index;
			Point_3 point = new Point_3(Pp.get(i, 0), Pp.get(i, 1), Pp.get(i, 2));
			pts[i] = point;
			
		}
		return pts;
		
	}
	
 	public Rotation_3 ComputeRotation(Vertex<Point_3> vi){
		Point_3 pip = vi.getPoint();
		Point_3 pi = this.originalPoints.get(vi);
		
		Halfedge<Point_3> h= vi.getHalfedge();
		Halfedge<Point_3> g = h.getOpposite();
		int count=0;
		
		Matrix S = new Matrix(3,3);
		
		while (count == 0 || !g.opposite.equals(h)){
			count++;
			
			Vertex<Point_3> vj = g.getVertex();
			Point_3 pjp = vj.getPoint();
			Point_3 pj = this.originalPoints.get(vj);
			double[] m = new double[]{pi.x-pj.x,pi.y-pj.y,pi.z-pj.z};
			Matrix Eij = new Matrix(m,1);
			
			double[] mp = new double[]{pip.x-pjp.x,pip.y-pjp.y,pip.z-pjp.z};
			Matrix Eijp = new Matrix(mp,1);
			
			S.plusEquals((Eij.transpose().times(Eijp)).times(this.weights.get(g)));
			
			g = g.opposite.next;
		}
		assert(count==this.P_prime.vertexDegree(vi));
		
		SingularValueDecomposition d = new SingularValueDecomposition(S);
		Matrix U = d.getU();
		Matrix V = d.getV();
		// to verify
		Matrix R = V.times(U.transpose());
		
		return (new Rotation_3(R));
	}

	
	public double Area(Face<Point_3> f){
		/*Computes the area of the triangle face f*/
		Halfedge<Point_3> h = f.getEdge();
		Point_3 p1 = h.getVertex().getPoint();
		Point_3 p2 = h.getNext().getVertex().getPoint();
		Point_3 p3 = h.getNext().getNext().getVertex().getPoint();
		Vector_3 a = (Vector_3)p2.minus(p1);
		Vector_3 b = (Vector_3)p3.minus(p1);
		Vector_ v = a.crossProduct(b);
		
		double res = 0.5* Math.sqrt(v.squaredLength().doubleValue());
		return res;
	}
	
	public double ComputeWeight2(Halfedge<Point_3> h){
		Point_3 A= h.getVertex().getPoint();
		Point_3 B= h.getOpposite().getVertex().getPoint();
		Point_3 C= h.getNext().getVertex().getPoint();
		Point_3 D= h.getOpposite().getNext().getVertex().getPoint();
		
		double d0c= A.squareDistance(B).doubleValue();
		double d1c= C.squareDistance(A).doubleValue();
		double d2c= C.squareDistance(B).doubleValue();
		double d3c= D.squareDistance(A).doubleValue();
		double d4c= D.squareDistance(B).doubleValue();
		
		double w;
		try {
			Face<Point_3> f1 = h.getFace(); Face<Point_3> f2 = h.opposite.face; //Checking if boundaries
			double cos1 = (d1c+d2c-d0c)/(2*Math.sqrt(d1c*d2c));
			double cot1 = cos1/Math.sqrt(1-cos1*cos1);
			double cos2 = (d3c+d4c-d0c)/(2*Math.sqrt(d3c*d4c));
			double cot2 = cos2/Math.sqrt(1-cos2*cos2);
			w = 1./2*(cot1+cot2);
		}
		catch(Exception ex) {
			try {
				Face<Point_3> f1 = h.getFace(); //Checking if boundaries
				double cos1 = (d1c+d2c-d0c)/(2*Math.sqrt(d1c*d2c));
				double cot1 = cos1/Math.sqrt(1-cos1*cos1);
				w = cot1;
				//System.out.println("boundaries here - 1");
			}
			catch(Exception ex2) {
				Face<Point_3> f2 = h.getOpposite().getFace(); //Checking if boundaries
				double cos2 = (d3c+d4c-d0c)/(2*Math.sqrt(d3c*d4c));
				double cot2 = cos2/Math.sqrt(1-cos2*cos2);
				w = cot2;
			}
			
		}
		
		//return (double) Math.round(w*100000d)/100000d;
		return w;
	}

	public double ComputeWeight(Halfedge<Point_3> h){
		Point_3 A= h.getVertex().getPoint();
		Point_3 B= h.getOpposite().getVertex().getPoint();
		Point_3 C= h.getNext().getVertex().getPoint();
		Point_3 D= h.getOpposite().getNext().getVertex().getPoint();
		
		double d0c= A.squareDistance(B).doubleValue();
		double d1c= C.squareDistance(A).doubleValue();
		double d2c= C.squareDistance(B).doubleValue();
		double d3c= D.squareDistance(A).doubleValue();
		double d4c= D.squareDistance(B).doubleValue();
		
		double A1, A2;
		double w;
		try {
			Face<Point_3> f1 = h.getFace(); Face<Point_3> f2 = h.getOpposite().getFace();
			A1 = Area(f1); 
			A2 = Area(f2);
			w = -1./8*((1./A1)*(d0c-d1c-d2c) + (1./A2)*(d0c-d3c-d4c));
		}
		catch(Exception ex) {
			try {
				Face<Point_3> f1 = h.getFace();
				A1 = Area(f1); w = -1./4*((1./A1)*(d0c-d1c-d2c));
				//System.out.println("boundaries here - 1");
			}
			catch(Exception ex2) {
				Face<Point_3> f2 = h.getOpposite().getFace();
				A2 = Area(f2); w = -1./4*((1./A2)*(d0c-d3c-d4c));
				//System.out.println("boundaries here - 2");
			}
			
		}
		
		//return (double) Math.round(w*100000d)/100000d;
		return w;
	}
		
	
}
