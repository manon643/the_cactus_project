
import java.util.LinkedList;

import Jcg.geometry.Point_3;
import Jcg.polyhedron.Polyhedron_3;
import Jcg.polyhedron.Vertex;



public abstract class SurfaceModeling {

	public  Polyhedron_3<Point_3> P_prime;
	
	public SurfaceModeling(Polyhedron_3<Point_3> polyhedron3D) {
		this.P_prime=polyhedron3D;
	}

	/**
	 * The main method performing the modeling process
	 * To be implemented
	 */
	public abstract void ComputeSorkineUntilThreshold(double epsilon);
	public abstract double ComputeSorkineIteration(int count);
	public abstract void startSorkine();
	public abstract void endSorkine(int count);

	public abstract void add_fixed(LinkedList<Vertex<Point_3>> fixed);

}
