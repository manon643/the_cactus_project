import java.util.LinkedList;

import Jcg.geometry.*;
import Jcg.polyhedron.Face;
import Jcg.polyhedron.Halfedge;
import Jcg.polyhedron.Polyhedron_3;
import Jcg.polyhedron.Vertex;

public class gui {
	
	public Polyhedron_3<Point_3> polyhedron3D;
	double pas = 0.1;
	
	public gui(Polyhedron_3<Point_3> polyhedron3D) {
		this.polyhedron3D=polyhedron3D;
	}
	
	void translate(int dir, LinkedList<Vertex<Point_3>> handle_vertices) {
		Vector_3 t = new Vector_3(); 
		if (Math.abs(dir)==1) {
			t = new Vector_3(dir*pas, 0, 0);
		}
		else if (Math.abs(dir)==2) {
			t = new Vector_3(0, dir/2*pas, 0);
		}
		else if (Math.abs(dir)==3) {
			t = new Vector_3(0, 0, dir/3*pas);
		}
		for (Vertex<Point_3> v : handle_vertices) {
			Point_3 point = v.getPoint().sum(t);
			v.setPoint(point);
		}
		return ;
	}
}
