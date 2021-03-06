import processing.core.*;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.LinkedList;

import Jcg.geometry.*;
import Jcg.polyhedron.*;

/**
 * A simple 3d viewer for visualizing surface meshes (based on Processing)
 * 
 * @author Luca Castelli Aleardi (INF555, 2012)
 *
 */
public class MeshViewer extends PApplet {

	SurfaceMesh mesh; // 3d surface mesh
	gui g;
	SurfaceModeling surface_modeling;
	
	int renderType=0; // choice of type of rendering
	int renderModes=3; // number of rendering modes
	int selectionMode = 0; //is selection mode activated
	int selectionModes = 3; //is selection mode activated
	boolean liveMode = false;
	double epsilon = 1e-3;
	int count = 0;
	LinkedList<Vertex<Point_3>> fixed_selection = new LinkedList<Vertex<Point_3>>(); //index of selected fixed vertex
	LinkedList<Vertex<Point_3>> handle_selection = new LinkedList<Vertex<Point_3>>(); //index of selected handle vertex
	
	String filename="OFF/high_genus.off";
	//String filename="OFF/sphere.off";
	//String filename="OFF/cube.off";
	//String filename="OFF/torus_33.off";
	//String filename="OFF/tore.off";
	//String filename="OFF/tri_hedra.off";
	//String filename="OFF/letter_a.off";
	//String filename="OFF/star.off";
	//String filename="OFF/tri_triceratops.off";
	//String filename="OFF/Meshes/bar1.off";//Vertices are too small 
	//String filename="OFF/Meshes/bar2.off";//Vertices are too small 
	//String filename="OFF/Meshes/bar3.off";
	//String filename="OFF/Meshes/cylinder_small.off";
	//String filename="OFF/Meshes/square_21.off";
	//String filename="OFF/Meshes/square_21_spikes.off"; //pb scale
	//String filename="OFF/Meshes/square_21_spikes_new.off"; //fixed pb scale
	//String filename="OFF/Meshes/dino.off"; //too large
	//String filename="OFF/Meshes/armadillo_1k.off";
	//String filename="OFF/Meshes/cactus_small.off";
	
	private PrintWriter writer;
	
	public void setup() {
		  size(800,600,P3D);
		  ArcBall arcball = new ArcBall(this);
		  
		  this.mesh=new SurfaceMesh(this, filename);
		  this.g = new gui(this.mesh.polyhedron3D);
		  LinkedList<Vertex<Point_3>> fixed = new LinkedList<Vertex<Point_3>>() ;
		  this.surface_modeling = new Sorkine(this.mesh.polyhedron3D);
		  
		  try {
			this.writer = new PrintWriter(new FileWriter("./running_time.txt", true));
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		  writer.println("NEW: "+filename);
		  writer.println(this.mesh.polyhedron3D.sizeOfVertices()+";"+this.mesh.polyhedron3D.sizeOfHalfedges()+";"+this.mesh.polyhedron3D.sizeOfFacets());
		  

	}
		 
		public void draw() {
		  background(255);
		  //this.lights();
		  directionalLight(101, 204, 255, -1, 0, 0);
		  directionalLight(51, 102, 126, 0, -1, 0);
		  directionalLight(51, 102, 126, 0, 0, -1);
		  directionalLight(102, 50, 126, 1, 0, 0);
		  directionalLight(51, 50, 102, 0, 1, 0);
		  directionalLight(51, 50, 102, 0, 0, 1);
		 
		  translate(width/2.f,height/2.f,-1*height/2.f);
		  this.strokeWeight(1);
		  stroke(150,150,150);
		  this.drawOptions();
		  
		  this.mesh.draw(renderType, fixed_selection, handle_selection);
		}
		public void drawOptions() {
			int hf = 13;
			fill(0);
			this.textMode(this.SCREEN);
			this.text("'r' for rendering modes", 10, hf);
			this.text("'c' for selection modes: fixed and handle", 10, 2*hf);
			
			this.text("'i' to run one iteration", 10, 4*hf);
			this.text("'s' for Surface Modeling", 10, 5*hf);
			this.text("'m'/'p' for (de)zooming", 10, 6*hf);
			this.text("arrows and 'f'/'b' to move the handle", 10, 7*hf);
			
			fill(255, 0, 0);
			if (count>0) this.text("Iteration "+count, 10, 9*hf);
			if (selectionMode==1) this.text("Selecting fixed points", 10, 10*hf);
			else if (selectionMode==2) this.text("Selecting handle points", 10, 10*hf);
			
		}
		
		public void keyPressed(){
	    		LinkedList<Vertex<Point_3>> fixed = new LinkedList<Vertex<Point_3>>() ;
	    		fixed.addAll(fixed_selection);
	    		fixed.addAll(handle_selection);
			int dir = 0; 
			  switch(key) {
			    case('r'): this.renderType=(this.renderType+1)%this.renderModes; break;
			    case('c'): case('C'): this.selectionMode=(this.selectionMode+1)%this.selectionModes; System.out.println("Selection Mode is "+this.selectionMode); break;
			    
			    //Zooming
			    case('m'): case('M'): this.mesh.scaleFactor *= 0.5; break;
			    case('p'): case('P'): this.mesh.scaleFactor *= 2; break;
			    
			    //Compute sorkine
			    case('s'):case('S'): 
			    		this.surface_modeling.add_fixed(fixed);
			    		this.surface_modeling.ComputeSorkineUntilThreshold(epsilon);
			    		this.surface_modeling = new Sorkine(this.mesh.polyhedron3D);
			    		break;
			    		
			    	//Iteration by iteration
			    case('i'):case('I'): 
			    		if (count==0){
			    			this.surface_modeling.add_fixed(fixed);
			    			this.surface_modeling.startSorkine();
			    		}
			    		long t = System.currentTimeMillis();
		    			this.surface_modeling.ComputeSorkineIteration(count);
		    			long d = (System.currentTimeMillis()-t);
		    			System.out.println("Iteration "+count+" took "+d+" ms.");
		    			writer.println(count+";"+d);
		    			count++;
		    		break;
		    		
			    case('e'):case('E'): 
		    			this.surface_modeling.endSorkine(count);
			    		this.surface_modeling = new Sorkine(this.mesh.polyhedron3D);
			    		writer.close();
		    		break;
			    		
			    	//Resets selection
			    case('x'): case('X'): this.fixed_selection = new LinkedList<Vertex<Point_3>>();this.handle_selection = new LinkedList<Vertex<Point_3>>();break; //Resets selections
			    
			    //Moves handle
			    case(CODED):		
			    		switch(keyCode) {
					    case(UP):  dir = 2; break;
					    case(DOWN): dir = -2; break;
					    case(LEFT):  dir = -1;	break;
					    case(RIGHT):  dir = 1; break;
			    		}
			    break;
			    case('f'): case('F'): dir = 3; break; //front
			    	case('b'): case('B'): dir = -3; break; //back
			  }
			  if (dir!=0) {
				  this.g.translate(dir, handle_selection);
				  //this.mesh.updateScaleFactor();
			  }
			  

		}
		
		public int select_closest_points(float mX, float mY) {
			float minDist2 = 1000000;
			int i_min = -1;
			//System.out.println(mouseX+","+ mouseY);
			float s = (float) this.mesh.scaleFactor;
			for (Vertex<Point_3> v :this.mesh.polyhedron3D.vertices) {
				Point_3 p = v.getPoint();
				float x = screenX(p.x.floatValue()*s, p.y.floatValue()*s, p.z.floatValue()*s);
				float y = screenY(p.x.floatValue()*s, p.y.floatValue()*s, p.z.floatValue()*s);
				float d = pow(x-mouseX,2) + pow(y-mouseY,2);
				if (d<minDist2){
					i_min = v.index;
					minDist2 = d;
				}
				//System.out.println(x+","+ y);
			}
			return i_min;
		}
		
		public void mouseClicked(){
			if (selectionMode>0) {
				int i_selection = this.select_closest_points(mouseX, mouseY);
				if (i_selection!=-1) {
					Vertex<Point_3> v = this.mesh.polyhedron3D.vertices.get(i_selection);
					if (fixed_selection.contains(v)) {
						fixed_selection.removeFirstOccurrence(v);
					}	
					else if (handle_selection.contains(v)) {
						handle_selection.removeFirstOccurrence(v);
					}
					else if (selectionMode == 1) {//selecting the fixed points
						fixed_selection.add(v);
					}
					else if (selectionMode == 2) {//selecting the handle points
						handle_selection.add(v);
					}
				}
			}
		}
		
		/**
		 * For running the PApplet as Java application
		 */
		public static void main(String args[]) {
			//PApplet pa=new MeshViewer();
			//pa.setSize(400, 400);
			PApplet.main(new String[] { "MeshViewer" });
		}

}
