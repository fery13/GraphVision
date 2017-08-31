import static org.lwjgl.opengl.GL11.glBlendFunc;
import static org.lwjgl.opengl.GL11.glDisable;
import static org.lwjgl.opengl.GL11.glEnable;

import java.awt.Font;
import java.awt.FontFormatException;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;

import org.lwjgl.opengl.GL11;
//import org.lwjgl.util.glu.Sphere;
import org.newdawn.slick.TrueTypeFont;
import org.newdawn.slick.util.ResourceLoader;

import GraphVision.Centrality;
import GraphVision.edgeCrossing;
import drawingAlgorithms.SyncBurst;
import graphGenerator.RandomTree;



//import com.sun.prism.paint.Color;



public class CFF {
	
	


	
public static double [][] FFF;
public static int chh =0; 
public static int vv =0; 
public static double scc =0; 
public static int [] loc;
public static double [][] sincos;
public static int [] check ;
public static double [] bound;
	
	public static double [][][] Cir_Force_Free(int adj [][], int v, boolean animation, double links, double time, boolean t3d) throws FileNotFoundException, InterruptedException, FontFormatException 
	{
		long end= System.nanoTime();
		
		//********************************
	RandomTree t = new RandomTree(v);
		
		
	int [][] edge = t.getEdgeList();
	//int [][] edge = adj.clone();
	
	
	SyncBurst graph = new SyncBurst(edge, v, false, 4, false);
	
	double [][] pos1 = graph.Cir_Force_Free();
	edgeCrossing cr = new edgeCrossing(pos1, edge);
	
	ArrayList<int []> edgeCr = cr.getCrossedEdges();
	edgeCr.forEach(e -> {
		System.out.println(e[0]+" "+e[1]+" "+e[2]+" "+e[3]);
	});
	
	System.out.println("Crs :   "+cr.number_of_edge_crossing());
	
	double [][][] pos = new double [1][3][v];
	for(int i=0;i<v;i++) {
		pos[0][0][i] = pos1[0][i];
		pos[0][1][i] = pos1[1][i];
		
	}
		
		//**********************************
		
		
		pos=main_fitting_in_center(v, pos, (int) time, !animation,t3d );
		pos=limit_in_screen(v, pos, (int) time, !animation, t3d);
		pos=covert_to_glortho(v, pos, (int) time, !animation, t3d);
		if(t3d)
			check_circular(v,  pos, 0);
		return pos;
	}

	public static double [][][] Cir_Force_Free22(int adj [][], int v, boolean animation, double links, double time, boolean t3d) throws FileNotFoundException, InterruptedException, FontFormatException 
	{
		
		int ed [][];
		
		//GraphMethods.preparing_for_graphviz(adj, v, (int) links);
		
		long start = System.nanoTime();
		bound = new double [4];
		bound[0]=10000;bound[1]=-10000;bound[2]=10000;bound[3]=-10000;
		ArrayList<double []> temp = new ArrayList<double []>();
		
		double pos [][][];
		//GraphMethods.matrix_change_to_distance_cells(adj, v, (int) links);
		
		if(animation)
			pos = new double [(int) time][3][v];
		else
			pos = new double [1][3][v];
		
		double n=v;
		double m=links;
	
		
		
		
		double radi = 1;
		double con =0;
		for(int i=0;i<v;i++)
		{
			pos[0][0][i]= (double) Math.random() * (1 - 0.0);
			pos[0][1][i]= (double) Math.random() * (1 - 0.0);
			
				pos[0][2][i]= (double) Math.random() * (1 - 0.0);
		}
		
		
		
		double FA=1;
		double FR=1;
	
		System.out.println("node: "+v+" links: "+links+" spars: "+((2*links)/(n*(n-1))));
		
		graph.messages = new String[(int) time];
		graph.messages[0]=" ";
		double veloc = 15;
		if(t3d)
			veloc=Math.sqrt(v);
		double signe = 1;
		FFF = new double [2][(int) time];
		int sync = (int) (time/2);
		//nf=(int) (15);
		double b=1;
		double V = (int) (Math.sqrt(n)*1);
		V=5;
		double disp = 0.0001;
		boolean hypo = true;
		boolean degree1=true;
		double burst = time-sync;
		double en1=0;
		double en2=0;
		double stress =0;
		en1=0;
		double summi=0;
		double rate = 1/(double)v;
		double M=10;
		double max = 0;
		double lo = Math.log1p((double)n);
		double po = (double)n/lo;
		po = Math.sqrt(po);
		System.out.println("injaaaaaaaaaaa Starts!");
		for(int q=1;q<time;q++)  // timmeeeeeeeeeeemmmmmmmmmmmmeeeeeeeeeeeeeeeemmmmmmmmmmmme
		{	
			
			stress =0;
			summi=0;
			//System.out.println("   "+q);
			boolean doQ = false;
			FR = (1.0/(double)v)*q; // very small for small graphs. 0.00000001
			//FR = 0.00000001*q;
			//FR = 0.000000001;
			//FR = 1/(double)v;
			//FR= (double)(v);
			max=0;
			if(((q+1)%V==0 )|| q>=burst)
			{
				
				doQ=true;
				//FA = FR/(double) 100000.000001;
				if(q>=burst)
					FR = (1.0/(double)v)*(q-burst+1);
				else
					FR = (1.0/(double)v)*q;
				
				FA=FR/(double)v;
				System.out.println(" dd  "+FR);
				//System.out.println(q+" QT "+n+" time: "+time+"  burst:  "+burst);
				
			}
			
			int t=0;
			if(animation)
				t=q;
			else
				t=1;
			//*********
			//System.out.println(q+" "+pos[t-1][0][5]+" "+pos[t-1][1][5]);
			//*********
			
			double FF [][] = new double [6][v];
			double cos [] = new double [v];
			double sin [] = new double [v];
			double cosz [] = new double [v];
			degree1=false;
		
			if(( q%V==0  && q<burst && !t3d) || ( q>=burst && t3d )) 
			{
					sincos = new double [2][v];
					//V = V + 15; 
					int ms =(int )po;
					ArrayList <ArrayList <ArrayList <Integer>>> mat = new ArrayList <ArrayList <ArrayList <Integer>>> (); 
					loc = new int [v]; 
					//mat = QuadTreeFixed.creating_the_matrix(bound, v, mat, ms, pos, t-1);
					//System.out.println("Quad Start .... ");
					mat = QuadTreeFixed.creating_the_matrix(v, mat, ms, pos, t-1);
					//System.out.println("created mat");
					
					for(int i=0;i<v;i++)
					{
						//System.out.println(i);
						double [] node = new double [2];
						node[0] = pos[t-1][0][i];
						node[1] = pos[t-1][1][i];
						double [] sincos = new double [2];
						//sincos = QuadTreeFixed.sin_cos(bound, v,  mat,  ms, pos, t-1, i, node);
						sincos = QuadTreeFixed.sin_cos(v,  mat,  ms, pos, t-1, i,  node);
						
						cos[i] = -sincos[0];
						sin[i] = -sincos[1];
						
						
						
					}
			}
			else
			if( t3d || ( q>=burst && !t3d ))
			{
				
			
				degree1=true;
				for(int i=0;i<v;i++)
				{
					for(int j=i+1;j<v;j++)
					{
						if(true)
						{
							double diss = 0;
							if(hypo)
							{
								if(!t3d)
									diss = Math.sqrt((Math.pow(pos[t-1][0][i]-pos[t-1][0][j], 2)) + ( Math.pow(pos[t-1][1][i]-pos[t-1][1][j], 2)));
								else
									diss = Math.sqrt((Math.pow(pos[t-1][0][i]-pos[t-1][0][j], 2)) + ( Math.pow(pos[t-1][1][i]-pos[t-1][1][j], 2))+
											(Math.pow(pos[t-1][2][i]-pos[t-1][2][j], 2)));
							}
							else
							{
								if(!t3d)
									diss = Math.abs(pos[t-1][0][i]-pos[t-1][0][j]) + Math.abs(pos[t-1][1][i]-pos[t-1][1][j]);
								else
									diss = Math.abs(pos[t-1][0][i]-pos[t-1][0][j]) + Math.abs(pos[t-1][1][i]-pos[t-1][1][j])
									+ Math.abs(pos[t-1][2][i]-pos[t-1][2][j]);
							}
						
							if(diss<0.00000000001)
							{	
								diss=0.00000000001;
							}
							double tsin = ((pos[t-1][1][j]-pos[t-1][1][i])/diss);
							double tcos = ((pos[t-1][0][j]-pos[t-1][0][i])/diss);
							
							sin[i] += -tsin;
							cos[i] += -tcos;
								
							sin[j] += tsin;
							cos[j] += tcos;
								
							if(t3d)
							{
								double tcosz = ((pos[t-1][2][j]-pos[t-1][2][i])/diss);
								cosz[i] += -tcosz;
								cosz[j] += tcosz;
							}	
						}
					}
				}
				
			}
			
			
			
			for(int i=0; i<links;i++)
			{
				int r=adj[0][i];
				int s=adj[1][i];
				
				if(q>=burst)
				{
						double dis=0;
						if(hypo)
						{
							if(!t3d)
								dis= Math.sqrt((Math.pow(pos[t-1][0][r]-pos[t-1][0][s], 2)) + ( Math.pow(pos[t-1][1][r]-pos[t-1][1][s], 2)));
							else 
								dis= Math.sqrt((Math.pow(pos[t-1][0][r]-pos[t-1][0][s], 2)) + ( Math.pow(pos[t-1][1][r]-pos[t-1][1][s], 2))+
										 ( Math.pow(pos[t-1][2][r]-pos[t-1][2][s], 2)));
						}
						else
						{
							if(!t3d)
								dis= Math.abs(pos[t-1][0][r]-pos[t-1][0][s]) + Math.abs(pos[t-1][1][r]-pos[t-1][1][s]);
							else
								dis= Math.abs(pos[t-1][0][r]-pos[t-1][0][s]) + Math.abs(pos[t-1][1][r]-pos[t-1][1][s])+
										  Math.abs(pos[t-1][2][r]-pos[t-1][2][s]);
						}
						if(dis<0.00000000001)
						{
							dis=0.00000000001;
						}
						double coss=(pos[t-1][0][r]-pos[t-1][0][s])/dis;
						double sinn=(pos[t-1][1][r]-pos[t-1][1][s])/dis;
								
					
						FF[3][s] += -FA*coss + FR*coss;
						FF[4][s] += -FA*sinn + FR*sinn;
						
						FF[3][r] += FA*coss - FR*coss;
						FF[4][r] += FA*sinn - FR*sinn;
						
						if(t3d)
						{
							double cos_z = ((pos[t-1][2][r]-pos[t-1][2][s])/dis);
							FF[5][s] += -FA*cos_z +FR*cos_z;
							FF[5][r] += FA*cos_z -FR*cos_z;
						}
					}
					else
					{
						
						FF[3][s] += pos[t-1][0][r];
						FF[4][s] += pos[t-1][1][r];
						
						
						FF[3][r] += pos[t-1][0][s];
						FF[4][r] += pos[t-1][1][s];
						
						//summi += Math.abs(FF[3][s]-FF[3][r])+Math.abs(FF[4][s]-FF[4][r]);
						
						
						if(t3d)
						{
							FF[5][s] += pos[t-1][2][r];
							FF[5][r] += pos[t-1][2][s];
						}
						if(!degree1)
						{
							if(graphvision.node_degree[s]==1)
							{
								FF[3][s] += Math.cos(Math.random()*Math.PI*2)*disp;
								FF[4][s] += Math.sin(Math.random()*Math.PI*2)*disp;
								//System.out.println("single  "+s);
							}
							if(graphvision.node_degree[r]==1)
							{
								FF[3][r] += Math.cos(Math.random()*Math.PI*2)*disp;
								FF[4][r] += Math.sin(Math.random()*Math.PI*2)*disp;
								//System.out.println("single  "+r);
							}
						}
					}
			}
			
			temp = new ArrayList<double []>();
			temp.clear();
			double big_cir = FR;
			en1 = 0;
			
			for(int i=0;i<v;i++)
			{	
				if(animation)
				{
					if(q>=burst)
					{
						pos[t][0][i] =((FR*cos[i])+(FF[3][i]));
						pos[t][1][i] =((FR*sin[i])+(FF[4][i]));
					}
					else
					{
						pos[t][0][i] =((big_cir*cos[i])+(FF[3][i])/(graphvision.node_degree[i]));
						pos[t][1][i] =((big_cir*sin[i])+(FF[4][i])/(graphvision.node_degree[i]));
						
					}
				}
				else
				{
					if(q>=burst)
					{
						pos[t-1][0][i] =((FR*cos[i])+(FF[3][i]));
						pos[t-1][1][i] =((FR*sin[i])+(FF[4][i]));
						en1 += pos[t-1][0][i]+pos[t-1][1][i];
					}
					else
						{
							double xi = (FF[3][i])/(graphvision.node_degree[i]);
							double yi = (FF[4][i])/(graphvision.node_degree[i]);
							
							pos[t-1][0][i] =((big_cir*cos[i])+xi);
							pos[t-1][1][i] =((big_cir*sin[i])+yi);
							//System.out.println(xi+"  "+yi);
							en1 += pos[t-1][0][i]+pos[t-1][1][i];
						}
						/*if(i==155)
						System.out.println(pos[t-1][0][155]+"  "+pos[t-1][1][155]);*/
						
					}
					
					if(t3d)
					{
						if(animation)
						{
							if(q>=burst)
								pos[t][2][i] =((FR*cosz[i])+(FF[5][i]));
							else
								pos[t][2][i] =((big_cir*cosz[i])+(FF[5][i])/(graphvision.node_degree[i]));
						}
						else
						{
							if(q>=burst)
								pos[t-1][2][i] =((FR*cosz[i])+(FF[5][i]));
							else
								pos[t-1][2][i] =((big_cir*cosz[i])+(FF[5][i])/(graphvision.node_degree[i]));
						}
					}  
			}	
			
			
			
			//System.out.println("s: "+ (summi/links)+"  "+q+" "+v+" "+links);
			if(q>10 )
			{
				//break;
				//stress = Math.abs((en2 - en1 )/en1);
				//stress = (en2/en1);
				/*if( Math.abs(en1-(summi/links))<0.1)
					break;*/
				//System.out.println( stress+"  "+en1+"  "+en2);
			}
			en2 = en1;
			//en1 = (summi/links);
		}
		long end= System.nanoTime();
		
		long run_time = end - start;
		System.out.println("Running time: "+run_time);
		
		 //writeToFile("cords.text", v, pos, (int)links, adj);
		
		pos=main_fitting_in_center(v, pos, (int) time, !animation,t3d );
		pos=limit_in_screen(v, pos, (int) time, !animation, t3d);
		pos=covert_to_glortho(v, pos, (int) time, !animation, t3d);
		if(t3d)
			check_circular(v,  pos, 0);
		return pos;
	}
	
	
	public static double [][][] Cir_Force_Free6(int adj [][], int v, boolean animation, double links, double time, boolean t3d) throws FileNotFoundException, InterruptedException, FontFormatException 
	{
		
		long start = System.nanoTime();
		double pos [][][];
		
		if(animation)
			pos = new double [(int) time][3][v];
		else
			pos = new double [1][3][v];
		
		double [] bound = new double [4];
		bound[0]=10000;bound[1]=-10000;bound[2]=10000;bound[3]=-10000;
		ArrayList<double []> temp = new ArrayList<double []>();
		
		
		double root_attr=0.9;
		double cof_attr=links;
		
		double root_repuls=1.0;
		double cof_repul=1.0;
		
		double n=v;
		double m=links;
		
		for(int i=0;i<v;i++)
		{
			pos[0][0][i]= (double) Math.random() * (1 - 0.0);
			pos[0][1][i]= (double) Math.random() * (1 - 0.0);
			pos[0][2][i]= (double) Math.random() * (1 - 0.0);
			double [] cord = new double [3];
			cord[0]=pos[0][0][i];
			cord[1]=pos[0][1][i];
			cord[2]=i;
			temp.add(cord);
			
			 if(bound[0]>pos[0][0][i])
	    		 bound[0]=pos[0][0][i];
			 
	    	 if(bound[1]<pos[0][0][i])
	    		 bound[1]=pos[0][0][i];
	    	 
	    	 if(bound[2]>pos[0][1][i])
	    		 bound[2]=pos[0][1][i];
	    	 
	    	 if(bound[3]<pos[0][1][i])
	    		 bound[3]=pos[0][1][i];
		}
		
		
		double C= (2*m*m)/(n*n*(n-1)*10);
		double D=10;
		for(int q=1;q<time;q++)  // timmeeeeeeeeeeemmmmmmmmmmmmeeeeeeeeeeeeeeeemmmmmmmmmmmme
		{	
			int t=0;
			if(animation )
				t=q;
			else
				t=1;
			//*********
			System.out.println(q);
			double tt= q;
			double force=0; 
			force = Math.pow(tt*C*D,10.0);
			D = 10 + (30*q)/time;
			//*********
			//System.out.println(force+"  "+pos[t-1][0][0]+"  "+pos[t-1][1][0]);
			double x_forces[]=new double [v];
			double y_forces[]=new double [v];
			double z_forces[]=new double [v];
			
			ArrayList<quadtree> chs = new ArrayList<quadtree>();
		    quadtree main_root = new quadtree(bound,null,temp,false,3,0,chs);
		     
		    quadtree.creat_childern(bound,main_root ,temp, 0 );
			
			for(int i=0;i<v;i++)
			{
				
				quadtree.xx=0;
				quadtree.yy=0;
				
			    quadtree.traversing_tree(main_root, pos[t-1][0][i],pos[t-1][1][i],i);
			    
			    double sin=quadtree.yy;
				double cos=quadtree.xx;
				
				double FR = cof_repul*(Math.pow(force,root_repuls));
				
				x_forces[i] += -(FR)*cos;
				y_forces[i] += -(FR)*sin;
			 
			}
			
			for(int i=0; i<links;i++)
			{
				int r=adj[0][i];
				int s=adj[1][i];
				double dis= Math.sqrt((Math.pow(pos[t-1][0][r]-pos[t-1][0][s], 2)) + ( Math.pow(pos[t-1][1][r]-pos[t-1][1][s], 2)));
								
				double cos=(pos[t-1][0][r]-pos[t-1][0][s])/dis;
				double sin=(pos[t-1][1][r]-pos[t-1][1][s])/dis;
							
				double FA = cof_attr*(Math.pow(force,root_attr));
				
				x_forces[s] += FA*cos;
				y_forces[s] += FA*sin;
				
				x_forces[r] += -FA*cos;
				y_forces[r] += -FA*sin;
				
				
				if(t3d)
				{
					double cos_z = ((pos[t-1][2][r]-pos[t-1][2][s])/dis);
					z_forces[s] += FA*cos_z;
					z_forces[r] += -FA*cos_z;
				}
			}
			temp = new ArrayList<double []>();
			temp.clear();
			for(int i=0;i<v;i++)
			{
				if(animation)
				{
				
					pos[t][0][i]=(x_forces[i]);
					pos[t][1][i]=(y_forces[i]);
					 
					if(bound[0]>pos[t][0][i])
			    		 bound[0]=pos[t][0][i];
					 
			    	if(bound[1]<pos[t][0][i])
			    		 bound[1]=pos[t][0][i];
			    	 
			    	if(bound[2]>pos[t][1][i])
			    		 bound[2]=pos[t][1][i];
			    	 
			    	if(bound[3]<pos[t][1][i])
			    		 bound[3]=pos[t][1][i];
			    	double [] cord = new double [3];
			    	cord[0]=pos[t][0][i];
					cord[1]=pos[t][1][i];
					cord[2]=i;
					temp.add(cord);
					
				}
				else
				{
					pos[t-1][0][i]=(x_forces[i]);
					pos[t-1][1][i]=(y_forces[i]);
					
				
					 
					if(bound[0]>pos[t-1][0][i])
			    		 bound[0]=pos[t-1][0][i];
					 
			    	if(bound[1]<pos[t-1][0][i])
			    		 bound[1]=pos[t-1][0][i];
			    	 
			    	if(bound[2]>pos[t-1][1][i])
			    		 bound[2]=pos[t-1][1][i];
			    	 
			    	if(bound[3]<pos[t-1][1][i])
			    		 bound[3]=pos[t-1][1][i];
			    	double [] cord = new double [3];
			    	cord[0]=pos[t-1][0][i];
					cord[1]=pos[t-1][1][i];
					cord[2]=i;
					temp.add(cord);
					
				}
				if(t3d)
				{
					if(animation)
					{
						pos[t][2][i]=(z_forces[i]);
					}
					else
					{
						pos[t-1][2][i]=(z_forces[i]);
					}
				}  
			}	
		}
		long end= System.nanoTime();
		
		long run_time = end - start;
		System.out.println("Running time: "+run_time);
		pos=main_fitting_in_center(v, pos, (int) time, !animation,t3d );
		pos=limit_in_screen(v, pos, (int) time, !animation, t3d);
		pos=covert_to_glortho(v, pos, (int) time, !animation, t3d);
		
	
		return pos;
	}
	
	
	public static double [][][] main_fitting_in_center(int v, double [][][] pos, int time, boolean stat, boolean t3d)
	{
		double max_x,max_y,min_x,min_y, ave_x,ave_y,max, max_z, min_z, ave_z;
		
		if(!t3d)
		{
			if(stat)
			{
				
					max_x=pos[0][0][0];
					max_y=pos[0][1][0];
					min_x=pos[0][0][0];
					min_y=pos[0][1][0];
					max=0;
					for(int i=0;i<v;i++)
					{
						if(max_x<pos[0][0][i])
							max_x=pos[0][0][i];
						
						if(max_y<pos[0][1][i])
							max_y=pos[0][1][i];
						
						if(min_x>pos[0][0][i])
							min_x=pos[0][0][i];
						
						if(min_y>pos[0][1][i])
							min_y=pos[0][1][i];
					}
					
					
					ave_x=(max_x+min_x)/2.0;
					ave_y=(max_y+min_y)/2.0;
					
					
					for(int i=0;i<v;i++)
					{
						pos[0][0][i]=pos[0][0][i]-ave_x;
						pos[0][1][i]=pos[0][1][i]-ave_y;
						
						double temp=Math.sqrt(Math.pow(pos[0][0][i], 2)+Math.pow(pos[0][1][i], 2));
						
						if(temp>max)
							max=temp;
					}
					
					
					for(int i=0;i<v;i++)
					{
						pos[0][0][i]=pos[0][0][i]/max;
						pos[0][1][i]=pos[0][1][i]/max;	
					}
				
			}
			else
			{
				for(int t=0;t<time;t++)
				{
					max_x=pos[t][0][0];
					max_y=pos[t][1][0];
					min_x=pos[t][0][0];
					min_y=pos[t][1][0];
					max=0;
					for(int i=0;i<v;i++)
					{
						if(max_x<pos[t][0][i])
							max_x=pos[t][0][i];
						
						if(max_y<pos[t][1][i])
							max_y=pos[t][1][i];
						
						if(min_x>pos[t][0][i])
							min_x=pos[t][0][i];
						
						if(min_y>pos[t][1][i])
							min_y=pos[t][1][i];
					}
					
					
					ave_x=(max_x+min_x)/2.0;
					ave_y=(max_y+min_y)/2.0;
					
					
					for(int i=0;i<v;i++)
					{
						pos[t][0][i]=pos[t][0][i]-ave_x;
						pos[t][1][i]=pos[t][1][i]-ave_y;
						
						double temp=Math.sqrt(Math.pow(pos[t][0][i], 2)+Math.pow(pos[t][1][i], 2));
						
						if(temp>max)
							max=temp;
					}
					
					
					for(int i=0;i<v;i++)
					{
						pos[t][0][i]=pos[t][0][i]/max;
						pos[t][1][i]=pos[t][1][i]/max;	
					}
				}
			}
		}
		else
		{
			if(stat)
			{
				
					max_x=pos[0][0][0];
					max_y=pos[0][1][0];
					min_x=pos[0][0][0];
					min_y=pos[0][1][0];
					
					min_z=pos[0][2][0];
					max_z=pos[0][2][0];
					max=0;
					for(int i=0;i<v;i++)
					{
						if(max_x<pos[0][0][i])
							max_x=pos[0][0][i];
						
						if(max_y<pos[0][1][i])
							max_y=pos[0][1][i];
						
						if(min_x>pos[0][0][i])
							min_x=pos[0][0][i];
						
						if(min_y>pos[0][1][i])
							min_y=pos[0][1][i];
						
						if(max_z<pos[0][2][i])
							max_z=pos[0][2][i];
						
						if(min_z>pos[0][2][i])
							min_z=pos[0][2][i];
						
					}
					
					
					ave_x=(max_x+min_x)/2.0;
					ave_y=(max_y+min_y)/2.0;
					ave_z=(max_z+min_z)/2.0;
					
					
					for(int i=0;i<v;i++)
					{
						pos[0][0][i]=pos[0][0][i]-ave_x;
						pos[0][1][i]=pos[0][1][i]-ave_y;
						pos[0][2][i]=pos[0][2][i]-ave_z;
						
						double temp=Math.sqrt(Math.pow(pos[0][0][i], 2)+Math.pow(pos[0][1][i], 2)+Math.pow(pos[0][2][i], 2));
						
						if(temp>max)
							max=temp;
					}
					
					
					for(int i=0;i<v;i++)
					{
						pos[0][0][i]=pos[0][0][i]/max;
						pos[0][1][i]=pos[0][1][i]/max;	
						pos[0][2][i]=pos[0][2][i]/max;	
					}
				
			}
			else
			{
				for(int t=0;t<time;t++)
				{
					max_x=pos[t][0][0];
					max_y=pos[t][1][0];
					min_x=pos[t][0][0];
					min_y=pos[t][1][0];
					max_z=pos[t][2][0];
					min_z=pos[t][2][0];
					
					max=0;
					for(int i=0;i<v;i++)
					{
						if(max_x<pos[t][0][i])
							max_x=pos[t][0][i];
						
						if(max_y<pos[t][1][i])
							max_y=pos[t][1][i];
						
						if(min_x>pos[t][0][i])
							min_x=pos[t][0][i];
						
						if(min_y>pos[t][1][i])
							min_y=pos[t][1][i];
						
						if(max_z<pos[t][2][i])
							max_z=pos[t][2][i];
						
						if(min_z>pos[t][2][i])
							min_z=pos[t][2][i];
					}
					
					
					ave_x=(max_x+min_x)/2.0;
					ave_y=(max_y+min_y)/2.0;
					ave_z=(max_z+min_z)/2.0;
					
					for(int i=0;i<v;i++)
					{
						pos[t][0][i]=pos[t][0][i]-ave_x;
						pos[t][1][i]=pos[t][1][i]-ave_y;
						pos[t][2][i]=pos[t][2][i]-ave_z;
						
						double temp=Math.sqrt(Math.pow(pos[t][0][i], 2)+Math.pow(pos[t][1][i], 2)+Math.pow(pos[t][2][i], 2));
						
						if(temp>max)
							max=temp;
					}
					
					
					for(int i=0;i<v;i++)
					{
						pos[t][0][i]=pos[t][0][i]/max;
						pos[t][1][i]=pos[t][1][i]/max;	
						pos[t][2][i]=pos[t][2][i]/max;	
					}
				}
			
			}
		}
		
		
		return pos;
		
	}
	
	
	public static double [][][] limit_in_screen(int v, double [][][] pos, int time, boolean stat, boolean t3d)
	{	
		double lm=1.0;
		int m_time=time;
		
		if(!t3d)
		{
			if(stat)
			{
				for (int i = 0; i < v; i++) 
					{
						if(pos[0][0][i]>lm || pos[0][0][i]<-lm)
						{
							double m=Math.abs(lm/pos[0][0][i]);
							for (int i1 = 0; i1 < v; i1++) 
							{
								pos[0][0][i1]=pos[0][0][i1]*m;
								pos[0][1][i1]=pos[0][1][i1]*m;
							}
						}
						
						if(pos[0][1][i]>lm || pos[0][1][i]<-lm)
						{
							double m=Math.abs(lm/pos[0][1][i]);
							for (int i1 = 0; i1 < v; i1++) 
							{
								pos[0][0][i1]=pos[0][0][i1]*m;
								pos[0][1][i1]=pos[0][1][i1]*m;
							}
						}
					}
				
			}
			else
			{
				for (int t = 0; t < m_time; t++) 
				{
					for (int i = 0; i < v; i++) 
					{
						if(pos[t][0][i]>lm || pos[t][0][i]<-lm)
						{
							double m=Math.abs(lm/pos[t][0][i]);
							for (int i1 = 0; i1 < v; i1++) 
							{
								pos[t][0][i1]=pos[t][0][i1]*m;
								pos[t][1][i1]=pos[t][1][i1]*m;
								pos[t][2][i1]=pos[t][2][i1]*m;								
							}
						}
						
						if(pos[t][1][i]>lm || pos[t][1][i]<-lm)
						{
							double m=Math.abs(lm/pos[t][1][i]);
							for (int i1 = 0; i1 < v; i1++) 
							{
								pos[t][0][i1]=pos[t][0][i1]*m;
								pos[t][1][i1]=pos[t][1][i1]*m;
								pos[t][2][i1]=pos[t][2][i1]*m;								
							}
						}
						if(pos[t][2][i]>lm || pos[t][2][i]<-lm)
						{
							double m=Math.abs(lm/pos[t][2][i]);
							for (int i1 = 0; i1 < v; i1++) 
							{
								pos[t][0][i1]=pos[t][0][i1]*m;
								pos[t][1][i1]=pos[t][1][i1]*m;
								pos[t][2][i1]=pos[t][2][i1]*m;								
							}
						}
					}
				}
			}
		}
		else
		{
			if(stat)
			{
				for (int i = 0; i < v; i++) 
					{
						if(pos[0][0][i]>lm || pos[0][0][i]<-lm)
						{
							double m=Math.abs(lm/pos[0][0][i]);
							for (int i1 = 0; i1 < v; i1++) 
							{
								pos[0][0][i1]=pos[0][0][i1]*m;
								pos[0][1][i1]=pos[0][1][i1]*m;
								pos[0][2][i1]=pos[0][2][i1]*m;
							}
						}
						
						if(pos[0][1][i]>lm || pos[0][1][i]<-lm)
						{
							double m=Math.abs(lm/pos[0][1][i]);
							for (int i1 = 0; i1 < v; i1++) 
							{
								pos[0][0][i1]=pos[0][0][i1]*m;
								pos[0][1][i1]=pos[0][1][i1]*m;
								pos[0][2][i1]=pos[0][2][i1]*m;
							}
						}
						if(pos[0][2][i]>lm || pos[0][2][i]<-lm)
						{
							double m=Math.abs(lm/pos[0][2][i]);
							for (int i1 = 0; i1 < v; i1++) 
							{
								pos[0][0][i1]=pos[0][0][i1]*m;
								pos[0][1][i1]=pos[0][1][i1]*m;
								pos[0][2][i1]=pos[0][2][i1]*m;
							}
						}
					}
				
			}
			else
			{
				for (int t = 0; t < m_time; t++) 
				{
					for (int i = 0; i < v; i++) 
					{
						if(pos[t][0][i]>lm || pos[t][0][i]<-lm)
						{
							double m=Math.abs(lm/pos[t][0][i]);
							for (int i1 = 0; i1 < v; i1++) 
							{
								pos[t][0][i1]=pos[t][0][i1]*m;
								pos[t][1][i1]=pos[t][1][i1]*m;
								pos[t][2][i1]=pos[t][2][i1]*m;
							}
						}
						
						if(pos[t][1][i]>lm || pos[t][1][i]<-lm)
						{
							double m=Math.abs(lm/pos[t][1][i]);
							for (int i1 = 0; i1 < v; i1++) 
							{
								pos[t][0][i1]=pos[t][0][i1]*m;
								pos[t][1][i1]=pos[t][1][i1]*m;
								pos[t][2][i1]=pos[t][2][i1]*m;
							}
						}
						
						if(pos[t][2][i]>lm || pos[t][2][i]<-lm)
						{
							double m=Math.abs(lm/pos[t][2][i]);
							for (int i1 = 0; i1 < v; i1++) 
							{
								pos[t][0][i1]=pos[t][0][i1]*m;
								pos[t][1][i1]=pos[t][1][i1]*m;
								pos[t][2][i1]=pos[t][2][i1]*m;
							}
						}
					}
				}
			}
		}
		
		double max_x=pos[0][0][0],max_y=pos[0][1][0], min_x=pos[0][0][0],min_y=pos[0][1][0];
		if(stat)
		{
			for(int i=0;i<v;i++)
			{
				if(max_x<pos[0][0][i])
					max_x=pos[0][0][i];
				if(max_y<pos[0][1][i])
					max_y=pos[0][1][i];
				
				if(min_x>pos[0][0][i])
					min_x=pos[0][0][i];
				if(min_y>pos[0][1][i])
					min_y=pos[0][1][i];
				
			}
		}
		scc  = Math.max( (Math.max(Math.abs(min_y), Math.abs(max_y))), (Math.max(Math.abs(min_x), Math.abs(max_x)))  );
		
		
		return pos;
	}
	

	public static double [][][] covert_to_glortho(int v, double [][][] pos, int time, boolean stat, boolean t3d)
	{
		double scale=2;
	
		int m_time=time;
		if(!t3d)
		{
			if(stat)
			{
				for(int j=0;j<v;j++)
				{
					pos[0][0][j]=((pos[0][0][j]*scale+1)*(graph.weidth/2));
					pos[0][1][j]=((pos[0][1][j]*scale+1)*(graph.hight/2));
				}
			}
			else
			{
				for(int i=0;i<m_time;i++)
				{
					for(int j=0;j<v;j++)
					{
						pos[i][0][j]=((pos[i][0][j]*scale+1)*(graph.weidth/2));
						pos[i][1][j]=((pos[i][1][j]*scale+1)*(graph.hight/2));
					}
					
				}
			}
		}
		else
		{
			if(stat)
			{
				for(int j=0;j<v;j++)
				{
					pos[0][0][j]=((pos[0][0][j]*scale+1)*(graph.weidth/2));
					pos[0][1][j]=((pos[0][1][j]*scale+1)*(graph.hight/2));
					pos[0][2][j]=((pos[0][2][j]*scale+1)*(graph.depth/2));
				}
			}
			else
			{
				for(int i=0;i<m_time;i++)
				{
					for(int j=0;j<v;j++)
					{
						pos[i][0][j]=((pos[i][0][j]*scale+1)*(graph.weidth/2));
						pos[i][1][j]=((pos[i][1][j]*scale+1)*(graph.hight/2));
						pos[i][2][j]=((pos[i][2][j]*scale+1)*(graph.depth/2));
					}
				}
			}
		}
		
		
		return pos;
	}


	public static void check_circular(int v, double [][][] pos, int t)
	{
		double maxes [] = new double [6];
		maxes[0]=pos[t][0][0];maxes[1]=pos[t][0][0];
		maxes[2]=pos[t][1][0];maxes[3]=pos[t][1][0];
		maxes[4]=pos[t][2][0];maxes[5]=pos[t][2][0];
		
		for(int i=0;i<v;i++)
		{
			if(maxes[0]>pos[t][0][i])
				maxes[0]=pos[t][0][i];
			if(maxes[1]<pos[t][0][i])
				maxes[1]=pos[t][0][i];
			
			if(maxes[2]>pos[t][1][i])
				maxes[2]=pos[t][1][i];
			if(maxes[3]<pos[t][1][i])
				maxes[3]=pos[t][1][i];
			
			if(maxes[4]>pos[t][2][i])
				maxes[4]=pos[t][2][i];
			if(maxes[5]<pos[t][2][i])
				maxes[5]=pos[t][2][i];
		}
		double cen [] = new double [3];
		cen[0] = (maxes[1]-maxes[0])/2+maxes[0];
		cen[1] = (maxes[3]-maxes[2])/2+maxes[2];
		cen[2] = (maxes[5]-maxes[4])/2+maxes[4];
		graphvision.cent = cen.clone();
		double dis [] = new double [v];
		double ave =0;
		double stnd =0;
		graphvision.radi =Math.sqrt(Math.pow(pos[t][0][0]-cen[0], 2)+Math.pow(pos[t][1][0]-cen[1], 2)+Math.pow(pos[t][2][0]-cen[2], 2));
		for(int i=0;i<v;i++)
		{
			dis[i] = Math.sqrt(Math.pow(pos[t][0][i]-cen[0], 2)+Math.pow(pos[t][1][i]-cen[1], 2)+Math.pow(pos[t][2][i]-cen[2], 2));
			if(graphvision.radi>dis[i])
				graphvision.radi=dis[i];
			//System.out.println("d: "+dis[i]);
			ave +=dis[i];
		}
		ave = ave/v;
		for(int i=0;i<v;i++)
		{
			stnd += Math.pow(dis[i]-ave,2);
		}
		stnd = stnd/v;
		stnd = Math.sqrt(stnd);
		System.out.println("std: "+ stnd);
		
		
	}
		
	public static void writeToFile(String addr, int v, double [][][] pos, int links, int adj[][]) throws FileNotFoundException {
		File simi = new File(addr);
		PrintWriter cords = new PrintWriter(simi);
		for (int i = 0; i < v; i++) {
			cords.println(pos[0][0][i]+"  "+pos[0][1][i]);
		}
		cords.println();cords.println();cords.println();cords.println("Edge List:");
		for (int i = 0; i < links; i++) {
			cords.println(adj[0][i]+" "+adj[1][i]);
		}
		cords.close();
		System.out.println("File has been written!");
		
	}
	
	
}
