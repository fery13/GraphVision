import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import edu.uci.ics.jung.algorithms.shortestpath.DijkstraShortestPath;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;




public class FR {
	
	static double largest_btw_cnt=0;
	static double smallest_btw_cnt=0;
	public static double [][][] force_based_layout(int adj [][], int v, boolean animation, double links, double time, boolean t3d) throws FileNotFoundException, InterruptedException
	{		
		//GraphMethods.cactus_creating(v, "", links);
		
		//btw_cnt(, v, (int) links);
		
		long start_point=0,end_point=0;
		long GA_run_start=0,GA_run_end=0;
		start_point=System.nanoTime();
		
		boolean GA=true;
		
		double eng0,eng1 = 0;
		boolean ch=false;
		eng0=0;
		
		double a,b;
	//	a=1.0;
		a=1.0;
		b=1.0;
		double pos [][][];
		
		if(animation)
			pos = new double [(int) time][3][v];
		else
			pos = new double [1][3][v];
		
		
		
		double interval =Math.PI/v;
		double vv=v;
		double k=Math.sqrt(1.0/vv);
		
		
		if(GA)
		{
			int adj_temp [][] = new int [v][v];
			for(int i=0;i<links;i++)
			{
				int r=adj[0][i];
				int s=adj[1][i];
				adj_temp[r][s]=1;
				adj_temp[s][r]=1;
			}
			double [][] posGA = new double [2][v];
			posGA=Geodesic_GA.GA_2015(adj_temp, v, 50, 4,  v,60, 80, v);
			
			for(int i=0;i<v;i++)
			{
				pos[0][0][i]=posGA[0][i];
				pos[0][1][i]=posGA[1][i];
			}
		}
		else
		{
			for(int i=0;i<v;i++)
			{
				pos[0][0][i]= (double) Math.random() * (2 - 0.0)-1;
				pos[0][1][i]= (double) Math.random() * (2 - 0.0)-1;
				pos[0][2][i]= (double) Math.random() * (2 - 0.0)-1;
			}
		}
		double c4=0.1;
		double aneal_rate=c4/time;
		
		double x_force [];
		double y_force [];
		double z_force [] = null;
	
		
		//start_point=System.nanoTime();
		//End metric --------------
		int t=0;
		for(int q=1;q<time;q++)
		{
			
			x_force = new double [v];
			y_force = new double [v];
			if(t3d)
				z_force = new double [v];
			
			if(animation && v<=100)
				t=q;
			else
				t=1;
			
			double max_displacemnet[]=new double [3];
						
			for(int i=0;i<v;i++)
			{
				x_force[i]=0;
				y_force[i]=0;
				if(t3d)
					z_force[i]=0;
				
					
				for(int j=0;j<v;j++)
				{
					if(i != j)
					{
						double diss=0;
						if(t3d)
							diss = Math.sqrt((Math.pow(pos[t-1][0][i]-pos[t-1][0][j], 2)) + ( Math.pow(pos[t-1][1][i]-pos[t-1][1][j], 2))+
									(Math.pow(pos[t-1][2][i]-pos[t-1][2][j], 2)));
						else
							diss= Math.sqrt((Math.pow(pos[t-1][0][i]-pos[t-1][0][j], 2)) + ( Math.pow(pos[t-1][1][i]-pos[t-1][1][j], 2)));
						
						
						
						/*dis= (Math.abs(variables.x_pos[t-1][i]-variables.x_pos[t-1][j]) + (Math.abs(variables.y_pos[t-1][i]-variables.y_pos[t-1][j])));
						double diss=1.0*dis;*/
						
						double cos=((pos[t-1][0][j]-pos[t-1][0][i]))/diss;
						double sin=((pos[t-1][1][j]-pos[t-1][1][i]))/diss;
						double cos_z=((pos[t-1][2][j]-pos[t-1][2][i]))/diss;
						
								
						x_force[i]+=cos*((k*k)/(-(diss)));
						y_force[i]+=sin*((k*k)/(-(diss)));
						if(t3d)
							z_force[i]+=cos_z*((k*k)/(-(diss)));
					}
				}
					
			if(max_displacemnet[0]<x_force[i])
				max_displacemnet[0]=x_force[i];
					
			if(max_displacemnet[1]<y_force[i])
				max_displacemnet[1]=y_force[i];
			if(t3d)				
				if(max_displacemnet[2]<z_force[i])
					max_displacemnet[2]=z_force[i];
					
		}
			
			
			
			
			for(int i=0; i<links;i++)
			{
				int r=adj[0][i];
				int s=adj[1][i];
				double diss= Math.sqrt((Math.pow(pos[t-1][0][r]-pos[t-1][0][s], 2)) + ( Math.pow(pos[t-1][1][r]-pos[t-1][1][s], 2)));
				
				if(t3d)
					diss = Math.sqrt((Math.pow(pos[t-1][0][r]-pos[t-1][0][s], 2)) + ( Math.pow(pos[t-1][1][r]-pos[t-1][1][s], 2))+
							(Math.pow(pos[t-1][2][r]-pos[t-1][2][s], 2)));
				else
					diss= Math.sqrt((Math.pow(pos[t-1][0][r]-pos[t-1][0][s], 2)) + ( Math.pow(pos[t-1][1][r]-pos[t-1][1][s], 2)));
				
				
				double cos=((pos[t-1][0][r]-pos[t-1][0][s]))/diss;
				double sin=((pos[t-1][1][r]-pos[t-1][1][s]))/diss;
				double cos_z=((pos[t-1][2][r]-pos[t-1][2][s]))/diss;
				
				x_force[s]+=cos*(((diss)*(diss))/k);
				y_force[s]+=sin*(((diss)*(diss))/k);
				if(t3d)	
				z_force[s]+=cos_z*(((diss)*(diss))/k);
				
				x_force[r]-=cos*(((diss)*(diss))/k);
				y_force[r]-=sin*(((diss)*(diss))/k);
				if(t3d)	
				z_force[r]-=cos_z*(((diss)*(diss))/k);
				
				
				
				
				if(max_displacemnet[0]<x_force[r])
					max_displacemnet[0]=x_force[r];
						
				if(max_displacemnet[1]<y_force[r])
					max_displacemnet[1]=y_force[r];
				
				if(max_displacemnet[0]<x_force[s])
					max_displacemnet[0]=x_force[s];
						
				if(max_displacemnet[1]<y_force[s])
					max_displacemnet[1]=y_force[s];
				
				if(t3d)	
				{
					if(max_displacemnet[2]<z_force[s])
						max_displacemnet[2]=z_force[s];
							
					if(max_displacemnet[2]<z_force[r])
						max_displacemnet[2]=z_force[r];
				}
			}
			
				
				for(int i=0;i<v;i++)
				{	
					if(!(i==1 || i==2 ))
					if(animation && v<=100)
					{
						pos[t][0][i]=(pos[t-1][0][i]+((x_force[i])*(c4))/max_displacemnet[0]);
						pos[t][1][i]=(pos[t-1][1][i]+((y_force[i])*(c4))/max_displacemnet[1]);
						if(t3d)	
						pos[t][2][i]=(pos[t-1][2][i]+((z_force[i])*(c4))/max_displacemnet[2]);
					}
					else
					{
						pos[t-1][0][i]=(pos[t-1][0][i]+((x_force[i])*(c4))/max_displacemnet[0]);
						pos[t-1][1][i]=(pos[t-1][1][i]+((y_force[i])*(c4))/max_displacemnet[1]);
						if(t3d)	
						pos[t-1][2][i]=(pos[t-1][2][i]+((z_force[i])*(c4))/max_displacemnet[2]);
					}
					
				}
				
					c4-=aneal_rate;
				//variables.forces_in_synchronization[1][t]=c4;
				//System.out.println(variables.forces_in_synchronization[0][t]);
		}
		
		end_point=System.nanoTime();
		
		
		
		
		
		pos=main_fitting_in_center(v, pos, (int) time, !animation,t3d );
		pos=limit_in_screen(v, pos, (int) time, !animation, t3d);
		pos=covert_to_glortho(v, pos, (int) time, !animation, t3d);
		
		File file = new File("FR paddy.txt");
		PrintWriter in = new PrintWriter(file);
	
		for(int i=0;i<v;i++){
			in.println(i+"\t"+pos[0][0][i]+"\t"+pos[0][1][i]);
			
		}
		in.close();
		
		
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
		double lm=0.9;
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
		
		return pos;
	}
	



	public static double [][][] covert_to_glortho(int v, double [][][] pos, int time, boolean stat, boolean t3d)
	{
		double scale=0.7;
		
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

	/*static double [] btw_cnt(int edge[][], int v, int link)
	{
		Graph g = new UndirectedSparseGraph();
		int lk=0;
		for(int i=0;i<v;i++)
		{	
			g.addVertex(i);
		}
		
		for(int i=0;i<link;i++)
		{
			g.addEdge(i+"",edge[0][i],edge[1][i]);
		}
		
				
		DijkstraShortestPath alg = new DijkstraShortestPath(g);
		
		ArrayList<Integer> [][] shortes  =  (ArrayList<Integer>[][]) new  ArrayList [v][v] ;
		
		for(int i=0;i<v;i++)
		{
			for(int j=0;j<v;j++)
			{
				List k  = alg.getPath(i, j);
				shortes[i][j] = new ArrayList<Integer>();
				if(i!=j)
				{
					for(int q=0;q<k.size();q++)
					{	
						int a1= Integer.parseInt((String) k.get(q));
						shortes [i][j].add(edge[0][a1]);
						shortes [i][j].add(edge[1][a1]);
					}
				}
			}
		}
		
		for(int i=0;i<v;i++)
			for(int j=0;j<v;j++)
				for(int q=0;q<shortes [i][j].size();q++)
					for(int p=q+1;p<shortes [i][j].size();p++)
						if(shortes [i][j].get(q) == shortes [i][j].get(p))
							shortes [i][j].remove(p);
		
		double centr [] = new double [v];
		largest_btw_cnt=0;
		smallest_btw_cnt=0;
		
		for(int i=0;i<v;i++)
			for(int j=i+1;j<v;j++)
				for(int q=0;q<shortes[i][j].size();q++)
					centr[shortes[i][j].get(q)]++;
		
		for(int i=0;i<v;i++)
			if(largest_btw_cnt<centr[i])
				largest_btw_cnt=centr[i];
		
		smallest_btw_cnt = largest_btw_cnt;
		
		for(int i=0;i<v;i++)
			if(smallest_btw_cnt>=centr[i] && centr[i]!=0)
				smallest_btw_cnt=centr[i];
		
		//**************STND
		double ave =0;
		for(int i=0;i<v;i++)
			 ave += centr[i]/largest_btw_cnt;
		
		ave = ave /(double)v;
		
		double var = 0;
		for(int i=0;i<v;i++)
			var += (ave-centr[i]/largest_btw_cnt)*(ave-centr[i]/largest_btw_cnt);
			
		var = var/ (double)v;
	
		
		
		for(int i=0;i<v;i++)
			System.out.println(centr[i]);
		return centr;
	}
	*/

}
