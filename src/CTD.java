import java.awt.FontFormatException;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;

import javax.swing.JOptionPane;

import Jama.Matrix;



public class CTD {

	static int longest_theoritical_disctance=0;
	static int initial=0;	
	public static double sumdis [];
	public static double doublemat [][];
	public static double logmat [][];
	public static int f1 =0;
	public static int f2 =0;
	private static ArrayList <Integer> floors  = new ArrayList <Integer>() ;
	
	public static int [][] met;
	
	public static double [][] shpth(double[][] adj, int v, int link) throws FileNotFoundException
	   {
		
		  double path [][] = new double [v][v];
		     sumdis = new double [v];
	           int totalNodes = v;
	           int totalEdges = link;

	           Map<Integer, List<Integer>> adjacencyList = new HashMap<Integer, List<Integer>>();
	           for (int i = 0; i < v; i++)
	           {
	           for (int j = i+1; j < v; j++)
	           {
	               if(adj[i][j]==1)
	               {
	        	   
	        	   int src = i;
	               int dest = j;

	               if (adjacencyList.get(src) == null)
	               {
	                   List<Integer> neighbours = new ArrayList<Integer>();
	                   neighbours.add(dest);
	                   adjacencyList.put(src, neighbours);
	               } else
	               {
	                   List<Integer> neighbours = adjacencyList.get(src);
	                   neighbours.add(dest);
	                   adjacencyList.put(src, neighbours);
	               }


	               if (adjacencyList.get(dest) == null)
	               {
	                   List<Integer> neighbours = new ArrayList<Integer>();
	                   neighbours.add(src);
	                   adjacencyList.put(dest, neighbours);
	               } else
	               {
	                   List<Integer> neighbours = adjacencyList.get(dest);
	                   neighbours.add(src);
	                   adjacencyList.put(dest, neighbours);
	               }
	               }
	           }
	           }
	       	longest_theoritical_disctance=(int)path[0][0];
	           for(int i=0;i<v;i++)
	           {
		           int start =i;
	
		           Queue<Integer> queue = new LinkedList<Integer>();
	
		           queue.add(start);
	
		           int[] costs = new int[totalNodes + 1];
	
		           Arrays.fill(costs, 0);
	
		           costs[start] = 0;
	
		           Map<String, Integer> visited = new HashMap<String, Integer>();
	
		           while (!queue.isEmpty())
		           {
		               int node = queue.remove();
	
		               if(visited.get(node +"") != null)
		               {
		                   continue;
		               }
	
		               visited.put(node + "", 1);
	
		               int nodeCost = costs[node];
	
		               List<Integer> children = adjacencyList.get(node);
	
		               if (children != null)
		               {
		                   for (Integer child : children)
		                   {
		                       int total = nodeCost + 6;
		                       String key = child + "";
	
		                       if (visited.get(key) == null)
		                       {
		                           queue.add(child);
	
		                           if (costs[child] == 0)
		                           {
		                               costs[child] = total;
		                           } else if (costs[child] > total)
		                           {
		                               costs[child] = total;
		                           }
		                       }
		                   }
		               }
		           }
	
		           for (int k = i+1; k < v; k++)
		           {
		              /* if (k == start)
		               {
		                   continue;
		               }*/
		        	   path[i][k]= costs[k]/6;
		        	   path[k][i]=path[i][k];
		        	   if(longest_theoritical_disctance<path[i][k])
		        	   {
		        		  longest_theoritical_disctance=(int) path[i][k];
		        		  f1 = i;
		        		  f2 = k;
		        	   }
		           }
		           
		           for (int k = 0; k < v; k++)
		           {
		        	   sumdis[i] += path[i][k];
		           }
	           
	           
	           }
	           System.out.println("long    "+longest_theoritical_disctance);
	          
	           
	           return path;
	   }
	
	 public static int [][] mult (int [][] lis, int v, double [][] mat, int [] visited) throws FileNotFoundException
		{
			
			
			int counter = 0;
			double lon =longest_theoritical_disctance;
			boolean exit = true;
			/*for(int i=0;i<v && exit;i++)
				for(int j=i+1;j<v && exit;j++)
					if((lon==mat[i][j] && visited[i]==0 && visited[j]==0 && initial==0) || 
							(lon==mat[i][j] && visited[i]==0 && visited[j]==0 && initial==2 && mat[i][lis[0][0]] == mat[i][lis[0][1]] ))
					{
						if(initial<3)
						{
							visited[i]=1;
							lis[0][counter]=i;
							lis[1][counter]=-1;
							initial++;
							counter++;
							System.out.println("i: "+i);
						}
						if(initial<3)
						{
							visited[j]=1;
							lis[0][counter]=j;
							lis[1][counter]=-1;
							initial++;
							counter++;
							System.out.println("j: "+j);
						}
						if(initial==3)
							exit=false;
					}*/
			
			initial = 2;
			lis[0][counter]=f1;
			lis[1][counter]=-1;
			counter++;
			lis[0][counter]=f2;
			lis[1][counter]=-1;
			counter++;
			visited[f1]=1;
			visited[f2]=1;
			
			
			if(initial ==2 )
			{
				for(int i=1;i<v;i++)
				{
					if((mat[i][lis[0][0]] == mat[i][lis[0][1]]) || (mat[i][lis[0][0]] == mat[i][lis[0][1]]-1) || (mat[i][lis[0][0]]-1 == mat[i][lis[0][1]]))
					{
						visited[i]=1;
						lis[0][counter]=i;
						lis[1][counter]=-1;
						initial++;
						counter++;
						System.out.println("i1: "+i);
						break;
					}
				}
			}
			floors.add(3);
			boolean check = true;
			int level =1;
			while(check)
			{
				int kk=0;
				for(int i=0;i<v;i++)
				{
					if(visited[i] == level)
					for(int j=0;j<v;j++)
					{
						if(visited[j]==0 && mat[i][j]==1)
						{
							//System.out.println(i +"  "+j);
							visited[j]=level+1;
							lis[0][counter]=j;
							lis[1][counter]=i;
							kk++; // no. of nodes at each floor;
							counter++;
						}
					}
				}
				int c = floors.get(floors.size()-1);
				floors.add(kk+c);
				level++;
				if(counter==v)
					check=false;
			}
			
			return lis;	
			
			
		}

		public static double [][][] Concentric_drawing_distance( int adj [][], int v, double time, boolean t3d, double links, boolean animation) throws FileNotFoundException
		{
	 		System.out.println("Start");
			long start_point,end_point;
			start_point= System.nanoTime();
			GraphMethods.preparing_for_graphviz(adj,  v, (int)links);
			
			long D1= System.nanoTime();
			
				double mat [][] = new double [v][v];
				
				for(int i=0;i<links;i++)
				{
					mat[(int)adj[0][i]][(int)adj[1][i]]=1;
					mat[(int)adj[1][i]][(int)adj[0][i]]=1;
				}
			long D2= System.nanoTime();
			
			
		
			/*System.out.println();System.out.println();
			for(int i=0;i<v;i++)
			{
				for(int j=0;j<v;j++)
				{
					System.out.print((int)mat[i][j]+" ");
				}System.out.println();
			}*/
			
			mat=shpth(mat, v, (int)links);
			
			
			
			
			////////////////
			/* File simi = new File("C:/Users/farshad.toosi.UL/Documents/matrices/for graphviz/tessttt.txt");
				PrintWriter similar = new PrintWriter(simi);
				for(int i=0;i<v;i++)
				{
					for(int j=0;j<v;j++)
					{
						similar.print(mat[i][j]+" ");
					}similar.println();
				}
			   similar.close();*/
			/////////////
			
			
			long disss= System.nanoTime();
			
			doublemat = new double [v][v];
			logmat = new double [v][v];
			
			boolean options = true;
			if(options)
			for(int i=0;i<v;i++)
			{
				for(int j=i+1;j<v;j++)
				{
					doublemat[i][j] = mat[i][j]*mat[i][j]*mat[i][j];
					//doublemat[i][j] = Math.sqrt(mat[i][j]);
					
					doublemat[j][i] = doublemat[i][j];
					logmat[i][j] = Math.log( mat[i][j]);
					logmat[j][i] = logmat[i][j];
					
					/*doublemat[i][j] = mat[i][j];
					doublemat[j][i] = doublemat[i][j];*/
					
				}
			}
			
			double n=v;
			int [] pp = new int [(int)n];
			int []visited = new int [(int)n];
			int [][] lis = new int [2][(int)n];
				
			int popu=0;
			
		
			
			//time = v/2 ;
			double lo = Math.log(v);
			
			
			long h1 = System.nanoTime();
			lis = mult (lis, v, mat, visited);
			long h2 = System.nanoTime();
			
			System.out.println("initial  "+initial+"   "+(h2-h1));

			double [][][] pos; 
			
			if(!animation) 
				pos= new double [1][3][(int)n];
			else
				pos= new double [(int) time][3][(int)n];
			
			
			
			double chunki = (Math.PI*2)/3;
			pos[0][0][lis[0][0]]=Math.cos(0);
			pos[0][1][lis[0][0]]=Math.sin(0);
			if(t3d)
				pos[0][1][lis[0][0]]=Math.sin(0);
			met = new int [(int) time][v];
			for(int i=0;i<time;i++){
			met[i][lis[0][0]]=1;
			met[i][lis[0][1]]=1;
			met[i][lis[0][2]]=1;
			}
			
			for(int i=1;i<initial;i++)
			{
				
				pos[0][0][lis[0][i]]=Math.cos(i*chunki);
				pos[0][1][lis[0][i]]=Math.sin(i*chunki);
				if(t3d)
					pos[0][2][lis[0][i]]=chunki/(Math.PI*2);
			}
			
			System.out.println("lll  : "+mat[lis[0][2]][lis[0][1]]+"  "+mat[lis[0][2]][lis[0][0]]+"  "+mat[lis[0][0]][lis[0][1]]);
			
			int fg=(int) initial;
			int gt=(int) initial +  10;//50;
			popu = (int) initial;
			System.out.println("initial "+initial+" log : "+ lo);
			double disp = 10.0/(v*v*v);
			boolean tri=false;
			int level = 0;
			int cc = 0;
			int t=0;
			int lt=0;
			int iter = 0;
			while(popu<(v+lo*5))  // timmeeeeeeeeeeemmmmmmmmmmmmeeeeeeeeeeeeeeeemmmmmmmmmmmme
			{
				iter++;
				//System.out.println(iter);
				if(animation) {
					t++;
					lt = t-1;
				}
				
				cc++;
				for(int i=fg;i<gt;i++)
				{
					
					pos[t][0][lis[0][i]]=pos[lt][0][lis[1][i]]+disp;
					pos[t][1][lis[0][i]]=pos[lt][1][lis[1][i]]-disp;
					for(int i1=t;i1<time;i1++){
						met[i1][lis[0][i]]=1;
					}
					
					if(t3d)
						pos[t][2][lis[0][i]]=pos[lt][2][lis[1][i]]+disp;
					popu++;
				}
				
				//******************************
				
				double x_forces[]=new double [v];
				double y_forces[]=new double [v];
				double z_forces[]=new double [v];
				
				double diss;
				
				for(int f1=0;f1<gt;f1++)
				{
					int i=lis[0][f1];
					for(int t1=f1+1;t1<gt;t1++)
					{
						int j=lis[0][t1];
						
						
						double p = pos[lt][0][j]-pos[lt][0][i];
						double q = pos[lt][1][j]-pos[lt][1][i];
						
						if(!tri)
						{
							if(!t3d)
								diss= Math.sqrt(p*p	+ q*q);
							else
								diss= Math.sqrt((pos[lt][0][i]-pos[lt][0][j])*(pos[lt][0][i]-pos[lt][0][j])
										+ (pos[lt][1][i]-pos[lt][1][j])*(pos[lt][1][i]-pos[lt][1][j])+
										(pos[lt][2][i]-pos[lt][2][j])*(pos[lt][2][i]-pos[lt][2][j]));
						}
						else
						{
							double p1 = p,q1 = q;
							
							if(p1<0)
								p1 = -p1;
							
							if(q1<0)
								q1 = -q1;
							
							if(!t3d)
								diss= (p1 + q1);
							else
								diss= (Math.abs(pos[lt][0][i]-pos[lt][0][j])+Math.abs(pos[lt][1][i]-pos[lt][1][j])+
									Math.abs(pos[lt][2][i]-pos[lt][2][j]));
						}
							if(diss<0.00000000002)
							{
								diss=0.00000000002;
								diss=0.00000000002;
							}
							//-----------
							double sin=(q/diss);
							double cos=(p/diss);
							
							double zSin=0;
							if(t3d)
								zSin=((pos[lt][2][j]-pos[lt][2][i])/diss);
					
						
							double x=0;
							double y=0;
							double z=0;
							
							if(popu>=v)
							{
								x=cos*((mat[i][j]));
								y=sin*((mat[i][j]));
								//System.out.print("ss");
								/*x=cos*(Math.log(mat[i][j]));
								y=sin*(Math.log(mat[i][j]));*/
								
								/*x=cos*(logmat[i][j]);
								y=sin*(logmat[i][j]);*/
							}
							else
							{
								x=cos*((mat[i][j]));
								y=sin*((mat[i][j]));
								//System.out.print("dds");
								/*x=cos*(Math.log(mat[i][j]));
								y=sin*(Math.log(mat[i][j]));*/
								x=cos*(doublemat[i][j]);
								y=sin*(doublemat[i][j]);
							}
							
							
							x_forces[i]+=-x;
							y_forces[i]+=-y;
							if(t3d)
								z_forces[i] +=-z;
							
							if(j % 2 ==0) {
							x_forces[j]+=x;
							y_forces[j]+=y;
							}
							else {
								x_forces[j]+=1*x;
								y_forces[j]+=1*y;
							}
							if(t3d)
								z_forces[j] +=z;
							
						}
				}
				
				double en2=0;
				for(int f1=0;f1<gt;f1++)
				{
					int i= lis[0][f1];
					pos[t][0][i]=(x_forces[i]);
					pos[t][1][i]=(y_forces[i]);
					if(t3d)
						pos[t][2][i]=(z_forces[i]);
					en2 += pos[t][0][i]+pos[t][1][i];
					
					
				}
				/*if(en1!=0)
					stress = Math.abs((en2-en1)/en1);
				System.out.println(stress+"   "+popu);
				if(stress<0.1 && popu>=v)
					break;
				en1=en2;*/
				//System.out.println(" grav:  "+en2);
				if(fg<v)
					fg = floors.get(level);
				if(gt<v)
					gt = floors.get(level+1);
				level++;
				
				if(popu>=v)
					popu++;
				
				if(gt>v)
					gt=(int) v;
				
			}
			
			end_point=System.nanoTime();
			long running_time = (end_point-start_point);

			long lll = System.nanoTime();
			long Final_time = (lll-disss)+(D2-D1);
			System.out.println("time only algorithm :  "+(Final_time)+"  node:  "+v+"  edges:  "+links);
			
			System.out.println("New time for this:    "+(running_time)+" itess:  "+cc);
			
			
			System.out.println("\nRunning time:  for a graph with:  "+v+"  and links: "+links+" is:   "+(running_time));
			
			
		
			
			pos=main_fitting_in_center(v, pos, (int) time, !animation,t3d );
			pos=limit_in_screen(v, pos, (int) time, !animation, t3d);
			pos=covert_to_glortho(v, pos, (int) time, !animation, t3d);
			
			return pos;
		}
		
	
	 
	 
	public static double [][][] Concentric_drawing_distance1( int adj [][], int v, double time, boolean t3d, double links, boolean animation) throws FileNotFoundException
		{
	 		System.out.println("Start");
			long start_point,end_point;
			start_point= System.nanoTime();
			//GraphMethods.preparing_for_graphviz(adj,  v, (int)links);
			
			long D1= System.nanoTime();
			
				double mat [][] = new double [v][v];
				
				for(int i=0;i<links;i++)
				{
					mat[(int)adj[0][i]][(int)adj[1][i]]=1;
					mat[(int)adj[1][i]][(int)adj[0][i]]=1;
				}
			long D2= System.nanoTime();
			
			
		
			/*System.out.println();System.out.println();
			for(int i=0;i<v;i++)
			{
				for(int j=0;j<v;j++)
				{
					System.out.print((int)mat[i][j]+" ");
				}System.out.println();
			}*/
			
			mat=shpth(mat, v, (int)links);
			
			
			
			
			////////////////
			/* File simi = new File("C:/Users/farshad.toosi.UL/Documents/matrices/for graphviz/tessttt.txt");
				PrintWriter similar = new PrintWriter(simi);
				for(int i=0;i<v;i++)
				{
					for(int j=0;j<v;j++)
					{
						similar.print(mat[i][j]+" ");
					}similar.println();
				}
			   similar.close();*/
			/////////////
			
			
			long disss= System.nanoTime();
			
			doublemat = new double [v][v];
			logmat = new double [v][v];
			
			boolean options = true;
			if(options)
			for(int i=0;i<v;i++)
			{
				for(int j=i+1;j<v;j++)
				{
					doublemat[i][j] = mat[i][j]*mat[i][j]*mat[i][j];
					//doublemat[i][j] = Math.sqrt(mat[i][j]);
					
					doublemat[j][i] = doublemat[i][j];
					logmat[i][j] = Math.log( mat[i][j]);
					logmat[j][i] = logmat[i][j];
					
					/*doublemat[i][j] = mat[i][j];
					doublemat[j][i] = doublemat[i][j];*/
					
				}
			}
			
			double n=v;
			int [] pp = new int [(int)n];
			int []visited = new int [(int)n];
			int [][] lis = new int [2][(int)n];
				
			int popu=0;
			
		
			
			//time = v/2 ;
			double lo = Math.log(v);
			
			
			long h1 = System.nanoTime();
			//lis = mult (lis, v, mat, visited);
			long h2 = System.nanoTime();
			
			System.out.println("initial  "+initial+"   "+(h2-h1));

			double [][][] pos; 
			
			if(!animation) 
				pos= new double [1][3][(int)n];
			else
				pos= new double [(int) time][3][(int)n];
			
			
			
			double chunki = (Math.PI*2)/v;
			pos[0][0][lis[0][0]]=Math.cos(0);
			pos[0][1][lis[0][0]]=Math.sin(0);
			if(t3d)
				pos[0][1][lis[0][0]]=Math.sin(0);
			met = new int [(int) time][v];
			for(int i=0;i<time;i++){
			met[i][lis[0][0]]=1;
			met[i][lis[0][1]]=1;
			met[i][lis[0][2]]=1;
			}
			
			for(int i=0;i<v;i++)
			{
				
				pos[0][0][i]=Math.random()*v;
				pos[0][1][i]=Math.random()*v;
				if(t3d)
					pos[0][2][lis[0][i]]=chunki/(Math.PI*2);
			}
			
			System.out.println("lll  : "+mat[lis[0][2]][lis[0][1]]+"  "+mat[lis[0][2]][lis[0][0]]+"  "+mat[lis[0][0]][lis[0][1]]);
			
			int fg=(int) initial;
			int gt=(int) initial +  10;//50;
			popu = 0;
			System.out.println("initial "+initial+" log : "+ lo);
			double disp = 10.0/(v*v*v);
			boolean tri=false;
			int level = 0;
			int cc = 0;
			int t=0;
			int lt=0;
			int iter = 50;
			while(popu<(iter))  // timmeeeeeeeeeeemmmmmmmmmmmmeeeeeeeeeeeeeeeemmmmmmmmmmmme
			{
				popu++;
				
				//System.out.println(iter);
				if(animation) {
					t++;
					lt = t-1;
				}
				
				cc++;
				
				
				//******************************
				
				double x_forces[]=new double [v];
				double y_forces[]=new double [v];
				double z_forces[]=new double [v];
				
				double diss;
				
				for(int i=0;i<v;i++)
				{
					for(int j=i+1;j<v;j++)
					{
						
						
						double p = pos[lt][0][j]-pos[lt][0][i];
						double q = pos[lt][1][j]-pos[lt][1][i];
						
						if(!tri)
						{
							if(!t3d)
								diss= Math.sqrt(p*p	+ q*q);
							else
								diss= Math.sqrt((pos[lt][0][i]-pos[lt][0][j])*(pos[lt][0][i]-pos[lt][0][j])
										+ (pos[lt][1][i]-pos[lt][1][j])*(pos[lt][1][i]-pos[lt][1][j])+
										(pos[lt][2][i]-pos[lt][2][j])*(pos[lt][2][i]-pos[lt][2][j]));
						}
						else
						{
							double p1 = p,q1 = q;
							
							if(p1<0)
								p1 = -p1;
							
							if(q1<0)
								q1 = -q1;
							
							if(!t3d)
								diss= (p1 + q1);
							else
								diss= (Math.abs(pos[lt][0][i]-pos[lt][0][j])+Math.abs(pos[lt][1][i]-pos[lt][1][j])+
									Math.abs(pos[lt][2][i]-pos[lt][2][j]));
						}
							if(diss<0.00000000002)
							{
								diss=0.00000000002;
								diss=0.00000000002;
							}
							//-----------
							double sin=(q/diss);
							double cos=(p/diss);
							
							double zSin=0;
							if(t3d)
								zSin=((pos[lt][2][j]-pos[lt][2][i])/diss);
					
						
							double x=0;
							double y=0;
							double z=0;
							
							if(popu>=iter/2)
							{
								x=cos*((mat[i][j]));
								y=sin*((mat[i][j]));
								//System.out.print("ss");
								/*x=cos*(Math.log(mat[i][j]));
								y=sin*(Math.log(mat[i][j]));*/
								
								/*x=cos*(logmat[i][j]);
								y=sin*(logmat[i][j]);*/
							}
							else
							{
								x=cos*(doublemat[i][j]);
								y=sin*(doublemat[i][j]);
							}
							
							
							x_forces[i]+=-x;
							y_forces[i]+=-y;
							if(t3d)
								z_forces[i] +=-z;
							
							x_forces[j]+=x;
							y_forces[j]+=y;
							if(t3d)
								z_forces[j] +=z;
							
						
							
						}
				}
				
				for(int i=0;i<v;i++)
				{
					pos[t][0][i]=(x_forces[i]);
					pos[t][1][i]=(y_forces[i]);
				}
				
				//System.out.println(x_forces[505]+" "+y_forces[505]);
				double en2=0;
				
				
			}
			
			end_point=System.nanoTime();
			long running_time = (end_point-start_point);

			long lll = System.nanoTime();
			long Final_time = (lll-disss)+(D2-D1);
			System.out.println("time only algorithm :  "+(Final_time)+"  node:  "+v+"  edges:  "+links);
			
			System.out.println("New time for this:    "+(running_time)+" itess:  "+cc);
			
			
			System.out.println("\nRunning time:  for a graph with:  "+v+"  and links: "+links+" is:   "+(running_time));
			
			
		
			
			pos=main_fitting_in_center(v, pos, (int) time, !animation,t3d );
			pos=limit_in_screen(v, pos, (int) time, !animation, t3d);
			pos=covert_to_glortho(v, pos, (int) time, !animation, t3d);
			
			return pos;
		}
		
	 	
	 	
	 	public static double [][][] Concentric_drawing_distance_adding( double mat [][], int v, double time, boolean t3d, double links, boolean animation, double [][][] pos, int newnodes) throws FileNotFoundException, FontFormatException
		{
			long start_point,end_point ;
			
			
			
			mat=shpth(mat, v, (int)links);
			
			start_point= System.nanoTime();
			System.out.println(start_point);
			int popu=0;
			System.out.println("initial  addinggg  "+initial);

			int upertime = 0;
			
			if(newnodes<5)
				upertime=10;
			else
				if(newnodes<10)
					upertime=15;
				else
					upertime=20;
				
			System.out.println("Time  :  "+upertime+"  "+ newnodes+"  "+v);
			
			
			while(popu<(upertime))  // timmeeeeeeeeeeemmmmmmmmmmmmeeeeeeeeeeeeeeeemmmmmmmmmmmme
			{
				
				//******************************
				
				double x_forces[]=new double [v];
				double y_forces[]=new double [v];
				double z_forces[]=new double [v];
				
				double diss;
				double [] cf =new double [v];
				double lo = Math.sqrt(longest_theoritical_disctance);
				for(int i=0;i<v;i++)
				{
					
					for(int j=i+1;j<v;j++)
					{
						if(!t3d)
							diss= Math.sqrt((pos[0][0][i]-pos[0][0][j])*(pos[0][0][i]-pos[0][0][j])
									+ (pos[0][1][i]-pos[0][1][j])*(pos[0][1][i]-pos[0][1][j]));
						else
							diss= Math.sqrt((pos[0][0][i]-pos[0][0][j])*(pos[0][0][i]-pos[0][0][j])
									+ (pos[0][1][i]-pos[0][1][j])*(pos[0][1][i]-pos[0][1][j])+
									(pos[0][2][i]-pos[0][2][j])*(pos[0][2][i]-pos[0][2][j]));
						
							if(diss<0.00000000002)
							{
								diss=0.00000000002;
								diss=0.00000000002;
							}
							//-----------
							double sin=((pos[0][1][j]-pos[0][1][i])/diss);
							double cos=((pos[0][0][j]-pos[0][0][i])/diss);
							double zSin=0;
							if(t3d)
								zSin=((pos[0][2][j]-pos[0][2][i])/diss);
					
						
							double x=0;
							double y=0;
							double z=0;
							
							if(mat[i][j]<lo)
							{
								cf[i] +=2;
								cf[j] +=2;
								x=cos*(mat[i][j]*2);
								y=sin*(mat[i][j]*2);
								z=zSin*(mat[i][j]*2);
								
								x_forces[i]+=-x;
								y_forces[i]+=-y;
								if(t3d)
									z_forces[i] +=-z;
								x_forces[j]+=x;
								y_forces[j]+=y;
								if(t3d)
									z_forces[j] +=z;
							}
							else
							{
								cf[i]++;
								cf[j]++;
							x=cos*(mat[i][j]);
							y=sin*(mat[i][j]);
							z=zSin*(mat[i][j]);
							
							x_forces[i]+=-x;
							y_forces[i]+=-y;
							if(t3d)
								z_forces[i] +=-z;
							x_forces[j]+=x;
							y_forces[j]+=y;
							if(t3d)
								z_forces[j] +=z;
							}
						}
				}
				
				
				for(int i=0;i<v;i++)
				{
					
					pos[0][0][i]=(x_forces[i]/cf[i]);
					pos[0][1][i]=(y_forces[i]/cf[i]);
					if(t3d)
						pos[0][2][i]=(z_forces[i]/cf[i]);
					
				}
				popu++;
				
			}
			
			
			
			
			end_point= System.nanoTime();
			System.out.println(end_point);
			long res =(end_point-start_point);
			System.out.println("Time for ne added:  "+res);
			
			
			pos=main_fitting_in_center(v, pos, (int) time, !animation,t3d );
			pos=limit_in_screen(v, pos, (int) time, !animation, t3d);
			pos=covert_to_glortho(v, pos, (int) time, !animation, t3d);
			graph CTD = new graph( pos,false, graph.adj,(int)((double)v/4), false, v,  links, false,5);
		
			
			
			CTD.visualization();
			return pos;
		}
		
	 	
	 	public static double [][][] adding_node(double [][][] oldpos, int v, ArrayList<ArrayList<Integer>> nodes, double [][] oldmat, int [][] oldadj, double links, boolean t3d) throws FileNotFoundException, FontFormatException
	 	{
	 		int newnodes = v+nodes.size();
	 		double [][] mat = new double [newnodes][newnodes];
	 		double [][][] pos = new double [1][3][newnodes];
	 		for(int i=0;i<v;i++)
	 		{
	 			pos[0][0][i]=oldpos[0][0][i];
	 			pos[0][1][i]=oldpos[0][1][i];
	 			pos[0][2][i]=oldpos[0][2][i];
	 			
	 			for(int j=i+1;j<v;j++)
		 		{
		 			mat[i][j]=oldmat[i][j];
		 			mat[j][i]=oldmat[j][i];
		 		}
	 		}
	 		double dis = 0.0000000001;
	 		for(int i=v;i<newnodes;i++)
	 		{
	 			for(int j=0;j<nodes.get(i-v).size();j++)
		 		{
		 			mat[i][nodes.get(i-v).get(j)]=1;
		 			mat[nodes.get(i-v).get(j)][i]=mat[i][nodes.get(i-v).get(j)];
		 			links++;
		 		}
	 			pos[0][0][i]=pos[0][0][nodes.get(i-v).get(0)]+Math.cos(Math.random()*Math.PI*2)*dis;
	 			pos[0][1][i]=pos[0][1][nodes.get(i-v).get(0)]+Math.sin(Math.random()*Math.PI*2)*dis;
	 			if(t3d)
	 				pos[0][2][i]=pos[0][2][nodes.get(i-v).get(0)]+dis;
	 		}
	 		int [][] adj = new int [2][(int) links];
	 		int cf=0;
	 		for(int i=0;i<newnodes;i++)
	 		{
	 			for(int j=i+1;j<newnodes;j++)
		 		{
	 				if(mat[i][j]==1)
	 				{
	 					adj[0][cf]=i;
	 					adj[1][cf]=j;
	 					cf++;
	 				}
		 		}
	 		}
	 		graph.adj=adj.clone();
	 		
	 		 		
	 		Concentric_drawing_distance_adding(  mat, newnodes,newnodes , t3d, links,false, pos, nodes.size());
	 		
	 		
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
			double scale=2.7;
			
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

		public static int check()
		{
			int addingnodes=0;
			int ans=addingnodes=JOptionPane.showConfirmDialog(null, "Do you want to add new nodes to the current layout?", "Adding new nodes", JOptionPane.YES_NO_OPTION);
        	if(ans==JOptionPane.YES_OPTION)
        		addingnodes=2;
        	else
        		addingnodes=1;
			return addingnodes;
		}
		public static void addings(double [][][] pos,  int v,  double [][] mat, int [][]adj, double links, boolean t3d) throws FileNotFoundException, FontFormatException
		{
        	ArrayList<ArrayList<Integer>> nodes = null;
        	
        		int q=-1;
        		while(q>50 || q<1)
        		{
        			try {
        				String p = JOptionPane.showInputDialog(null,"", "How many new nodes?",	JOptionPane.INFORMATION_MESSAGE);
            			q=Integer.parseInt(p);
        		    } 
        			catch (NumberFormatException e) 
        			{
        		        q = 0;
        		        JOptionPane.showMessageDialog(null, "eneter a value in range 1-10");
        		    }
        		}
        		nodes = new ArrayList<ArrayList<Integer>> ();
        		for(int i=0;i<q;i++)
        		{
        				String p = JOptionPane.showInputDialog(null,"", "Enter Nighburs for node  "+i+"  , E.g. 4-102-83",	JOptionPane.INFORMATION_MESSAGE);
            			String [] l =p.split("-");
            			ArrayList<Integer> n = new ArrayList<Integer>();
            			boolean check=true;
            			for(int j=0;j<l.length;j++)
            			{
            				try 
            				{	
            					n.add(Integer.parseInt(l[j]));
            				}
            				catch (NumberFormatException e) 
                			{
                		        JOptionPane.showMessageDialog(null, "eneter values separated by '-'");
                		        i--;
                		        check = false;
                		        break;
                		    }
            			}
            			if(check)
            				nodes.add(n);
        		    } 
        		
        		
        		
        		
        	adding_node(pos, v,  nodes, mat,  adj, links, t3d);
        	
        	}
        
		
}
