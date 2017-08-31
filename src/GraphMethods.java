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

import edu.uci.ics.jung.algorithms.importance.BetweennessCentrality;
import edu.uci.ics.jung.algorithms.shortestpath.DijkstraShortestPath;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;
import edu.uci.ics.jung.graph.event.GraphEvent.Edge;
import edu.uci.ics.jung.visualization.renderers.Renderer.Vertex;


public class GraphMethods {
	static double longest_theoritical_disctance=0;
	static double largest_btw_cnt=0;
	static double smallest_btw_cnt=0;
	static double [][] matrix_change_to_distance_cells(int edge[][], int v, int link)
	{
		double[][] temp = new double [v][v];
		Graph g = new UndirectedSparseGraph();
		
		for(int i=0;i<v;i++)
		{	
			g.addVertex(i);
		}
		for(int i=0;i<link;i++)
		{
			g.addEdge(i+"",edge[0][i],edge[1][i]);
		}
		
				
		BetweennessCentrality btw = new BetweennessCentrality(g);
		DijkstraShortestPath alg = new DijkstraShortestPath(g);
		
		ArrayList<String> name[] = new ArrayList[9];
		
		
		
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
						//System.out.print(k.get(q)+" : "+ edge[0][a1]+"  "+edge[1][a1]);
						
						
						shortes [i][j].add(edge[0][a1]);
						shortes [i][j].add(edge[1][a1]);
					}
				}
				//System.out.print(shortes[i][j]);
			}
			//System.out.println();
		}
		for(int i=0;i<v;i++)
		{
			for(int j=0;j<v;j++)
			{
				for(int q=0;q<shortes [i][j].size();q++)
				{
					for(int p=q+1;p<shortes [i][j].size();p++)
					{
						if(shortes [i][j].get(q) == shortes [i][j].get(p))
						{
							shortes [i][j].remove(p);
						}
					}
				}
			}
		}
		double longest =0;
		for(int i=0;i<v;i++)
		{
			for(int j=i+1;j<v;j++)
			{
				if(longest<shortes[i][j].size())
					longest=shortes[i][j].size();
			}
		}
		
		
		double centr [] = new double [v];
		for(int i=0;i<v;i++)
		{
			for(int j=i+1;j<v;j++)
			{
				if(i!=j)
				{
					for(int q=0;q<shortes[i][j].size();q++)
					{
						centr[shortes[i][j].get(q)]++;
					}
				}
			}
		}
		
		double aveg =0;
		double maxx=0;
		for(int i=0;i<v;i++)
		{
			centr[i] -= (v-1);
		}
		
		for(int i=0;i<v;i++)
		{
			if(maxx<centr[i])
				maxx=centr[i];
		}
		
		double slops [][] = new double [2][v];
		
		for(int i=0;i<v;i++)
		{
			//System.out.println(centr[i]);
			centr[i] /= maxx;
			aveg += centr[i];
			//System.out.println(centr[i]);
		}
		double pivot=0;
		for(int q=0;q<v;q++)
		{
			double min =v*v;
			int w=0;
			for(int i=0;i<v;i++)
			{
				if(min>centr[i] && slops[1][i]==0)
				{
					slops[0][q]=centr[i];
					w=i;
					min = centr[i];
				}
			}
			slops[1][w]=1;
		}
		
		
		
		
		//System.out.println();
		
		aveg = aveg / (double)v;
		double stnd = 0;
		double common = 0;
		double te=0;
		
		for(int i=v-1;i>=0;i--)
		{
			//common += slops[0][i+1]-slops[0][i];
			System.out.println("s"+"\t"+slops[0][i]);
		}
		
		for(int i=0;i<v;i++)
		{
			
			if(centr[i]>aveg)
			{
				//common +=centr[i];
				te++;
			}
				//common += centr[i];
			stnd += Math.pow((centr[i]-aveg), 2);
		}
		stnd /= (double) v;
		stnd = Math.sqrt(stnd);
		//System.out.println();System.out.println();
	 	//System.out.println("   Stnd cent:  "+(stnd)+"  lin: "+ (stnd/link)+"  :nod  "+(stnd/(double)v)+"  :long  "+(stnd/(double)longest));
	 	//System.out.println(common);
		////////////////
		int [][] adj = new int [v][v];
		
		for(int i=0;i<link;i++)
		{
			int a1= edge[0][i];
			int a2= edge[1][i];
			adj[a1][a2]=1;
			adj[a2][a1]=1;
		}
		
		
		temp= shpth( adj,  v,  link);
		
		
		double nodes_centerality [] = new double [v];
		
		
		double min= v*v*v*v;
		double longest_theoritical_disctance=0;
		double total_dist =0;
		double st=0;
		for(int i=0;i<v;i++)
		{
			for(int j=i+1;j<v;j++)
			{	
				total_dist += temp[i][j];
				if(temp[i][j]> longest_theoritical_disctance)
					 longest_theoritical_disctance=temp[i][j];
			}
		}
		total_dist /= (double)(v*(v-1)/2);
		
		for(int i=0;i<v;i++)
		{
			for(int j=i+1;j<v;j++)
			{
				st += Math.pow(total_dist-temp[i][j], 2);
			}
		}
		st /= (double)(v*(v-1)/2);
	
		//////////
		for(int i=0;i<v;i++)
		{
			double c=0;
			for(int j=0;j<v;j++)
			{
				c += temp[i][j]/longest_theoritical_disctance;
			}
			nodes_centerality[i] = c;
			//System.out.println(nodes_centerality[i]);
			total_dist += nodes_centerality[i];
		}
		
		total_dist /= (double) v;
		
		
		double maximum_value_for_centerality=0;
		double average_node_centrality =0;
		for(int i=0;i<v;i++)
		{
			nodes_centerality[i]-=(min-1);
			if( maximum_value_for_centerality< nodes_centerality[i])
				 maximum_value_for_centerality= nodes_centerality[i];
				
			average_node_centrality+= nodes_centerality[i];
			
		}
		
		double v_temp=v;
		average_node_centrality= average_node_centrality/(v_temp);
		
			int tem [] = new int [v];
			int res [] = new int [v];
			int c=1;
			int last_value=0;
			for(int i=0;i<v;i++)
			{
				int ind=0;
				int value=v*v*v;
				
				for(int j=0;j<v;j++)
				{
					if(value> nodes_centerality[j] && tem[j]==0)
					{
						ind=j;
						value= (int) nodes_centerality[j];
					}
				}
				
				if(value==last_value)
				{
					res[ind]=c-1;
				}
				else
				{
					res[ind]=c;
					c++;				
				}
				tem[ind]=-1;
				last_value=value;
			}
			 maximum_value_for_centerality=0;
			for(int i=0;i<v;i++)
			{
				 nodes_centerality[i]=res[i];
				if( maximum_value_for_centerality< nodes_centerality[i])
					 maximum_value_for_centerality= nodes_centerality[i];
				
			}
			
		
		
		return temp;
	}
	
	public static int [][] star_shaped_creating(int v, double deg)
	{
	
		int star [][] = new int [v][v];
		
		double sht [][] = new double [v][v];
		double btw [] = new double [v];
		
		double visited [] = new double [v];
		
		int a1 = (int)(Math.random()*v);
		int a2 = (int)(Math.random()*v);
		int a3 = (int)(Math.random()*v);
		int popu = 0;
		
		while(a1==a2 || a1==a2 || a2==a3)
		{
			a1 = (int)(Math.random()*v);
			a2 = (int)(Math.random()*v);
			a3 = (int)(Math.random()*v);
		}
		
		visited[a1]=1;visited[a2]=1; visited[a3]=1; popu++;popu++;popu++;
		star[a1][a2]=1;star[a2][a1]=1;
		star[a3][a1]=1;star[a1][a3]=1;
		star[a3][a2]=1;star[a2][a3]=1;
		int lnk =3;
		while(popu<v)
		{
			int temp_edge [][] = new int [2][lnk];
			int y=0;
			for(int i=0;i<v;i++)
			{
				for(int j=i+1;j<v;j++)
				{
					if(star[i][j]==1)
					{
						temp_edge[0][y]=i;
						temp_edge[1][y]=j;
						y++;
					}
				}
			}
			
			
			btw = btw_cnt(temp_edge, v, lnk);
			
			//System.out.println("1: "+lnk);
			
			double degree = (int)(Math.random()*deg);
			int n1 = (int)(Math.random()*v);
			
			while(visited[n1]==1)
				n1 = (int)(Math.random()*v);
			
			//System.out.println("2: "+lnk);
			
			boolean ch = true;
			double max = (largest_btw_cnt+1)*(largest_btw_cnt+1-1)/2;
			
			//System.out.println("max: "+max);
			
			while(ch)
			{
				
				int nigh = (int)(Math.random()*v);
				double chance = Math.random();
				//System.out.println(((btw[nigh]-largest_btw_cnt+1)/max)+"  "+chance+"  "+visited[nigh]+" "+max+" "+largest_btw_cnt);
				if(visited[nigh]==1)
				{
					if((btw[nigh]-largest_btw_cnt)/max <= chance)
					{
						if(star[n1][nigh]==0)
							lnk++;
						
						star[n1][nigh]=1;
						star[nigh][n1]=1;
						popu++;
						ch=false;
					}
				}
			}
			sht = shpth(star, v, 0);
			for(int i=0;i<degree;i++)
			{
				int ij=0;
				int j = (int)(Math.random()*v);
				while(visited[j]==0 && ij<3*v && sht[n1][j]>2)
				{
					j = (int)(Math.random()*v);
					ij++;
				}
				if(visited[j]==1 && sht[n1][j]<=2 && n1!=j)
				{
					if(star[n1][j]==0)
						lnk++;
					
					star[n1][j]=1;
					star[j][n1]=1;
					//popu++;
					break;
				}
			}
			visited[n1]=1;
		}
		
		int link=0;
		
		
		for(int i=0;i<v;i++)
		{
			for(int j=i+1;j<v;j++)
			{
				if(star[i][j]==1)
					link++;
			}
		}
		int edge [][] = new int [2][link];
		link = 0;
		for(int i=0;i<v;i++)
		{
			for(int j=i+1;j<v;j++)
			{
				if(star[i][j]==1)
				{
					edge[0][link]=i;
					edge[1][link]=j;
					link++;
				}
			}
		}
		
		return edge;
	}
	
	
	static double [] btw_cnt(int edge[][], int v, int link)
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
		
		
		return centr;
	}
	
	public static double [][] shpth(int[][] adj, int v, int link)
	   {
		     double path [][] = new double [v][v];
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
	        		  longest_theoritical_disctance=(int) path[i][k];
	           }
	           }
	         //  System.out.println("long    "+longest_theoritical_disctance);
	          
	           
	           return path;
	   }
	
	
	static double [][] shtbfs(double[][] adj, int v, int links)
	{
		double sh [][] = new double [v][v];
		Graph g = new UndirectedSparseGraph();
		int lk=0;
		int edge [][] = new int [2][links];
		for(int i=0;i<v;i++)
		{	
			g.addVertex(i);
		}
		int u=0;
		for(int i=0;i<v;i++)
		{
			for(int j=i+1;j<v;j++)
			{
				if(adj[i][j]==1)
				{
					g.addEdge(u+"",i,j);
					edge[0][u]=i;
					edge[1][u]=j;
					u++;
				}
			}
		}
		DijkstraShortestPath alg = new DijkstraShortestPath(g);
		
		int ch =0;
		for(int i=0;i<v;i++)
		{
			for(int j=v-1;j>i;j--)
			{
				if(sh[i][j]==0)
				{
					ch++;
					List k  = alg.getPath(i, j);
					for(int q=0;q<k.size();q++)
					{
						for(int p=q+1;p<k.size();p++)
						{
							int a1 = edge[0][Integer.parseInt(k.get(q).toString())];
							int a2 = edge[1][Integer.parseInt(k.get(q).toString())];
							
							int b1 = edge[0][Integer.parseInt(k.get(p).toString())];
							int b2 = edge[1][Integer.parseInt(k.get(p).toString())];
							
							sh[a1][b2]=(p-q);
							sh[b2][a1]=(p-q);
							
							sh[a1][b1]=(p-q);
							sh[b1][a1]=(p-q);
							
							sh[a2][b2]=(p-q);
							sh[b2][a2]=(p-q);
							
							
						}
						int a1 = edge[0][Integer.parseInt(k.get(q).toString())];
						int a2 = edge[1][Integer.parseInt(k.get(q).toString())];
						
						sh[i][a1]=q+1;
						sh[a1][i]=q+1;
						
						sh[i][a2]=q+2;
						sh[a2][i]=q+2;
					}
				}
			}
		}
		System.out.println("ch: "+ch);
		/*for(int i=0;i<v;i++)
		{
			for(int j=0;j<v;j++)
			{
				System.out.print(sh[i][j]+" ");
			}System.out.println();
		}*/
		
		return sh;
	}


	public static void preparing_for_graphviz(int adj[][], int v, int links)throws FileNotFoundException
	{
		//adj=matrix_manipulation.randomly_reordering(adj, v);
		File simi = new File("C:/Users/farshad.toosi/Documents/matrices/for graphviz/graphviz.txt");
		PrintWriter similar = new PrintWriter(simi);
		
		similar.println("graph g");
		similar.println("{");
		similar.println("node [shape=point];");
		similar.println("edge [style=bold];");
		similar.println("start =\"rand\";");
		
		for(int i=0;i<links;i++)
		{
			similar.println(""+adj[0][i]+" -- "+adj[1][i]+";");
		}
		similar.println("}");
		similar.close();
		
		System.out.println("finish writing");
	}

	public static int [][] cactus_creating(int v, String tit, double links) throws FileNotFoundException
	{
		//String add="C:/Users/farshad.toosi/workspace/open/GA clusters mats/";
		String add="";
		//File res = new File(add+tit+"cactus.txt");
		//PrintWriter wres = new PrintWriter(res);
		links=0;
		int adj[][] = new int [v][v];
		int nodes=v;
		v=v-(int) Math.sqrt(v)-1;
		System.out.println(nodes+"  "+v);
		int y=0;
		ArrayList <Integer> temp = new ArrayList <Integer>();
		while(y<v)
		{
			int x= (int) (Math.random() * (Math.sqrt(v))+3);
			y+=x;
			
				if(y<v)
				{
					adj[y-x][y-1]=1;
					adj[y-1][y-x]=1;
					adj[y-x][y-x+1]=1;
					adj[y-x+1][y-x]=1;
					temp.add(y-x+1);
					for(int j=y-x+1;j<(y-1);j++)
					{
						adj[j][j+1]=1;
						adj[j+1][j]=1;
						
						adj[j][j-1]=1;
						adj[j-1][j]=1;
					}
					
				}
			//-----------
				int p=(int) (Math.random() * (100));
				if(p<50)
				{
				x= (int) (Math.random() * (Math.sqrt(v))+3);
				y+=x;
				
					if(y<v)
					{
						adj[y-x-1][y-1]=1;
						adj[y-1][y-x-1]=1;
						adj[y-x][y-x]=1;
						adj[y-x][y-x]=1;
						
						for(int j=y-x;j<(y-1);j++)
						{
							adj[j][j+1]=1;
							adj[j+1][j]=1;
							
							adj[j][j-1]=1;
							adj[j-1][j]=1;
						}
					}
				}
		}
		
		
		adj[v][v+3]=1;
		adj[v+3][v]=1;
		
		adj[v][v+1]=1;
		adj[v+1][v]=1;
		
		adj[v+2][v+1]=1;
		adj[v+1][v+2]=1;
		
		adj[v+2][v+3]=1;
		adj[v+3][v+2]=1;
		
		int h=v;
		for(int i=0;i<temp.size();i++)
		{
			adj[h][temp.get(i)]=1;
			adj[temp.get(i)][h]=1;
			if(h==v+3)
				h=v;
			else
				h++;
		}
		
		
		for(int i=0;i<nodes;i++)
		{
			int kl=0;
			for(int j=0;j<nodes;j++)
			{
				if(adj[i][j]==1)
					kl=1;
					
			}
			if(kl==0)
			{
				int x= (int) (Math.random() * v);
				adj[x][i]=1;
				adj[i][x]=1;
			}
		}
		
	
		
		for(int i=0;i<nodes;i++)
		{
			for(int j=i+1;j<nodes;j++)
			{
				if(adj[i][j]==1)
					links++;
			}
		}
		int [][] edge = new int [2][nodes];
		int count =0;
		for(int i=0;i<nodes;i++)
		{
			for(int j=i+1;j<nodes;j++)
			{
				if(adj[i][j]==1)
				{
					edge[0][count]=i;
					edge[1][count]=j;
					System.out.println(edge[0][count]+" "+edge[1][count]);
					count++;
				}
			}
		}
		return edge;
		
	}
	
	


}
