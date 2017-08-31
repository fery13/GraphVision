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

public class Geodesic_GA {

	


		static int [][] ordering;
		static int [][] up ;
		static double longest_theoritical_disctance=0;
		static double largest_btw_cnt=0;
		static double smallest_btw_cnt=0;
		
		
		
		static double [][] GA_2015(int adj[][], int v, int population, int mating_pool,  int gins,int cp, int crs, int dists) throws FileNotFoundException
		{
			
			double [][] pos = new double [2][v];
			
			ordering = new  int [population][v+1];
			up = new int [population][v+1];
			int link=0;
			
			for(int i=0;i<v;i++)
				for(int j=i+1;j<v;j++)
					if(adj[i][j]==1)
						link++;
			
			
			adj = shpth(adj, v);
			
			ArrayList <Integer> temp = new ArrayList <Integer>();
			
			
			for(int i=0;i<population;i++)
			{
				for(int i1=0;i1<v;i1++)
					temp.add(i1);
				
				
				int q=v;
				for(int j=0;j<v;j++)
				{
					int p=(int) ((double) Math.random()*(q-0));
					ordering[i][j]=temp.get(p);
					temp.remove(p);
					q--;
				}

				//***********checking none-repeated orderings ***********
				if(i>0)
				for(int k=0;k<(i-1);k++)
				{
					int t1[] = new int [v];
					int t2[] = new int [v];
					
					for(int j=0;j<v;j++)
					{
						t1[j]=ordering[k][j];
						t2[j]=ordering[i][j];
					}
					
					if(Arrays.equals(t1, t2))
					{
						k=population+10;
						i--;
					}
				}
				//**********End checking none-repeated orderings ********
			}
			
			
			
			
			
			
			
			for(int i=0;i<population;i++)
			{
				int mat [] = new int [v];
				for(int j=0;j<v;j++)
				{
					mat[j]=ordering[i][j];
					
				}
				ordering[i][v]=fitness(adj,mat,v, dists);
				
			}
			
			
			
			//*************Above code has been called only once at beginning***************
		
			
			up=mating_pool( v, population, mating_pool, ordering);
			
			int abs_min[]= new int [v+1];
			
			
			for(int genet=0;genet<gins;genet++)
			{		
				ordering=probablity_operations(up, v, population, crs , cp);
				
				for(int i=0;i<population;i++)
				{
					int mat [] = new int [v];
					
					for(int j=0;j<v;j++)
					{
						mat[j]=ordering[i][j];
					}
					
					//ordering[i][v]=fitness_calculation(mat,v,adj);
					ordering[i][v]=fitness(adj,mat,v, dists);
				}
				
				up=mating_pool( v, population, mating_pool, ordering);
				
				abs_min=minimum_fitness(up, v, population);
			}
		
			
			
			//extracting the ordering with a minimum fittness
			
			abs_min=minimum_fitness(up, v, population);
		    
			double chunk = Math.PI/(double)v;
			
			for(double i=0;i<v;i++)
			{
				pos[0][abs_min[(int)i]]=Math.cos(i*chunk);
				pos[1][abs_min[(int)i]]=Math.sin(i*chunk);
			}
			
			return pos;
			
		}
		
		static int [] minimum_fitness_temp(int ordering[][], int v, int population)
		{
			int max=-1;
			int order []= new int[v+1];
			for(int i=0;i<population;i++)
			{
				if(max<ordering[i][v])
				{
					for(int j=0;j<v+1;j++)
					{
						order[j]=ordering[i][j];
					}
					max=ordering[i][v];
				}
			}
			return order;
		}
		
		
		static int [] minimum_fitness(int ordering[][], int v, int population)
		{
			int order []= new int[v+1];
			
			for(int j=0;j<v+1;j++)
				order[j]=ordering[0][j];
			
			for(int i=0;i<population;i++)
			{
				if(order[v]>ordering[i][v])
				{
					for(int j=0;j<v+1;j++)
					{
						order[j]=ordering[i][j];
					}
					order[v]=ordering[i][v];
				}
			}

			return order;
		}
		
		
		static int maximum_fitness(int ordering[][], int v, int population)
		{
			int max=0;
			
			for(int i=0;i<population;i++)
			{
				if(max<ordering[i][v])
					max=ordering[i][v];
			}
			
			return max;
		}
		
		
		static double average_fitness(int ordering[][], int v, int population)
		{
			double ave=0;
			
			for(int i=0;i<population;i++)
			{
				ave+=ordering[i][v];
			}
			
			return (ave/population);
		}
		
		
		static int [][] cross_over(int order1[], int order2[], int v)
		{
			int [][]cross_over =new int [2][v];
			int []cross1 =new int [v];
			ArrayList <Integer> par1 = new ArrayList <Integer>();
			int []cross2 =new int [v];
			ArrayList <Integer> par2 = new ArrayList <Integer>();
			//***** New Version ******
			int [] all = new int [2*v];
			int j=0;
			
			int p=(int) ((double) Math.random()*(v-2)+1);
			int cross_type=1;
			//cross_type=1 POS
			//cross_type=3 voting crossover
			//cross_type=2 partially
			//cross_type=4 alternating
			if(cross_type==2)
			{
				
				for(int i=0;i<v;i++)
				{
					cross_over[0][i]=-1;
					cross_over[1][i]=-1;
				}
				int cr1=0;
				
				for(int i=0;i<v;i++)
				{
					cross1[cr1]=order1[i];
					cr1++;
				}
				cr1=0;
				for(int i=0;i<v;i++)
				{
					cross2[cr1]=order2[i];
					cr1++;
				}
				
				for(int i=0;i<v;i++)
				{
					if(cross1[i]==cross2[i])
					{
						cross_over[0][i]=cross1[i];
						cross_over[1][i]=cross1[i];
						cross1[i]=-1;
						cross2[i]=-1;
					}
				}
				for(int i=v-1;i>=0;i--)
				{
					if(cross1[i]!=-1)
						for(int j1=0;j1<v;j1++)
						{
							if(cross_over[0][j1]==-1)
							{
								cross_over[0][j1]=cross1[i];
								break;
							}
						}
				}
				
				for(int i=v-1;i>=0;i--)
				{
					if(cross2[i]!=-1)
						for(int j1=0;j1<v;j1++)
						{
							if(cross_over[1][j1]==-1)
							{
								cross_over[1][j1]=cross2[i];
								break;
							}
						}
				}
				
				
				
				
			}
			
			
			
			
			if(cross_type==1)
			{
				int n=(int) ((double) Math.random()*(v-1));
				//n--;
				for(int i=0;i<v;i++)
				{
					cross_over[0][i]=-1;
					cross_over[1][i]=-1;
				}
				
				
				for(int i=0;i<v;i++)
				{
					par1.add(order1[i]);
				}
			
				for(int i=0;i<v;i++)
				{
					par2.add(order2[i]);
				}
				
				for(int i=0;i<n;i++)
				{
					int m=(int) ((double) Math.random()*(v-0));
					if(par1.get(m)==-1)
						i--;
					else
					{
						cross_over[0][m]=par1.get(m);
						par1.set(m, -1);
						for(int j1=0;j1<v;j1++)
						{
							if(par2.get(j1)==cross_over[0][m])
							{
								par2.remove(j1);
								break;
							}
						}
					}
					
				}
				int c=0;
				for(int i=0;i<v;i++)
				{
					if(cross_over[0][i]==-1)
					{
						cross_over[0][i]=par2.get(c);
						c++;
					}
				}
				
				//**second child
				
				ArrayList <Integer> par11 = new ArrayList <Integer>();
				ArrayList <Integer> par21 = new ArrayList <Integer>();
				n=(int) ((double) Math.random()*(v-1));
				//n--;
				
				
				for(int i=0;i<v;i++)
				{
					par11.add(order1[i]);
				}
			
				for(int i=0;i<v;i++)
				{
					par21.add(order2[i]);
				}
				
				for(int i=0;i<n;i++)
				{
					int m=(int) ((double) Math.random()*(v-0));
					if(par21.get(m)==-1)
						i--;
					else
					{
						cross_over[1][m]=par21.get(m);
						par21.set(m, -1);
						for(int j1=0;j1<v;j1++)
						{
							if(par11.get(j1)==cross_over[1][m])
							{
								par11.remove(j1);
								break;
							}
						}
					}
					
				}
				c=0;
				for(int i=0;i<v;i++)
				{
					if(cross_over[1][i]==-1)
					{
						cross_over[1][i]=par11.get(c);
						c++;
					}
				}
				
				//*********
			}
			
			if(cross_type==2)
			{
			//***new version
			int cr1=0;
			for(int i=0;i<p;i++)
			{
				cross1[cr1]=order1[i];
				cr1++;
			}
			for(int i=0;i<v-p;i++)
			{
				cross1[cr1]=order2[i];
				cr1++;
			}
			
			cr1=0;
			for(int i=p;i<v;i++)
			{
				cross2[cr1]=order1[i];
				cr1++;
			}
			for(int i=v-p;i<v;i++)
			{
				cross2[cr1]=order2[i];
				cr1++;
			}
			}
			
			//***end new version
			
			
			if(cross_type==4)
			{
			for(int i=0;i<v;i++)
			{
				all[j]=order1[i];
				j++;
				all[j]=order2[i];
				j++;
			}
			
			
			for(int i=0;i<v;i++)
			{
				cross1[i]=all[i];
			}
			j=0;
			for(int i=v;i<2*v;i++)
			{
				cross2[j]=all[i];
				j++;
			}
			}
			
			if(cross_type==2 || cross_type==4)
			{
			boolean check=true;
			while(check)
			{
				
				check=false;
				int c1=0,c2=0;
				for(int i=0;i<v;i++)
				{
					for(j=i+1;j<v;j++)
					{
						if(cross1[i]==cross1[j])
						{
							c1=i;
							j=v*10;
							i=v*10;
							check=true;
						}
					}
				}
				
				for(int i=0;i<v;i++)
				{
					for(j=i+1;j<v;j++)
					{
						if(cross2[i]==cross2[j])
						{
							c2=i;
							j=v*10;
							i=v*10;
							check=true;
						}
					}
				}
				if(check)
				{
					int temp=cross1[c1];
					cross1[c1]=cross2[c2];
					cross2[c2]=temp;
				}
				
				
			}
			
			for(int i=0;i<v;i++)
			{
				cross_over[0][i]=cross1[i];
				cross_over[1][i]=cross2[i];
			}
			//****End New Version ****
			
			}
			
			
			
			
			return cross_over;
		}
		
		
		static int [][] probablity_operations(int ordering[][], int v, int population, int cross,  int copy)
		{
			int [][] updated_mating = new int [population][v+1];
			
			ArrayList <Integer> temp = new ArrayList <Integer>();
			
			for(int i=0;i<population;i++)
				temp.add(i);
			
			//***********
			int pop=population;
			int mating=0;
			while(mating<population)
			{
				int p=(int) ((double) Math.random()*(100-0));
				
				if(p<copy)
				{
						int q=0;
						for(int i=0;i<v;i++)
						{
							updated_mating[mating][i]=ordering[temp.get(q)][i];
						}
						temp.remove(q);
						pop--;
						mating++;
				}
				if(p>=copy && p<(cross+copy) && mating <(population-1))
				{
					//**************CrossOver******
					int crossing[][] = new int [2][v];
					int order1[] = new int [v];
					int order2[] = new int [v];
					
					int h1=0;
					int h2=0;
					while(h1==h2)
					{
						h1=(int) ((double) Math.random()*(pop-0));
						h2=(int) ((double) Math.random()*(pop-0));
					}
					for(int i=0;i<v;i++)
					{
						order1[i]=ordering[temp.get(h1)][i];
						order2[i]=ordering[temp.get(h2)][i];
					}
					
					crossing=cross_over( order1, order2, v);
					
					for(int i=0;i<v;i++)
					{
						updated_mating[mating][i]=crossing[0][i];
						updated_mating[mating+1][i]=crossing[1][i];
					}
					
					//*********************
					if(h1>h2)
					{
						temp.remove(h1);
						temp.remove(h2);
					}
					else
					{
						temp.remove(h2);
						temp.remove(h1);
					}
					pop--;
					pop--;
					mating+=2;
				}
				if(p>=(cross+copy))
				{
					//************Mutation*****************
					int q=(int) ((double) Math.random()*(pop-0));
					int perm1=0;
					int perm2=0;
					while(perm1==perm2)
					{
						perm1=(int) ((double) Math.random()*(v-1)+1);
						perm2=(int) ((double) Math.random()*(v-0));
					}
					for(int i=0;i<v;i++)
					{
						updated_mating[mating][i]=ordering[temp.get(q)][i];
					}
					
					int sw=updated_mating[temp.get(q)][perm1];
					updated_mating[temp.get(q)][perm1]=updated_mating[temp.get(q)][perm2];
					updated_mating[temp.get(q)][perm2]=sw;
		
					
					
					perm1=0;
					perm2=0;
					while(perm1==perm2)
					{
						perm1=(int) ((double) Math.random()*(v-0));
						perm2=(int) ((double) Math.random()*(v-0));
					}
					
					
					sw=updated_mating[temp.get(q)][perm1];
					updated_mating[temp.get(q)][perm1]=updated_mating[temp.get(q)][perm2];
					updated_mating[temp.get(q)][perm2]=sw;
					
			
					
					//****************************
					temp.remove(q);
					pop--;			
					
					mating++;
				}
				
				
			}
			temp.clear();
			return updated_mating;
		}
		
		
		static int [][] mating_pool(int v, int population, int mating_pool, int ordering[][])
			{	
				int [][] mating = new int [population][v+1];
			
				for(int i=0;i<population;i++)
				{
						int [][] temp = new int [mating_pool][v+1];
							
						ArrayList <Integer> pop = new ArrayList<Integer> ();
						
						for(int j=0;j<population;j++)
							pop.add(j);
						
						int g=population;
						for(int j=0;j<mating_pool;j++)
						{
								
							int p =(int) ((double) Math.random()*(g-0));
							for(int k=0;k<v+1;k++)
								temp[j][k]=ordering[pop.get(p)][k];
							pop.remove(p);
							g--;
						}
						pop.clear();
						
						int [] min = new int [v+1];
						
						for(int q=0;q<v+1;q++)
							min[q]=temp[0][q];
						
						/*if(variables.GA_factor)
						{
							for(int j=0;j<mating_pool;j++)
							{
									if(min[v]<temp[j][v])
										for(int q=0;q<v+1;q++)
											min[q]=temp[j][q];
							}
						}
						else
						{*/
							for(int j=0;j<mating_pool;j++)
							{
									if(min[v]>temp[j][v])
										for(int q=0;q<v+1;q++)
											min[q]=temp[j][q];
							}
						//}
						for(int j=0;j<v+1;j++)
						{
								mating[i][j]=min[j];
						}
				}
				
				return mating;
		}
		
		
		static int fitness (int adj[][], int order[], int v, int dists)
		{
			int fit=0;
			double chunk = Math.PI/(double)v;
			for(int i=0;i<v;i++)
			{
				for(int j=i+1;j<v;j++)
				{	
					if(adj[order[i]][order[j]]<=dists)
						fit+= Math.abs(adj[order[i]][order[j]]-Math.abs(i*chunk-j*chunk));
				}
			}
			return fit;
		}
				
		
		
		public static int numbero_edge_crossing(double [][] poss, int t, int v, int [][] adj)
		{
			int ans=0;
			for(int i=0;i<v;i++)
			{
				for(int j=i+1;j<v;j++)
				{
					if(adj[i][j]==1)
					{
						for(int p=0;p<v;p++)
						{
							for(int q=p+1;q<v;q++)
							{
								if(adj[p][q]==1 && (p!=i && p!=j && q!=i && q!=j))
								{
									double x1=poss[0][i];
									double x2=poss[0][j];
									double x3=poss[0][p];
									double x4=poss[0][q];
									
									double y1=poss[1][i];
									double y2=poss[1][j];
									double y3=poss[1][p];
									double y4=poss[1][q];
									
									double d = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4);
									if (d != 0)
									{
										
										double xi = ((x3-x4)*(x1*y2-y1*x2)-(x1-x2)*(x3*y4-y3*x4))/d;
										double yi = ((y3-y4)*(x1*y2-y1*x2)-(y1-y2)*(x3*y4-y3*x4))/d;
										boolean check1=false;
										boolean check2=false;
										
										if(x1==x2)
										{
											if(y1>y2)
											{
												if(yi>=y2 && y1>=yi)
													check1=true;
											}
											else
												if(y1<y2)
													if(yi<=y2 && y1<=yi)
														check1=true;
										}
										else
											if(x1>x2)
											{
												if(xi>=x2 && x1>=xi)
													check1=true;
											}
											else
												if(x1<x2)
													if(xi<=x2 && x1<=xi)
														check1=true;
										
											
										if(x3==x4)
										{
											if(y3>y4)
											{
												if(yi>=y4 && y3>=yi)
													check2=true;
											}
											else
												if(y3<y4)
													if(yi<=y4 && y3<=yi)
														check2=true;
										}
										else
											if(x3>x4)
											{
												if(xi>=x4 && x3>=xi)
													check2=true;
											}
											else
												if(x3<x4)
													if(xi<=x4 && x3<=xi)
														check2=true;
										
										
										if(check1 && check2)
										{
											ans++;
											
										}
									}				
								}
							}
						}
					}
				}
			}
				
			
			return (ans/2);
		}


		public static int [][] barabasi_albert_graph_generating(int v, int m, String tit) throws FileNotFoundException
		{
			//String add="C:/Users/farshad.toosi/Dropbox/PhD/publications/GD 2015/GA/barabashi 10 - 50//matrics/";
			String add="";
			File res = new File(add+tit+"barabashi.txt");
			PrintWriter wres = new PrintWriter(res);
			int adj [][] = new int [v][v];
			double p=0;
			for(int i=0;i<m;i++)
			{
				p = Math.random();
				if(p<0.95)
				{
					adj[m][i]=1;
					adj[i][m]=1;
				}
			}
			
			
			for(int i=m+1;i<v;i++)
			{
				double k=0;
				double s=0;
				for(int j=0;j<v;j++)
				{
					for(int q=0;q<v;q++)
					{
						if(adj[j][q]==1)
							s++;
					}
				}
				double r= Math.random();
				int n=0;
				boolean ch=true;
				int deg=(int) (Math.random()*2)+1;
				while(ch)
				{
					int b=(int)Math.round(Math.random() * (i-1));
					
					if(adj[b][i]==0)
					{
						for(int j=0;j<v;j++)
						{
							if(adj[j][b]==1)
								k++;
						}
						
						p=k/s;
						
						if(r<p && n<deg)
						{
							adj[b][i]=1;
							adj[i][b]=1;
							n++;
						}
						if(n==deg)
							ch=false;
					}
				}
			}
			
			int link = 0;
			for(int i=0;i<v;i++)
			{
				for(int j=i+1;j<v;j++)
				{
					if(adj[i][j]==1)
						link++;
				}
			}
			int edge [][] = new int [2][link];
			int h=0;
			for(int i=0;i<v;i++)
			{
				for(int j=i+1;j<v;j++)
				{
					if(adj[i][j]==1)
					{
						edge[0][h]=i;
						edge[1][h]=j;
						h++;
					}
				}
			}
			
			return edge;
			
		}
		
		public static int [][] star_shaped_creating(int v, double deg)
		{
		
			int star [][] = new int [v][v];
			
			int sht [][] = new int [v][v];
			double btw [] = new double [v];
			
			double visited [] = new double [v];
			
			int a1 = (int)(Math.random()*v);
			int a2 = (int)(Math.random()*v);
			int popu = 0;
			
			while(a1==a2)
			{
				a1 = (int)(Math.random()*v);
				a2 = (int)(Math.random()*v);
			}
			
			visited[a1]=1;visited[a2]=1; popu++;popu++;
			star[a1][a2]=1;star[a2][a1]=1;
			int lnk =2;
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
				sht = shpth(star, v);
				for(int i=0;i<degree;i++)
				{
					for(int j=0;j<v;j++)
					{
						if(visited[j]==1 && sht [n1][j]<=2 && n1!=j)
						{
							if(star[n1][j]==0)
								lnk++;
							
							star[n1][j]=1;
							star[j][n1]=1;
							//popu++;
							break;
						}
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
					if(i!=j)
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
		
		public static int [][] shpth(int [][] adj, int v)
		   {
			     int path [][] = new int [v][v];
		           int totalNodes = v;
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
		          // System.out.println("long    "+longest_theoritical_disctance);
		          
		           
		           return path;
		   }

		
}
