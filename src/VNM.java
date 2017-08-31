	import java.io.File;
	import java.io.FileNotFoundException;
	import java.io.IOException;
	import java.io.PrintWriter;
	import java.util.ArrayList;
	import java.util.Scanner;

	import javax.swing.JOptionPane;

public class VNM {
	
		
		public static int [][] matx;
		public static int [][] lists;
		public static ArrayList<int []> ll = new ArrayList<int []>();
		public static ArrayList<Integer> degree = new ArrayList<Integer>();
		public static int [][] lists1;
		public static int [][] lists2;
		public static int vertices []; 
		private static double force_x[];
		private static double force_y[];
		private static boolean ev = true;
		public static int level=1;
		
		public static boolean [][] mat;
		public static int [][] mat1;
		private static ArrayList<String> Alists = new ArrayList<String>();
		private static int com=0;
		private static double hypotenuse=0;
		private static double energy0=10E100;
		private static double energy1=0;
		private static double k=0;
		private static double popu=0;
		private static boolean [][] sh2;
		private static double c = 0.5;
		private static double progress=0; 
		private static double int_pr=0; 
		
		private static double x_cent=0;
		private static double y_cent=0;
		
		public static int[] degree_met;
		public static double force=0;
		public static double tol=0;
		public static double deg_met=0;
		public static long main_run_time=0;
		public static ArrayList<int []> seen_links = new ArrayList<int []> ();
		public static double ch1 = 2;
		public static double ch2 = 10;
		private static double pos [][][];
	
			
		static boolean [][] matrix_change_to_distance_cells_boolean(boolean adj[][], int v)
		{
			boolean temp [][]=new boolean [v][v];
			
			for(int i=0;i<v;i++)
			{
				for(int j=0;j<v;j++)
				{
					if(adj[i][j])
					{
						for(int z=j;z<v;z++)
						{
							if(z!=j && adj[i][z] && !adj[j][z] && !adj[z][j])
							{
								temp[j][z]=true;
								temp[z][j]=true;
							}	
						}
					}
				}
			}
			
		
			
			return temp;
		}
		
		
		public static double [][][] Multilevel_centric   (int adj [][], int v, boolean animation, double links, double time, boolean t3d) throws FileNotFoundException
		{
		
			long start_time=0;
			start_time=System.nanoTime();
			lists =adj.clone();
			degree_met = new int [v];
			double node_degree [] = new double [v]; 
			for(int j=0;j< links;j++)
			{
				int a1= lists[0][j];
				int a2= lists[1][j];
				node_degree[a1]++;
				node_degree[a2]++;
			}
			if(animation)
				pos = new double [(int) time][3][v];
			else
				pos = new double [1][3][v];
			
			long end_time=0;
			boolean bigraph=true;
			popu=0;
			
			int lin = (int) (Math.random()* links);
			lin=0;
			System.out.println("lin: "+lin);
			vertices= new int [v];
			
			double vv= v;
			double max_x0=-1000;
			double min_x0=1000;
			double max_y0=-1000;
			double min_y0=1000;
			double dir= (double) Math.PI/vv;
			int lev=0;
			
			level = 1;
			double disp=1.0/vv;
			tol=Math.sqrt( v);
			tol=0.1;
			
			boolean check=false;
			
			//*************************************
			
			int a1=0, a2=0;
			int ini [] = new int [2];
							
			a1= lists[0][lin];
			a2= lists[1][lin];
			
			/*a1= lists.get(lin)[0];
			a2= lists.get(lin)[1];*/
			
			/*a1=matx[0][lin];
			a2=matx[1][lin];*/


			//System.out.println( node_degree[a1]+"  "+ node_degree[a2]+"  ");
			ini[0]=a1;
			ini[1]=a2;
			seen_links.add(ini);
			System.out.println(lin+"  "+a1+"   "+a2);
					
			 pos[0][0][a1]= (double)Math.cos((double)0)* 0.025;
			 pos[0][1][a1]= (double)Math.sin((double)0)* 0.025;
			popu++;
			vertices[a1]=1;
			degree_met[a1]--;
			
			 pos[0][0][a2]= (double)Math.cos((double)1.0)* 0.025;
			 pos[0][1][a2]= (double)Math.sin((double)1.0)* 0.025;
			popu++;
			vertices[a2]=1;
			degree_met[a2]--;
			
			//k=0.5/popu;
			//k=di/popu;
			
			for(int i=0;i<2;i++)
			{
				if( pos[0][0][ini[i]]>max_x0)
					max_x0= pos[0][0][ini[i]];
				if( pos[0][0][ini[i]]<min_x0)
					min_x0= pos[0][0][ini[i]];
				
				if( pos[0][1][ini[i]]>max_y0)
					max_y0= pos[0][1][ini[i]];
				if( pos[0][1][ini[i]]<min_y0)
					min_y0= pos[0][1][ini[i]];
			}
			//*************************************
			//hypotenuse = Math.sqrt(Math.pow(max_x0-min_x0, 2)+Math.pow(max_y0-min_y0, 2));
			hypotenuse = (Math.abs(max_x0-min_x0)+Math.abs(max_y0-min_y0));
			
			
			//*************************
			//*************************
			System.out.println("Starts algorithm parallel: ... ");
			
			
			System.out.println( v+ " Vertices and  "+  links+" links  exist!");
			
			//for(t=1;t< time;t++)
			//double di=1.0/Math.sqrt(popu);
			//double	ti = 0.1/Math.sqrt(vv); // DO NOT TOUCH IT :( YOU UNDERSTAND ??????????
			
			double di=0.5;
			double ti=0.0001;
			
			//di=1.75;
			
			
			double ti_lev=ti;
			
			double ti_rate=0.9;
			
			double k_rate=0.99;
			
			double roun = Math.round(Math.sqrt(vv));
			roun = 2;
			double last=0;
			//k=disp;
			// pos=initial_placement ((int) vv,  lists);
			int q=0;
			int p=0;
			//ti = ti*0.9; // mesh
			//ti = ti*0.99; // sipenski
			//ti= ti*0.999; // small
			//-----------energy checking
			int t=0;
			int y=0;
			while( y<time-1)
			{	
				y++;
				if(animation)
				{
					t=y;
					p=y-1;
				}
				else
				{
					t=0;
					p=1;
				}
				
				
				
				force_x = new double [v];
				force_y = new double [v];
				

				//if(y%2==0)
				if(animation)
					forc(p,(int) links,0, v); /////////////////
				else
					forc(t,(int) links,0, v); /////////////////
				
				/*else
					forc_less(p,(int) links,0, v);*/
				
				
				boolean maxmin=true;
				for(int i=0;i< v;i++)
				{	
					if(vertices[i]!=0)
					{	
						/*if(y%2==0)
						{
							pos[t][0][i] = pos[p-1][0][i] + (force_x[i]*ti);
							pos[t][1][i] = pos[p-1][1][i] + (force_y[i]*ti);
						}
						else
						{*/
							pos[t][0][i] = (force_x[i]/node_degree[i]);
							pos[t][1][i] = (force_y[i]/node_degree[i]);
						//}
						
					}
				}


				hypotenuse = (Math.abs(max_x0-min_x0)+Math.abs(max_y0-min_y0));
				
					
					
					//**********************
				
				if(/*Math.abs(energy1-energy0) < tol*disp*/  true || (y % roun)==0)
				{
					
					di =di*k_rate; // mesh
					System.out.println(y+" "+popu);	
					for(int i=0;i< v;i++)
					{	
						if(vertices[i]==level)
						{		
							
							double s=0;
							//double chunk = (Math.PI*2.0)/(double) node_degree[i];
							double chunk = (Math.PI*2.0)/(double)degree_met[i];
							for(int j=0;j< links;j++)
							{
								a1= lists[0][j];
								a2= lists[1][j];
							//	disp =0.1;
								if((a1==i && vertices[a2]==0))
								{
								
									vertices[a2]=level+1;
									
									 pos[t][0][a2]= pos[t][0][i]+ (Math.cos(chunk*s)*disp);
									 pos[t][1][a2]= pos[t][1][i]+ (Math.sin(chunk*s)*disp);								
									s +=1;
									popu++;
									degree_met[a2]--;
								}
								else
									if(a2==i && vertices[a1]==0)
									{
										
										vertices[a1]=level+1;
										
										 pos[t][0][a1]= pos[t][0][i]+ (Math.cos(chunk*s)*disp);
										 pos[t][1][a1]= pos[t][1][i]+ (Math.sin(chunk*s)*disp);
										s +=1;
										popu++;
										degree_met[a1]--;
									}
							}
						} // if vertice -> level
					}  // for loop
					level++;
					//k=di/Math.sqrt(popu);
					//System.out.println(popu+" "+k+" "+ti+" "+hypotenuse);
					ti = ti_lev;
					//ti_lev = ti_lev - (20.0/vv);
					lev++;
					//tol=tol*1.02;
					int_pr=0;
				}
				int_pr++;
				energy0=energy1;
				//tol = tol*0.9;
				if(popu==vv)
					last++;
				if(check)
					break;
			}
			
			//out.close();
				
			main_run_time += (end_time-start_time);
			
			//3d --------
			
			
			//3d ----------
			
			
			int j=0;
			System.out.println("Here");
			pos=main_fitting_in_center(v, pos, (int) time, !animation,t3d );
			pos=limit_in_screen(v, pos, (int) time, !animation, t3d);
			pos=covert_to_glortho(v, pos, (int) time, !animation, t3d);
			
			return pos;
			
		}
		
		public static void forc1(int t, int e, int s, int v)
		{
			
			ch1=Math.sqrt((double)  v);
			ch2=ch1*5;
			
			if(hypotenuse<0.1)
			{
				ch1=ch1*(0.1);
				ch2=ch2*(0.1);
			}
			
			
			for(int j=s;j<e;j++)
			{
				
				int v1 = lists[0][j];
				int v2 = lists[1][j];
			
				if(vertices[v1]!=0 && vertices[v2]!=0   )
				{
					double rep=0;
					double att=0;
					force=0;
					com++;
					//double dis=Math.abs( pos[t][0][v1]- pos[t][0][v2])+Math.abs( pos[t][1][v1]- pos[t][1][v2]);
					
					double dis=Math.sqrt(Math.pow(pos[t][0][v1]- pos[t][0][v2],2)+Math.pow( pos[t][1][v1]- pos[t][1][v2],2));
					
					
					if(dis<0.00000000002)
						dis=0.00000000002;
					double cos=(( pos[t][0][v2]- pos[t][0][v1])/dis);
					double sin=(( pos[t][1][v2]- pos[t][1][v1])/dis);
					
							
					rep = ((k*k)/-dis);
					att = dis*dis/k;
				
					/*if(popu== v)
						System.out.println(dis+"  "+rep+" "+att);*/
					
					
					
					force= rep+att;
					
					if(force>0 && force>hypotenuse/ch1)
						force=hypotenuse/ch1;
					if(force<0 && Math.abs(force)>hypotenuse/ch2)
						force=-hypotenuse/ch2;
					
					energy1 +=force*2;
									
					double x_t=cos*force;
					double y_t=sin*force;
					
					force_x[v1]+=x_t;
					force_y[v1]+=y_t;

				
					force_x[v2]-=x_t;
					force_y[v2]-=y_t;
					
					/*force_x[v1]+=cos*force;
					force_y[v1]+=sin*force;

				
					force_x[v2]-=cos*force;
					force_y[v2]-=sin*force;*/
					
									
					//*********replusion for nighbures
				}
			}
			
		/*	for(int i=0;i< v;i++)
			{
				//for(int j=i+1;j< v;j++)
				for(int b=0;b<progress;b++)
				{
					int j= (int) Math.random()* v;
					if(vertices[i]!=0 && vertices[j]!=0 && !mat[i][j] && i!=j)
					{
						
						double dis=Math.abs( pos[0][i]- pos[0][j])+Math.abs( pos[1][i]- pos[1][j]);
						double cos=(( pos[0][j]- pos[0][i])/dis);
						double sin=(( pos[1][j]- pos[1][i])/dis);
						
						double force = (k*k)/(-dis);
						
						if(Math.abs(force)>hypotenuse/1)
							force=-(double) hypotenuse/1;
												
						force_x[i]+=cos*force;
						force_y[i]+=sin*force;
										
						force_x[j]-=cos*force;
						force_y[j]-=sin*force;
												
												
						energy1 += force*2;
					}
				}
			}*/
			
			
		}
		
		public static void forc(int t, int e, int s, int v)
		{
			
			ch1=Math.sqrt((double)  v);
			ch2=ch1*5;
			
			if(hypotenuse<0.1)
			{
				ch1=ch1*(0.1);
				ch2=ch2*(0.1);
			}
			
			
			for(int j=s;j<e;j++)
			{
				
				int v1 = lists[0][j];
				int v2 = lists[1][j];
			
				if(vertices[v1]!=0 && vertices[v2]!=0   )
				{
					double rep=0;
					double att=0;
					force=0;
					com++;
					//double dis=Math.abs( pos[t][0][v1]- pos[t][0][v2])+Math.abs( pos[t][1][v1]- pos[t][1][v2]);
					
					double dis=Math.sqrt(Math.pow(pos[t][0][v1]- pos[t][0][v2],2)+Math.pow( pos[t][1][v1]- pos[t][1][v2],2));
					
					
					if(dis<0.00000000002)
						dis=0.00000000002;
					double cos=(( pos[t][0][v2]- pos[t][0][v1])/dis);
					double sin=(( pos[t][1][v2]- pos[t][1][v1])/dis);
					
					k=10/(double)v;
					k= 2.1;		
					//rep = ((k*k)/-dis);
					att = dis/k;
				
					
					force= rep+att;
					
					/*if(force>0 && force>hypotenuse/ch1)
						force=hypotenuse/ch1;
					if(force<0 && Math.abs(force)>hypotenuse/ch2)
						force=-hypotenuse/ch2;
					
					energy1 +=force*2;*/
					double kh = 1/(double)v;	
					double x_t=cos*force*kh;
					double y_t=sin*force*kh;
					
					force_x[v1]+= x_t;
					force_y[v1]+= y_t;

				
					force_x[v2]+= -x_t;
					force_y[v2]+= -y_t;
			
									
					//*********replusion for nighbures
					
					force_x[v1] += pos[t][0][v2];
					force_y[v1] += pos[t][1][v2];

				
					force_x[v2] += pos[t][0][v1];
					force_y[v2] += pos[t][1][v1];
			
					
				}
			}
			
		
			
			
		}
		
		public static void forc_less(int t, int e, int s, int v)
		{
			
		
			
			
			for(int j=s;j<e;j++)
			{
				
				int v1 = lists[0][j];
				int v2 = lists[1][j];
			
				if(vertices[v1]!=0 && vertices[v2]!=0   )
				{
					force_x[v1] += pos[t][0][v2];
					force_y[v1] += pos[t][1][v2];

				
					force_x[v2] += pos[t][0][v1];
					force_y[v2] += pos[t][1][v1];
			
									
					//*********replusion for nighbures
				}
			}
			
		
			
			
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
			double scale=2.0;
			
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


}
