import java.io.FileNotFoundException;
import java.util.ArrayList;

import Jama.*; 



public class BiStress {
	
	
	public static double [][][] Bi_Stress(int adj [][], int v, boolean animation, double links, double time, boolean t3d) throws FileNotFoundException, InterruptedException
	{		
		
		int[][] dij = new int [v][v];
		for(int i=0;i<links;i++)
		{
			dij[adj[1][i]][adj[0][i]]=1;
			dij[adj[0][i]][adj[1][i]]=1;
		}
		
		
		long start =0;
		
		////////////////////////
		///////////////////
		double pos [][][];
		
		if(animation)
			pos = new double [(int) time][3][v];
		else
			pos = new double [1][3][v];
		
		
		
		
		for(int i=0;i<v;i++)
		{
			pos[0][0][i]= (double) Math.random() * (2 - 0.0)-1;
			pos[0][1][i]= (double) Math.random() * (2 - 0.0)-1;
			pos[0][2][i]= (double) Math.random() * (2 - 0.0)-1;
		}
		
		
		
		double a=100;
		double c=100;
		double degree=0;
		double M [][] = new double [v][v];
		double L [][] = new double [v][v];
		Matrix AA= new Matrix(v,v);
		for(int i=0;i<v;i++)
		{
			degree=0;
			for(int j=0;j<v;j++)
			{	
		
				if(i!=j)
					M[i][j]=-1;
				if(i==j)
					M[i][j]=(v-1);
				
				
				if(dij[i][j]==1)
				{
					L[i][j]=-1;
					degree++;
				}
				else
					if(i!=j)
						L[i][j]=0;
			}
			L[i][i]=degree;
		}
	
		for(int i=0;i<v;i++)
		{
			for(int j=0;j<v;j++)
			{
				AA.set(i, j, (M[i][j]+a*L[i][j]));
			}	
		}
		
		
		for(int t=1;t<time;t++)
		{
			
			System.out.println(t);
			
			if(t>time/2)
			{
				//a=100;
				for(int i=0;i<v;i++)
				{
					for(int j=0;j<v;j++)
					{
						AA.set(i, j, (M[i][j]+a*L[i][j]));
					}	
				}
			}
			
			
			double cos [] = new double [v];
			double sin [] = new double [v];
			double sin_z [] = new double [v];
			
			
			
			Matrix x= new Matrix(v,1);
			Matrix y= new Matrix(v,1);
			Matrix z= new Matrix(v,1);
			Matrix bx= new Matrix(v,1);
			Matrix by= new Matrix(v,1);
			Matrix bz= new Matrix(v,1);
		
			
			
			for(int i=0;i<v;i++)
			{
				for(int j=0;j<v;j++)
				{	
					double dist=0;
					if(animation)
					{
						if(t3d)
							dist= Math.sqrt((Math.pow(pos[t-1][0][i]-pos[t-1][0][j], 2)) + ( Math.pow(pos[t-1][1][i]-pos[t-1][1][j], 2))+ (Math.pow(pos[t-1][2][i]-pos[t-1][2][j], 2)));
						else
							dist= Math.sqrt((Math.pow(pos[t-1][0][i]-pos[t-1][0][j], 2)) + ( Math.pow(pos[t-1][1][i]-pos[t-1][1][j], 2)));
					}
					else
					{
						if(t3d)
							dist= Math.sqrt((Math.pow(pos[0][0][i]-pos[0][0][j], 2)) + ( Math.pow(pos[0][1][i]-pos[0][1][j], 2))+ (Math.pow(pos[0][2][i]-pos[0][2][j], 2)));
						else
							dist= Math.sqrt((Math.pow(pos[0][0][i]-pos[0][0][j], 2)) + ( Math.pow(pos[0][1][i]-pos[0][1][j], 2)));
					}
					
					
					
					//***********************Cos and Sin
					if(i!=j)
					{
						if(animation)
						{
							cos[i]+=(pos[t-1][0][i]-pos[t-1][0][j])/dist;
							sin[i]+=(pos[t-1][1][i]-pos[t-1][1][j])/dist;
							if(t3d)
								sin_z[i]+=(pos[t-1][2][i]-pos[t-1][2][j])/dist;
						}
						else
						{
							cos[i]+=(pos[0][0][i]-pos[0][0][j])/dist;
							sin[i]+=(pos[0][1][i]-pos[0][1][j])/dist;
							if(t3d)
								sin_z[i]+=(pos[0][2][i]-pos[0][2][j])/dist;
						}
						bx.set(i, 0, cos[i]);
						by.set(i, 0, sin[i]);
						if(t3d)
							bz.set(i, 0, sin_z[i]);
					}
					
				}
				
			}
			
			
		
			x=AA.solve(bx);
			y=AA.solve(by);
			if(t3d)
				z=AA.solve(bz);
			
			
			for(int i=0;i<v;i++)
			{
				
				if(animation)
				{
					pos[t][0][i]=x.get(i, 0);
					pos[t][1][i]=y.get(i, 0);
					if(t3d)
						pos[t][2][i]=z.get(i, 0);
				}
				else
				{
					pos[0][0][i]=x.get(i, 0);
					pos[0][1][i]=y.get(i, 0);
					if(t3d)
						pos[0][2][i]=z.get(i, 0);
				}
			}
		}
		
		long end=0;
		System.out.println("Time:  "+(end-start));
		pos=main_fitting_in_center(v, pos, (int) time, !animation,t3d );
		pos=limit_in_screen(v, pos, (int) time, !animation, t3d);
		pos=covert_to_glortho(v, pos, (int) time, !animation, t3d);
		
		/*for(int i=0;i<v;i++)
		{
			System.out.println(pos[0][0][i]+"  "+pos[0][1][i]+"  ");
		}*/
		
		check_circular(v,  pos, 0);
			return pos;
		
		
		
	}
	
	//Qtree
	public static double [][][] Bi_Stress1(int adj [][], int v, boolean animation, double links, double time, boolean t3d) throws FileNotFoundException, InterruptedException
	{		
		
		int[][] dij = new int [v][v];
		for(int i=0;i<links;i++)
		{
			dij[adj[1][i]][adj[0][i]]=1;
			dij[adj[0][i]][adj[1][i]]=1;
		}
		
		
		long start =0;
		
		////////////////////////
		///////////////////
		double pos [][][];
		
		if(animation)
			pos = new double [(int) time][3][v];
		else
			pos = new double [1][3][v];
		
		double [] bound = new double [4];
		bound[0]=10000;bound[1]=-10000;bound[2]=10000;bound[3]=-10000;
		ArrayList<double []> temp = new ArrayList<double []>();
		
		
		for(int i=0;i<v;i++)
		{
			pos[0][0][i]= (double) Math.random() * (2 - 0.0)-1;
			pos[0][1][i]= (double) Math.random() * (2 - 0.0)-1;
			pos[0][2][i]= (double) Math.random() * (2 - 0.0)-1;
			
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
		
		System.out.println(bound[0]+"  "+bound[1]+"  "+bound[2]+"  "+bound[3]);
		
		double a=1;
		double c=100;
		double degree=0;
		double M [][] = new double [v][v];
		double L [][] = new double [v][v];
		Matrix AA= new Matrix(v,v);
		
		for(int i=0;i<v;i++)
		{
			degree=0;
			for(int j=0;j<v;j++)
			{	
		
				if(i!=j)
					M[i][j]=-1;
				if(i==j)
					M[i][j]=(v-1);
				
				
				if(dij[i][j]==1)
				{
					L[i][j]=-1;
					degree++;
				}
				else
					if(i!=j)
						L[i][j]=0;
			}
			L[i][i]=degree;
		}
	
		for(int i=0;i<v;i++)
		{
			for(int j=0;j<v;j++)
			{
				AA.set(i, j, (M[i][j]+a*L[i][j]));
			}	
		}
		
		
		for(int t=1;t<time;t++)
		{
			
			System.out.println(t);
			
			//if(t>time/2)
				//a=1;
			Matrix x= new Matrix(v,1);
			Matrix y= new Matrix(v,1);
			Matrix bx= new Matrix(v,1);
			Matrix by= new Matrix(v,1);
		
			
			ArrayList<quadtree> chs = new ArrayList<quadtree>();
			quadtree main_root = new quadtree(bound,null,temp,false,3,0,chs);
			    
			quadtree.creat_childern(bound,main_root ,temp, 0 );
				
			    
			for(int i=0;i<v;i++)
			{	
				quadtree.xx=0;
				quadtree.yy=0;
				if(animation)
					quadtree.traversing_tree(main_root, pos[t-1][0][i],pos[t-1][1][i],i);
				else
					quadtree.traversing_tree(main_root, pos[0][0][i],pos[0][1][i],i);
				
						    
				double tsin=quadtree.yy;
				double tcos=quadtree.xx;
				
				bx.set(i, 0, tcos);
				by.set(i, 0, tsin);
			}
		
		
			x=AA.solve(bx);
			y=AA.solve(by);
		
			
			
			double maxx=0;
			for(int i=0;i<v;i++)
			{
				
				if(animation)
				{
					pos[t][0][i]=x.get(i, 0);
					pos[t][1][i]=y.get(i, 0);
					
					if(maxx<Math.abs(pos[t][0][i]))
			    		 maxx=Math.abs(pos[t][0][i]);
					 					
					if(maxx<Math.abs(pos[t][1][i]))
			    		 maxx=Math.abs(pos[t][1][i]);
				}
				else
				{
					pos[0][0][i]=x.get(i, 0);
					pos[0][1][i]=y.get(i, 0);
				}
			}
			
			temp = new ArrayList<double []>();
			temp.clear();
			
			double disp = 10.1;
			for(int i=0;i<v;i++)
			{
				for(int j=i+1;j<v;j++)
				{
					if(animation)
					{
						if(pos[t][0][i]==pos[t][0][j] && pos[t][1][i]==pos[t][1][j])
							pos[t][0][i] +=disp;
					}
					else
					{
						if(pos[0][0][i]==pos[0][0][j] && pos[0][1][i]==pos[0][1][j])
						{
							pos[0][0][i] +=disp;
							System.out.println("dsafdsaf ");
						}
					}
				}
			}
			
			if(animation)
			{
				bound[0]=bound[1]=(pos[t][0][0]);
				bound[2]=bound[3]=(pos[t][1][0]);
			}
			else
			{
				bound[0]=bound[1]=(pos[0][0][0]);
				bound[2]=bound[3]=(pos[0][1][0]);
			}
			for(int i=0;i<v;i++)
			{
				
				if(animation)
				{
					
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
					
					if(bound[0]>pos[0][0][i])
			    		 bound[0]=pos[0][0][i];
					 
			    	if(bound[1]<pos[0][0][i])
			    		 bound[1]=pos[0][0][i];
			    	 
			    	if(bound[2]>pos[0][1][i])
			    		 bound[2]=pos[0][1][i];
			    	 
			    	if(bound[3]<pos[0][1][i])
			    		 bound[3]=pos[0][1][i];
			    	double [] cord = new double [3];
			    	cord[0]=pos[0][0][i];
					cord[1]=pos[0][1][i];
					cord[2]=i;
					temp.add(cord);
					
				}
			}
		}
		
		long end=0;
		System.out.println("Time:  "+(end-start));
		pos=main_fitting_in_center(v, pos, (int) time, !animation,t3d );
		pos=limit_in_screen(v, pos, (int) time, !animation, t3d);
		pos=covert_to_glortho(v, pos, (int) time, !animation, t3d);
		
		/*for(int i=0;i<v;i++)
		{
			System.out.println(pos[0][0][i]+"  "+pos[0][1][i]+"  ");
		}*/
		
		
			return pos;
		
		
		
	}
	
	
	public static double [][][] main_fitting_in_center_mean(int v, double [][][] pos, int t)
	{
		double max_x,max_y,min_x,min_y, ave_x,ave_y,max, max_z, min_z, ave_z;
		double posi [][][] = new double [1][2][v]; 
			max_x=pos[0][0][0];
			max_y=pos[0][1][0];
			min_x=pos[0][0][0];
			min_y=pos[0][1][0];
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
						
						double temp=Math.sqrt(Math.pow(pos[0][0][i], 2)+Math.pow(pos[0][1][i], 2));
						
						if(temp>max)
							max=temp;
					}
					
					
					for(int i=0;i<v;i++)
					{
						posi[0][0][i]=pos[0][0][i]/max;
						posi[0][1][i]=pos[0][1][i]/max;	
					}
		
		return posi;
		
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
		

}
