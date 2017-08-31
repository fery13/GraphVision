import java.io.FileNotFoundException;

import javax.swing.Icon;
import javax.swing.JOptionPane;
import javax.swing.UIManager;

import Jama.Matrix;

public class CMD {

	
	
	public static double [][][] Multidimentional_forcing(double sim [],int col [], int v, double time, boolean t3d) throws FileNotFoundException
	{
			
		Object[] possibilities = {1, 2, 3};
	    
	    Integer iii = (Integer) JOptionPane.showInputDialog(null,
	            "1: Normal - 2:Medium - 3:Strong", "Groups Emphasising",
	            JOptionPane.PLAIN_MESSAGE, null, possibilities, "Numbers");

		//////////////////////
		graph.colorings();
		////////////////////////////////////
		/////////////////////////////////////////////////
		
		
		
		long ss = System.nanoTime();
		
		double n=v;
		double [][][] pos; 
		
		pos= new double [1][3][(int)n];
		int [] visited = new int [v];
		
		int popu=0;
		popu++;
			
		
		
		double minn=0;
		double maxx=0;
		
		double min=1001E216;
		double max=0;
		int inde=0;
		
		for(int i=0;i<n-1;i++)
		{
			minn += sim[i];
		}
		double [][] dists = new double [2][v];
		int ord [] = new int [v];
		for(int i=0;i<n;i++)
		{
			double temp=0;
			for(int j=0;j<n;j++)
			{
				if(i!=j)
				{
					int index = 0;
					
					if(i<j)
						index = i*(v-1)+j-(i*i+i)/2-1;
					else
						index = j*(v-1)+i-(j*j+j)/2-1;
					
					temp += sim[index];
					dists[0][i] +=sim[index];
					if(max<sim[index])
						max=sim[index];
					if(min>sim[index])
						min=sim[index];
				}
			}
			
			if(temp<minn)
			{
				minn=temp;
				inde=i;
			}
			if(temp>maxx)
				maxx=temp;
		}
		ord[0]=inde;
		dists[1][inde]=1;
		int cf=1;
		int tmin=0;
		
		while(cf<v)
		{
			double tm = maxx*10;
			for(int i=0;i<n;i++)
			{
				if(dists[1][i]==0)
				{	
					if(dists[0][i]<tm)
					{
						tmin=i;
						tm=dists[0][i];
					}
				}
			}
			ord[cf]=tmin;
			dists[1][tmin]=1;
			cf++;
		}
		
		
		visited[ord[0]]=1;
		pos[0][0][ord[0]]=0;
		pos[0][1][ord[0]]=0;
		
		if(t3d)
			pos[0][2][ord[0]]=0;
			
			
		
		double interations=n/2;
		double rad=1.0/((double)v*(double)v);
		double iter=0;
		
		double b= Math.log(n);
		//double pow = b/interations;
		double pp=popu;
		
		double step = (maxx-minn)/b;
		int fg=1;
		int gt=6;
		int rgt=gt;
		if(v<1000)
		{
			gt=2;
			rgt=2;
			interations=n/1.1;
		}
		else
		if(v<2000)
		{
			gt=3;
			rgt=3;
			interations=n/1.5;
		}
	//	interations = n*2;
		while( iter<interations || pp<v+b*2 )  // timmeeeeeeeeeeemmmmmmmmmmmmeeeeeeeeeeeeeeeemmmmmmmmmmmme
		{
			
			 if(popu==v)
				 pp++;
			 else
				 pp=popu;
			//System.out.println(iter+"  "+pp+" "+popu);
			minn += step;
			//******************************
			double x_forces[]=new double [v];
			double y_forces[]=new double [v];
			double z_forces[]=new double [v];
			
			double diss;
			
			if(popu<v)
			{
				for(int i=fg;i<gt;i++)
				{
					double s = Math.random()*0.025;
					pos[0][0][ord[i]]=pos[0][0][ord[i-1]]+Math.cos(s)*rad;
					pos[0][1][ord[i]]=pos[0][1][ord[i-1]]+Math.sin(s)*rad;
							
					if(t3d)
						pos[0][2][ord[i]]=pos[0][2][ord[i-1]]+(s/6.29)*rad;
							
					visited[ord[i]]=1;
					popu++;
				}
				fg=gt;
				gt +=rgt;
				if(gt>v)
					gt=v;
			}
			
			
				for(int i=0;i<v;i++)
				{
					for(int j=i+1;j<v;j++)
					{
						if( visited[i]!=0 && visited[j]!=0)
						{
							//diss[i][j]= Math.sqrt((Math.pow(pos[0][i]-pos[0][j], 2)) + ( Math.pow(pos[1][i]-pos[1][j], 2)));
							if(t3d)
								diss= ((Math.abs(pos[0][0][i]-pos[0][0][j])) + ( Math.abs(pos[0][1][i]-pos[0][1][j]))+( Math.abs(pos[0][2][i]-pos[0][2][j])));
							else
								diss= ((Math.abs(pos[0][0][i]-pos[0][0][j])) + ( Math.abs(pos[0][1][i]-pos[0][1][j])));
							
							
							if(diss<0.000000000000002)
							{
								diss=0.000000000000002;
								diss=0.000000000000002;
							}
							//-----------
							double sin=((pos[0][1][j]-pos[0][1][i])/diss);
							double cos=((pos[0][0][j]-pos[0][0][i])/diss);
							
							double z_hend =0;
							if(t3d)
								z_hend = ((pos[0][2][j]-pos[0][2][i])/diss);
							
							
							int index = i*(v-1)+j-(i*i+i)/2-1;
							double x=0;
							double y=0;
							double zz=0;
							
							double force = 0;
							force = (sim[index]);
							
							if(iii==3)
								force=force*force*force*force;
							else
								if(iii==2)
									force=force*force*force;
								else
								if(iii==1 && pp<v+10)
									force=force*force;
								else
									force=force;
							
									
							x =cos*(force);
							y =sin*(force);
							
							if(t3d)
								zz =z_hend*force;	
							
							
							x_forces[i]+=-x;
							y_forces[i]+=-y;
							
							x_forces[j]+=x;
							y_forces[j]+=y;
							
							if(t3d)
							{
								z_forces[i] += -zz;
								z_forces[j] += zz;
							}
							
							//-----------
						}
					}
				}
				for(int i=0;i<v;i++)
				{
					if(visited[i]!=0)
					{
						
					    pos[0][0][i] =(x_forces[i]/popu);
					    pos[0][1][i] =(y_forces[i]/popu);
					    if(t3d)
					    	pos[0][2][i] =(z_forces[i]/popu);
					     			   
					}
				}
			
			iter++;
			//popu++;
		
		}
		long ee = System.nanoTime();
		
		System.out.println("time: " +(ee-ss));
		double stress=0;
		double sorat=0;
		double makh =0;
		for(int i=0;i<v;i++)
		{
			for(int j=i+1;j<v;j++)
			{
				double dis = Math.sqrt(Math.pow(pos[0][0][i]-pos[0][0][j], 2)+Math.pow(pos[0][1][i]-pos[0][1][j], 2));
				int index = i*(v-1)+j-(i*i+i)/2-1;
				sorat += (Math.pow(sim[index]-dis,2));
				makh += sim[index]*sim[index];
			}
		}
		stress = Math.sqrt(sorat/makh);
		System.out.println("stress:   "+stress);
			
			pos=main_fitting_in_center(v, pos, (int) time,true,t3d );
			pos=limit_in_screen(v, pos, (int) time, true,t3d);
			pos=covert_to_glortho(v, pos, (int) time,  true,t3d);
			
			
			
				
			
			
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
		double lm=4.9;
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
		double scale=2.5;
		
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
