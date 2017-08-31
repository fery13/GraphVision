
	import static org.lwjgl.opengl.GL11.GL_COLOR_BUFFER_BIT;
	import static org.lwjgl.opengl.GL11.GL_LINES;
	import static org.lwjgl.opengl.GL11.glBegin;
	import static org.lwjgl.opengl.GL11.glClear;
	import static org.lwjgl.opengl.GL11.glClearColor;
	import static org.lwjgl.opengl.GL11.glColor3d;
	import static org.lwjgl.opengl.GL11.glEnd;
	import static org.lwjgl.opengl.GL11.glFlush;
	import static org.lwjgl.opengl.GL11.glOrtho;
	import static org.lwjgl.opengl.GL11.glVertex2f;
	import static org.lwjgl.opengl.GL11.glViewport;

	import java.awt.Color;
	import java.awt.Font;
	import java.awt.FontFormatException;
	import java.awt.image.BufferedImage;
	import java.io.File;
	import java.io.FileNotFoundException;
	import java.io.IOException;
	import java.io.PrintWriter;
	import java.nio.ByteBuffer;
	import java.util.ArrayList;
	import java.util.Scanner;

	import javax.imageio.ImageIO;
	import javax.swing.JOptionPane;

	import org.lwjgl.BufferUtils;
	import org.lwjgl.LWJGLException;
	import org.lwjgl.input.Keyboard;
	import org.lwjgl.input.Mouse;
	import org.lwjgl.opengl.Display;
	import org.lwjgl.opengl.GL11;
	//import org.lwjgl.util.glu.Sphere;
	import org.omg.CORBA.portable.InputStream;

	//import sun.font.TrueTypeFont;

	import Jama.EigenvalueDecomposition;
	import Jama.Matrix;

	//import com.sun.xml.internal.ws.api.ResourceLoader;
	import java.awt.Font;
	public class MDS {
	
		public static int [][] lists;
		public static int [][] lists1;
		public static int [][] lists2;
	//	public static double [][] adjj ;
	
		
		 
		
		
		
	public static double [][][] Multidimentional_forcing(double sim [],int col [], int v, double time, boolean t3d) throws FileNotFoundException
	{
		
		
		graph.colorings();
		
			double pos [][][];
			
			
			double n=v;
			
			
			pos= new double [1][3][(int)n];
			
			int [] visited = new int [v];
			
			int popu=0;
			popu++;
			
			double maxxx=0;
			double min=1001E216;
			double max=0;
			int inde=0;
			int total = (v*v-v)/2;
			
			for(int i=1;i<v;i++)
			{
				int index = i*(v-1)+0-(i*i+i)/2-1;
				maxxx += sim[index];
			}
			
			for(int i=0;i<n;i++)
			{
				double temp=0;
				
				for(int j=0;j<n;j++)
				{
					if(i!=j)
					{
						int index = i*(v-1)+j-(i*i+i)/2-1;
						temp += sim[index];
						
						if(max<sim[index])
							max=sim[index];
						if(min>sim[index])
							min=sim[index];
					}
				}
				
				if(temp<maxxx)
				{
					maxxx=temp;
					inde=i;
				}
			}
			
			double step = (max-min)/Math.sqrt(v);
			
			visited[inde]=1;
			for(int i=0;i<v;i++)
			{
				pos[0][0][i]=Math.random();
				pos[0][1][i]=Math.random();
				if(t3d)
					pos[0][2][i]=Math.random();
			}
			int t=1;
			int s_point=1;
			
			double k = 1.0;
			double visit=0;
			double interations=n/8;
			int level=0;
			double rad=1.0/(double)v;
			double iter=0;
			double circle =0;
			double circle_step = (Math.PI*2)/(double)v;
			
			double cooling = 1.0/interations;
			double b= Math.sqrt(n);
			double pow = b/interations;
			Matrix A= new Matrix(v,v);
			
			Matrix XX= new Matrix(v,1);
			Matrix YY= new Matrix(v,1);
			Matrix ZZ= new Matrix(v,1);
			for(int i=0;i<v;i++)
			{
				A.set(i, i, 0);
				for(int j=i+1;j<v;j++)
				{
					int index = i*(v-1)+j-(i*i+i)/2-1;
					A.set(i, j, sim[index]);
					A.set(j, i, sim[index]);
				}
			}
			
			while( iter<interations)  // timmeeeeeeeeeeemmmmmmmmmmmmeeeeeeeeeeeeeeeemmmmmmmmmmmme
			{
				level++;
				//System.out.println(iter+"  "+popu);
				min += step;
				//******************************
				double x_forces[]=new double [v];
				double y_forces[]=new double [v];
				double z_forces[]=new double [v];
				
				Matrix Coss= new Matrix(v,1);
				Matrix Sinn= new Matrix(v,1);
				Matrix Z_Sinn= new Matrix(v,1);
				
				
				double diss = 0;
				
				
				
				for(int i=0;i<v;i++)
				{
					double sin=0;
					double cos=0;
					double z_hend=0;
					for(int j=0;j<v;j++)
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
						sin += ((pos[0][1][j]-pos[0][1][i])/diss);
						cos += ((pos[0][0][j]-pos[0][0][i])/diss);
								
						if(t3d)
							z_hend += ((pos[0][2][j]-pos[0][2][i])/diss);
					}
					Coss.set(i, 0, cos);
					Sinn.set(i, 0, sin);
					if(t3d)
					{
						Z_Sinn.set(i, 0, z_hend);
					}
								
								//-----------
						
					}
					
					XX=A.times(Coss);
					YY=A.times(Sinn);
					
				    if(t3d)
				    	ZZ=A.times(Z_Sinn);
				    
					for(int i=0;i<v;i++)
					{
					    pos[0][0][i] = XX.get(i, 0);
					    pos[0][1][i] = YY.get(i, 0);
					    if(t3d)
					    	pos[0][2][i] = ZZ.get(i, 0);
					}
			
				//System.out.println((b+1-iter*pow)+"  "+Math.floor(b+1-iter*pow));
				k += 0.05;
				iter++;
				//popu++;
			
			}
			
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
		
			
			//colors
			double r1=0;
			double r2=0;
			double r3=0;
			
			
				
			
			
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


	 
	
	}

	

