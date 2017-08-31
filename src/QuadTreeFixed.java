import java.util.ArrayList;

public class QuadTreeFixed {
	
	private double [] bounds = new double [4]  ;
	private ArrayList<Integer> mat[][] ;
	private int v;
	private int mat_size;
	private double [][][] nodes;
	private int time;
	
	public QuadTreeFixed(double [] bounds ,ArrayList <Integer> mat [][], int v, int mat_size, double [][][] nodes, int time)
	{
		this.bounds = bounds.clone();
		this.v=v;
		this.mat= (ArrayList<Integer> [][]) mat.clone();	
		this.mat_size = mat_size;
		this.nodes= nodes.clone();
		this.time = time;
	}
	public double[] getCentre()
	{
		double xy []= new double [2];
		xy[0]=(this.bounds[1]-this.bounds[0])/2+this.bounds[0];
		xy[1]=(this.bounds[3]-this.bounds[2])/2+this.bounds[2];
		
		return xy;
	}
	
	public int getNoChildCell(int x, int y)
	{
		return mat[x][y].size();
	}
	
	
	public static ArrayList <ArrayList <ArrayList <Integer>>> creating_the_matrix(int v,  ArrayList <ArrayList <ArrayList <Integer>>> mat, int mat_size, double [][][] pos, int time)
	{
		CFF.bound = new double [4];
		
		CFF.bound[0]=pos[time][0][0];
		CFF.bound[1]=pos[time][0][0];
		CFF.bound[2]=pos[time][1][0];
		CFF.bound[3]=pos[time][1][0];
	
		for(int i=0;i<v;i++)
		{
			 if(CFF.bound[0]>pos[time][0][i])
				 CFF.bound[0]=pos[time][0][i];
			 
	    	 if(CFF.bound[1]<pos[time][0][i])
	    		 CFF.bound[1]=pos[time][0][i];
	    	 
	    	 if(CFF.bound[2]>pos[time][1][i])
	    		 CFF.bound[2]=pos[time][1][i];
	    	 
	    	 if(CFF.bound[3]<pos[time][1][i])
	    		 CFF.bound[3]=pos[time][1][i];
		}

		for(int i=0;i<mat_size;i++)
		{
			ArrayList <ArrayList <Integer>> temp = new ArrayList <ArrayList <Integer>> ();
			for(int j=0;j<mat_size;j++)
			{
				ArrayList <Integer> t1 = new ArrayList <Integer> ();
				temp.add(t1);
			}
			mat.add(temp);
		}
		
		for(int i=0;i<v;i++)
		{
			double x,y;
			x=CFF.bound[1]-CFF.bound[0];
			y=CFF.bound[3]-CFF.bound[2];
			x = x/(double) mat_size;
			y = y/(double) mat_size;
			int X,Y;
			X = (int) Math.floor((pos[time][0][i]-CFF.bound[0])/x);
			Y = (int) Math.floor((pos[time][1][i]-CFF.bound[2])/y);
			
			if(pos[time][0][i]==CFF.bound[1] || X>=mat_size)
				X=mat_size-1;
			if(pos[time][1][i]==CFF.bound[3] || Y>=mat_size)
				Y=mat_size-1;
			
			
				//System.out.println("herrrrrrr"+X+"  "+Y+CFF.bound[0]+" "+CFF.bound[1]+" "+CFF.bound[2]+" "+CFF.bound[3]);
			
			mat.get(X).get(Y).add(i);
			CFF.loc[i] = Y*mat_size+X;
		}
		return mat;
		
	}
	

	
	public static double [] sin_cos(int v,   ArrayList <ArrayList <ArrayList <Integer>>> mat, int ms, double [][][] pos, int time, int n,  double[] node)
	{
		double [] SinCos = new double [2];
		int Y= (int) (CFF.loc[n]/ms);
		int X= (int) (CFF.loc[n]%ms);
		
		
		for(int i=0;i<ms;i++)
		{
			for(int j=0;j<ms;j++)
			{
				if(i==X && j==Y)
				{
					int size = mat.get(i).get(j).size();
					for(int q=0;q<size;q++)
					{
						int u = mat.get(i).get(j).get(q);
						if(n!=u)
						{
							double d = Math.sqrt(Math.pow(node[0]-pos[time][0][u], 2)+Math.pow(node[1]-pos[time][1][u], 2));
							if(d>0)
							{
								SinCos[0] += (pos[time][0][u]-node[0])/d;
								SinCos[1] += (pos[time][1][u]-node[1])/d;
							}
						}
					}
				}
				else
				{
					double step_x = CFF.bound[1]-CFF.bound[0];
					double step_y = CFF.bound[3]-CFF.bound[2];
					step_x = step_x/ ms;
					step_y = step_y/ ms;
					
					double tx = CFF.bound[0]+ i*step_x + step_x/2;
					double ty = CFF.bound[2]+ j*step_y + step_y/2;
					
					
					double m = mat.get(i).get(j).size();
					if(m>0)
					{
						double d = Math.sqrt(Math.pow(node[0]-tx, 2)+Math.pow(node[1]-ty, 2));
						if(d>0)
						{
							SinCos[0] += m*((tx-node[0])/d);
							SinCos[1] += m*((ty-node[1])/d);
						}
					}
				}
			}
			
		}
		
		
		
		
		return SinCos;
	}
	
}
