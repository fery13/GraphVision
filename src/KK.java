import java.io.FileNotFoundException;


public class KK {
	
	
	public static double longest_theoritical_disctance;
	
	public static double [][][] kamada_kawai(int adjlink [][], int v, boolean animation, double links, double time) throws FileNotFoundException, InterruptedException
	{
		
		
		
		double pos [][][];
		
		if(animation)
			pos = new double [(int) time][3][v];
		else
			pos = new double [1][3][v];
		
		int[][] dij = new int [v][v];
		for(int i=0;i<links;i++)
		{
			dij[adjlink[1][i]][adjlink[0][i]]=1;
			dij[adjlink[0][i]][adjlink[1][i]]=1;
		}
		
		dij=  matrix_change_to_distance_cells(dij, v);
	
	
		for(int i=0;i<v;i++)
		{
			pos[0][0][i]= (double) Math.random() * (2 - 0.0)-1;
			pos[0][1][i]= (double) Math.random() * (2 - 0.0)-1;
			
		}
		
		
		
		double vv=v;
		double L = 2.0/(double)longest_theoritical_disctance;
		double K=2.0;
		
		double landa = 0.0001;
		double max_delta=0;
		int max_ind = 0;
		double cof = 1.0;
		int t=1;
		int q=1;
		///--------------	
		double deltax =0;
		double deltay =0;
		double delta=0;
		for(int i=0;i<vv;i++)
		{
			deltax =0;
			deltay =0;
			delta=0;
			for(int j=0;j<vv;j++)
			{
				if(i!=j)
				{
					double xpp=(pos[t-1][0][i]-pos[t-1][0][j]);
					double ypp=(pos[t-1][1][i]-pos[t-1][1][j]);
					deltax += (K/Math.pow((double)dij[i][j],2.0))*(xpp-((L*(double)dij[i][j]*xpp)/(Math.sqrt(Math.pow(xpp, 2)+Math.pow(ypp, 2)))));
					deltay += (K/Math.pow((double)dij[i][j],2.0))*(ypp-((L*(double)dij[i][j]*ypp)/(Math.sqrt(Math.pow(xpp, 2)+Math.pow(ypp, 2)))));
				}
			}
			delta = Math.sqrt(Math.pow(deltax, 2)+Math.pow(deltay, 2));
			if(delta > max_delta)
			{
				max_delta = delta;
				max_ind=i;
			}
		}
			
		//---------------	
			int y=max_ind;
			//time = v*50;
			while(max_delta>landa && q<time)
			{
				
				y=max_ind;
				delta = max_delta;
				int x=0;
					while(delta>landa && q<time && x<5)
					{
						
						x++;
						double d2Edx2=0;
						double d2Edxdy=0;
						double d2Edy2=0;
						
						double dx=0;
						double dy=0;
						deltax=0;
						deltay=0;
						
						for(int j=0;j<vv;j++)
						{
							if(max_ind!=j)
							{
							double xpp=(pos[t-1][0][max_ind]-pos[t-1][0][j]);
							double ypp=(pos[t-1][1][max_ind]-pos[t-1][1][j]);
							d2Edx2 += (K/Math.pow((double)dij[max_ind][j],2.0))*(1-((L*(double)dij[max_ind][j]*ypp*ypp)/(Math.pow(Math.pow(xpp, 2)+Math.pow(ypp, 2), 1.5))));
							d2Edxdy +=(K/Math.pow((double)dij[max_ind][j],2.0))*(((L*(double)dij[max_ind][j]*xpp*ypp)  /(Math.pow(Math.pow(xpp, 2)+Math.pow(ypp, 2), 1.5))));
							d2Edy2 += (K/Math.pow((double)dij[max_ind][j],2.0))*(1-((L*(double)dij[max_ind][j]*xpp*xpp)/(Math.pow(Math.pow(xpp, 2)+Math.pow(ypp, 2), 1.5))));
							
							deltax += (K/Math.pow((double)dij[max_ind][j],2.0))*(xpp-((L*(double)dij[max_ind][j]*xpp)/(Math.sqrt(Math.pow(xpp, 2)+Math.pow(ypp, 2)))));
							deltay += (K/Math.pow((double)dij[max_ind][j],2.0))*(ypp-((L*(double)dij[max_ind][j]*ypp)/(Math.sqrt(Math.pow(xpp, 2)+Math.pow(ypp, 2)))));
							}
						}
						deltax=-deltax;
						deltay=-deltay;
						
						dy = (deltay-(deltax*d2Edxdy)/d2Edx2)/(d2Edy2-(d2Edxdy*d2Edxdy)/d2Edx2);
						dx = (deltax-(d2Edxdy*dy))/d2Edx2;
						
						if(animation)
						{
							for(int i=0;i<vv;i++)
							{
								pos[t][0][i]=pos[t-1][0][i];
								pos[t][1][i]=pos[t-1][1][i];
							}
						}
						else
						{
							for(int i=0;i<vv;i++)
							{
								pos[t-1][0][i]=pos[t-1][0][i];
								pos[t-1][1][i]=pos[t-1][1][i];
							}
						}
						if(animation)
						{
							pos[t][0][max_ind]=pos[t-1][0][max_ind]+dx*cof;
							pos[t][1][max_ind]=pos[t-1][1][max_ind]+dy*cof;
						}
						else
						{
							pos[t-1][0][max_ind]=pos[t-1][0][max_ind]+dx*cof;
							pos[t-1][1][max_ind]=pos[t-1][1][max_ind]+dy*cof;							
						}
						if(animation)
							t++;   
						q++;
						delta = Math.sqrt(Math.pow(deltax, 2)+Math.pow(deltay, 2));
					}
					
					///--------------	
					deltax =0;
					deltay =0;
					delta=0;
					max_delta=0;
					for(int i=0;i<vv;i++)
					{
						if(i!=y)
						{
							deltax =0;
							deltay =0;
							delta=0;
							for(int j=0;j<vv;j++)
							{
								if(i!=j)
								{
								double xpp=(pos[t-1][0][i]-pos[t-1][0][j]);
								double ypp=(pos[t-1][1][i]-pos[t-1][1][j]);
								deltax += (K/Math.pow((double)dij[i][j],2.0))*(xpp-((L*(double)dij[i][j]*xpp)/(Math.sqrt(Math.pow(xpp, 2)+Math.pow(ypp, 2)))));
								deltay += (K/Math.pow((double)dij[i][j],2.0))*(ypp-((L*(double)dij[i][j]*ypp)/(Math.sqrt(Math.pow(xpp, 2)+Math.pow(ypp, 2)))));
								}
							}
							delta = Math.sqrt(Math.pow(deltax, 2)+Math.pow(deltay, 2));
							if(delta > max_delta  )
							{
								max_delta = delta;
								max_ind=i;
							}
						}
					}
						
					//---------------	
					
				}
			//System.out.println("t: "+t);
			if(animation)
			{
				q--;
				if(q<time)
				{
					for(int b=t;b<time;b++)
					{
						for(int i=0;i<vv;i++)
						{
							pos[b][0][i]=pos[b-1][0][i];
							pos[b][1][i]=pos[b-1][1][i];
						}
					}
				}
			}	
			pos=main_fitting_in_center(v, pos, (int) time, !animation);
			pos=limit_in_screen(v, pos, (int) time, !animation);
			pos=covert_to_glortho(v, pos, (int) time, !animation);
			return pos;
		
		}
		
	
	
	public static double [][][] main_fitting_in_center(int v, double [][][] pos, int time, boolean stat)
	{
		double max_x,max_y,min_x,min_y, ave_x,ave_y,max;
		
		
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
		return pos;
		
	}
	
	
	public static double [][][] limit_in_screen(int v, double [][][] pos, int time, boolean stat)
	{	
		double lm=0.9;
		int m_time=time;
		
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
						}
					}
					
					if(pos[t][1][i]>lm || pos[t][1][i]<-lm)
					{
						double m=Math.abs(lm/pos[t][1][i]);
						for (int i1 = 0; i1 < v; i1++) 
						{
							pos[t][0][i1]=pos[t][0][i1]*m;
							pos[t][1][i1]=pos[t][1][i1]*m;
						}
					}
				}
			}
		}
		
		return pos;
	}
	



	public static double [][][] covert_to_glortho(int v, double [][][] pos, int time, boolean stat)
	{
		double scale=0.7;
		
		int m_time=time;
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
		
		return pos;
	}


	
	
	
	
	
	
	static int [][] matrix_change_to_distance_cells(int adj[][], int v)
	{
		int temp [][]=new int[v][v];
		
		for(int i=0;i<v;i++)
			for(int j=0;j<v;j++)
				if(adj[i][j]==1)
					temp[i][j]=1;
				else
					temp[i][j]=0;
		
		for(int dis=2;dis<v-1;dis++)
		{
			for(int i=0;i<v;i++)
			{
				for(int j=0;j<v;j++)
				{
					if(temp[i][j]!=0)
					{
						for(int z=j;z<v;z++)
						{
							if(z!=j && temp[i][z]!=0 && (temp[i][j]+temp[i][z])==dis && temp[j][z]==0 && temp[z][j]==0)
							{
								temp[j][z]=dis;
								temp[z][j]=dis;
							}	
						}
					}
				}
			}
		}
		longest_theoritical_disctance=0;
		for(int i=0;i<v;i++)
		{
			int c=0;
			for(int j=0;j<v;j++)
			{	
				c+=temp[i][j];
				if(temp[i][j]>longest_theoritical_disctance)
					longest_theoritical_disctance=temp[i][j];
			}
		}
		
		
		return temp;
	}
	
	
	
}
