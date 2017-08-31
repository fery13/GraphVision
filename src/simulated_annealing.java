
import java.util.ArrayList;

public class simulated_annealing {
	
	static double fitness;
	static int check [][];
	public static double [][] initial_placement (int v, int [][] adj, double gen)
	{
		
		
		int ordering [] = new int [v];
		
		double initial_placement [][] = new double [3][v];
		
		ArrayList <Integer> temp = new ArrayList <Integer>();
		for(int j=0;j<v;j++)
			temp.add(j);
		
		int q=v;
		
		
		for(int i=0;i<v;i++)
		{
			int b=(int) (Math.random()*(q-0));
			ordering[i]=temp.get(b);
			temp.remove(b);
			q--;
		}


		fitness=fitness(v, ordering, adj, adj);
		
		int generation = (int) (v*gen);
		for(int i=generation;i>0;i--)
		{
			double t=(double) i / gen;
			//t=Math.sqrt(v);
			ordering=mutation (v,ordering, t, fitness, adj);
			
		}
		
		int[] temp1 = new int [v];
		temp1= ordering.clone();
		double min=0;
		for(int i=0;i<v;i++)
		{
			for(int j=0;j<v;j++)
			{
				if(temp1[j]==min)
				{
					ordering[i]=j;
					min++;
					break;
				}
			}
		}
		for(int i=0;i<v;i++)
		{
			initial_placement[0][i]=Math.cos(((Math.PI*2.0)/((double) v))* ((double)ordering[i]));
			initial_placement[1][i]=Math.sin(((Math.PI*2.0)/((double) v))* ((double)ordering[i]));
			initial_placement[2][i]=(double)i/v;
		}
		
		return initial_placement;
	}
	
	public static int [] mutation (int v, int [] ordering, double t, double fittt, int adj [][])
	{
		int ordered [] = new int [v];
		
		ordered=ordering.clone();
		t=t/2.0;
		for(int i=0;i<t;i++)
		{
			int k1=(int) (Math.random()*(v-0));
			int k2=(int) (Math.random()*(v-0));
			
			while(k1==k2)
			{
				k1=(int) (Math.random()*(v-0));
				k2=(int) (Math.random()*(v-0));
			}
			 int temp = ordered[k1];
			 ordered[k1]=ordered[k2];
			 ordered[k2]=temp;
			 
			 //*****************
			
			 //*****************
		}
		
		double fit=fitness(v, ordered, adj, adj);
		if(fit<fittt)
		{
			fitness=fit;
			
			//System.out.println(fitness);
			return ordered;
		}
		else
		{
			//System.out.println(fitness);
			return ordering;
		}
		
		
		
	}
	
	
 	public static double fitness(int v, int pos[], int dis [][], int adj [][])
	{
		double fitness=0;
		double chunk = (Math.PI*2)/(double) v;
		
		for(int i=0;i<v;i++)
		{
			for(int j=0;j<v;j++)
			{
				if(i!=j )
				{
					
					
					double physical =0;
				
					if(dis[i][j]==1)
					{
						int p=0;
						int q=0;
						
						for(p=0;p<v;p++)
							if(pos[p]==i)
								break;
						
						for(q=0;q<v;q++)
							if(pos[q]==j)
								break;
						
						
						if(Math.abs(p-q)> ((double)(v/2) ))
							physical=v-Math.abs(p-q);
						else
							physical=Math.abs(p-q);
						
						//physical=Math.abs(pos[i]-pos[j])*chunk;
						fitness+=Math.abs(   (physical)-dis[i][j]  );
					}
				}
			}
		}
		
		return fitness;
	}

 	//*******************************************
 	//*********************************
 	//**********************************
 	
	public static double [][] initial_placement_w (int v, int [][] adj)
	{
		
		
		
		int ordering [] = new int [v];
		
		double initial_placement [][] = new double [2][v];
		
		ArrayList <Integer> temp = new ArrayList <Integer>();
		for(int j=0;j<v;j++)
			temp.add(j);
		
		int q=v;
		
		
		for(int i=0;i<v;i++)
		{
			int b=(int) (Math.random()*(q-0));
			ordering[i]=temp.get(b);
			temp.remove(b);
			q--;
		}

		double gen =100;
		int generation = (int) (v*gen);
		double pp= v/5;
		for(int i=generation;i>0;i--)
		{
			double t=(double) i / gen;
			//t=Math.sqrt(v);
			
			ordering=mutation_w (v,ordering, pp,adj);
			pp = pp-0.005;
		}
		
		for(int i=0;i<v;i++)
		{
			initial_placement[0][i]=Math.cos((Math.PI*2.0)/((double) v)* ((double)ordering[i]));
			initial_placement[1][i]=Math.sin((Math.PI*2.0)/((double) v)* ((double)ordering[i]));
		}
		
		return initial_placement;
	}
	
	public static int [] mutation_w (int v, int [] ordering, double t, int dis [][])
	{
		int ordered [] = new int [v];
		
		ordered=ordering.clone();
		t=t/2.0;
		for(int i=0;i<t;i++)
		{
			int k1=(int) (Math.random()*(v-0));
			int k2=(int) (Math.random()*(v-0));
			
			while(k1==k2)
			{
				k1=(int) (Math.random()*(v-0));
				k2=(int) (Math.random()*(v-0));
			}
			 int temp = ordered[k1];
			 ordered[k1]=ordered[k2];
			 ordered[k2]=temp;
			 
			 //*****************
			double fit1=0;double fit2=0;
			for(int j=0;j<v;j++)
			{
				if(k1!=j)
				{
					if(dis[k1][j]==1)
					{
						double physical =0;
						
						if(Math.abs(ordering[k1]-ordering[j])> ((double)(v/2) ))
							physical=v-Math.abs(ordering[k1]-ordering[j]);
						else
							physical=Math.abs(ordering[k1]-ordering[j]);
						
						fit1+=Math.abs(  ( dis[k1][j])- (physical)  );
					}
				}
				if(k2!=j)
				{
					if(dis[k2][j]==1)
					{
					double physical =0;
					
					if(Math.abs(ordering[k2]-ordering[j])> ((double)(v/2) ))
						physical=v-Math.abs(ordering[k2]-ordering[j]);
					else
						physical=Math.abs(ordering[k2]-ordering[j]);
					
					fit1+=Math.abs(  ( dis[k2][j])- (physical)  );
					}
				}
				if(k1!=j)
				{
					if(dis[k1][j]==1)
					{
					double physical =0;
					
					if(Math.abs(ordered[k1]-ordered[j])> ((double)(v/2) ))
						physical=v-Math.abs(ordered[k1]-ordered[j]);
					else
						physical=Math.abs(ordered[k1]-ordered[j]);
					
					fit2+=Math.abs(  ( dis[k1][j])- (physical)  );
					}
				}
				if(k2!=j)
				{
					if(dis[k2][j]==1)
					{
					double physical =0;
					
					if(Math.abs(ordered[k2]-ordered[j])> ((double)(v/2) ))
						physical=v-Math.abs(ordered[k2]-ordered[j]);
					else
						physical=Math.abs(ordered[k2]-ordered[j]);
					
					fit2+=Math.abs(  ( dis[k2][j])- (physical)  );
					}
				}
			}
			if(fit1<fit2)
				ordered=ordering.clone();
			else
				ordering=ordered.clone();
			 //*****************
		}
		

		return ordered;
	}
	
}
