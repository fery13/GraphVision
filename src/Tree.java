import java.util.ArrayList;



public class Tree {
	
	private static int root1 =-1;
	private static int root2 =-1;
	private static int [][] tree_feature;
	private static double [][][] sync;
	private static double max_layer=0;
	public static double [][][] tree_angular_radial1(int adj [][], int v, boolean animation, double links, double time, boolean t3d)
	{
		
		tree_feature_extraction(adj , v,links);
		
		double intervale=4.0/v;
		sync = new double [1][1][v];
		for(int i=0;i<v;i++)
		{
			sync[0][0][i]=(double) (Math.random() * ((Math.PI*2 - 0.0)));
		}
		double interval [] = new double [v];
		
		double pos [][][];
		if(animation)
			pos = new double [(int) time][3][v];
		else
			pos = new double [1][3][v];
		
		
		//**********************************
		//tree_feature[0][i]  ==> degree of each node
		//tree_feature[1][i]  ==> index of its parent
		//tree_feature[2][i]  ==> number of sibling of the node
		//tree_feature[3][i]  ==> layer index of the node
		//tree_feature[4][i]  ==> index of the node among its sibling
		//tree_feature[5][i]  ==> number of nodes in the layer which the node is placed
		//tree_feature[6][i]  ==> weight of the nodes 
		//**********************************
		
		if(root2!=-1)
		{
			interval[root2]=Math.PI;
			interval[root1]=Math.PI;
			
			sync[0][0][root1]=0;
			sync[0][0][root2]=Math.PI;
		}
		else
		{
			interval[root1]=Math.PI*2;
			
			sync[0][0][root1]=0;
		}
		
		double l=0;
		for(int t=1;t<time;t++)
		{
			l++;
			for(int i=0;i<v;i++)
			{
				if(tree_feature[3][i]==l)
				{
					
					sync[0][0][i]=(interval[tree_feature[1][i]]/(tree_feature[2][i]) )*((tree_feature[4][i]-1))+sync[0][0][tree_feature[1][i]];
					interval[i]=interval[tree_feature[1][i]]/((tree_feature[2][i]) );
					
				}
				if(animation)
				{
					pos[t][0][i]=Math.cos(sync[0][0][i])*(tree_feature[3][i]/max_layer);
					pos[t][1][i]=Math.sin(sync[0][0][i])*(tree_feature[3][i]/max_layer);
				}
				else
				{
					pos[0][0][i]=Math.cos(sync[0][0][i])*(tree_feature[3][i]/max_layer);
					pos[0][1][i]=Math.sin(sync[0][0][i])*(tree_feature[3][i]/max_layer);
				}
			}
			
			
			
		}
		
		
		
		pos=main_fitting_in_center(v, pos, (int) time, !animation,t3d );
		pos=limit_in_screen(v, pos, (int) time, !animation, t3d);
		pos=covert_to_glortho(v, pos, (int) time, !animation, t3d);
		
		
	return pos;
		
	}

	public static double [][][] tree_angular_radial2(int adj [][], int v, boolean animation, double links, double time, boolean t3d, boolean initial)
	{
		
		tree_feature_extraction(adj , v,links);
		
		double intervale=4.0/v;
		sync = new double [1][1][v];
		for(int i=0;i<v;i++)
		{
			sync[0][0][i]=(double) (Math.random() * ((Math.PI*2 - 0.0)));
		}
		double interval [] = new double [v];
		
		double pos [][][];
		if(animation)
			pos = new double [(int) time][3][v];
		else
			pos = new double [1][3][v];
		
		//**********************************
		//tree_feature[0][i]  ==> degree of each node
		//tree_feature[1][i]  ==> index of its parent
		//tree_feature[2][i]  ==> number of sibling of the node
		//tree_feature[3][i]  ==> layer index of the node
		//tree_feature[4][i]  ==> index of the node among its sibling
		//tree_feature[5][i]  ==> number of nodes in the layer which the node is placed
		//tree_feature[6][i]  ==> weight of the nodes 
		//**********************************
		
	
		if(root2!=-1)
		{
			interval[root2]=Math.PI;
			interval[root1]=Math.PI;
			
			sync[0][0][root1]=0;
			sync[0][0][root2]=Math.PI;
		}
		else
		{
			interval[root1]=Math.PI*2;
			
			sync[0][0][root1]=0;
		}
		tree_feature[6][root1]=v;
		double l=0;
		for(int t=1;t<time;t++)
		{
			l++;
			for(int i=0;i<v;i++)
			{
				if(tree_feature[3][i]==l)
				{
					
					sync[0][0][i]=(interval[tree_feature[1][i]]/(tree_feature[6][tree_feature[1][i]]) )*((tree_feature[7][i]))+sync[0][0][tree_feature[1][i]];
					interval[i]=(interval[tree_feature[1][i]]/((tree_feature[6][tree_feature[1][i]])  ))*tree_feature[6][i];
					
				}
				
				if(animation)
				{
					pos[t][0][i]=Math.cos(sync[0][0][i])*(tree_feature[3][i]/max_layer);
					pos[t][1][i]=Math.sin(sync[0][0][i])*(tree_feature[3][i]/max_layer);
				}
				else
				{
					pos[0][0][i]=Math.cos(sync[0][0][i])*(tree_feature[3][i]/max_layer);
					pos[0][1][i]=Math.sin(sync[0][0][i])*(tree_feature[3][i]/max_layer);
				}
				
				
			}
			
			
			
		}
		
		if(!initial)
		{
			pos=main_fitting_in_center(v, pos, (int) time, !animation,t3d );
			pos=limit_in_screen(v, pos, (int) time, !animation, t3d);
			pos=covert_to_glortho(v, pos, (int) time, !animation, t3d);
		}
		
		
	return pos;
		
	}


	

	public static void tree_feature_extraction(int adj1 [][], int v,double links)
	{
	int temp[][]=new int[v][v];
	int adj[][]=new int[v][v];
	tree_feature = new int [10][v];
	
	for(int i=0;i<links;i++)
	{
		int r=adj1[0][i];
		int s=adj1[1][i];
		adj[r][s]=1;
		adj[s][r]=1;
		temp[r][s]=1;
		temp[s][r]=1;
	}
	
	 
	// ROOT identification 
	ArrayList<Integer> temp_adj=new ArrayList<Integer>();
	for(int i =0; i<v;i++)
	{
		temp_adj.add(i);
	}
	//*************
	
	while(temp_adj.size()>2)
	{
		for(int i=0;i<v;i++)
		{
			int p=0;
			for(int j=0;j<v;j++)
			{
				if(temp[i][j]==1)
					p++;			
			}		
			if(p==1 )
			{
				for(int j=0;j<v;j++)
				{
					temp[i][j]=0;
					temp[j][i]=0;
				}
				for(int j=0;j<temp_adj.size();j++)
				{
					if(temp_adj.get(j)==i)
						temp_adj.remove(j);
				}
			}
		}
	}
	 
	if(temp_adj.size()==1)
		root1=(int) temp_adj.get(0);
	
	if(temp_adj.size()==2)
	{
		root1=(int) temp_adj.get(0);
		root2=(int) temp_adj.get(1);		
	}
	
	//*********************
	/*root1=0;
	root2=-1;*/
	
	  
   // END ROOT identification 
      
     
      //DEGREE assigning
      for(int i=0;i<v;i++)
      {
    	  int s=0;
    	  for(int j=0;j<v;j++)
    	  {
    		  if(adj[i][j]==1)
    			  s++;
    	  }
    	  tree_feature[0][i]=s;
    	  tree_feature[1][i]=-1;
    	  tree_feature[2][i]=-1;
    	  tree_feature[3][i]=-1;
    	  tree_feature[4][i]=-1;
      }
      //END DEGREE ASSIGNING
      
      
      //Parent finding
      if(root2==-1)
    	 tree_feature[1][root1]=0;
      else
      {
    	  tree_feature[1][root1]=0;
    	  tree_feature[1][root2]=0;
      }
    	int u=0;
      while(u==0)
    	  {
    	  for(int i=0;i<v;i++)
    		  for(int j=0;j<v;j++)
    			  if(adj[i][j]==1 && tree_feature[1][j]!=-1 && tree_feature[1][i]==-1)
    				 tree_feature[1][i]=j;
    		u=1;
    	  for(int i=0;i<v;i++)
    		  if(tree_feature[1][i]==-1)
    			  u=0;
    	  }
    	
     
    //END Parent finding
      
      
      //number of sibling
      for(int i=0;i<v;i++)
	  {
    	  if(i==root1 || i==root2)
    		  tree_feature[2][i]=-1;
    	  else
    		  {
    		  if(tree_feature[1][i]==root1 ||tree_feature[1][i]==root2)
    			  tree_feature[2][i]=tree_feature[0][tree_feature[1][i]];
    	      else
    	    	  tree_feature[2][i]=tree_feature[0][tree_feature[1][i]]-1;
    	      }
	  }
      //End number of sibling
     
      
      //Finding index of the layer
      u=0;
      while(u==0)
      {
      for(int i=0;i<v;i++)
	  {
    	  if(i==root1 || i==root2)
    		  tree_feature[3][i]=0;
    	  else
    		  for(int j=0;j<v;j++)
    		    if(adj[i][j]==1 && tree_feature[3][j]!=-1 && tree_feature[3][i]==-1)
    		    	tree_feature[3][i]=tree_feature[3][j]+1;
	  }
      
      u=1;
	  for(int i=0;i<v;i++)
		  if(tree_feature[3][i]==-1)
			  u=0;
      
      }
      max_layer=0;
      for(int i=0;i<v;i++)
	   if(tree_feature[3][i]>max_layer)
		  max_layer=tree_feature[3][i];
      
     
	  //End of Finding index of layer
      
     
      //Finding the index for each node among its siblings
      if(root2!=-1)
      {
    	  tree_feature[4][root1]=1;
    	  tree_feature[4][root2]=1;
      }
      else
    	  tree_feature[4][root1]=1;
      
      for(int i=0;i<v;i++)
    	  {
    	  int ch=1;
    		  for(int j=0;j<v;j++)
    		  {
    			  if(adj[i][j]==1 && tree_feature[3][i]<tree_feature[3][j])
    			  {
    				 tree_feature[4][j]=ch;
    				  ch++;
    			  }
    		  }
    	  }
       //END Finding the index for each node among its siblings
      
      
      
      //Number of nodes in the level that i is placed
      for(int i=0;i<v;i++)
	   {
		int c=0;
		for(int j=0;j<v;j++)
		   {
			if(tree_feature[3][j]==tree_feature[3][i])
			 c++;
		   }
		tree_feature[5][i]=c; 
	   }
      //END Number of nodes in the level that i is placed
     
      
     
      
      
      //Finding the weight for each node
      for(int j=0;j<v;j++)
    	  tree_feature[6][j]=1;
       
		for(int i=(int)max_layer-1;i>0;i--)
		 for(int j=0;j<v;j++)
		  if(tree_feature[3][j]==i)
		   for(int q=0;q<v;q++)
			if(adj[j][q]==1 && tree_feature[3][j]<tree_feature[3][q])
				tree_feature[6][j]=tree_feature[6][q]+tree_feature[6][j];
      //End of finding the weight for nodes
		
		
		 
		//Finding the accumlative of nodes' weight among its siblings
		
		for(int j=0;j<v;j++)
		{
			 tree_feature[7][j]=0;
			 tree_feature[8][j]=-1;
		}
		
		for(int i=(int) (max_layer);i>=0;i--)
		{
			 for(int j=0;j<v;j++)
			 {
			  if(tree_feature[3][j]==i && tree_feature[8][j]==-1)
			  {
				  int acc=tree_feature[6][j];
				  tree_feature[7][j]=0;
				  tree_feature[8][j]=0;
				  
				  for(int q=0;q<v;q++)
				  {
					  if( j!=q && tree_feature[1][j] == tree_feature[1][q] && tree_feature[8][q]==-1)
					  {
						   tree_feature[7][q]=acc;
						   acc+=tree_feature[6][q];
						   tree_feature[8][q]=0;
						  
					  }
				  }
			  }
			}
		}
		
		
		
		
		//End of finding the accumlative of nodes' weight among its siblings

		
	
		System.out.println("Finishing tree feature extracting ...  !!  ");
      
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
