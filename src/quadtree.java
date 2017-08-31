
import java.util.ArrayList;

public class quadtree implements Cloneable {

		private double [] bounds = new double [4]  ;
		private int indTree;
		private quadtree parent;
		private ArrayList <double []> nodes;
		private boolean leaf;
		private int ind;  // index amoung childern
		private int level;
		private ArrayList<quadtree> childeren = new ArrayList<quadtree>();
		static quadtree myTree ;
		//static int temp_tree_index=0;
		static int test =0;
		static double dist =0 ;
		static double counter =0 ;
		static double whole_force=0; 
		static double [][] simi;
		static double xx ;
		static double yy ;
		
		
		public quadtree(double [] bounds ,quadtree parent, ArrayList <double []> nodes, boolean leaf, int ind, int level, ArrayList<quadtree> childeren)
		{
			this.bounds = bounds.clone();
			this.parent=parent;
			this.nodes= (ArrayList<double []>) nodes.clone();
			this.leaf=leaf;
			this.ind=ind;
			this.level=level;
			this.indTree = indTree;
			//this.childeren= (ArrayList<Tree>) childeren.clone();
				
		}
		
	
		public double[] getCentre()
		{
			double xy []= new double [2];
			xy[0]=(this.bounds[1]-this.bounds[0])/2+this.bounds[0];
			xy[1]=(this.bounds[3]-this.bounds[2])/2+this.bounds[2];
			
			return xy;
		}
		
		public int getLevel()
		{
			return this.level;
		}
		
		public boolean getLeaf()
		{
			return this.leaf;
		}
		
		public ArrayList<double []> getChil()
		{
			return this.nodes;
		}
		
		public int noChild()
		{
			return this.nodes.size();
		}
		
		public ArrayList <quadtree> getChilds_Tree()
		{
			return this.childeren;
		}
		
		
		public int inexTree()
		{
			return this.indTree;
		}
		
		public void setleaf(boolean leaf)
		{
			this.leaf=leaf;
		}
		public void setChild(quadtree c)
		{
			this.childeren.add(c);
		}
		
		
		public static int creat_childern(double [] bounds,quadtree parent, ArrayList <double[]> nodes, int level )
		{
			CFF.chh++;
			quadtree t0,t1,t2,t3;
			double [] bound = new double [4];
			//System.out.println(1);
			/// first child
			bound[0]=bounds[0];bound[1]=(bounds[1]-bounds[0])/2+bounds[0];
			bound[2]=bounds[2];bound[3]=(bounds[3]-bounds[2])/2+bounds[2];
			ArrayList  <double []> temp = new ArrayList<double []>();
			temp.clear();
			for(int i=0;i<nodes.size();i++)
			{
				if(nodes.get(i)[0]>=bound[0] &&  nodes.get(i)[0]<=bound[1] &&
						nodes.get(i)[1]>=bound[2] &&  nodes.get(i)[1]<=bound[3])
					temp.add(nodes.get(i));
			}
			if(temp.size()>0)
			{
				ArrayList<quadtree> c0 = new ArrayList<quadtree>();
				t0 = new quadtree(bound,parent,temp,false,0,level+1,c0);
				parent.setChild(t0);
				
				if(temp.size()==1)
					t0.setleaf(true);
				
				if(temp.size()>1 && CFF.chh<CFF.vv)
					creat_childern(bound,t0 ,temp, level+1 );
			}
			
			
			/// second child
			bound[0]=(bounds[1]-bounds[0])/2+bounds[0];bound[1]=bounds[1];
			bound[2]=bounds[2];bound[3]=(bounds[3]-bounds[2])/2+bounds[2];
			temp.clear();
			for(int i=0;i<nodes.size();i++)
			{
				if(nodes.get(i)[0]>bound[0] &&  nodes.get(i)[0]<=bound[1] &&
						nodes.get(i)[1]>=bound[2] &&  nodes.get(i)[1]<=bound[3])
					temp.add(nodes.get(i));
			}
			
			if(temp.size()>0)
			{
				ArrayList<quadtree> c1 = new ArrayList<quadtree>();;
				t1 = new quadtree(bound,parent,temp,false,1,level+1,c1);
				parent.setChild(t1);
				
				if(temp.size()==1)
					t1.setleaf(true);
				
				if(temp.size()>1 && CFF.chh<CFF.vv)
					creat_childern(bound,t1 ,temp, level+1 );
			}
			
			
			/// third child
			bound[0]=bounds[0];bound[1]=(bounds[1]-bounds[0])/2+bounds[0];
			bound[2]=(bounds[3]-bounds[2])/2+bounds[2];bound[3]=bounds[3];
			temp.clear();
			for(int i=0;i<nodes.size();i++)
			{
				if(nodes.get(i)[0]>=bound[0] &&  nodes.get(i)[0]<=bound[1] &&
						nodes.get(i)[1]>bound[2] &&  nodes.get(i)[1]<=bound[3])
					temp.add(nodes.get(i));
			}
			if(temp.size()>0)
			{
				ArrayList<quadtree> c2 = new ArrayList<quadtree>();;
				t2 = new quadtree(bound,parent,temp,false,2,level+1,c2);
				parent.setChild(t2);
				
				if(temp.size()==1)
					t2.setleaf(true);
				
				if(temp.size()>1 && CFF.chh<CFF.vv)
					creat_childern(bound,t2 ,temp, level+1 );
			}
			
			
			/// forth child
			bound[0]=(bounds[1]-bounds[0])/2+bounds[0];bound[1]=bounds[1];
			bound[2]=(bounds[3]-bounds[2])/2+bounds[2];bound[3]=bounds[3];
			temp.clear();
			for(int i=0;i<nodes.size();i++)
			{
				if(nodes.get(i)[0]>bound[0] &&  nodes.get(i)[0]<=bound[1] &&
						nodes.get(i)[1]>bound[2] &&  nodes.get(i)[1]<=bound[3])
					temp.add(nodes.get(i));
			}
			if(temp.size()>0)
			{
				ArrayList<quadtree> c3 = new ArrayList<quadtree>();;
				t3 = new quadtree(bound,parent,temp,false,3,level+1,c3);
				parent.setChild(t3);
				
				if(temp.size()==1)
					t3.setleaf(true);
				
				if(temp.size()>1 && CFF.chh<CFF.vv)
					creat_childern(bound,t3 ,temp, level+1 );
			}
			
			return 0;
		}
		
	
		public static void print_tree(quadtree res)
		{
			 System.out.println("Leaf:  "+res.getLeaf()+"  level:  "+res.getLevel()+" No.Childs  "+ res.noChild());
			 if(!res.getLeaf())
			 {
				 ArrayList<quadtree> children = new ArrayList<quadtree>();
				 children=res.getChilds_Tree();
				 for(int i=0;i<children.size();i++)
					 print_tree(children.get(i));
			 }
		}
		
		
		public  static void traversing_tree(quadtree tree, double x, double y, int vx)
		{
			
			ArrayList<quadtree> childs = new ArrayList<quadtree>();
			childs = tree.getChilds_Tree(); 
			int c0=0,c1=0,c2=0,c3=0;
			for(int i=childs.size()-1;i>=0;i--)
		    {
				double xy [] = new double [2];
				double L=0;
				xy = childs.get(i).getCentre();
				L = (childs.get(i).bounds[1]-childs.get(i).bounds[0]);
				double d = Math.sqrt((xy[0]-x)*(xy[0]-x)+(xy[1]-y)*(xy[1]-y));
				//double d = Math.abs(xy[0]-x)+Math.abs(xy[1]-y);
				if(d!=0)
					if((L/d)<0.5 || (childs.get(i).getChil().size()==1 && childs.get(i).getChilds_Tree().size()==0 ) )
						if(childs.get(i).getChil().get(0)[2]!=vx)
						{
							System.out.println((L/d));
							
							double x1 = (xy[0]-x);
							double y1 = (xy[1]-y);
							double n = childs.get(i).getChil().size();
							
							xx += (x1/d)*n;
							yy += (y1/d)*n;
						
							//System.out.println(n+"  "+d+"  "+xx+" "+yy+" ");
							if(i==0) c0=1;
							if(i==1) c1=1;
							if(i==2) c2=1;
							if(i==3) c3=1;
							
							//childs.remove(i);
						}
				counter++;
		    }
			//childs = tree.getChilds_Tree(); 
			for(int i=childs.size()-1;i>=0;i--)
		    {
				if(i==0 && c0==0)
					traversing_tree(childs.get(i), x, y,vx);
				if(i==1 && c1==0)
					traversing_tree(childs.get(i), x, y,vx);
				if(i==2 && c2==0)
					traversing_tree(childs.get(i), x, y,vx);
				if(i==3 && c3==0)
					traversing_tree(childs.get(i), x, y,vx);
		    }
			
			
		}
		
		
		public  static  void XYtravers_tree(quadtree tree, double x, double y, int vx)
		{
			
			ArrayList<quadtree> childs = new ArrayList<quadtree>();
			childs = tree.getChilds_Tree(); 
			int c0=0,c1=0,c2=0,c3=0;
			for(int i=childs.size()-1;i>=0;i--)
		    {
				
				
				if (childs.get(i).getChil().size()==1 && childs.get(i).getChilds_Tree().size()==0 ) 
				{
					double xy [] = new double [2];
					double L=0;
					xy = childs.get(i).getCentre();
					L = childs.get(i).bounds[1]-childs.get(i).bounds[0];
					double d = Math.sqrt((xy[0]-x)*(xy[0]-x)+(xy[1]-y)*(xy[1]-y));
					
					
					if(d!=0)
					if((L/d)>=0.1 /*|| (childs.get(i).getChil().size()==1 && childs.get(i).getChilds_Tree().size()==0 )*/ )
					if(childs.get(i).getChil().get(0)[2]!=vx)
					{
						double x1 = (xy[0]-x);
						double y1 = (xy[1]-y);
						double n = childs.get(i).getChil().size();
						xx=0;
						yy=0;
						xx = (x1/d)*n;
						yy = (y1/d)*n;
						System.out.println(xx+" "+yy+"  "+x1+"  "+y1+" "+d+" "+n);
						//dist += n*d;
						if(i==0) c0=1;
						if(i==1) c1=1;
						if(i==2) c2=1;
						if(i==3) c3=1;
						whole_force=0;
						/*for(int b=0;b<n;b++)
						{
							whole_force += simi[(int) childs.get(i).getChil().get(b)[2]][vx];
						}
						//whole_force = whole_force*dist;
						xx = whole_force*xx;
						yy = whole_force*yy;*/
						
						
						//System.out.println("     "+d+"     "+n);
						//childs.remove(i);
					}
				}
				counter++;
		    }
			//childs = tree.getChilds_Tree(); 
			for(int i=childs.size()-1;i>=0;i--)
		    {
				if(i==0 && c0==0)
					XYtravers_tree(childs.get(i), x, y,vx);
				if(i==1 && c1==0)
					XYtravers_tree(childs.get(i), x, y,vx);
				if(i==2 && c2==0)
					XYtravers_tree(childs.get(i), x, y,vx);
				if(i==3 && c3==0)
					XYtravers_tree(childs.get(i), x, y,vx);
		    }
			
			
		}
		
		
		
	
		
		public static void main(String[] args) throws CloneNotSupportedException  
		{
		     int v = 225;
		     double [][] vert = new double [2][v];
		     double [] bound = new double [4];
		     ArrayList<double []> temp = new ArrayList<double []>();
		     bound[0]=1000;bound[1]=-10000;bound[2]=10000;bound[3]=-10000;
		     simi = new double [v][v];
		   
		     
		     for(int i=0;i<v;i++)
		     {
		    	 for(int j=i+1;j<v;j++)
		    	 {
		    		 simi[i][j] =Math.random()*v;
		    		 simi[j][i] = simi[i][j];
		    	 }
		    	 
		    	 
		    	 double [] cord = new double [3];
		    	 vert[0][i] = Math.random()* (1 - 0.0);;
		    	 vert[1][i] = Math.random()* (1 - 0.0);;
		    	 
		    	 
		    	 cord[0] = vert[0][i];
		    	 cord[1] = vert[1][i];
		    	 cord[2] = i;
		    	 temp.add(cord);
		    	 if(bound[0]>vert[0][i])
		    		 bound[0]=vert[0][i];
		    	 
		    	 if(bound[1]<vert[0][i])
		    		 bound[1]=vert[0][i];
		    	 
		    	 if(bound[2]>vert[1][i])
		    		 bound[2]=vert[1][i];
		    	 
		    	 if(bound[3]<vert[1][i])
		    		 bound[3]=vert[1][i];
		    	 
		    	 //System.out.println("X:  "+vert[0][i]+"  Y:  "+vert[1][i]);
		     }
		     long start = System.nanoTime();
		     ArrayList<quadtree> chs = new ArrayList<quadtree>();
		     quadtree main_root = new quadtree(bound,null,temp,false,3,0,chs);
		     
		     creat_childern(bound,main_root ,temp, 0 );
		     for(int i=0;i<v;i++)
		     {
		    	dist=0;
				counter=0;
				xx=0;
				yy=0;
				whole_force=0;
				traversing_tree(main_root, vert[0][i],vert[1][i],i);
				System.out.println(xx+"  "+yy);
				 //*******************************************
		     }
		     
			 
			
		     long end = System.nanoTime();
		    // print_tree(main_root);
		     System.out.println("Time:  "+(end-start));
		}
		
		
				
	}
