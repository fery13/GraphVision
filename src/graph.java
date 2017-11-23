import static org.lwjgl.opengl.GL11.GL_COLOR_BUFFER_BIT;
import static org.lwjgl.opengl.GL11.GL_LINES;
import static org.lwjgl.opengl.GL11.glBegin;
import static org.lwjgl.opengl.GL11.glClear;
import static org.lwjgl.opengl.GL11.glClearColor;
import static org.lwjgl.opengl.GL11.glColor3d;
import static org.lwjgl.opengl.GL11.glEnd;
import static org.lwjgl.opengl.GL11.glFlush;
import static org.lwjgl.opengl.GL11.glVertex2f;
import static org.lwjgl.opengl.GL11.glViewport;
import static org.lwjgl.opengl.GL11.glOrtho;
import static org.lwjgl.opengl.GL11.glBlendFunc;
import static org.lwjgl.opengl.GL11.glEnable;
import static org.lwjgl.opengl.GL11.glDisable;

import java.awt.EventQueue;
import java.awt.image.BufferedImage;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.ByteBuffer;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;

import javax.imageio.ImageIO;
import javax.swing.JOptionPane;

import org.lwjgl.BufferUtils;
import org.lwjgl.LWJGLException;
import org.lwjgl.input.Keyboard;
import org.lwjgl.input.Mouse;
import org.lwjgl.opengl.Display;
import org.lwjgl.opengl.GL11;
//import org.lwjgl.util.glu.Sphere;
import org.newdawn.slick.Color;
import org.newdawn.slick.TrueTypeFont;
import org.newdawn.slick.util.ResourceLoader;

import java.awt.Font;
import java.awt.FontFormatException;
import java.io.InputStream;
 
import org.lwjgl.opengl.DisplayMode;





public class graph implements in1{
	
	private static int v;
	
	private static boolean antiAlias = true;
	static int colors;
	static int iii;
	private static double [][][] pos;
	
	private static boolean edge_show = false;
	
	public static int [][] adj;
	private static int [][] edge;
	public static double [][] col;
	public static int cs [];
	private static boolean sORa;
	private static boolean tree;
	private static boolean adj_edge;
	private static boolean big_data;
	private static int groups;
	private static boolean stop=false;
	private static double time;
	
	public static int weidth=1000;
	public static int hight=1000;
	public static int depth=1000;
	private static int keyy =0;
	private static double zoom=1;
	
	
	private static double links;
	
	private static int counter;
	 public static boolean clicked = false;
	
	 private static float initial_x =0;
	 
	 private static float initial_y =0;
	 
	 private static float initial_z =0;
	 
	 public static String [] messages ;
	 
	
	public graph(double [][][] pos_ani, boolean tree,  int adj[][], double time, boolean sORa, int v, double links, boolean big_data, int groups )
	{
		this.adj=adj;
		this.tree=tree;
		this.pos=pos_ani;
		this.time=time;
		this.sORa=sORa;
		this.v=v;
		this.links=links;
		this.big_data=big_data;
		this.groups = groups;
	}
	
	
	public static void visualization() throws FileNotFoundException, FontFormatException
	{
		try 
        {	
    		counter=0;
        	Display.setDisplayMode(new org.lwjgl.opengl.DisplayMode(weidth, hight));
        	Display.setLocation(0, 0);
        	Display.setTitle("Graph");
            Display.create();
            glClearColor(1.0f,1.0f,1.0f,1.0f);
        	glOrtho(-weidth, weidth*2, hight*2, -hight,  -depth*5,depth*5);
            glFlush();
            
        } 
        
        catch (LWJGLException e) 
        {
        	 e.printStackTrace();
           // System.err.println("Display wasn't initialized correctly.");
            System.exit(1);
        }
		 int h=0;
		 boolean snapp = true;
		 DateFormat df = new SimpleDateFormat("dd/MM/yy HH:mm:ss");
	       Date dateobj = new Date();
	       //System.out.println(df.format(dateobj));
		while (!Display.isCloseRequested() /*&& counter !=variables.time-5*/ ) 
        {	
    		update();
    		Display.update();
            Display.sync(60);
            glFlush();
            if(graphvision.SaveOne.isSelected() && sORa && snapp)
            {
            	snapp=false;
            	//snapshot(df.format(dateobj)+"");
            	snapshot(df.format(dateobj)+"");
            	
            }
            if(graphvision.SaveAll.isSelected() && !sORa && snapp)
            {
            	//snapshot(df.format(dateobj)+" "+counter);
            	snapshot(" "+counter);
            }
            	
            
            
            /*if(graphvision.choice.getSelectedItem()=="CTD" && h==0)
            {
            	
            	h=CTD.check();
            	if(h==2)
            		break;
            }*/
        }
    	Display.destroy();
    	
    	/*if(h==2)
    		CTD.addings(pos,  v,  graphvision.mat, adj, links, false);*/
    	
    	
    	
    	
    	
    	
    	
	}
	
	  public static void update() throws FileNotFoundException, FontFormatException
	    {
	    	 glClear(GL_COLOR_BUFFER_BIT);
	    	 glClearColor(1.0f,1.0f,1.0f,1.0f);
	    	 glViewport(0,0,weidth,hight);
	    	 
	    	 /*if(graphvision.choice.getSelectedItem()=="Circular-Force-Free")
	    		 text(50, messages[counter],-992, -992,Color.red );*/
	    	
	    	 
	    	 if(Keyboard.isKeyDown(Keyboard.KEY_1))
	    	 {
	    		 keyy=1;
	    	 }
	    	 if(Keyboard.isKeyDown(Keyboard.KEY_2))
	    	 {
	    		 keyy=2;
	    	 }
	    	 if(Keyboard.isKeyDown(Keyboard.KEY_3))
	    	 {
	    		 keyy=3;
	    	 }
	    	 if(Keyboard.isKeyDown(Keyboard.KEY_O))
	    	 {
	    		 zoom = 1.05;
	    		 GL11.glScaled(zoom,zoom,zoom);
	    		 
	    	 }
	    	 if(Keyboard.isKeyDown(Keyboard.KEY_I))
	    	 {
	    		 zoom = 0.95;
	    		 GL11.glScaled(zoom,zoom,zoom);
	    	 }
	    	 
	    	 double movement=1;
	    	 if(Keyboard.isKeyDown(Keyboard.KEY_UP))
	    	 {
	    		 GL11.glTranslated(0,-movement, 0);
	    	 }
	    	 if(Keyboard.isKeyDown(Keyboard.KEY_DOWN))
	    	 {
	    		 GL11.glTranslated(0,movement, 0);
	    	 }
	    	 if(Keyboard.isKeyDown(Keyboard.KEY_RIGHT))
	    	 {
	    		 GL11.glTranslated(movement, 0,0);
	    	 }
	    	 if(Keyboard.isKeyDown(Keyboard.KEY_LEFT))
	    	 {
	    		 GL11.glTranslated(-movement, 0,0);
	    	 }
	    	 if(Keyboard.isKeyDown(Keyboard.KEY_S))
	    	 	 	stop=true;
	    	 
	    	 if(Keyboard.isKeyDown(Keyboard.KEY_R))
	    	 	 	stop=false;
	    	 
	    	 // 3d interactive 
	    	 while (Mouse.next())
	    	 {
	    		    if (Mouse.getEventButtonState()) 
	    		    {
	    		        if (Mouse.getEventButton() == 0) 
	    		        {
	    		        	initial_x=Mouse.getX();
	    		        	initial_y=Mouse.getY();
	    		        	clicked=true;	    		        
	    		        }
	    		    }
	    		    else 
	    		    {
	    		        if (Mouse.getEventButton() == 0) 
	    		        {
	    		        	clicked=false;
	    		        }
	    		    }
	    		    
	    		    if(clicked)
	    		    {
	    		    	float endx = Mouse.getX();
    		        	float endy = Mouse.getY();
	    		    	
	    		    	float zavx = (endx-initial_x);
	    		    	float zavy = (endy-initial_y);
	    		    	
	    		    	
	    		    	
	    		    	
	    		    	float ang = (float) zavx/20;
	    		    
	    		    	if(keyy==1)
	    		    		GL11.glRotated(ang, 0,1,0);
	    		    	if(keyy==2)
	    		    		GL11.glRotated(ang, 1,0,0);
	    		    	if(keyy==3)
	    		    		GL11.glRotated(ang, 0,0,1);
	    		    	
	    		    	
	    		    	initial_x=endx;
	    		    	initial_y=endy;
	    		 	    		    	
	    		    }
	    		    
	    		}
	    	 
	    	 	
	    	 // sphere
	  if(graphvision.choice.getSelectedItem()=="Circular-Force-Free" || graphvision.choice.getSelectedItem()=="Bi-Stress")
	  	 if(graphvision.t3d.isSelected() && graphvision.sphere.isSelected())
	   		 {
	    		 renderSphere((float)graphvision.cent[0], (float)graphvision.cent[1], (float)graphvision.cent[2], (float)graphvision.radi, (float) 0.05, (float) 0.07);
	    	 }
	    	 //end
	    	
	    	
	    	 
    		 if(!graphvision.show_edges.isSelected() && !(graphvision.choice.getSelectedItem()=="MDS") && !(graphvision.choice.getSelectedItem()=="CMD") )
		    		 {
		    			 GL11.glLineWidth(0.005f);
		        		 glBegin(GL_LINES);
		        		 glColor3d(34.0/255, 177.0/255.0, 76.0/255.0);
		        		 glColor3d(0, 0, 0);
			    		 if(sORa && adj_edge)
			    		 {
			    			 for(int i=0;i<links;i++)
				    		 {
			    				 GL11.glVertex3f( (float)pos[0][0][adj[0][i]], (float)pos[0][1][adj[0][i]], (float)pos[0][2][adj[0][i]]);
			    				 GL11.glVertex3f( (float)pos[0][0][adj[1][i]], (float)pos[0][1][adj[1][i]], (float)pos[0][2][adj[1][i]]);
			    			
				    		 }
			    			 
			    		 }
			    		 
			    		 if(sORa && !adj_edge)
			    		 {
			    			 for(int i=0;i<links;i++)
				    		 {
			    				 int v1 = adj[0][i];
					    		 int v2 = adj[1][i];
			    				 
			    				 GL11.glVertex3f( (float)pos[0][0][v1], (float)pos[0][1][v1], (float)pos[0][2][v1]);
			    				 GL11.glVertex3f( (float)pos[0][0][v2], (float)pos[0][1][v2], (float)pos[0][2][v2]);
				    		 }
			    			 
			    		 }
			    		 
			    		 
			    		 if(!sORa && adj_edge)
			    		 {
			    			 for(int i=0;i<links;i++)
				    		 {
			    				 GL11.glVertex3f( (float)pos[0][0][adj[0][i]], (float)pos[0][1][adj[0][i]], (float)pos[0][2][adj[0][i]]);
			    				 GL11.glVertex3f( (float)pos[0][0][adj[1][i]], (float)pos[0][1][adj[1][i]], (float)pos[0][2][adj[1][i]]);
				    		 }
			    		 }
			    		 
			    		 if(!sORa && !adj_edge)
			    		 {
			    			 
			    			/* System.out.println("dsds");
			    			 glColor3d(1, 0, 0);
			    			 for(int i=0;i<v;i++)
				    		 {
			    				 for(int j=0;j<v;j++)
					    		 {
			    				 GL11.glVertex3f( (float)pos[counter][0][i], (float)pos[counter][1][i], 0);
			    				 GL11.glVertex3f( (float)pos[counter][0][j], (float)pos[counter][1][j], 0);
					    		 }
				    		 }
			    			 glColor3d(1, 1, 1);*/
			    			 for(int i=0;i<links;i++)
				    		 {
			    				 int v1 = adj[0][i];
					    		 int v2 = adj[1][i];
			    				/*if(CTD.met[counter][v1]==1 && CTD.met[counter][v2]==1)
			    				{*/
			    					
	    						 GL11.glVertex3f( (float)pos[counter][0][v1], (float)pos[counter][1][v1], (float)pos[counter][2][v1]);
			    				 GL11.glVertex3f( (float)pos[counter][0][v2], (float)pos[counter][1][v2], (float)pos[counter][2][v2]);	
			    				//}
			    				 //graph.text(12f, "test ", 10f,20f, org.newdawn.slick.Color.blue);
	    						 
				    		 }
			    			
			    		 }
			    		 glEnd();	
		    		 }
		    		 
		    		
		    		 
					if (sORa) {
						for (int i = 0; i < v; i++)
							circle((float) pos[0][0][i], (float) pos[0][1][i], (float) pos[0][2][i], 0.38f, i, col);
					} else {
						for (int i = 0; i < v /*&& CTD.met[counter][i]==1*/ ; i++)
							circle((float) pos[counter][0][i], (float) pos[counter][1][i], (float) pos[counter][2][i], 0.38f, i,
									col);
						org.newdawn.slick.Color col = new org.newdawn.slick.Color(10,10,250) ;
						//text( (float) 120, ""+counter, (float) 10, (float) 10, col);
					}
			    	 
					glFlush();
	    	 	 	
					if (counter < time - 1 && !stop)
						counter++;
	    	 	 	 
	    	 	 	try 
	    	         {
	    	        	 Thread.sleep(graphvision.layout_rate_frame.getValue());
	    	         }
	    	         catch (InterruptedException e) 
	    	         {
	    				 e.printStackTrace();
	    			 }
	    }

	  public static void colorings()
		{
			colors=	Integer.parseInt(graphvision.col_choice.getSelectedItem());
			
			
			if(colors>1)
			{
				String ranges = JOptionPane.showInputDialog(null, "Enter the ranges of the groups, e.g. 1:200-201:400-401:600");
				String r [] = ranges.split("-");
				while(colors!=r.length)
				{
					JOptionPane.showMessageDialog(null, "Wrong ranges!");
					ranges = JOptionPane.showInputDialog(null, "Enter the ranges of the groups, e.g. 1:200-201:400-401:600");
					r = ranges.split("-");
				}
				graph.cs = new int [r.length];
				for(int i=0;i<r.length;i++)
				{
					String h [] = r[i].split(":");
					graph.cs[i]=Integer.parseInt(h[0])-1;
				}
			}
		}

		public static void circle(float x, float y,float z, float size, int n, double [][] col)
		{
			
	    	//****************
			//****************
			glBegin(GL_LINES);
	    	
			
	    	double c1=0,c2=0,c3=0;
	    	
	    	
	    	double cols [][] = new double [3][13];
	    	
	    	cols[0][0]=1;   	cols[1][0]=0;		cols[2][0]=0; //red
	    	cols[0][1]=0;   	cols[1][1]=1;		cols[2][1]=0; //lime
	    	cols[0][2]=0;   	cols[1][2]=0;		cols[2][2]=1; //blue
	    	cols[0][3]=1;   	cols[1][3]=0;		cols[2][3]=1; //magenta
	    	cols[0][4]=0;   	cols[1][4]=1;		cols[2][4]=1; //cyan
	    	
	    	cols[0][5]=0;   	cols[1][5]=0;		cols[2][5]=0; //black
	    	cols[0][6]=0.5; 	cols[1][6]=0;		cols[2][6]=0.5; //purple
	    	cols[0][7]=0.5; 	cols[1][7]=0;		cols[2][7]=0; //maroon
	    	cols[0][8]=0.5; 	cols[1][8]=0.5;		cols[2][8]=0; //olive
	    	cols[0][9]=1;   	cols[1][9]=1;		cols[2][9]=0; //yellow
	    	cols[0][10]=0;  	cols[1][10]=0.5;	cols[2][10]=0.5; //teal
	    	cols[0][11]=0;  	cols[1][11]=0;		cols[2][11]=0.5; //navy
	    	cols[0][12]=0.5;  	cols[1][12]=0.5;	cols[2][12]=0.5; //black
	    	
	    	
	    	if(colors>1)
	    	{
	    		for(int i=0;i<graph.cs.length;i++)
	    		{
	    			if(i<graph.cs.length-1)
		    			if(n>=graph.cs[i] && n<graph.cs[i+1])
		    			{
		    				c1=cols[0][i];
			    	    	c2=cols[1][i];
			    	    	c3=cols[2][i];
			    	    	break;
		    			}
	    			if(i==graph.cs.length-1)
	    			{
	    				if(n>=graph.cs[i])
		    			{
		    				c1=cols[0][i];
			    	    	c2=cols[1][i];
			    	    	c3=cols[2][i];
			    	    	break;
		    			}
	    			}
	    		}
	    	}
	    	else
	    	{
		    	c1=1.0;
	    	   	c2=0.0;
	    	   	c3=0.0;
		    }
		    	if(n==111 || n==93 || n==103 || n==66){
	    	c1=0.0;
    	   	c2=0.0;
    	   	c3=1.0;
		    	}
	    	//----------------
	        	//glColor3d(col[0][n],col[1][n],col[2][n]);
		    	c1=0.0;
	    	   	c2=1.0;
	    	   	c3=0.0;
		    	glColor3d(c1,c2,c3);
	        	float interval =(float) (10.510/(float)Math.sqrt(v));
	        	//interval =0.1f;
	        	if(!graphvision.t3d.isSelected())
	        	{
		        	for(double i=0;i<628;i++)
		        	{
			        	GL11.glVertex3f( x, (float)(y), z);	        	
			        	GL11.glVertex3f( x+(float)Math.cos(i/100)*interval*5, y+(float)(Math.sin(i/100)*interval*5),z);
		        	}
	        	}
	        	else
	        	{
		        	for(double i=0;i<12;i++)
		        	{
		        		for(double j=0;j<12;j++)
			        	{
			        		GL11.glVertex3f( x, (float)(y), z);
				        	GL11.glVertex3f(x + interval*5*(float)Math.cos(i/2)*(float)Math.sin(j/2),
			        				y +	interval*50*(float)Math.sin(i/2)*(float)Math.sin(j/2),
			        				z + interval*50*(float)Math.cos(j/2));
			        	}
		        	}
	        	}
	        	
	        	glEnd();
	        
	        glFlush();
		}
		
		private static void renderSphere(float x, float y, float z, float rad, float slic, float smooth) 
		{
			for(float i=0;i<Math.PI*2;i += slic)
			{
				for(float j=0;j<Math.PI*1;j += smooth)
				{
					GL11.glLineWidth(1.5f);
	       		 	glBegin(GL_LINES);
	       		 	glColor3d(0, 0, 0);
	       		 	
		    		GL11.glVertex3f((float)(rad*Math.cos(j)*Math.cos(i)+x),
		    				(float) (rad*Math.sin(i)+z),
		    				(float) (rad*Math.sin(j)*Math.cos(i)+y)
		    				
		    				);
		    		
		    		GL11.glVertex3f((float)(rad*Math.cos(j)*Math.cos(i+slic)+x),
		    				(float) (rad*Math.sin(i+slic)+z),
		    				(float) (rad*Math.sin(j)*Math.cos(i+slic)+y)
		    				);
		    		
		    		GL11.glEnd();
				}	
			}
		}
		
		static void snapshot(String name)
		{
			System.out.println(name);
			GL11.glReadBuffer(GL11.GL_FRONT);
			int width = Display.getDisplayMode().getWidth();
			int height= Display.getDisplayMode().getHeight();
			int bpp = 4; // Assuming a 32-bit display with a byte each for red, green, blue, and alpha.
			ByteBuffer buffer = BufferUtils.createByteBuffer(width * height * bpp);
			GL11.glReadPixels(0, 0, width, height, GL11.GL_RGBA, GL11.GL_UNSIGNED_BYTE, buffer );
			
			//*********
			/*DateFormat df = new SimpleDateFormat("dd/MM/yy HH:mm:ss");
			String fname= df.toString();*/
			 File dir1 = new File("Last_reults");
			    dir1.mkdir();
			 File dir2 = new File(dir1,"pics");
			    dir2.mkdir();
			    File file = new File (dir2,name+".jpg");
				   
			//*********
			
			//File file = new File (name+".jpg");
			String format = "JPG";
			BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
			for(int x = 0; x < width; x++)
			for(int y = 0; y < height; y++)
			{
				int i = (x + (width * y)) * bpp;
				int r = buffer.get(i) & 0xFF;
				int g = buffer.get(i + 1) & 0xFF;
				int b = buffer.get(i + 2) & 0xFF;
				image.setRGB(x, height - (y + 1), (0xFF << 24) | (r << 16) | (g << 8) | b);
			}
			
			try
			{
				ImageIO.write(image, format, file);
			} 
			catch (IOException e) 
			{ 
				e.printStackTrace(); 
			}
			
		}
		
		 public static void text(float size, String text, float x, float y, org.newdawn.slick.Color col) throws FontFormatException
		    {
		    	
		    	TrueTypeFont font;
		    	TrueTypeFont font2 = null;
		    	
		    	  Font awtFont = new Font("Times New Roman", Font.BOLD, (int)size);
		          font = new TrueTypeFont(awtFont, antiAlias);
		           
		          // load font from file
		          try {
		              InputStream inputStream = ResourceLoader.getResourceAsStream("myfont.ttf");
		               
		              Font awtFont2 = Font.createFont(Font.TRUETYPE_FONT, inputStream);
		              awtFont2 = awtFont2.deriveFont(size); // set font size
		              font2 = new TrueTypeFont(awtFont2, antiAlias);
		               
		          } catch (Exception e) {
		              e.printStackTrace();
		          }
		          
		         //font2.drawString(x, y,text, col);
		          glEnable(GL11.GL_BLEND);
		    	  	
		          glBlendFunc(GL11.GL_SRC_ALPHA, GL11.GL_ONE_MINUS_SRC_ALPHA);
		          font.drawString(x, y, text, col);
		          glDisable(GL11.GL_BLEND);
		    }


		@Override
		public void method1(int a) {
			// TODO Auto-generated method stub
			
		}
}
