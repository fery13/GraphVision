import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Scanner;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.JTextField;


import java.awt.*;

public class graphvision extends JFrame{
	
	private static int groups=5;
	
	static JPanel contentPane;
	static JFileChooser browse = new JFileChooser();
	
	static JButton loadbtm 			 = new JButton("Load");
	static JButton btnbrws 			 = new JButton("Brows");
	
	static JLabel lable1 = new JLabel();
	static JLabel lable2 = new JLabel();
	static JLabel info1 = new JLabel("GraphVision. 	Computer Science and Information Systems Department, 	University of Limerick");
	static JLabel info2 = new JLabel("By: Farshad Ghassemi Toosi");
	static JLabel info3 = new JLabel("Email: farshad.toosi@ul.ie");
	
	static JCheckBox animation = new JCheckBox("Animation");
	
	static JCheckBox t3d = new JCheckBox("3D");
	
	public static JCheckBox show_edges = new JCheckBox("Do not show edges");
	public static JCheckBox MorL = new JCheckBox("Read as Matrix!");
	public static JCheckBox SaveOne = new JCheckBox("Save the last image!");
	public static JCheckBox SaveAll = new JCheckBox("Save all images if animation!");
	public static JCheckBox sphere = new JCheckBox("Surroundign by Sphere !");
	
	static public Choice choice = new Choice();
	
	static public Choice col_choice = new Choice();
		
	static public JSlider layout_rate_frame = new JSlider(JSlider.VERTICAL,1, 1000,1);
	
	static String ff;
	static private String file_add="";
	static private double links=0;
	static private int v=0;
	static public int adj [][];
	static public double mat [][];
	
	static private double sim [];
	static private int col [];
	static public double node_degree [];
	
	public static double radi=0;
	public static double cent [] = new double [3];
	
	static private double pos [][][];
	
	static private boolean anim=false;
	static private boolean tree=false;
	static private double round =150;
	static public boolean correct=false;
	static public boolean threeD = false;
	static public boolean initial = false;
	
	static public boolean input_read = false;
	static public long S1 = 0;
	static public long E1 = 0;
	
	// main form to be created here.
	public graphvision() {
		
		
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setBounds(40, 40, 1100, 800);
		contentPane = new JPanel();
		contentPane.setLayout(null);
		setContentPane(contentPane);
		
		col_choice.addItem("1");
		col_choice.addItem("2");
		col_choice.addItem("3");
		col_choice.addItem("4");
		col_choice.addItem("5");
		col_choice.addItem("6");
		col_choice.addItem("7");
		col_choice.addItem("8");
		col_choice.addItem("9");
		col_choice.addItem("10");
		col_choice.addItem("11");
		col_choice.addItem("12");
		col_choice.addItem("13");
		
		
		choice.addItem("fruchterman-reingold");
		choice.addItem("kamada-kawai");
		choice.addItem("Circular-Force-Free");
		choice.addItem("Tree1");
		choice.addItem("Tree2");
		choice.addItem("MDS");
		choice.addItem("Bi-Stress");
		choice.addItem("CMD");
		choice.addItem("CTD");
		choice.addItem("VNM");
		
		
		
		
		
		
		
		//BROWS Bottom
		loadbtm.setBounds(335, 250, 89, 23);
		contentPane.add(loadbtm);
		
		
		btnbrws.setBounds(335, 228, 89, 23);
		contentPane.add(btnbrws);
		
		btnbrws.setEnabled(true);
		loadbtm.setEnabled(false);
		
		
		
		
		lable1.setBounds(63, 187, 244, 14);
		contentPane.add(lable1);
		
		
		lable2.setBounds(10, 175, 284, 14);
		contentPane.add(lable2);
		
		info1.setBounds(10, 675, 800, 20);
		contentPane.add(info1);
		info1.setForeground(Color.blue);
		
		info2.setBounds(10, 695, 800, 20);
		contentPane.add(info2);
		info2.setForeground(Color.blue);
		
		info3.setBounds(10, 715, 800, 20);
		contentPane.add(info3);
		info3.setForeground(Color.blue);
		
		
		
		choice.setBounds(318, 184, 156, 20);
		contentPane.add(choice);
		
	
		col_choice.setBounds(588, 184, 56, 20);
		contentPane.add(col_choice);
		
		JLabel colLable = new JLabel("No.Colors");
		colLable.setBounds(508, 184, 156, 20);
		contentPane.add(colLable);
		
		
		show_edges.setBounds(266, 57, 224, 23);
		contentPane.add(show_edges);
		
		
		
		
		
		
		animation.setBounds(266, 38, 97, 23);
		contentPane.add(animation);
		
		t3d.setBounds(266, 78, 97, 23);
		contentPane.add(t3d);
		
		MorL.setBounds(266, 118, 197, 23);
		contentPane.add(MorL);
		
		SaveOne.setBounds(6, 158, 197, 23);
		contentPane.add(SaveOne);
		
		SaveAll.setBounds(6, 198, 197, 23);
		contentPane.add(SaveAll);
		
		sphere.setBounds(6, 258, 197, 23);
		contentPane.add(sphere);
		
		//JSlider slider = new JSlider();
		layout_rate_frame.setBounds(10, 42, 33, 86);
		contentPane.add(layout_rate_frame);
		
		JLabel lblNewLabel = new JLabel("Animation Speed");
		lblNewLabel.setBounds(10, 1, 175, 54);
		contentPane.add(lblNewLabel);
		
		
	
		
		//button thingssssssssssss

		btnbrws.addActionListener(new ActionListener() 
		{
			public void actionPerformed(ActionEvent arg0) 
			{
				File adjacency_matrix = null;
				browse.showOpenDialog(getComponent(0));
				adjacency_matrix=browse.getSelectedFile();
				file_add= adjacency_matrix.getAbsolutePath();
				System.out.println(file_add+"     KK");
				lable1.setText("Adjacency Matrix:   "+adjacency_matrix.getName());
				
				//******loading
				
				correct=true;
				try 
				{	
					if(graphvision.choice.getSelectedItem().equals("MDS") || graphvision.choice.getSelectedItem().equals("CMD"))
						sim= read_MDS(file_add);
					else 
					{
						S1 = System.nanoTime();
					
							adj=read_adj(file_add);
						E1 = System.nanoTime();
					}
				} 
				catch (FileNotFoundException e) 
				{
					e.printStackTrace();
					correct=false;
				}
				if(!correct)
				{
					JOptionPane.showMessageDialog(null, "Fix the input file and try again ... !");
				}
			
				if(correct)
					loadbtm.setEnabled(true);
				
				
				
			}
		}
		);
		
		//END  LOAD BOTTOM
		
		
	

		loadbtm.addActionListener(new ActionListener() 
		{
			
			public void actionPerformed(ActionEvent arg0) 
			{	
				
				if(animation.isSelected())
					anim=true;
				else
					anim=false;
				
				if(t3d.isSelected())
					threeD=true;
			
				
				
			    if(graphvision.choice.getSelectedItem()=="fruchterman-reingold")
				try {
						round = 10;
						pos=FR.force_based_layout(adj, v, anim,links, v*round,threeD);
						btnbrws.setEnabled(true);
						graph FRgraph = new graph( pos,tree, adj, v*round, !anim, v,  links, false,1);
						FRgraph.visualization();
					} 
			     catch (FileNotFoundException e5) 
			     {
							// TODO Auto-generated catch block
							e5.printStackTrace();
						} catch (InterruptedException e5) {
							// TODO Auto-generated catch block
							e5.printStackTrace();
				} catch (FontFormatException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
						
				if(graphvision.choice.getSelectedItem()=="kamada-kawai")
				try {
						pos=KK.kamada_kawai(adj, v, anim,links, v*round);
						btnbrws.setEnabled(true);
						graph KKgraph = new graph( pos,tree, adj,v*round, !anim, v,  links, false,1);
						KKgraph.visualization();
					} 
				catch (FileNotFoundException e5) 
				{
					// TODO Auto-generated catch block
					e5.printStackTrace();
					} catch (InterruptedException e5) {
					// TODO Auto-generated catch block
					e5.printStackTrace();
				} catch (FontFormatException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					
				if(graphvision.choice.getSelectedItem()=="Circular-Force-Free")
					try {
							round = 35; ///// temprorary
							if(v<100)
								round=10;
							if(v>=100 && v<200)
								round=20;
							
							int time=(int) (v*round);
							//time = 8000;
							pos=CFF.Cir_Force_Free(adj, v, anim,links, time, threeD);
							btnbrws.setEnabled(true);
							graph CFFgraph = new graph( pos,tree, adj,time, !anim, v,  links, false,1);
							CFFgraph.visualization();
						} 
					catch (FileNotFoundException e5) 
					{
						// TODO Auto-generated catch block
						e5.printStackTrace();
						} catch (InterruptedException e5) {
						// TODO Auto-generated catch block
						e5.printStackTrace();
					} catch (FontFormatException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
					
				
				if(graphvision.choice.getSelectedItem()=="VNM")
					try {
							round =5; ///// temprorary
							int time=(int) (v*round);
							pos=VNM.Multilevel_centric(adj, v, anim,links, time, threeD);
							btnbrws.setEnabled(true);
							graph VNMgraph = new graph( pos,tree, adj,time, !anim, v,  links, false,1);
							VNMgraph.visualization();
						} 
					catch (FileNotFoundException e5) 
					{
						// TODO Auto-generated catch block
						e5.printStackTrace();
						} catch (FontFormatException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}
				
				
					
				if(graphvision.choice.getSelectedItem()=="Tree1")
					try {
							round=1;
							double time = v*round;
							time = 5;
							long alg1=System.nanoTime();
							pos= Tree.tree_angular_radial1(adj, v, anim,links, time, threeD);
							long alg2=System.nanoTime();
							System.out.println("Running time for:  "+graphvision.choice.getSelectedItem()+":  "+(alg2-alg1));
							btnbrws.setEnabled(true);
							graph Three = new graph( pos,tree, adj,v*round, !anim, v,  links, false,1);
							Three.visualization();
						} 
					catch (FileNotFoundException e5) 
					{
						// TODO Auto-generated catch block
						e5.printStackTrace();
					} catch (FontFormatException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				
				if(graphvision.choice.getSelectedItem()=="Tree2")
					try {
							round=1;
							double time = v*round;
							time = 15;
							long alg1=System.nanoTime();
							pos= Tree.tree_angular_radial2(adj, v, anim,links, time, threeD, initial);
							long alg2=System.nanoTime();
							System.out.println("Running time for:  "+graphvision.choice.getSelectedItem()+":  "+(alg2-alg1));
							btnbrws.setEnabled(true);
							graph Three = new graph( pos,tree, adj,v*round, !anim, v,  links, false,1);
							Three.visualization();
						} 
					catch (FileNotFoundException e5) 
					{
						// TODO Auto-generated catch block
						e5.printStackTrace();
					} catch (FontFormatException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				
				if(graphvision.choice.getSelectedItem()=="MDS")
					try {
							pos= MDS.Multidimentional_forcing(sim, col, v, (double)v/8, threeD);
							btnbrws.setEnabled(true);
							graph MDS = new graph( pos,tree, adj,(int)((double)v/8), !anim, v,  links, false,groups);
							MDS.visualization();
						} 
					catch (FileNotFoundException e5) 
					{
						// TODO Auto-generated catch block
						e5.printStackTrace();
					} catch (FontFormatException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				
				if(graphvision.choice.getSelectedItem()=="CMD")
					try {
							pos= CMD.Multidimentional_forcing(sim, col, v, (double)v/8, threeD);
							btnbrws.setEnabled(true);
							graph CMD = new graph( pos,tree, adj,(int)((double)v/4), !anim, v,  links, false,groups);
							CMD.visualization();
						} 
					catch (FileNotFoundException e5) 
					{
						// TODO Auto-generated catch block
						e5.printStackTrace();
					} catch (FontFormatException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				
				if(graphvision.choice.getSelectedItem()=="CTD")
					try {
						
							int time = v*2;
							pos= CTD.Concentric_drawing_distance(adj,  v,(double)time, threeD, links, anim);
							btnbrws.setEnabled(true);
							graph CTD = new graph( pos,tree, adj,(int)((double)v/4), !anim, v,  links, false,groups);
							CTD.visualization();
						} 
					catch (FileNotFoundException e5) 
					{
						// TODO Auto-generated catch block
						e5.printStackTrace();
					} catch (FontFormatException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				
				
				if(graphvision.choice.getSelectedItem()=="Bi-Stress")
					try {
						
							round=40;
							int time = (int) (round*v);
							time=150;
							pos= BiStress.Bi_Stress(adj, v, anim,links, time,threeD);
							btnbrws.setEnabled(true);
							graph BiStress = new graph( pos,tree, adj,time, !anim, v,  links, false,groups);
							BiStress.visualization();
						} 
					catch (FileNotFoundException e5) 
					{
						// TODO Auto-generated catch block
						e5.printStackTrace();
					} catch (FontFormatException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					} catch (InterruptedException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				
				
				
				
				
					
					
			}
		}
		);
		
	}
		
	
	
		
	public static double [][] read_Matrix(String filename) throws FileNotFoundException
	{
		File aFile = new File(filename);
		Scanner in = new Scanner (aFile);
		ArrayList <String> lines = new ArrayList <String>();
		
		lines.clear();
		
		while(in.hasNextLine())
		{
			String line = in.nextLine();
			lines.add(line);
		}
		v=lines.size();
		node_degree= new double [v];
		mat= new double[v][v];
		links=0;
		
		for(int i=0;i<(lines.size());i++)
		{
				String [] line = lines.get(i).split("\\s+");
				if(line.length!=lines.size())
				{
					JOptionPane.showMessageDialog(null, "Matrix is wrong  ... !");
					break;
				}
				String [] line_deg = (lines.get(i)+" ").split("1");
				node_degree[i]=line_deg.length-1;
				for(int j=i+1;j<(lines.size());j++)
				{
					if(line[j].equals("1"))
					{
						
						mat[i][j]=1;
						mat[j][i]=1;
						
						links++;
						
					}
					else
					{
						mat[i][j]=0;
						mat[j][i]=0;
					}
				}
			}
		int cf=0;
		adj = new int [2][(int) links];
		for(int i=0;i<v;i++)
		{
			for(int j=i+1;j<v;j++)
			{
				if(mat[i][j]==1)
				{
					adj[0][cf]=i;
					adj[1][cf]=j;
					cf++;
				}
			}
		}
		
		in.close();
		return mat;
	}
	
	public static double [][] read_Matrix_pure(String filename) throws FileNotFoundException
	{
		File aFile = new File(filename);
		Scanner in = new Scanner (aFile);
		ArrayList <String> lines = new ArrayList <String>();
		
		lines.clear();
		
		while(in.hasNextLine())
		{
			String line = in.nextLine();
			lines.add(line);
		}
		v=lines.size();
		node_degree= new double [v];
		mat= new double[v][v];
		links=0;
		
		for(int i=0;i<(lines.size());i++)
		{
				String [] line = lines.get(i).split("\\s+");
				if(line.length!=lines.size())
				{
					JOptionPane.showMessageDialog(null, "Matrix is wrong  ... !");
					break;
				}
				
				for(int j=i+1;j<(lines.size());j++)
				{
						mat[i][j]=Double.parseDouble(line[j]);
						mat[j][i]=Double.parseDouble(line[j]);
						
					
				}
			}
		
		
		in.close();
		return mat;
	}
	
	
	public static int [][] read_adj(String filename) throws FileNotFoundException
	{
		links=0;
		correct=true;
		File aFile = new File(filename);
		int [][] edges = null;
		Scanner in = null;
		ArrayList<String> lines = new ArrayList<String>();
		
		/*File real = new File("Nw.txt");
		PrintWriter out = new PrintWriter(real);*/
		
		if (!aFile.exists())
			JOptionPane.showMessageDialog(null,filename + " was not found");
		else
		{
			in = new Scanner(aFile);
			String lineFromFile="";
			
			while(in.hasNext())
			{
				lineFromFile=in.nextLine();
				lines.add(lineFromFile);
				links++;
			}

			edges = new int [2][(int) links];
			int min=0;
			
			
			
			for(int i=0;i<lines.size();i++)
			{
				String [] temp = lines.get(i).split(" ");
				if(temp.length!=2)
				{
					correct=false;
					break;
				}
				
				try
				{
					if(Integer.parseInt(temp[0])!=Integer.parseInt(temp[1]))
					{
						edges[0][i]=Integer.parseInt(temp[0]);
						edges[1][i]=Integer.parseInt(temp[1]);
						//out.println(edges[0][i]+" "+edges[1][i]);
					}
				}
				catch (NumberFormatException ex )
			    {
					correct=false;
					break;
			    }
				if(min<edges[0][i])
					min=edges[0][i];
				
				if(min<edges[1][i])
					min=edges[1][i];
				
			}	
			//out.close();
			v=min+1;
			node_degree = new double [v];
			for(int i=0;i<links;i++)
			{
				node_degree[edges[0][i]]++;
				node_degree[edges[1][i]]++;
			}
			
		}

		in.close();
		return edges;
	}
	
	public static double [][] read_adjD(String filename) throws FileNotFoundException
	{
		links=0;
		correct=true;
		File aFile = new File(filename);
		double [][] edges = null;
		Scanner in = null;
		ArrayList<String> lines = new ArrayList<String>();
		if (!aFile.exists())
			JOptionPane.showMessageDialog(null,filename + " was not found");
		else
		{
			in = new Scanner(aFile);
			String lineFromFile="";
			
			while(in.hasNext())
			{
				lineFromFile=in.nextLine();
				lines.add(lineFromFile);
				links++;
			}

			edges = new double [2][(int) links];
			double min=0;
			
			for(int i=0;i<lines.size();i++)
			{
				String [] temp = lines.get(i).split(" ");
				if(temp.length!=2)
				{
					correct=false;
					break;
				}
				
				try
				{
					edges[0][i]=Integer.parseInt(temp[0]);
					edges[1][i]=Integer.parseInt(temp[1]);		
					
					
					
					
				}
				catch (NumberFormatException ex )
			    {
					correct=false;
					break;
			    }
				if(min<edges[0][i])
					min=edges[0][i];
				
				if(min<edges[1][i])
					min=edges[1][i];
				
			}	
			
			
			
			v=(int) (min+1);
			node_degree = new double [v];
			for(int i=0;i<links;i++)
			{
				node_degree[(int) edges[0][i]]++;
				node_degree[(int) edges[1][i]]++;
			}
		}

		in.close();
		return edges;
	}
	
	
	public static double [] read_MDS(String filename) throws FileNotFoundException
	{
		File aFile = new File(filename);
		Scanner in = new Scanner (aFile);
		ArrayList <String> lines = new ArrayList <String>();
		
		lines.clear();
		
		while(in.hasNextLine())
		{
			String line = in.nextLine();
			lines.add(line);
		}
		
		
		v=lines.size();
		
		double [] similarities = new double [(v*v-v)/2];
		
		if(lines.get(lines.size()-1).equals("-1"))
			input_read=false;
		else
			input_read=true;
		
		
		
		if(input_read)
		{
			for(int i=0;i<(lines.size());i++)
			{
				String [] line = lines.get(i).split("\t");
				
				
				
				if(line.length!=lines.size())
				{
					JOptionPane.showMessageDialog(null, "Matrix is wrong  ... !");
					System.out.println(i);
					input_read=false;
					break;
				}
				for(int j=i+1;j<line.length;j++)
				{
					int index = i*(v-1)+j-(i*i+i)/2-1;
					double mn = Math.abs(Double.parseDouble(line[j]));
					similarities[index]=mn;
				}
			}
		}
		correct = input_read;
		
		return similarities;
	}
	
	
	
	
	
	public static void main(String [] args)
	{
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				try {
		
		graphvision frame = new graphvision();
		frame.setVisible(true);
		
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
		});
	}

}