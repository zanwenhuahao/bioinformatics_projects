package seqalignment;

import java.io.*;
//import java.util.ArrayList;

public class SeqAlign {
	 //Inputs
	
	static String seqOne;
	static String seqTwo;
	static String [][]scorematrix = new String[22][22];
	
	 //Algorithm needs
	static int [][]globalTable;
	static int [][]localTable;
	
	static int local_max_row;
	static int local_max_col;
	
	//Outputs
	static String globalOne;
	static String globalTwo;
	static int globalScore;
	
	static String localOne;
	static String localTwo;
	static int localScore;
	
	//Helper function to get scorematrix score
	static int getscore(String hor, String vert)
	{
		int theRow = 0;
		int theCol = 0;
		
		//System.out.println("DEBUG: " + hor + vert);
		
		for(int i = 0; i<scorematrix.length;++i)
		{
			if (scorematrix[i][0].equals(hor))
			{
				theRow = i;
				break;
			}
		}
		for(int j = 0; j<scorematrix[0].length; ++j)
		{
			if (scorematrix[0][j].equals(vert))
			{
				theCol = j;
				break;
			}
		}
		/*System.out.println("DEBUG: " + theRow);
		System.out.println("DEBUG: " + theCol);
		System.out.println("DEBUG: Scorematrix: "+ scorematrix[theRow][theCol-1]);*/
		int scoreval = Integer.parseInt(scorematrix[theRow][theCol]);
		//System.out.println("DEBUG: Score: "+ scoreval);
		return scoreval;
	}
	
	//Algorithm - Global Alignment
	static void initialization_global()
	{
		int numRows = seqOne.length();
		int numCols = seqTwo.length();
		globalTable = new int [numRows+1][numCols+1];
		globalTable[0][0] = 0;
		int score = 0;
		
		for(int i=1; i<=numRows; i++)
		{
			//System.out.println("DEBUG: " + ""+seqOne.charAt(i-1));
			score = getscore(""+seqOne.charAt(i-1),"*");
			globalTable[i][0] = i*score;
		}
		for(int j=1; j<=numCols; j++)
		{
			//System.out.println("DEBUG: " + ""+seqTwo.charAt(j-1));
			score = getscore("*", ""+seqTwo.charAt(j-1));
			globalTable[0][j] = j*score;
		}
	}
	
	static void alignment_global()
	{
		int numRows = seqOne.length();
		int numCols = seqTwo.length();
		
		for(int i=1;i<=numRows;i++)
		{
			for(int j=1;j<=numCols;j++)
			{
				//Note that I use i-1 and j-1 to the function getscore to mean the character at seqOne/seqTwo
				//'s index position, but they are on position i and j on the dp table.
				int diag = globalTable[i-1][j-1]+getscore(""+seqOne.charAt(i-1), ""+seqTwo.charAt(j-1));
				int up = globalTable[i-1][j]+getscore(""+seqOne.charAt(i-1),"*");
				int left = globalTable[i][j-1]+getscore("*", ""+seqTwo.charAt(j-1));
				int maxij = Math.max(Math.max(diag, left), up); //Gets the maximum of the 3
				globalTable[i][j] = maxij;
			}
		}
	}
	
	static void backtracking_global()
	{
		int row = seqOne.length();
		int col = seqTwo.length();
		globalOne = "";
		globalTwo = "";
		globalScore = globalTable[row][col];
		
		while (row>0 || col>0)
		{
			if(row == 0)
			{
				globalTwo = ""+seqTwo.charAt(col-1)+globalTwo;
				globalOne = "-" + globalOne;
				col--;
			}
			else if(col == 0)
			{
				globalOne = ""+seqOne.charAt(row-1) + globalOne;
				globalTwo = "-" + globalTwo;
				row--;
			}
			else
			{
				if(globalTable[row][col] == globalTable[row-1][col] + getscore(""+seqOne.charAt(row-1),"*"))
				{
					globalOne = ""+seqOne.charAt(row-1) + globalOne;
					globalTwo = "-" + globalTwo;
					row--;
				}
				else if(globalTable[row][col] == globalTable[row][col-1] + getscore("*",""+seqTwo.charAt(col-1)))
				{
					globalTwo = ""+seqTwo.charAt(col-1) + globalTwo;
					globalOne = "-" + globalOne;
					col--;
				}
				else
				{
					globalOne = ""+seqOne.charAt(row-1) + globalOne;
					globalTwo = ""+seqTwo.charAt(col-1) + globalTwo;
					row--;
					col--;
				}
			}
				
		}

	}
	
	//Algorithm - local alignment
	static void initialization_local()
	{
		int numRows = seqOne.length();
		int numCols = seqTwo.length();
		localTable = new int [numRows+1][numCols+1];
		localTable[0][0] = 0;
		int score = 0;
		
		for(int i=1; i<=numRows; i++)
		{
			//System.out.println("DEBUG: " + ""+seqOne.charAt(i-1));
			//score = getscore(""+seqOne.charAt(i-1),"*");
			localTable[i][0] = score;
		}
		for(int j=1; j<=numCols; j++)
		{
			//System.out.println("DEBUG: " + ""+seqTwo.charAt(j-1));
			//score = getscore("*", ""+seqTwo.charAt(j-1));
			localTable[0][j] = score;
		}
	}
	
	static void alignment_local()
	{
		int numRows = seqOne.length();
		int numCols = seqTwo.length();
		
		for(int i=1;i<=numRows;i++)
		{
			for(int j=1;j<=numCols;j++)
			{
				//Note that I use i-1 and j-1 to the function getscore to mean the character at seqOne/seqTwo
				//'s index position, but they are on position i and j on the dp table.
				int diag = localTable[i-1][j-1]+getscore(""+seqOne.charAt(i-1), ""+seqTwo.charAt(j-1));
				int up = localTable[i-1][j]+getscore(""+seqOne.charAt(i-1),"*");
				int left = localTable[i][j-1]+getscore("*", ""+seqTwo.charAt(j-1));
				int maxij = Math.max(Math.max(Math.max(diag, left), up), 0); //Gets the maximum of the 3
				localTable[i][j] = maxij;
			}
		}
	}
	
	static int findmax(int[][] dptable)
	{
		int max_local_score = 0;
		for(int i=0;i<dptable.length;i++)
		{
			for(int j=0;j<dptable[0].length;j++)
			{
				if(dptable[i][j]>max_local_score)
				{
					max_local_score=dptable[i][j];
					local_max_row = i;
					local_max_col = j;
				}
			}
		}
		//System.out.println("DEBUG: maximum local score: " + max_local_score);
		//System.out.println("DEBUG:LocalRow: "+ local_max_row);
		//System.out.println("DEBUG:LocalCol: "+ local_max_col);
		return max_local_score;
	}
	
	static void backtracking_local()
	{
		localOne = "";
		localTwo = "";
		localScore = findmax(localTable);
		int row = local_max_row;
		int col = local_max_col;
		//System.out.println("DEBUG:LocalRow: "+ row);
		//System.out.println("DEBUG:LocalCol: "+ col);
		
		while (row>0 || col>0)
		{
			if (localTable[row][col] == 0)
			{
				break;
			}
			else
			{
				if(localTable[row][col] == localTable[row-1][col] + getscore(""+seqOne.charAt(row-1),"*"))
				{
					localOne = ""+seqOne.charAt(row-1) + localOne;
					localTwo = "-" + localTwo;
					row--;
				}
				else if(localTable[row][col] == localTable[row][col-1] + getscore("*",""+seqTwo.charAt(col-1)))
				{
					localTwo = ""+seqTwo.charAt(col-1) + localTwo;
					localOne = "-" + localOne;
					col--;
				}
				else
				{
					localOne = ""+seqOne.charAt(row-1) + localOne;
					localTwo = ""+seqTwo.charAt(col-1) + localTwo;
					row--;
					col--;
				}
			}
				
		}

	}
	
	//Debug functions
	//Print scoring matrix to check it's read in properly
	static void printscorematrix(String[][] mat)
	{
		for(int i=0;i<mat.length;i++)
		{
			for(int j=0;j<mat[0].length;j++)
			{
				System.out.print(mat[i][j]+ " ");
			}
			System.out.println();
		}
	}
	static void printdptable(int [][] dptable)
	{
		for(int i=0;i<dptable.length;i++)
		{
			for(int j=0;j<dptable[0].length;j++)
			{
				System.out.print(dptable[i][j]+" ");
			}
			System.out.println();
		}
	}
	
	static void printsequences(String seqlist)
	{
			System.out.println(seqlist);
	}
	
	static void printscore(String which)
	{
		if(which.equals("global"))
		{
			System.out.println(globalScore);
		}
		else
		{
			System.out.println(localScore);
		}
	}
	
	//Reading in files
	public static void main(String[] args) throws IOException 
	{
		File FASTAseq = new File("longfasta.txt");
		File ScoringMatrix = new File("A2_ScoringMatrixSample.txt");
		
		BufferedReader seqbr = new BufferedReader(new FileReader(FASTAseq));
		BufferedReader matrixbr = new BufferedReader(new FileReader(ScoringMatrix));
		
		//Reading the two protein sequences
		String line;
		int seqnum = 0;
		seqOne = "";
		seqTwo = "";
		
		while((line = seqbr.readLine())!=null)
		{
			if((line.trim().length() > 0) && (line.trim().substring(0, 1).equals(">"))) //If it's not a header line
			{
				seqnum++;
			}
			else
			{	
				if(seqnum == 1)
				
				{
					seqOne = seqOne + line.trim();
				}
				else if(seqnum == 2)
				{
					seqTwo = seqTwo + line.trim();
				}
				else
				{
					continue;
				}
			}
		}
		
		boolean headrowflag = true;
		//Reading the Scoring Matrix
		int linenum = 1;
		while((line = matrixbr.readLine())!=null)
		{
			if((line.trim().length()>0)
					&& !(line.substring(0, 1).equals("#")))
			{
				if(headrowflag)
				{	
					String[] linesep = line.trim().split("\\s+");
					scorematrix[0][0] = " ";
					for(int i=0;i<linesep.length;i++)
					{
						scorematrix[0][i+1] = linesep[i];	
					}
					headrowflag = false;
				}
				else
				{
					String[]linesep = line.trim().split("\\s+");
					for(int j=0;j<linesep.length;j++)
					{
						scorematrix[linenum][j] = linesep[j];
					}
					linenum++;
				}
			}
			//System.out.println("DEBUG: Horizontal row: " + scorematrix[0][1]);
		}
		//Testing purpose
		//printscorematrix(scorematrix);
		//printsequences(seqOne);
		//printsequences(seqTwo);
		initialization_global();
		//printdptable(globalTable);
		alignment_global();
		//printdptable(globalTable);
		backtracking_global();
		printscore("global");
		printsequences(globalOne);
		printsequences(globalTwo);
		initialization_local();
		//printdptable(localTable);
		alignment_local();
		//printdptable(localTable);
		backtracking_local();
		printscore("local");
		printsequences(localOne);
		printsequences(localTwo);
		
		
	}
}
