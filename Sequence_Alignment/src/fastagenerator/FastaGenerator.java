package fastagenerator;

public class FastaGenerator {
	
	static void matrixgenerate()
	{
		char[] aminoAcids = {'A', 'C', 'D','E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y' , '*'};
		System.out.print("  ");
		for(int i=0;i<22;i++)
		{
			for(int j=0;j<21;j++)
			{
				if(i==0)
				{
					System.out.print(aminoAcids[j]+"  ");
				}
				else
				{
					if(j==0)
					{
						System.out.print(aminoAcids[i-1] + "  ");
					}
					if(i==(j+1))
					{
						System.out.print(1+"  ");
					}
					else
					{
						System.out.print(-1 + "  ");
					}
				}
			}
			System.out.println();
		}
	}
	public static void main(String []args)
	{
		/*
		char[] aminoAcids = {'A', 'R', 'N', 'D', 'C', 
			'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'};
		System.out.println(">seq 1");
		for(int i=0;i<1325;i++)
		{
			int randAA = (int)(Math.random()*20+0);
			System.out.print(aminoAcids[randAA]);
		}
		System.out.println();
		System.out.println(">seq 2");
		for(int j=0;j<1991;j++)
		{
			int randAA = (int)(Math.random()*20+0);
			System.out.print(aminoAcids[randAA]);
		}
		
		System.out.println();
		*/
		matrixgenerate();
	}
}
